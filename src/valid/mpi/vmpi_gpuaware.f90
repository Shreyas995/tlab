! vmpi_gpuaware.f90
! ============================================================================
! Phase-0 validation harness for "Option 2": GPU-aware MPI (GPUDirect RDMA)
! for the INTER-NODE leg of the K-transpose.
!
! It mirrors production exactly so a PASS ports straight into fabricdirect-K:
!   - decomposition NPRO_I=6 x NPRO_K=8 = 48 ranks, 2 nodes x 24 ranks/node
!   - world rank = pro_k*NPRO_I + pro_i  (Cartesian dims=(NPRO_K,NPRO_I))
!   - K traffic on fabric_comm_k = MPI_Comm_split(WORLD, color=pro_i, key=pro_k)
!   - the strided GPU pack a(pk*nlines + i*npage + j) -> contiguous, identical
!     to TLabMPI_Trp_ExecK_Forward_Real
!   - source buffers are plain `allocate` under !$omp requires unified_shared_memory,
!     exactly like wrk_mpi_dp / the flow arrays in the main code.
!
! Correctness encode: a(idx) = world_rank*ENC + idx. Sender packs the chunk it
! owes destination Q at source offset Q*nlines. Receiver Q (from sender world pw)
! therefore expects  recv = pw*ENC + (Q*nlines + i*npage + j + 1)  -- exact in dp.
! Mismatches are split inter-node vs intra-node and reduced to rank 0.
!
! ===> THE BLOCK THAT PORTS TO PRODUCTION is the METH_PACK case below:
!      GPU !$omp target strided pack into sendbuf, optional hipDeviceSynchronize,
!      then MPI_Isend(sendbuf, ...). That replaces the CPU pack + Isend currently
!      in fabricdirect-K's inter-node peers.
!
! Modes (argv1); the PBS runs each as a separate, timeout-guarded mpirun:
!   1  PACK  unified send buf  + flush     -> GATE: does GPU-aware Isend deliver inter-node?
!   2  PACK  unified send buf  NO flush    -> is the hipDeviceSynchronize flush required?
!   3  PACK  device send buf   + flush     -> is omp_target_alloc needed, or is unified enough?
!   4  STRIDED MPI_TYPE_VECTOR of GPU `a`  -> can we skip the pack (send strided GPU buffer)?
!   5  BANDWIDTH: CPU-pack baseline vs GPU-pack GPU-aware, inter-node peers, production size.
! ============================================================================
module gpuaware_mod
    use mpi_f08
    use iso_c_binding
    use omp_lib
    implicit none
    !$omp requires unified_shared_memory

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: NPRO_I = 6, NPRO_K = 8, RANKS_PER_NODE = 24
    real(dp), parameter :: ENC = 1.0d9            ! value = world*ENC + src_index (exact: < 2^53)

    integer, parameter :: METH_PACK = 1, METH_STRIDED = 2, METH_CPUPACK = 3

    interface
        function hipDeviceSynchronize() bind(C, name='hipDeviceSynchronize') result(r)
            use iso_c_binding
            integer(c_int) :: r
        end function
    end interface

    integer :: ims_rank, ims_nprocs, ims_pro_i, ims_pro_k, my_node
    type(MPI_Comm) :: fabric_comm_k

contains

    subroutine setup_topology()
        type(MPI_Comm) :: comm_xz
        integer :: dims(2), coord(2), ierr
        logical :: period(2)
        dims(1) = NPRO_K; dims(2) = NPRO_I; period = .true.
        call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, .false., comm_xz, ierr)
        call MPI_Cart_coords(comm_xz, ims_rank, 2, coord, ierr)
        ims_pro_k = coord(1); ims_pro_i = coord(2)
        ! Same MPI_COMM_WORLD split production uses (color=pro_i,key=pro_k): no Cartesian taint,
        ! local rank in fabric_comm_k = ims_pro_k.
        call MPI_Comm_split(MPI_COMM_WORLD, ims_pro_i, ims_pro_k, fabric_comm_k, ierr)
        my_node = ims_rank / RANKS_PER_NODE
    end subroutine

    integer function peer_world(pk)        ! local K-rank pk -> world rank
        integer, intent(in) :: pk
        peer_world = pk*NPRO_I + ims_pro_i
    end function
    integer function peer_node(pk)
        integer, intent(in) :: pk
        peer_node = (pk*NPRO_I + ims_pro_i) / RANKS_PER_NODE
    end function

    subroutine report_env()
        ! Confirm the !$omp target regions actually run on the GPU (not a host fallback),
        ! so an M1 PASS cannot be a false positive from host-resident buffers.
        logical :: on_host
        on_host = .true.
        !$omp target map(tofrom: on_host)
        on_host = omp_is_initial_device()
        !$omp end target
        if (ims_rank == 0) &
            write(*,'(a,i0,a,l1,a)') '[ENV] omp_get_num_devices=', omp_get_num_devices(), &
                '   target_on_host=', on_host, '   (target_on_host=F means real GPU offload)'
    end subroutine report_env

    subroutine run_test(label, method, alloc_dev, do_flush, inter_only, nmax_p, nlines_p, niter)
        character(*), intent(in) :: label
        integer, intent(in) :: method, nmax_p, nlines_p, niter
        logical, intent(in) :: alloc_dev, do_flush, inter_only

        integer :: mas, npage, total, i, j, k, pk, it, ierr, nreq, funit
        integer :: n_inter_err, n_intra_err, g_inter, g_intra, n_peers, pw
        integer(c_int) :: hip_err, c_dev
        integer(c_size_t) :: nbytes
        real(dp), allocatable, target :: a(:), recvbuf(:)
        real(dp), pointer :: sendbuf(:) => null()
        type(c_ptr) :: sptr
        type(MPI_Request), allocatable :: req(:)
        type(MPI_Datatype) :: vtype
        real(dp) :: t0, t1, exp_val
        logical :: is_inter

        mas   = nmax_p * nlines_p
        npage = nlines_p * NPRO_K
        total = mas * NPRO_K
        funit = 300 + ims_rank

        allocate(a(total), recvbuf(total), req(2*NPRO_K))
        c_dev = omp_get_default_device()    ! defined unconditionally so the cleanup branch is clean
        sptr  = c_null_ptr
        if (alloc_dev) then
            nbytes = int(total, c_size_t) * 8_c_size_t
            sptr   = omp_target_alloc(nbytes, c_dev)
            call c_f_pointer(sptr, sendbuf, [total])
        else
            allocate(sendbuf(total))
        end if
        if (method == METH_STRIDED) then
            call MPI_Type_vector(nmax_p, nlines_p, npage, MPI_DOUBLE_PRECISION, vtype, ierr)
            call MPI_Type_commit(vtype, ierr)
        end if

        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        t0 = MPI_Wtime()
        do it = 1, niter
            ! Fresh GPU-produced source data (physics kernels write `a` each transpose).
            !$omp target teams distribute parallel do
            do k = 1, total
                a(k) = real(ims_rank, dp)*ENC + real(k, dp)
            end do
            !$omp end target teams distribute parallel do

            nreq = 0
            do pk = 0, NPRO_K - 1                          ! post IRECVs (always into unified recvbuf)
                if (pk == ims_pro_k) cycle
                if (inter_only .and. peer_node(pk) == my_node) cycle
                nreq = nreq + 1
                call MPI_Irecv(recvbuf(pk*mas + 1), mas, MPI_DOUBLE_PRECISION, pk, 0, fabric_comm_k, req(nreq), ierr)
            end do

            select case (method)
            case (METH_PACK)                               ! <== PRODUCTION-BOUND PATH
                do pk = 0, NPRO_K - 1
                    if (pk == ims_pro_k) cycle
                    if (inter_only .and. peer_node(pk) == my_node) cycle
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            sendbuf(pk*mas + i*nlines_p + j + 1) = a(pk*nlines_p + i*npage + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end do
                if (do_flush) hip_err = hipDeviceSynchronize()   ! GPU L2 -> HBM before NIC reads sendbuf
                do pk = 0, NPRO_K - 1
                    if (pk == ims_pro_k) cycle
                    if (inter_only .and. peer_node(pk) == my_node) cycle
                    nreq = nreq + 1
                    call MPI_Isend(sendbuf(pk*mas + 1), mas, MPI_DOUBLE_PRECISION, pk, 0, fabric_comm_k, req(nreq), ierr)
                end do

            case (METH_STRIDED)                            ! send strided GPU `a` directly, no pack
                if (do_flush) hip_err = hipDeviceSynchronize()
                do pk = 0, NPRO_K - 1
                    if (pk == ims_pro_k) cycle
                    if (inter_only .and. peer_node(pk) == my_node) cycle
                    nreq = nreq + 1
                    call MPI_Isend(a(pk*nlines_p + 1), 1, vtype, pk, 0, fabric_comm_k, req(nreq), ierr)
                end do

            case (METH_CPUPACK)                            ! current fabricdirect baseline
                hip_err = hipDeviceSynchronize()           ! mandatory: CPU must read GPU-written `a`
                do pk = 0, NPRO_K - 1
                    if (pk == ims_pro_k) cycle
                    if (inter_only .and. peer_node(pk) == my_node) cycle
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            sendbuf(pk*mas + i*nlines_p + j + 1) = a(pk*nlines_p + i*npage + j + 1)
                        end do
                    end do
                    nreq = nreq + 1
                    call MPI_Isend(sendbuf(pk*mas + 1), mas, MPI_DOUBLE_PRECISION, pk, 0, fabric_comm_k, req(nreq), ierr)
                end do
            end select

            call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
        end do
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        t1 = MPI_Wtime()

        ! Verify last iteration's received data.
        n_inter_err = 0; n_intra_err = 0; n_peers = 0
        do pk = 0, NPRO_K - 1
            if (pk == ims_pro_k) cycle
            is_inter = (peer_node(pk) /= my_node)
            if (inter_only .and. .not. is_inter) cycle
            n_peers = n_peers + 1
            pw = peer_world(pk)
            do i = 0, nmax_p - 1
                do j = 0, nlines_p - 1
                    exp_val = real(pw, dp)*ENC + real(ims_pro_k*nlines_p + i*npage + j + 1, dp)
                    if (recvbuf(pk*mas + i*nlines_p + j + 1) /= exp_val) then
                        if (is_inter) then
                            n_inter_err = n_inter_err + 1
                        else
                            n_intra_err = n_intra_err + 1
                        end if
                    end if
                end do
            end do
        end do

        write(funit,'(a,a,a,i3,a,i2,a,i2,a,i2,a,i3,a,i12,a,i12)') &
            '[', trim(label), '] rank', ims_rank, ' node', my_node, ' pro_i', ims_pro_i, &
            ' pro_k', ims_pro_k, ' npeers', n_peers, ' inter_err', n_inter_err, ' intra_err', n_intra_err
        flush(funit)

        call MPI_Reduce(n_inter_err, g_inter, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_Reduce(n_intra_err, g_intra, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        if (ims_rank == 0) then
            write(*,'(a)') '=================================================================='
            write(*,'(2a)')     '[RESULT] ', trim(label)
            write(*,'(a,i0,a,i0)') '   global mismatches: inter-node = ', g_inter, '   intra-node = ', g_intra
            if (g_inter == 0 .and. g_intra == 0) then
                write(*,'(a)') '   ==> PASS (all bytes delivered correctly)'
            else
                write(*,'(a)') '   ==> FAIL'
            end if
            write(*,'(a,es12.4,a,f9.2,a,i0)') '   time/iter = ', (t1-t0)/real(niter,dp), ' s   per-peer = ', &
                real(mas,dp)*8.0d0/(1024.0d0**2), ' MB   niter = ', niter
            write(*,'(a)') '=================================================================='
        end if

        if (method == METH_STRIDED) call MPI_Type_free(vtype, ierr)
        if (alloc_dev) then
            call omp_target_free(sptr, c_dev)
        else
            deallocate(sendbuf)
        end if
        deallocate(a, recvbuf, req)
    end subroutine run_test

end module gpuaware_mod


program vmpi_gpuaware
    use gpuaware_mod
    implicit none
    integer :: ierr, mode, nmax_p, nlines_p, niter, nargs
    character(len=32) :: arg

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, ims_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, ims_nprocs, ierr)
    if (ims_nprocs /= NPRO_I*NPRO_K) then
        if (ims_rank == 0) write(*,*) 'ERROR: need ', NPRO_I*NPRO_K, ' ranks, got ', ims_nprocs
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    call setup_topology()
    call report_env()

    mode = 1; nmax_p = 128; nlines_p = 4096; niter = 5          ! defaults (correctness size)
    nargs = command_argument_count()
    if (nargs >= 1) then; call get_command_argument(1, arg); read(arg,*) mode;     end if
    if (nargs >= 2) then; call get_command_argument(2, arg); read(arg,*) nmax_p;   end if
    if (nargs >= 3) then; call get_command_argument(3, arg); read(arg,*) nlines_p; end if
    if (nargs >= 4) then; call get_command_argument(4, arg); read(arg,*) niter;    end if

    select case (mode)
    case (1); call run_test('M1 PACK unified flush',  METH_PACK,    .false., .true.,  .false., nmax_p, nlines_p, niter)
    case (2); call run_test('M2 PACK unified NOflush', METH_PACK,   .false., .false., .false., nmax_p, nlines_p, niter)
    case (3); call run_test('M3 PACK device flush',   METH_PACK,    .true.,  .true.,  .false., nmax_p, nlines_p, niter)
    case (4); call run_test('M4 STRIDED type flush',  METH_STRIDED, .false., .true.,  .false., nmax_p, nlines_p, niter)
    case (5)                                                        ! bandwidth: inter-node peers, prod size
        call run_test('M5a CPU-pack baseline', METH_CPUPACK, .false., .true., .true., nmax_p, nlines_p, niter)
        call run_test('M5b GPU-pack GPU-aware', METH_PACK,   .false., .true., .true., nmax_p, nlines_p, niter)
    case default
        if (ims_rank == 0) write(*,*) 'unknown mode ', mode
    end select

    call MPI_Finalize(ierr)
end program vmpi_gpuaware
