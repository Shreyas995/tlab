! vmpi_nodewin.f90
! ============================================================================
! Node-window + GPU-aware-MPI COMBINATION validation (Option-2 de-risk). On a 2-node job, one
! K-transpose runs the intra-node peers via node-local shared-window GPU writes (+ MPI_Win_fence)
! AND the inter-node peers via GPU-aware MPI_Isend of a GPU-resident buffer (no flush) -- the exact
! production fabricdirect-K pattern once Option 2 lands. The window mechanism alone passed earlier
! (CPU-pack inter leg; backup: vmpi_nodewin.f90.bak.cpumix); this version swaps the inter leg to
! GPU-aware MPI to prove the two coexist (the one untested interaction before editing production).
!
! It answers:
!   Q1  Does MPI_Win_allocate_shared on a 24-rank NODE comm (split by node, with the
!       contiguous-segment info key) give one contiguous cross-XCD mapping?  -> [INIT] delta==expect
!   Q2  Can a GPU kernel write ACROSS XCDs into a peer's segment coherently under
!       MPI_Win_fence?  (the 4 intra-node K-peers sit in 4 different XCDs)
!   Q3  Do BOTH legs deliver correctly when GPU-aware inter-node MPI_Isend runs alongside the
!       shared-window writes + fence in the same transpose?  -> [VERIFY] errors==0
!
! Gotchas applied (see CLAUDE.md "Hunter harness/build gotchas"):
!   - use mpi (NOT mpi_f08): mpi_f08 predefined datatype constants miscompile under -fopenmp offload.
!   - NO module load craype-accel-amd-gfx942 in the PBS (inherit .bashrc HLRS/APU/2026.1).
!   - literal args, no block-local array constructors (Cray zeroes them under offload+USM).
!   - ONE fused c_f_pointer over the whole contiguous window (apudirect style); per-peer c_f_pointer
!     inside a target region caused HSA aperture violations on cross-XCD writes (Bug C).
!
! Encode a(idx)=world*ENC+idx. Sender S (K-rank sk) writes into dest D's segment at slot sk,
! source a_S(D_krank*nlines + i*npage + j). Receiver D verifies its slot sk == w_S*ENC +
! (D_krank*nlines + i*npage + j + 1).  Topology: world = pro_k*NPRO_I + pro_i; node = world/24.
! ============================================================================
module nodewin_mod
    use mpi
    use iso_c_binding
    use omp_lib
    implicit none
    !$omp requires unified_shared_memory
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: NPRO_I = 6, NPRO_K = 8, RANKS_PER_NODE = 24
    real(dp), parameter :: ENC = 1.0d9
    integer :: ims_rank, ims_nprocs, ims_pro_i, ims_pro_k, my_node
contains
    subroutine report_env()
        logical :: on_host
        on_host = .true.
        !$omp target map(tofrom: on_host)
        on_host = omp_is_initial_device()
        !$omp end target
        if (ims_rank == 0) write(*,'(a,i0,a,l1)') '[ENV] omp_get_num_devices=', omp_get_num_devices(), &
            '   target_on_host=', on_host
    end subroutine
end module

program vmpi_nodewin
    use nodewin_mod
    implicit none
    !$omp requires unified_shared_memory
    integer :: ierr, comm_xz, node_comm, node_rank, node_sz, fabric_comm_k, nreq
    integer :: dims(2), coord(2), lp, dk, sk, i, j, funit, disp_unit, win, win_info
    integer :: nmax_p, nlines_p, niter, chunk, npage, segsize, it, wd, ws
    integer :: n_err, g_err, n_intra, n_inter, my_noncontig, g_noncontig
    logical :: period(2)
    integer(MPI_ADDRESS_KIND) :: segbytes, qsize, va0, va
    integer(8) :: base
    type(c_ptr) :: baseptr
    real(dp), pointer :: my_recv(:) => null(), all_win(:) => null()
    real(dp), allocatable, target :: a(:)
    real(dp), allocatable :: inter_recv(:), c_send(:)      ! inter-node MPI recv / pack-send buffers
    integer, allocatable :: req(:)
    type(c_ptr), allocatable :: peer_cptr(:)
    real(dp) :: t0, t1, exp_val

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, ims_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, ims_nprocs, ierr)
    if (ims_nprocs /= NPRO_I*NPRO_K) then
        if (ims_rank == 0) write(*,*) 'ERROR need ', NPRO_I*NPRO_K, ' ranks, got ', ims_nprocs
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    dims(1) = NPRO_K; dims(2) = NPRO_I; period = .true.
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, .false., comm_xz, ierr)
    call MPI_Cart_coords(comm_xz, ims_rank, 2, coord, ierr)
    ims_pro_k = coord(1); ims_pro_i = coord(2)
    my_node = ims_rank / RANKS_PER_NODE
    funit = 300 + ims_rank
    call report_env()

    ! NODE-local comm: split by node (NOT MPI_Comm_split_type, which fragments per-XCD).
    call MPI_Comm_split(MPI_COMM_WORLD, my_node, ims_rank, node_comm, ierr)
    call MPI_Comm_rank(node_comm, node_rank, ierr)
    call MPI_Comm_size(node_comm, node_sz, ierr)
    ! K-comm for the INTER-node MPI leg (same split fabricdirect uses): color=pro_i, key=pro_k.
    ! local rank in fabric_comm_k = ims_pro_k; K-peer dk = local rank dk.
    call MPI_Comm_split(MPI_COMM_WORLD, ims_pro_i, ims_pro_k, fabric_comm_k, ierr)

    nmax_p = 64; nlines_p = 1024; niter = 10           ! small = correctness; per-seg ~4 MB
    chunk = nmax_p*nlines_p; npage = nlines_p*NPRO_K; segsize = NPRO_K*chunk

    call MPI_Info_create(win_info, ierr)
    call MPI_Info_set(win_info, 'alloc_shared_noncontig', 'false', ierr)   ! force contiguous segments
    segbytes = int(segsize, MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
    call MPI_Win_allocate_shared(segbytes, 8, win_info, node_comm, baseptr, win, ierr)
    call MPI_Info_free(win_info, ierr)
    call c_f_pointer(baseptr, my_recv, [segsize])
    allocate(peer_cptr(0:node_sz-1))
    do lp = 0, node_sz-1
        call MPI_Win_shared_query(win, lp, qsize, disp_unit, peer_cptr(lp), ierr)
    end do

    ! Q1: contiguity. delta(lp) must == lp*segbytes for the fused span to be valid cross-XCD.
    va0 = transfer(peer_cptr(0), va0)
    my_noncontig = 0
    do lp = 0, node_sz-1
        va = transfer(peer_cptr(lp), va)
        if (va - va0 /= int(lp, MPI_ADDRESS_KIND)*segbytes) my_noncontig = my_noncontig + 1
        write(funit,'(a,i3,a,i3,a,i20,a,i20)') '[INIT] noderank', node_rank, ' peer', lp, &
            ' delta=', va - va0, ' expect=', int(lp, MPI_ADDRESS_KIND)*segbytes
    end do
    flush(funit)
    call MPI_Reduce(my_noncontig, g_noncontig, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! ONE fused pointer over the whole contiguous node window (peer lp at offset lp*segsize).
    call c_f_pointer(peer_cptr(0), all_win, [int(segsize,8)*int(node_sz,8)])
    allocate(a(segsize), inter_recv(segsize), c_send(segsize), req(2*NPRO_K))
    n_intra = 0; n_inter = 0
    do dk = 0, NPRO_K - 1
        if ((dk*NPRO_I + ims_pro_i)/RANKS_PER_NODE == my_node) then
            n_intra = n_intra + 1
        else
            n_inter = n_inter + 1
        end if
    end do

    n_err = 0
    t0 = MPI_Wtime(); t1 = t0
    if (my_noncontig == 0) then          ! only write if the window is one contiguous mapping
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        t0 = MPI_Wtime()
        do it = 1, niter
            !$omp target teams distribute parallel do
            do i = 1, segsize
                a(i) = real(ims_rank, dp)*ENC + real(i, dp)
            end do
            !$omp end target teams distribute parallel do
            ! THE MIX (exact Option-1 K-Forward pattern): inter-node K-peers via MPI, intra via window.
            ! 1. post IRECVs for INTER-node K-peers into inter_recv slot dk.
            nreq = 0
            do dk = 0, NPRO_K - 1
                if ((dk*NPRO_I + ims_pro_i)/RANKS_PER_NODE == my_node) cycle   ! intra -> window
                nreq = nreq + 1
                call MPI_Irecv(inter_recv(dk*chunk + 1), chunk, MPI_DOUBLE_PRECISION, dk, 0, fabric_comm_k, req(nreq), ierr)
            end do
            ! 2. open window epoch.
            call MPI_Win_fence(0, win, ierr)
            ! 3. INTRA-node peers: ONE fused cross-XCD GPU write over ALL peers (apudirect-style collapse(3)
            !    over dk,i,j); inter peers are skipped by the in-kernel mask. The base offset is inlined into
            !    the index (no shared scalars), exactly the production fusion pattern being de-risked here.
            !$omp target teams distribute parallel do collapse(3)
            do dk = 0, NPRO_K - 1
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        if ((dk*NPRO_I + ims_pro_i)/RANKS_PER_NODE == my_node) then
                            all_win(int(mod(dk*NPRO_I + ims_pro_i, RANKS_PER_NODE),8)*int(segsize,8) &
                                    + ims_pro_k*chunk + i*nlines_p + j + 1) = a(dk*nlines_p + i*npage + j + 1)
                        end if
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            ! 4. INTER-node peers: GPU pack a -> c_send, then GPU-AWARE ISEND of the GPU-resident c_send
            !    (NO hipDeviceSynchronize flush; MPICH orders the GPU stream — vmpi_gpuaware M2). This is the
            !    production Option-2 path, tested here COEXISTING with the node-window writes + MPI_Win_fence
            !    (the one previously-untested interaction before porting Option 2 into the K transpose).
            do dk = 0, NPRO_K - 1
                if ((dk*NPRO_I + ims_pro_i)/RANKS_PER_NODE == my_node) cycle
                !$omp target teams distribute parallel do collapse(2)
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        c_send(dk*chunk + i*nlines_p + j + 1) = a(dk*nlines_p + i*npage + j + 1)
                    end do
                end do
                !$omp end target teams distribute parallel do
                nreq = nreq + 1
                call MPI_Isend(c_send(dk*chunk + 1), chunk, MPI_DOUBLE_PRECISION, dk, 0, fabric_comm_k, req(nreq), ierr)
            end do
            ! 5. close epoch (intra writes committed) then wait for inter MPI.
            call MPI_Win_fence(0, win, ierr)
            if (nreq > 0) call MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE, ierr)
        end do
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        t1 = MPI_Wtime()

        ! verify ALL 8 K-peers: intra from my window segment, inter from the MPI recv buffer.
        do sk = 0, NPRO_K - 1
            ws = sk*NPRO_I + ims_pro_i
            do i = 0, nmax_p - 1
                do j = 0, nlines_p - 1
                    exp_val = real(ws, dp)*ENC + real(ims_pro_k*nlines_p + i*npage + j + 1, dp)
                    if (ws/RANKS_PER_NODE == my_node) then           ! intra -> window segment
                        if (my_recv(sk*chunk + i*nlines_p + j + 1) /= exp_val) n_err = n_err + 1
                    else                                              ! inter -> MPI recv buffer
                        if (inter_recv(sk*chunk + i*nlines_p + j + 1) /= exp_val) n_err = n_err + 1
                    end if
                end do
            end do
        end do
    end if
    write(funit,'(a,i3,a,i2,a,i2,a,i2,a,i2,a,i2,a,i12,a,i6)') '[VERIFY] rank', ims_rank, ' node', my_node, &
        ' pro_i', ims_pro_i, ' pro_k', ims_pro_k, ' intra', n_intra, ' inter', n_inter, ' errors', n_err, &
        ' noncontig', my_noncontig
    flush(funit)
    call MPI_Reduce(n_err, g_err, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (ims_rank == 0) then
        write(*,'(a)') '=================================================================='
        write(*,'(a,i0,a)') '[Q1 contiguity]  ranks with NON-contiguous window = ', g_noncontig, &
            '   (0 = every node window is one contiguous cross-XCD mapping)'
        write(*,'(a,i0,a)') '[correctness]    global errors = ', g_err, &
            '   (all 8 K-peers: intra via window + inter via GPU-AWARE MPI, mixed in one transpose)'
        if (g_noncontig == 0 .and. g_err == 0) then
            write(*,'(a)') '   ==> PASS: node-window + GPU-aware MPI coexist safely on 2 nodes (no hang, correct)'
        else
            write(*,'(a)') '   ==> FAIL (see fort.300..347 [INIT]/[VERIFY])'
        end if
        write(*,'(a,es12.4,a,f8.2,a)') '   time/iter = ', (t1 - t0)/real(niter,dp), ' s   per-seg = ', &
            real(segsize, dp)*8.0d0/(1024.0d0**2), ' MB'
        write(*,'(a)') '=================================================================='
    end if

    call MPI_Win_free(win, ierr)
    deallocate(a, inter_recv, c_send, req, peer_cptr)
    call MPI_Finalize(ierr)
end program vmpi_nodewin
