! vmpi_xcd_write.f90
!
! Validates the FABRIC_DIRECT I-direction setup on Hunter (MI300A, 2 nodes).
!
! Key finding from 2026-05-22 run: the I-direction IS the inter-node direction.
! Node 0 holds pro_i=0..2, node 1 holds pro_i=3..5 (for npro_i=6 on 2 nodes).
! MPI_Win_allocate_shared on the FULL ims_comm_x (cross-node) returns NULL from
! MPI_Win_shared_query for cross-node peers -> job hung at MPI_Barrier.
!
! Correct design (mirrors tlab_mpi_transpose.f90 post-fix):
!   - MPI_Comm_split_type(ims_comm_x, MPI_COMM_TYPE_SHARED) -> intra-node shmem sub-comm
!   - MPI_Win_allocate_shared on the intra-node sub-comm (3 members per node)
!   - hip_write_with_fence for intra-node I-peers (GPU writes)
!   - MPI ISEND/IRECV on dup'd mpi_comm_i for cross-node I-peers
!
! Build: make -f Makefile.xcd_write
! Run:   ./vmpi_xcd_write.x <npro_k> <npro_i> [chunk]
! Example: mpirun -np 48 ./vmpi_xcd_write.x 8 6 512

program vmpi_xcd_write
    use mpi_f08
    use iso_c_binding
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    ! -------------------------------------------------------------------
    ! HIP C interface (hip_write_fence.hip)
    ! -------------------------------------------------------------------
    interface
        subroutine hip_write_with_fence(src, dst, n) bind(C, name='hip_write_with_fence')
            use iso_c_binding
            real(c_double), intent(in)  :: src(*)
            real(c_double), intent(out) :: dst(*)
            integer(c_int), value       :: n
        end subroutine

        function hipHostRegister(ptr, sz, flags) bind(C, name='hipHostRegister') result(ierr)
            use iso_c_binding
            integer(c_int) :: ierr
            type(c_ptr), value       :: ptr
            integer(c_size_t), value :: sz
            integer(c_int), value    :: flags
        end function
    end interface

    ! -------------------------------------------------------------------
    ! MPI state
    ! -------------------------------------------------------------------
    type(MPI_Comm)  :: ims_comm_xz, ims_comm_x, ims_comm_z
    type(MPI_Comm)  :: mpi_comm_i, mpi_comm_k           ! dup'd, untainted
    type(MPI_Comm)  :: shmem_comm_i, shmem_comm_k       ! node-local sub-comms
    type(MPI_Win)   :: win_i, win_k
    type(MPI_Group) :: dir_group_i, world_group

    integer :: ims_pro, ims_npro
    integer :: npro_i, npro_k
    integer :: ims_pro_i, ims_pro_k
    integer :: shmem_size_i, shmem_size_k
    integer :: ims_err, hip_reg_err

    integer, allocatable :: shmem_to_dir_i(:), shmem_to_dir_k(:)

    ! -------------------------------------------------------------------
    ! Window / peer state
    ! -------------------------------------------------------------------
    type(c_ptr)  :: win_baseptr
    integer(MPI_ADDRESS_KIND) :: seg_size_abi, win_size_abi
    integer :: disp_unit_i

    real(dp), pointer, contiguous :: recv_i(:)  => null()  ! own I recv buffer
    real(dp), pointer, contiguous :: pfptr_m(:) => null()  ! per-peer scratch
    type(c_ptr),  allocatable :: peer_i(:)                 ! c_ptr per I-peer
    logical,      allocatable :: is_local_i(:)             ! true if same-node peer

    ! -------------------------------------------------------------------
    ! Test data
    ! -------------------------------------------------------------------
    real(dp), allocatable :: send_buf(:)
    integer :: chunk, seg_size_elems, flat_off

    ! MPI requests / status for cross-node peers
    type(MPI_Request), allocatable :: request(:)
    type(MPI_Status),  allocatable :: status(:)

    ! Verification
    integer, allocatable :: i_global_ranks(:), local_ranks(:)
    real(dp) :: expected, actual
    integer :: g_m, m, j, ip
    integer :: errors, total_errors_global
    character(len=64) :: arg

    ! -------------------------------------------------------------------
    ! MPI init
    ! -------------------------------------------------------------------
    call MPI_Init(ims_err)
    call MPI_Comm_rank(MPI_COMM_WORLD, ims_pro, ims_err)
    call MPI_Comm_size(MPI_COMM_WORLD, ims_npro, ims_err)

    write(500+ims_pro,'(a,i4,a,i4)') '[INIT] MPI init done, rank=', ims_pro, ' npro=', ims_npro
    flush(500+ims_pro)

    if (command_argument_count() < 2) then
        if (ims_pro == 0) write(*,*) 'Usage: vmpi_xcd_write.x <npro_k> <npro_i> [chunk]'
        call MPI_Finalize(ims_err); stop
    end if
    call get_command_argument(1, arg); read(arg,*) npro_k
    call get_command_argument(2, arg); read(arg,*) npro_i
    chunk = 512
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read(arg,*) chunk
    end if

    if (npro_i * npro_k /= ims_npro) then
        if (ims_pro == 0) write(*,'(a,3i6)') 'ERROR: npro_k*npro_i /= nproc:', npro_k, npro_i, ims_npro
        call MPI_Finalize(ims_err); stop
    end if

    ! -------------------------------------------------------------------
    ! Directional comms (Cartesian: dim0=K, dim1=I)
    ! -------------------------------------------------------------------
    block
        integer :: dims(2)
        logical :: period(2), remain_dims(2), reorder
        dims(1) = npro_k; dims(2) = npro_i
        period = .true.; reorder = .false.
        call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xz, ims_err)
        remain_dims(1) = .false.; remain_dims(2) = .true.
        call MPI_Cart_sub(ims_comm_xz, remain_dims, ims_comm_x, ims_err)
        remain_dims(1) = .true.;  remain_dims(2) = .false.
        call MPI_Cart_sub(ims_comm_xz, remain_dims, ims_comm_z, ims_err)
    end block
    call MPI_Comm_rank(ims_comm_x, ims_pro_i, ims_err)
    call MPI_Comm_rank(ims_comm_z, ims_pro_k, ims_err)
    write(500+ims_pro,'(a,i4,a,i4,a,i4)') '[INIT] rank=', ims_pro, ' pro_i=', ims_pro_i, ' pro_k=', ims_pro_k
    flush(500+ims_pro)

    ! -------------------------------------------------------------------
    ! Comm/window order (mirrors production tlab_mpi_transpose.f90):
    !   Step 1: K-dup  (before any window)
    !   Step 2: I-dup  (before any window)
    !   Step 3: K-shmem split + K-window
    !   Step 4: I-shmem split + I-window
    ! -------------------------------------------------------------------

    ! Step 1: K-dup
    call MPI_Comm_dup(ims_comm_z, mpi_comm_k, ims_err)

    ! Step 2: I-dup (BEFORE any window — avoids Cray taint on ims_comm_x)
    call MPI_Comm_dup(ims_comm_x, mpi_comm_i, ims_err)
    write(500+ims_pro,'(a)') '[INIT] K-dup and I-dup done'
    flush(500+ims_pro)

    ! Step 3: K-shmem split + K-window (K is intra-node: all 8 K-peers on same node)
    call MPI_Comm_split_type(ims_comm_z, MPI_COMM_TYPE_SHARED, ims_pro_k, MPI_INFO_NULL, shmem_comm_k, ims_err)
    call MPI_Comm_size(shmem_comm_k, shmem_size_k, ims_err)
    win_size_abi = int(npro_k * chunk, MPI_ADDRESS_KIND) * int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND)
    call MPI_Win_allocate_shared(win_size_abi, int(c_sizeof(1.0_dp)), MPI_INFO_NULL, &
                                  shmem_comm_k, win_baseptr, win_k, ims_err)
    write(500+ims_pro,'(a,i4,a,i4)') '[INIT] K-window done, shmem_size_k=', shmem_size_k, ' pro_k=', ims_pro_k
    flush(500+ims_pro)

    ! Step 4: I-shmem split + I-window (I is inter-node: only intra-node I-peers share memory)
    seg_size_elems = npro_i * chunk
    win_size_abi   = int(seg_size_elems, MPI_ADDRESS_KIND) * int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND)

    call MPI_Comm_split_type(ims_comm_x, MPI_COMM_TYPE_SHARED, ims_pro_i, MPI_INFO_NULL, shmem_comm_i, ims_err)
    call MPI_Comm_size(shmem_comm_i, shmem_size_i, ims_err)
    write(500+ims_pro,'(a,i4,a,i4)') '[INIT] I-shmem split done, shmem_size_i=', shmem_size_i, ' pro_i=', ims_pro_i
    flush(500+ims_pro)

    ! Gather shmem-rank → dir-rank mapping
    allocate(shmem_to_dir_i(0:shmem_size_i-1))
    call MPI_Allgather(ims_pro_i, 1, MPI_INTEGER, shmem_to_dir_i, 1, MPI_INTEGER, shmem_comm_i, ims_err)

    call MPI_Win_allocate_shared(win_size_abi, int(c_sizeof(1.0_dp)), MPI_INFO_NULL, &
                                  shmem_comm_i, win_baseptr, win_i, ims_err)
    write(500+ims_pro,'(a,i4)') '[INIT] I-window on shmem_comm_i done, err=', ims_err
    flush(500+ims_pro)

    ! Populate peer_i and is_local_i
    allocate(peer_i(0:npro_i-1), is_local_i(0:npro_i-1))
    peer_i    = c_null_ptr
    is_local_i = .false.
    do ip = 0, shmem_size_i - 1
        is_local_i(shmem_to_dir_i(ip)) = .true.
        call MPI_Win_shared_query(win_i, ip, seg_size_abi, disp_unit_i, &
                                  peer_i(shmem_to_dir_i(ip)), ims_err)
    end do
    deallocate(shmem_to_dir_i)

    ! Bug A fix: bind recv_i via the own-rank query stored in peer_i(ims_pro_i).
    ! (win_baseptr on Cray may point 1 segsize before the actual segment.)
    call c_f_pointer(peer_i(ims_pro_i), recv_i, [seg_size_elems])

    ! hipHostRegister for intra-node peers (enables GPU writes to their HBM pages)
    do ip = 0, npro_i - 1
        if (is_local_i(ip)) then
            hip_reg_err = hipHostRegister(peer_i(ip), win_size_abi, 0_c_int)
            write(500+ims_pro,'(a,i4,a,i4,a,i0)') &
                '[HIREG] PE', ims_pro, '  peer_pro_i=', ip, '  err=', hip_reg_err
        end if
    end do
    flush(500+ims_pro)

    ! Log local vs remote classification and peer VAs
    block
        integer(MPI_ADDRESS_KIND) :: vap
        do ip = 0, npro_i - 1
            if (is_local_i(ip)) then
                vap = transfer(peer_i(ip), vap)
                write(500+ims_pro,'(a,i4,a,i4,a,i22)') &
                    '[PEER_LOCAL] PE', ims_pro, '  pro_i=', ip, '  VA=', vap
            else
                write(500+ims_pro,'(a,i4,a,i4)') '[PEER_REMOTE] PE', ims_pro, '  pro_i=', ip
            end if
        end do
        flush(500+ims_pro)
    end block

    ! Global rank map for I-comm (needed for verification)
    allocate(i_global_ranks(0:npro_i-1), local_ranks(0:npro_i-1))
    call MPI_Comm_group(ims_comm_x, dir_group_i, ims_err)
    call MPI_Comm_group(MPI_COMM_WORLD, world_group, ims_err)
    do m = 0, npro_i-1; local_ranks(m) = m; end do
    call MPI_Group_translate_ranks(dir_group_i, npro_i, local_ranks, world_group, i_global_ranks, ims_err)

    allocate(send_buf(chunk))
    do j = 1, chunk
        send_buf(j) = dble(ims_pro) * 1000.0d0 + dble(j)
    end do

    ! flat_off: our slot within any peer's recv segment
    flat_off = ims_pro_i * chunk

    ! MPI requests: 2 per cross-node peer (1 IRECV + 1 ISEND)
    allocate(request(2*npro_i), status(2*npro_i))

    ! ================================================================
    ! TEST: hip_write_with_fence (intra-node) + MPI ISEND/IRECV (cross-node)
    ! This mirrors the production FABRIC_DIRECT I-Forward S1..S5 exactly.
    ! ================================================================
    recv_i = 0.0_dp
    call MPI_Barrier(MPI_COMM_WORLD, ims_err)

    ! S1: post IRECVs for cross-node peers
    !   recv into our recv buffer at the slot for peer m
    !   (recv_i(m*chunk+1 : (m+1)*chunk) is our expected location for peer m's data)
    !   But FABRIC_DIRECT uses a separate staging buffer for MPI recv — here for
    !   simplicity we recv directly into recv_i since we own the buffer.
    !   NOTE: in production, MPI recv goes into c_wrk_dp_recv to avoid overwriting
    !   GPU-written data. Here recv_i is zero so we can use it directly.
    block
        integer :: l
        l = 0
        do m = 0, npro_i - 1
            if (.not. is_local_i(m)) then
                l = l + 1
                call MPI_IRECV(recv_i(m*chunk + 1), chunk, MPI_DOUBLE_PRECISION, &
                               m, 99, mpi_comm_i, request(l), ims_err)
            end if
        end do
        write(500+ims_pro,'(a,i4,a,i4)') '[S1] PE', ims_pro, ' IRECVs posted, count=', l
        flush(500+ims_pro)

        ! S2: barrier — all IRECVs posted before any ISEND
        call MPI_Barrier(mpi_comm_i, ims_err)
        write(500+ims_pro,'(a,i4)') '[S2] PE', ims_pro, ' barrier passed'
        flush(500+ims_pro)

        ! S3: post ISENDs for cross-node peers
        do m = 0, npro_i - 1
            if (.not. is_local_i(m)) then
                l = l + 1
                call MPI_ISEND(send_buf(1), chunk, MPI_DOUBLE_PRECISION, &
                               m, 99, mpi_comm_i, request(l), ims_err)
            end if
        end do
        write(500+ims_pro,'(a,i4,a,i4)') '[S3] PE', ims_pro, ' ISENDs posted, total_requests=', l
        flush(500+ims_pro)

        ! S4: GPU writes for intra-node peers via hip_write_with_fence
        do m = 0, npro_i - 1
            if (is_local_i(m)) then
                call c_f_pointer(peer_i(m), pfptr_m, [seg_size_elems])
                call hip_write_with_fence(send_buf(1:chunk), &
                    pfptr_m(flat_off + 1 : flat_off + chunk), int(chunk, c_int))
                nullify(pfptr_m)
            end if
        end do
        write(500+ims_pro,'(a,i4)') '[S4] PE', ims_pro, ' GPU writes done'
        flush(500+ims_pro)

        ! S4b: barrier — GPU writes visible to all before unpack
        call MPI_Barrier(mpi_comm_i, ims_err)
        write(500+ims_pro,'(a,i4)') '[S4b] PE', ims_pro, ' post-write barrier passed'
        flush(500+ims_pro)

        ! S5: WAITALL for cross-node MPI
        if (l > 0) call MPI_WAITALL(l, request(1:l), status(1:l), ims_err)
        write(500+ims_pro,'(a,i4)') '[S5] PE', ims_pro, ' WAITALL done'
        flush(500+ims_pro)
    end block

    ! ================================================================
    ! Verify: recv_i(m*chunk + j) == i_global_ranks(m)*1000 + j
    ! ================================================================
    errors = 0
    do m = 0, npro_i - 1
        g_m = i_global_ranks(m)
        do j = 1, chunk
            expected = dble(g_m) * 1000.0d0 + dble(j)
            actual   = recv_i(m * chunk + j)
            if (abs(actual - expected) > 1.0d-6) then
                errors = errors + 1
                if (errors <= 3) then
                    write(500+ims_pro,'(a,i4,a,i3,a,i5,a,g20.12,a,g20.12,a,l1)') &
                        '[ERR] PE', ims_pro, ' peer=', m, ' j=', j, &
                        '  exp=', expected, '  got=', actual, '  local=', is_local_i(m)
                    flush(500+ims_pro)
                end if
            end if
        end do
    end do
    if (errors == 0) then
        write(500+ims_pro,'(a,i4)') '[PASS] PE', ims_pro
    else
        write(500+ims_pro,'(a,i4,a,i6)') '[FAIL] PE', ims_pro, '  errors=', errors
    end if
    flush(500+ims_pro)

    ! ================================================================
    ! Global report
    ! ================================================================
    call MPI_Reduce(errors, total_errors_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
    if (ims_pro == 0) then
        write(*,'(a)') '================================================'
        write(*,'(a,i4,a,i4,a,i5)') ' npro_k=', npro_k, '  npro_i=', npro_i, '  chunk=', chunk
        write(*,'(a,i4)') ' intra-node I-peers per rank: ', shmem_size_i
        if (total_errors_global == 0) then
            write(*,'(a)') ' PASS: fabricdirect I-direction (shmem GPU + MPI)'
        else
            write(*,'(a,i8,a)') ' FAIL: ', total_errors_global, ' total mismatches'
        end if
        write(*,'(a)') '================================================'
    end if

    ! -------------------------------------------------------------------
    ! Cleanup
    ! -------------------------------------------------------------------
    call MPI_Win_free(win_i, ims_err)
    call MPI_Win_free(win_k, ims_err)
    call MPI_Comm_free(shmem_comm_i, ims_err)
    call MPI_Comm_free(shmem_comm_k, ims_err)
    call MPI_Comm_free(mpi_comm_i, ims_err)
    call MPI_Comm_free(mpi_comm_k, ims_err)
    nullify(recv_i)
    deallocate(send_buf, peer_i, is_local_i, i_global_ranks, local_ranks, request, status)
    call MPI_Finalize(ims_err)

end program vmpi_xcd_write
