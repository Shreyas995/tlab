! vmpi_hip_shmwrite.f90
!
! Validation test for Approach A (HIP kernel + system-scope fence) as a
! replacement for MPI_Win_fence + !$omp target in the fabricdirect transport.
!
! What it tests
! -------------
! Each rank holds a "send chunk" filled with a rank-specific value.
! It delivers that chunk to every peer's recv window:
!   - intra-node peers  -> hip_write_with_fence() (system-scope flush, no Win_fence)
!   - inter-node peers  -> MPI_ISEND/IRECV on MPI_COMM_WORLD
! Process-level ordering uses only MPI_Barrier(MPI_COMM_WORLD), which is
! guaranteed outside the tainted shmem-window lineage.
! Both I-direction and K-direction topologies are tested in sequence.
!
! Expected output: PASS for both I and K with zero mismatches.
! If I passes but K fails (or vice versa), that narrows the topology that hangs.
!
! Build on Hunter (Cray / ROCm):
!   hipcc -c hip_write_fence.hip -o hip_write_fence.o
!   ftn -DUSE_APU vmpi_hip_shmwrite.f90 hip_write_fence.o -lamdhip64 -o vmpi_hip_shmwrite.x
!
! Run:
!   srun -N 2 -n 96 ./vmpi_hip_shmwrite.x <npro_k> <npro_i> [chunk]
!   e.g. srun -N 2 -n 96 ./vmpi_hip_shmwrite.x 8 12 512
!
! Reader-side cache invalidation (optional):
!   If the test reports wrong values with the write side reporting success,
!   uncomment the two hip_invalidate_recv() calls below.  This would indicate
!   that MPI_Win_allocate_shared on this system is COARSE_GRAINED memory and
!   the reader's L2 is not invalidated by the writer's __threadfence_system().

program vmpi_hip_shmwrite
    use mpi_f08
    use iso_c_binding
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    ! -------------------------------------------------------------------
    ! HIP C wrapper interface (hip_write_fence.hip)
    ! -------------------------------------------------------------------
    interface
        subroutine hip_write_with_fence(src, dst, n) bind(C, name='hip_write_with_fence')
            use iso_c_binding
            real(c_double), intent(in)  :: src(*)
            real(c_double), intent(out) :: dst(*)
            integer(c_int), value       :: n
        end subroutine

        subroutine hip_invalidate_recv(buf, n) bind(C, name='hip_invalidate_recv')
            use iso_c_binding
            real(c_double), intent(inout) :: buf(*)
            integer(c_int), value         :: n
        end subroutine

        subroutine hip_write_with_fence_diag(rank, src, dst, n) &
                bind(C, name='hip_write_with_fence_diag')
            use iso_c_binding
            integer(c_int), value       :: rank
            real(c_double), intent(in)  :: src(*)
            real(c_double), intent(out) :: dst(*)
            integer(c_int), value       :: n
        end subroutine

        ! Register a host pointer with the HIP runtime so the GPU MMU has
        ! a mapping for it.  Required for MPI_Win_allocate_shared memory
        ! when HSA_XNACK is not effective.  Returns hipError_t (0 on success).
        function hipHostRegister(ptr, sz, flags) bind(C, name='hipHostRegister') result(ierr)
            use iso_c_binding
            integer(c_int) :: ierr
            type(c_ptr), value :: ptr
            integer(c_size_t), value :: sz
            integer(c_int), value :: flags
        end function
    end interface

    integer :: hip_reg_err

    ! -------------------------------------------------------------------
    ! MPI state
    ! -------------------------------------------------------------------
    type(MPI_Comm)    :: ims_comm_xz, ims_comm_x, ims_comm_z
    type(MPI_Comm)    :: mpi_comm_i, mpi_comm_k  ! dup'd BEFORE any Win_allocate_shared
    type(MPI_Comm)    :: shmem_comm_i, shmem_comm_k
    type(MPI_Win)     :: win_i, win_k
    type(MPI_Group)   :: dir_group_i, dir_group_k, shmem_group_i, shmem_group_k
    type(MPI_Request), allocatable :: req_i(:), req_k(:)
    type(MPI_Status),  allocatable :: sta_i(:), sta_k(:)

    integer :: ims_npro, ims_pro
    integer :: npro_i, npro_k
    integer :: ims_pro_i, ims_pro_k     ! my dir rank in each comm (0-indexed)
    integer :: shmem_size_i, shmem_rank_i, shmem_size_k, shmem_rank_k
    integer :: ims_err
    integer :: chunk                     ! elements per peer-pair (cmd arg, default 512)

    ! -------------------------------------------------------------------
    ! Shared windows and peer pointers
    ! -------------------------------------------------------------------
    real(dp), pointer, contiguous :: recv_i(:) => null()   ! my I-recv window
    real(dp), pointer, contiguous :: recv_k(:) => null()   ! my K-recv window
    real(dp), pointer, contiguous :: pfptr(:) => null()    ! scratch pointer into peer window

    ! peer_win_i(m): c_ptr to base of dir-rank-m peer's shared window (I direction)
    type(c_ptr), allocatable :: peer_win_i(:)  ! 0:npro_i-1
    type(c_ptr), allocatable :: peer_win_k(:)  ! 0:npro_k-1
    logical,     allocatable :: is_local_i(:)
    logical,     allocatable :: is_local_k(:)
    integer,     allocatable :: shmem_to_dir_i(:)  ! shmem-rank -> I dir-rank
    integer,     allocatable :: shmem_to_dir_k(:)  ! shmem-rank -> K dir-rank

    ! -------------------------------------------------------------------
    ! Send buffer (shared between both directions)
    ! -------------------------------------------------------------------
    real(dp), allocatable :: send_buf(:)

    ! Separate inter-shmem MPI recv buffers — plain heap allocations, NOT
    ! shared-window memory. Cray MPICH cannot do RDMA into a buffer that is
    ! both an MPI shared window AND hipHostRegister'd; mixing them deadlocks
    ! MPI_WAITALL. The fabricdirect production code uses the same pattern.
    real(dp), allocatable :: mpi_recv_i(:)   ! npro_i * chunk
    real(dp), allocatable :: mpi_recv_k(:)   ! npro_k * chunk

    ! -------------------------------------------------------------------
    ! Misc
    ! -------------------------------------------------------------------
    integer(MPI_ADDRESS_KIND) :: seg_size, win_size
    integer :: disp_unit
    type(c_ptr) :: win_baseptr, seg_cptr
    integer :: dims(2)
    logical :: period(2), remain_dims(2), reorder
    integer :: m, j, ip, l
    integer :: g_m                           ! global rank of peer m
    real(dp) :: expected, actual
    integer :: errors_i, errors_k, total_errors
    integer :: total_errors_global
    character(len=32) :: arg

    ! -------------------------------------------------------------------
    ! MPI init
    ! -------------------------------------------------------------------
    ! ===== LINE-BY-LINE TRACE: every step writes to fort.1000+rank =====
    ! Use unit 999 for pre-MPI_Init output (single file, all ranks contend but
    ! at least we'll see *something* if any rank reaches this).
    open(unit=999, file='trace_pre_init.log', status='replace', action='write')
    write(999,*) 'L114: about to print ALIVE-PRE'; flush(999)
    write(*,'(a)') '[L114-stdout] alive before MPI_Init'; flush(6)
    print *, '[L114-print] alive before MPI_Init'
    write(999,*) 'L115: about to call MPI_Init'; flush(999)
    close(999)

    call MPI_Init(ims_err)

    call MPI_Comm_rank(MPI_COMM_WORLD, ims_pro, ims_err)
    call MPI_Comm_size(MPI_COMM_WORLD, ims_npro, ims_err)

    ! From here on, every rank writes to its own fort.<1000+rank> file.
    ! The 'fort.1xxx' files will exist on disk even if the job is killed.
    write(1000+ims_pro,*) 'L120: MPI_Init done, rank=', ims_pro, ' npro=', ims_npro; flush(1000+ims_pro)
    write(*,'(a,i4,a,i4)') '[L120-stdout] PE', ims_pro, ' MPI_Init OK, npro=', ims_npro; flush(6)

    write(1000+ims_pro,*) 'L123: about to read args'; flush(1000+ims_pro)
    if (command_argument_count() < 2) then
        write(1000+ims_pro,*) 'L124: not enough args, exiting'; flush(1000+ims_pro)
        if (ims_pro == 0) write(*,*) 'Usage: vmpi_hip_shmwrite.x <npro_k> <npro_i> [chunk]'
        call MPI_Finalize(ims_err); stop
    end if
    call get_command_argument(1, arg); read(arg,*) npro_k
    write(1000+ims_pro,*) 'L129: read npro_k=', npro_k; flush(1000+ims_pro)
    call get_command_argument(2, arg); read(arg,*) npro_i
    write(1000+ims_pro,*) 'L131: read npro_i=', npro_i; flush(1000+ims_pro)
    chunk = 512
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read(arg,*) chunk
    end if
    write(1000+ims_pro,*) 'L136: chunk=', chunk; flush(1000+ims_pro)

    if (npro_i * npro_k /= ims_npro) then
        write(1000+ims_pro,*) 'L139: rank mismatch, exiting'; flush(1000+ims_pro)
        if (ims_pro == 0) write(*,'(a,3i6)') 'ERROR: npro_k*npro_i /= nproc:', npro_k, npro_i, ims_npro
        call MPI_Finalize(ims_err); stop
    end if
    write(1000+ims_pro,*) 'L143: rank check OK'; flush(1000+ims_pro)

    dims(1) = npro_k; dims(2) = npro_i; period = .true.; reorder = .false.
    write(1000+ims_pro,*) 'L146: dims set, calling Cart_create'; flush(1000+ims_pro)
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xz, ims_err)
    write(1000+ims_pro,*) 'L148: Cart_create done, err=', ims_err; flush(1000+ims_pro)

    remain_dims(1) = .false.; remain_dims(2) = .true.
    call MPI_Cart_sub(ims_comm_xz, remain_dims, ims_comm_x, ims_err)
    write(1000+ims_pro,*) 'L152: Cart_sub I done, err=', ims_err; flush(1000+ims_pro)
    remain_dims(1) = .true.;  remain_dims(2) = .false.
    call MPI_Cart_sub(ims_comm_xz, remain_dims, ims_comm_z, ims_err)
    write(1000+ims_pro,*) 'L155: Cart_sub K done, err=', ims_err; flush(1000+ims_pro)

    call MPI_Comm_rank(ims_comm_x, ims_pro_i, ims_err)
    write(1000+ims_pro,*) 'L158: my pro_i=', ims_pro_i; flush(1000+ims_pro)
    call MPI_Comm_rank(ims_comm_z, ims_pro_k, ims_err)
    write(1000+ims_pro,*) 'L160: my pro_k=', ims_pro_k; flush(1000+ims_pro)

    call MPI_Comm_dup(ims_comm_x, mpi_comm_i, ims_err)
    write(1000+ims_pro,*) 'L163: Comm_dup I done, err=', ims_err; flush(1000+ims_pro)
    call MPI_Comm_dup(ims_comm_z, mpi_comm_k, ims_err)
    write(1000+ims_pro,*) 'L165: Comm_dup K done, err=', ims_err; flush(1000+ims_pro)

    allocate(send_buf(chunk))
    write(1000+ims_pro,*) 'L168: send_buf allocated'; flush(1000+ims_pro)
    do j = 1, chunk
        send_buf(j) = dble(ims_pro) * 1000.0d0 + dble(j)
    end do
    write(1000+ims_pro,*) 'L172: send_buf filled'; flush(1000+ims_pro)

    allocate(req_i(2*npro_i), sta_i(2*npro_i))
    allocate(req_k(2*npro_k), sta_k(2*npro_k))
    allocate(mpi_recv_i(npro_i * chunk), mpi_recv_k(npro_k * chunk))
    write(1000+ims_pro,*) 'L176: req/sta arrays + mpi_recv buffers allocated'; flush(1000+ims_pro)

    allocate(is_local_i(0:npro_i-1), peer_win_i(0:npro_i-1))
    is_local_i = .false.;  peer_win_i = c_null_ptr
    write(1000+ims_pro,*) 'L180: I peer arrays allocated'; flush(1000+ims_pro)

    write(1000+ims_pro,*) 'L182: calling Comm_split_type I'; flush(1000+ims_pro)
    call MPI_Comm_split_type(ims_comm_x, MPI_COMM_TYPE_SHARED, ims_pro_i, &
                             MPI_INFO_NULL, shmem_comm_i, ims_err)
    write(1000+ims_pro,*) 'L185: Comm_split_type I done, err=', ims_err; flush(1000+ims_pro)
    call MPI_Comm_rank(shmem_comm_i, shmem_rank_i, ims_err)
    write(1000+ims_pro,*) 'L187: shmem_rank_i=', shmem_rank_i; flush(1000+ims_pro)
    call MPI_Comm_size(shmem_comm_i, shmem_size_i, ims_err)
    write(1000+ims_pro,*) 'L189: shmem_size_i=', shmem_size_i; flush(1000+ims_pro)

    allocate(shmem_to_dir_i(0:shmem_size_i-1))
    write(1000+ims_pro,*) 'L192: calling Allgather I'; flush(1000+ims_pro)
    call MPI_Allgather(ims_pro_i, 1, MPI_INTEGER, shmem_to_dir_i, 1, MPI_INTEGER, &
                       shmem_comm_i, ims_err)
    write(1000+ims_pro,*) 'L195: Allgather I done, err=', ims_err; flush(1000+ims_pro)

    win_size = int(npro_i, MPI_ADDRESS_KIND) * int(chunk, MPI_ADDRESS_KIND) &
             * int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND)
    write(1000+ims_pro,*) 'L199: calling Win_allocate_shared I, size=', win_size; flush(1000+ims_pro)
    call MPI_Win_allocate_shared(win_size, int(c_sizeof(1.0_dp)), MPI_INFO_NULL, &
                                  shmem_comm_i, win_baseptr, win_i, ims_err)
    write(1000+ims_pro,*) 'L202: Win_allocate_shared I done, err=', ims_err; flush(1000+ims_pro)

    call MPI_Comm_group(ims_comm_x, dir_group_i, ims_err)
    call MPI_Comm_group(shmem_comm_i, shmem_group_i, ims_err)
    write(1000+ims_pro,*) 'L206: Comm_group calls done'; flush(1000+ims_pro)
    do ip = 0, shmem_size_i - 1
        is_local_i(shmem_to_dir_i(ip)) = .true.
        call MPI_Win_shared_query(win_i, ip, seg_size, disp_unit, &
                                  peer_win_i(shmem_to_dir_i(ip)), ims_err)
    end do
    write(1000+ims_pro,*) 'L212: Win_shared_query loop done'; flush(1000+ims_pro)
    call c_f_pointer(peer_win_i(ims_pro_i), recv_i, [npro_i * chunk])
    write(1000+ims_pro,*) 'L214: I recv_i bound'; flush(1000+ims_pro)
    ! Register each intra-shmem peer's segment with HIP so the GPU MMU has a
    ! mapping for it.  Without this, hipDeviceSynchronize blocks forever waiting
    ! on a kernel that silently faulted on an unmapped address.  Flag 0 = default.
    do ip = 0, npro_i - 1
        if (is_local_i(ip)) then
            hip_reg_err = hipHostRegister(peer_win_i(ip), &
                int(npro_i * chunk * 8, c_size_t), 0_c_int)
            write(1000+ims_pro,*) 'L222: hipHostRegister I peer ', ip, ' err=', hip_reg_err; flush(1000+ims_pro)
        end if
    end do
    deallocate(shmem_to_dir_i)
    write(1000+ims_pro,*) 'L216: I setup complete'; flush(1000+ims_pro)

    ! ================================================================
    ! K-DIRECTION SETUP
    ! ================================================================
    allocate(is_local_k(0:npro_k-1), peer_win_k(0:npro_k-1))
    is_local_k = .false.;  peer_win_k = c_null_ptr
    write(1000+ims_pro,*) 'L220: K peer arrays allocated'; flush(1000+ims_pro)

    write(1000+ims_pro,*) 'L222: calling Comm_split_type K'; flush(1000+ims_pro)
    call MPI_Comm_split_type(ims_comm_z, MPI_COMM_TYPE_SHARED, ims_pro_k, &
                             MPI_INFO_NULL, shmem_comm_k, ims_err)
    write(1000+ims_pro,*) 'L225: Comm_split_type K done, err=', ims_err; flush(1000+ims_pro)
    call MPI_Comm_rank(shmem_comm_k, shmem_rank_k, ims_err)
    write(1000+ims_pro,*) 'L227: shmem_rank_k=', shmem_rank_k; flush(1000+ims_pro)
    call MPI_Comm_size(shmem_comm_k, shmem_size_k, ims_err)
    write(1000+ims_pro,*) 'L229: shmem_size_k=', shmem_size_k; flush(1000+ims_pro)

    allocate(shmem_to_dir_k(0:shmem_size_k-1))
    write(1000+ims_pro,*) 'L232: calling Allgather K'; flush(1000+ims_pro)
    call MPI_Allgather(ims_pro_k, 1, MPI_INTEGER, shmem_to_dir_k, 1, MPI_INTEGER, &
                       shmem_comm_k, ims_err)
    write(1000+ims_pro,*) 'L235: Allgather K done, err=', ims_err; flush(1000+ims_pro)

    win_size = int(npro_k, MPI_ADDRESS_KIND) * int(chunk, MPI_ADDRESS_KIND) &
             * int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND)
    write(1000+ims_pro,*) 'L239: calling Win_allocate_shared K, size=', win_size; flush(1000+ims_pro)
    call MPI_Win_allocate_shared(win_size, int(c_sizeof(1.0_dp)), MPI_INFO_NULL, &
                                  shmem_comm_k, win_baseptr, win_k, ims_err)
    write(1000+ims_pro,*) 'L242: Win_allocate_shared K done, err=', ims_err; flush(1000+ims_pro)

    call MPI_Comm_group(ims_comm_z, dir_group_k, ims_err)
    call MPI_Comm_group(shmem_comm_k, shmem_group_k, ims_err)
    write(1000+ims_pro,*) 'L246: Comm_group K calls done'; flush(1000+ims_pro)
    do ip = 0, shmem_size_k - 1
        is_local_k(shmem_to_dir_k(ip)) = .true.
        call MPI_Win_shared_query(win_k, ip, seg_size, disp_unit, &
                                  peer_win_k(shmem_to_dir_k(ip)), ims_err)
    end do
    write(1000+ims_pro,*) 'L252: K Win_shared_query loop done'; flush(1000+ims_pro)
    call c_f_pointer(peer_win_k(ims_pro_k), recv_k, [npro_k * chunk])
    write(1000+ims_pro,*) 'L254: K recv_k bound'; flush(1000+ims_pro)
    do ip = 0, npro_k - 1
        if (is_local_k(ip)) then
            hip_reg_err = hipHostRegister(peer_win_k(ip), &
                int(npro_k * chunk * 8, c_size_t), 0_c_int)
            write(1000+ims_pro,*) 'L258: hipHostRegister K peer ', ip, ' err=', hip_reg_err; flush(1000+ims_pro)
        end if
    end do
    deallocate(shmem_to_dir_k)
    write(1000+ims_pro,*) 'L256: ALL SETUP COMPLETE'; flush(1000+ims_pro)

    ! Diagnostic: print peer topology for PE 0
    if (ims_pro == 0) then
        write(*,'(a,i4,a,i4,a,i4,a,i4)') '[INIT] npro_k=', npro_k, ' npro_i=', npro_i, &
            ' chunk=', chunk, ' total_ranks=', ims_npro
        write(*,'(a,i4,a,i4,a,i4)') '[INIT] PE 0: pro_k=', ims_pro_k, &
            ' pro_i=', ims_pro_i, ' shmem_size_i=', shmem_size_i
        do m = 0, npro_i - 1
            if (is_local_i(m)) then
                write(*,'(a,i3,a)') '[INIT] I peer ', m, '  -> INTRA (HIP write)'
            else
                write(*,'(a,i3,a)') '[INIT] I peer ', m, '  -> INTER (MPI send)'
            end if
        end do
    end if

    ! Global sync: all init done, recv windows allocated, query complete
    write(1000+ims_pro,*) 'L260: entering first Barrier WORLD'; flush(1000+ims_pro)
    call MPI_Barrier(MPI_COMM_WORLD, ims_err)
    write(1000+ims_pro,*) 'L262: passed first Barrier WORLD, err=', ims_err; flush(1000+ims_pro)

    ! ================================================================
    ! TEST I-DIRECTION
    ! ================================================================
    recv_i = 0.0_dp
    write(1000+ims_pro,*) 'L267: recv_i zeroed'; flush(1000+ims_pro)

    ! Step 1: post IRECVs for inter-shmem I-peers into the PLAIN buffer
    ! (not into recv_i, which is shared-window + HIP-registered memory).
    mpi_recv_i = 0.0_dp
    l = 0
    do m = 0, npro_i - 1
        if (.not. is_local_i(m)) then
            l = l + 1
            g_m = ims_pro_k * npro_i + m
            write(1000+ims_pro,*) 'L274: posting IRECV from peer dir=', m, ' g_m=', g_m; flush(1000+ims_pro)
            call MPI_IRECV(mpi_recv_i(m*chunk + 1), chunk, MPI_DOUBLE_PRECISION, &
                           g_m, 1001, MPI_COMM_WORLD, req_i(l), ims_err)
            write(1000+ims_pro,*) 'L277: IRECV posted, err=', ims_err; flush(1000+ims_pro)
        end if
    end do
    write(1000+ims_pro,*) 'L280: all IRECVs done, l=', l; flush(1000+ims_pro)

    write(1000+ims_pro,*) 'L282: entering Barrier (post-IRECV)'; flush(1000+ims_pro)
    call MPI_Barrier(MPI_COMM_WORLD, ims_err)
    write(1000+ims_pro,*) 'L284: passed Barrier (post-IRECV)'; flush(1000+ims_pro)

    ! Step 3: post ISENDs for inter-node I-peers
    do m = 0, npro_i - 1
        if (.not. is_local_i(m)) then
            l = l + 1
            g_m = ims_pro_k * npro_i + m
            write(1000+ims_pro,*) 'L291: posting ISEND to peer dir=', m, ' g_m=', g_m; flush(1000+ims_pro)
            call MPI_ISEND(send_buf(1), chunk, MPI_DOUBLE_PRECISION, &
                           g_m, 1001, MPI_COMM_WORLD, req_i(l), ims_err)
            write(1000+ims_pro,*) 'L294: ISEND posted, err=', ims_err; flush(1000+ims_pro)
        end if
    end do
    write(1000+ims_pro,*) 'L297: all ISENDs done, l=', l; flush(1000+ims_pro)

    ! Step 4: HIP write to intra-shmem I-peers (diagnostic version)
    do m = 0, npro_i - 1
        if (is_local_i(m)) then
            write(1000+ims_pro,*) 'L302: HIP write to intra peer dir=', m; flush(1000+ims_pro)
            call c_f_pointer(peer_win_i(m), pfptr, [npro_i * chunk])
            ! Use the DIAGNOSTIC wrapper — writes hip_trace_<rank>.log so we can
            ! see exactly where inside the HIP layer the call wedges.
            call hip_write_with_fence_diag(int(ims_pro, c_int), send_buf(1), &
                pfptr(ims_pro_i * chunk + 1), int(chunk, c_int))
            nullify(pfptr)
            write(1000+ims_pro,*) 'L306: HIP write done for peer dir=', m; flush(1000+ims_pro)
        end if
    end do
    write(1000+ims_pro,*) 'L309: all HIP writes done'; flush(1000+ims_pro)

    if (l > 0) then
        write(1000+ims_pro,*) 'L312: entering WAITALL, l=', l; flush(1000+ims_pro)
        call MPI_WAITALL(l, req_i(1:l), sta_i(1:l), ims_err)
        write(1000+ims_pro,*) 'L314: WAITALL done, err=', ims_err; flush(1000+ims_pro)
    else
        write(1000+ims_pro,*) 'L316: no inter-shmem peers, skipping WAITALL'; flush(1000+ims_pro)
    end if

    write(1000+ims_pro,*) 'L319: entering final I-Barrier'; flush(1000+ims_pro)
    call MPI_Barrier(MPI_COMM_WORLD, ims_err)
    write(1000+ims_pro,*) 'L321: passed final I-Barrier'; flush(1000+ims_pro)

    ! Copy inter-shmem MPI recv data into recv_i so verification sees it
    do m = 0, npro_i - 1
        if (.not. is_local_i(m)) then
            do j = 1, chunk
                recv_i(m*chunk + j) = mpi_recv_i(m*chunk + j)
            end do
        end if
    end do
    write(1000+ims_pro,*) 'L329: inter-shmem data copied to recv_i'; flush(1000+ims_pro)

    ! Optional reader-side L2 invalidation — uncomment if mismatches are seen
    ! despite the writer reporting success (indicates COARSE_GRAINED memory):
    ! call hip_invalidate_recv(recv_i(1), int(npro_i * chunk, c_int))
    ! call MPI_Barrier(MPI_COMM_WORLD, ims_err)

    ! Verify: recv_i(m*chunk + j) must equal (ims_pro_k*npro_i + m)*1000.0 + j
    errors_i = 0
    do m = 0, npro_i - 1
        g_m = ims_pro_k * npro_i + m
        do j = 1, chunk
            expected = dble(g_m) * 1000.0d0 + dble(j)
            actual   = recv_i(m * chunk + j)
            if (abs(actual - expected) > 1.0d-6) then
                errors_i = errors_i + 1
                if (errors_i <= 5) then
                    write(500+ims_pro,'(a,i4,a,i3,a,i5,a,g20.12,a,g20.12)') &
                        '[I-ERR] PE', ims_pro, ' peer_dir=', m, ' j=', j, &
                        ' expected=', expected, ' got=', actual
                    flush(500+ims_pro)
                end if
            end if
        end do
    end do

    ! ================================================================
    ! TEST K-DIRECTION
    ! ================================================================
    recv_k = 0.0_dp
    mpi_recv_k = 0.0_dp
    write(1000+ims_pro,*) 'L335: K test, recv_k zeroed'; flush(1000+ims_pro)
    l = 0
    do m = 0, npro_k - 1
        if (.not. is_local_k(m)) then
            l = l + 1
            g_m = m * npro_i + ims_pro_i
            call MPI_IRECV(mpi_recv_k(m*chunk + 1), chunk, MPI_DOUBLE_PRECISION, &
                           g_m, 1002, MPI_COMM_WORLD, req_k(l), ims_err)
        end if
    end do
    write(1000+ims_pro,*) 'L344: K IRECVs done, l=', l; flush(1000+ims_pro)
    call MPI_Barrier(MPI_COMM_WORLD, ims_err)
    write(1000+ims_pro,*) 'L346: K post-IRECV Barrier passed'; flush(1000+ims_pro)

    do m = 0, npro_k - 1
        if (.not. is_local_k(m)) then
            l = l + 1
            g_m = m * npro_i + ims_pro_i
            call MPI_ISEND(send_buf(1), chunk, MPI_DOUBLE_PRECISION, &
                           g_m, 1002, MPI_COMM_WORLD, req_k(l), ims_err)
        end if
    end do
    write(1000+ims_pro,*) 'L355: K ISENDs done, l=', l; flush(1000+ims_pro)

    do m = 0, npro_k - 1
        if (is_local_k(m)) then
            call c_f_pointer(peer_win_k(m), pfptr, [npro_k * chunk])
            call hip_write_with_fence(send_buf(1), pfptr(ims_pro_k * chunk + 1), int(chunk, c_int))
            nullify(pfptr)
        end if
    end do
    write(1000+ims_pro,*) 'L364: K HIP writes done'; flush(1000+ims_pro)

    if (l > 0) then
        write(1000+ims_pro,*) 'L367: K entering WAITALL, l=', l; flush(1000+ims_pro)
        call MPI_WAITALL(l, req_k(1:l), sta_k(1:l), ims_err)
        write(1000+ims_pro,*) 'L369: K WAITALL done'; flush(1000+ims_pro)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ims_err)
    write(1000+ims_pro,*) 'L372: K final Barrier passed'; flush(1000+ims_pro)

    do m = 0, npro_k - 1
        if (.not. is_local_k(m)) then
            do j = 1, chunk
                recv_k(m*chunk + j) = mpi_recv_k(m*chunk + j)
            end do
        end if
    end do
    write(1000+ims_pro,*) 'L381: K inter-shmem data copied'; flush(1000+ims_pro)

    ! Optional: call hip_invalidate_recv(recv_k(1), int(npro_k * chunk, c_int))
    ! call MPI_Barrier(MPI_COMM_WORLD, ims_err)

    ! Verify K: recv_k(m*chunk + j) must equal (m*npro_i + ims_pro_i)*1000.0 + j
    errors_k = 0
    do m = 0, npro_k - 1
        g_m = m * npro_i + ims_pro_i
        do j = 1, chunk
            expected = dble(g_m) * 1000.0d0 + dble(j)
            actual   = recv_k(m * chunk + j)
            if (abs(actual - expected) > 1.0d-6) then
                errors_k = errors_k + 1
                if (errors_k <= 5) then
                    write(500+ims_pro,'(a,i4,a,i3,a,i5,a,g20.12,a,g20.12)') &
                        '[K-ERR] PE', ims_pro, ' peer_dir=', m, ' j=', j, &
                        ' expected=', expected, ' got=', actual
                    flush(500+ims_pro)
                end if
            end if
        end do
    end do

    ! ================================================================
    ! REPORT
    ! ================================================================
    total_errors = errors_i + errors_k
    call MPI_Reduce(total_errors, total_errors_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)

    ! Per-rank summary to fort.500+rank
    write(500+ims_pro,'(a,i4,a,i5,a,i5)') &
        '[RESULT] PE', ims_pro, ' errors_i=', errors_i, ' errors_k=', errors_k
    flush(500+ims_pro)

    if (ims_pro == 0) then
        write(*,'(a)') '============================================'
        if (total_errors_global == 0) then
            write(*,'(a)') 'PASS: HIP shmwrite + system fence CORRECT'
            write(*,'(a)') '  No mismatches in I or K direction.'
        else
            write(*,'(a,i8,a)') 'FAIL: ', total_errors_global, ' total mismatches'
            write(*,'(a)') '  Check fort.500+rank for per-rank [I-ERR]/[K-ERR] lines.'
            write(*,'(a)') '  If errors only in I/K: that is the direction that still hangs.'
            write(*,'(a)') '  If COARSE_GRAINED memory is suspected, uncomment the'
            write(*,'(a)') '  hip_invalidate_recv() calls in vmpi_hip_shmwrite.f90.'
        end if
        write(*,'(a)') '============================================'
    end if

    ! Clean up
    call MPI_Win_free(win_i, ims_err)
    call MPI_Win_free(win_k, ims_err)
    call MPI_Comm_free(shmem_comm_i, ims_err)
    call MPI_Comm_free(shmem_comm_k, ims_err)
    call MPI_Comm_free(mpi_comm_i, ims_err)
    call MPI_Comm_free(mpi_comm_k, ims_err)
    nullify(recv_i, recv_k)
    deallocate(send_buf, is_local_i, is_local_k, peer_win_i, peer_win_k)
    deallocate(req_i, sta_i, req_k, sta_k)
    deallocate(mpi_recv_i, mpi_recv_k)

    call MPI_Finalize(ims_err)

end program vmpi_hip_shmwrite
