! vmpi_fabricdirect_commtest.f90
!
! Diagnoses the FABRIC_DIRECT I-direction IFR hang and tests communicator fixes.
!
! Decomposition: NPRO_K x NPRO_I  (default 8x6 = 48 ranks, matching production)
!   global rank r -> ims_pro_k = r/NPRO_I,  ims_pro_i = r mod NPRO_I
!
! The test runs four sequential phases. Each phase mirrors the FABRIC_DIRECT IFR
! block in tlab_mpi_transpose.f90: post IRECVs, call hipDeviceSynchronize, post
! ISENDs, then WAITALL. The four phases vary the communicator choice and the
! position of hipDeviceSynchronize:
!
!   Phase 1 — dup(ims_comm_x),   hipDevSync AFTER  IRECVs  [current code — hangs]
!   Phase 2 — dup(ims_comm_x),   hipDevSync BEFORE IRECVs
!   Phase 3 — MPI_Comm_split,    hipDevSync AFTER  IRECVs
!   Phase 4 — MPI_Comm_split,    hipDevSync BEFORE IRECVs  [proposed fix]
!
! BEFORE the phases, this test allocates a K-direction MPI_Win_allocate_shared
! on a MPI_Comm_split_type sub-comm of ims_comm_z — the exact taint source in
! the real FABRIC_DIRECT init. This ensures the test sees the same tainted state.
!
! What to look for:
!   - Phases that print "[Px_S9]" for all ranks: PASS (WAITALL completed)
!   - Phases that hang between [Px_S4] and [Px_Ssync]: hipDevSync is the blocker
!   - Phases that hang between [Px_S8] and [Px_S9]: WAITALL is the blocker
!   - "[Px] rank N errors=0 PASS" lines: communicator routes to correct peers
!
! Build on Hunter (Cray / ROCm):
!   hipcc -c hip_write_fence.hip -o hip_write_fence.o
!   ftn -DUSE_APU vmpi_fabricdirect_commtest.f90 hip_write_fence.o -lamdhip64 \
!       -o vmpi_fabricdirect_commtest.x
!
! Run (2 nodes, 24 ranks each = 48 total, decomp 8x6):
!   srun -N 2 --ntasks-per-node=24 ./vmpi_fabricdirect_commtest.x 8 6 500

program vmpi_fabricdirect_commtest
    use mpi_f08
    use iso_c_binding
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    ! ---------------------------------------------------------------
    ! HIP interface
    ! ---------------------------------------------------------------
    interface
#ifdef USE_APU
        function hipDeviceSynchronize() bind(C, name='hipDeviceSynchronize') result(res)
            use iso_c_binding
            integer(c_int) :: res
        end function
#endif
    end interface

    ! ---------------------------------------------------------------
    ! Decomposition
    ! ---------------------------------------------------------------
    integer :: NPRO_I, NPRO_K, CHUNK
    integer :: ims_pro_i, ims_pro_k, ims_rank, ims_nprocs
    integer :: coord(2), dims(2), ierr_mpi
    logical :: period(2), reorder, remain_dims(2)

    ! Communicators
    type(MPI_Comm) :: ims_comm_xz, ims_comm_x, ims_comm_z
    type(MPI_Comm) :: shmem_comm_k, test_comm

    ! K-direction window (the taint source)
    type(MPI_Win) :: k_win
    type(c_ptr)   :: k_win_baseptr_c
    integer(MPI_ADDRESS_KIND) :: win_seg_size
    integer :: win_disp_unit_int

    ! Send/recv maps (same logic as tlab_mpi_transpose.f90 lines 303-310)
    integer, allocatable :: maps_send_i(:), maps_recv_i(:)

    ! Staging buffer layout: [1..CHUNK*NPRO_I] = recv, [CHUNK*NPRO_I+1..2*CHUNK*NPRO_I] = send
    real(dp), allocatable :: staging(:)

    ! Requests and status
    type(MPI_Request), allocatable :: req(:)
    type(MPI_Status),  allocatable :: stat(:)

    ! Misc
    integer :: ip, m, l, nr, ns, ipr, ips, n_errors, phase, i
    integer :: expected_rank
    integer :: phase_lo, phase_hi   ! run phases phase_lo..phase_hi only
#ifdef USE_APU
    integer(c_int) :: hip_err
#endif
    character(len=80) :: arg

    ! ---------------------------------------------------------------
    ! Init MPI
    ! ---------------------------------------------------------------
    call MPI_Init(ierr_mpi)
    call MPI_Comm_rank(MPI_COMM_WORLD, ims_rank, ierr_mpi)
    call MPI_Comm_size(MPI_COMM_WORLD, ims_nprocs, ierr_mpi)

    ! Argument order: <NPRO_K> <NPRO_I> <CHUNK> [PHASE]
    ! PHASE: single digit 1-4 runs only that phase; 0 or absent runs all 4.
    ! Run each phase as a separate job because a hang blocks subsequent phases.
    NPRO_K = 8; NPRO_I = 6; CHUNK = 500; phase_lo = 1; phase_hi = 4
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg); read(arg,*) NPRO_K
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg); read(arg,*) NPRO_I
    end if
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read(arg,*) CHUNK
    end if
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg); read(arg,*) phase_lo
        if (phase_lo >= 1 .and. phase_lo <= 4) phase_hi = phase_lo
    end if

    if (ims_nprocs /= NPRO_I * NPRO_K) then
        if (ims_rank == 0) &
            write(*,'(a,i0,a,i0)') 'ERROR: need ', NPRO_I*NPRO_K, ' ranks, got ', ims_nprocs
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr_mpi)
    end if

    ! ---------------------------------------------------------------
    ! Cartesian topology — identical to tlab_mpi_procs.f90 lines 75-95
    ! dims(1)=NPRO_K (first coord = k), dims(2)=NPRO_I (second coord = i)
    ! rank r -> coord(1)=k=r/NPRO_I, coord(2)=i=r mod NPRO_I
    ! ---------------------------------------------------------------
    dims(1) = NPRO_K; dims(2) = NPRO_I
    period = .true.; reorder = .false.
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xz, ierr_mpi)
    call MPI_Cart_coords(ims_comm_xz, ims_rank, 2, coord, ierr_mpi)
    ims_pro_k = coord(1); ims_pro_i = coord(2)

    ! ims_comm_x: fix dim[0]=K, keep dim[1]=I free
    remain_dims(1) = .false.; remain_dims(2) = .true.
    call MPI_Cart_sub(ims_comm_xz, remain_dims, ims_comm_x, ierr_mpi)

    ! ims_comm_z: keep dim[0]=K free, fix dim[1]=I
    remain_dims(1) = .true.; remain_dims(2) = .false.
    call MPI_Cart_sub(ims_comm_xz, remain_dims, ims_comm_z, ierr_mpi)

    ! ---------------------------------------------------------------
    ! Send/recv maps — identical to tlab_mpi_transpose.f90 lines 303-310
    ! Values in 0..NPRO_I-1 = local ranks in the I-comm = ims_pro_i of peers
    ! ---------------------------------------------------------------
    allocate(maps_send_i(NPRO_I), maps_recv_i(NPRO_I))
    do ip = 0, NPRO_I - 1
        maps_send_i(ip + 1) = ip
        maps_recv_i(ip + 1) = mod(NPRO_I - ip, NPRO_I)
    end do
    maps_send_i = cshift(maps_send_i, ims_pro_i)
    maps_recv_i = cshift(maps_recv_i, -ims_pro_i)

    ! ---------------------------------------------------------------
    ! K-direction taint source — replicates FABRIC_DIRECT init exactly:
    !   MPI_Comm_split_type(ims_comm_z, SHARED) -> MPI_Win_allocate_shared
    ! On Cray MPICH this taints ims_comm_z and ims_comm_x.
    ! ---------------------------------------------------------------
    call MPI_Comm_split_type(ims_comm_z, MPI_COMM_TYPE_SHARED, ims_pro_k, &
                             MPI_INFO_NULL, shmem_comm_k, ierr_mpi)
    call MPI_Win_allocate_shared( &
        int(CHUNK, MPI_ADDRESS_KIND) * int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
        int(c_sizeof(1.0_dp)), MPI_INFO_NULL, shmem_comm_k, k_win_baseptr_c, k_win, ierr_mpi)

    if (ims_rank == 0) then
        write(*,'(a,i0,a,i0,a,i0)') &
            '[INIT] decomp=', NPRO_K, 'x', NPRO_I, '  chunk=', CHUNK
        write(*,'(a)') '[INIT] K-direction MPI_Win_allocate_shared done (taint now active)'
        flush(6)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)

    ! ---------------------------------------------------------------
    ! Staging and request arrays
    ! ---------------------------------------------------------------
    allocate(staging(2 * CHUNK * NPRO_I))
    allocate(req(2 * NPRO_I), stat(2 * NPRO_I))

    ! ---------------------------------------------------------------
    ! Four test phases
    ! ---------------------------------------------------------------
    do phase = phase_lo, phase_hi

        ! --- Select communicator for this phase ---
        select case (phase)
        case (1, 2)
            ! dup(ims_comm_x) — same as apu_async_mpi_comm_i in the real code
            call MPI_Comm_dup(ims_comm_x, test_comm, ierr_mpi)
        case (3, 4)
            ! Fresh MPI_Comm_split from MPI_COMM_WORLD: color=k, key=i
            ! Local rank in this comm = ims_pro_i (key ordering is ascending)
            call MPI_Comm_split(MPI_COMM_WORLD, ims_pro_k, ims_pro_i, test_comm, ierr_mpi)
        end select

        call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)
        if (ims_rank == 0) then
            write(*,*)
            write(*,'(a,i0,a)') '=== Phase ', phase, ' ==='
            select case (phase)
            case (1)
                write(*,'(a)') '    comm=dup(ims_comm_x)   hipDevSync AFTER  IRECVs [current code]'
            case (2)
                write(*,'(a)') '    comm=dup(ims_comm_x)   hipDevSync BEFORE IRECVs'
            case (3)
                write(*,'(a)') '    comm=MPI_Comm_split     hipDevSync AFTER  IRECVs'
            case (4)
                write(*,'(a)') '    comm=MPI_Comm_split     hipDevSync BEFORE IRECVs [proposed fix]'
            end select
            flush(6)
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)

        ! --- Fill send slots in second half of staging ---
        ! Each slot m: send real(ims_rank) so receiver can verify global rank
        do m = 0, NPRO_I - 1
            staging(CHUNK * NPRO_I + CHUNK * m + 1 : CHUNK * NPRO_I + CHUNK * (m + 1)) &
                = real(ims_rank, dp)
        end do
        staging(1 : CHUNK * NPRO_I) = -999.0_dp   ! poison recv area

        ! --- Phase 2 / 4: flush GPU L2 -> HBM BEFORE IRECVs ---
        if (phase == 2 .or. phase == 4) then
#ifdef USE_APU
            hip_err = hipDeviceSynchronize()
#endif
            write(*,'(a,i0,a,i3,a,i2,a,i2,a)') &
                '[P', phase, '_pre] rank ', ims_rank, &
                ' (k=', ims_pro_k, ',i=', ims_pro_i, ') hipDevSync done before IRECVs'
            flush(6)
        end if

        ! --- Post IRECVs (recv into first half of staging, slot = local rank of sender) ---
        l = 0
        do m = 1, NPRO_I
            nr  = maps_recv_i(m) + 1   ! local rank in test_comm we expect to recv from
            ipr = nr - 1
            l = l + 1
            call MPI_IRECV(staging((nr - 1) * CHUNK + 1), CHUNK, &
                           MPI_DOUBLE_PRECISION, ipr, 77, test_comm, req(l), ierr_mpi)
        end do
        write(*,'(a,i0,a,i3,a,i2,a,i2,a)') &
            '[P', phase, '_S4] rank ', ims_rank, &
            ' (k=', ims_pro_k, ',i=', ims_pro_i, ') IRECVs posted'
        flush(6)

        ! --- Phase 1 / 3: flush GPU L2 -> HBM AFTER IRECVs (current code order) ---
        if (phase == 1 .or. phase == 3) then
#ifdef USE_APU
            hip_err = hipDeviceSynchronize()
#endif
            write(*,'(a,i0,a,i3,a,i2,a,i2,a)') &
                '[P', phase, '_Ssync] rank ', ims_rank, &
                ' (k=', ims_pro_k, ',i=', ims_pro_i, ') hipDevSync done after IRECVs'
            flush(6)
        end if

        ! --- Post ISENDs from send half of staging ---
        do m = 1, NPRO_I
            ns  = maps_send_i(m) + 1   ! local rank in test_comm to send to
            ips = ns - 1
            l = l + 1
            call MPI_ISEND(staging(CHUNK * NPRO_I + (ns - 1) * CHUNK + 1), CHUNK, &
                           MPI_DOUBLE_PRECISION, ips, 77, test_comm, req(l), ierr_mpi)
        end do
        write(*,'(a,i0,a,i3,a,i2,a,i2,a)') &
            '[P', phase, '_S8] rank ', ims_rank, &
            ' (k=', ims_pro_k, ',i=', ims_pro_i, ') ISENDs posted'
        flush(6)

        call MPI_WAITALL(l, req, stat, ierr_mpi)
        write(*,'(a,i0,a,i3,a,i2,a,i2,a)') &
            '[P', phase, '_S9] rank ', ims_rank, &
            ' (k=', ims_pro_k, ',i=', ims_pro_i, ') WAITALL done'
        flush(6)

        ! --- Verify received data ---
        ! Slot nr-1 (local rank nr-1) holds data from the peer whose ims_pro_i = nr-1.
        ! That peer's global rank = ims_pro_k * NPRO_I + (nr-1).
        n_errors = 0
        do m = 1, NPRO_I
            nr = maps_recv_i(m) + 1
            expected_rank = ims_pro_k * NPRO_I + (nr - 1)
            do i = 1, CHUNK
                if (abs(staging((nr - 1) * CHUNK + i) - real(expected_rank, dp)) > 0.5_dp) &
                    n_errors = n_errors + 1
            end do
        end do

        if (n_errors == 0) then
            write(*,'(a,i0,a,i3,a,i2,a,i2,a)') &
                '[P', phase, '] rank ', ims_rank, &
                ' (k=', ims_pro_k, ',i=', ims_pro_i, ') PASS'
        else
            write(*,'(a,i0,a,i3,a,i2,a,i2,a,i0,a)') &
                '[P', phase, '] rank ', ims_rank, &
                ' (k=', ims_pro_k, ',i=', ims_pro_i, ') FAIL  errors=', n_errors
        end if
        flush(6)

        call MPI_Comm_free(test_comm, ierr_mpi)
        call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)

    end do   ! phase

    ! ---------------------------------------------------------------
    ! Summary from rank 0
    ! ---------------------------------------------------------------
    call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)
    if (ims_rank == 0) then
        write(*,*)
        write(*,'(a)') '=== All phases complete ==='
        write(*,'(a)') 'Search output for: hangs (missing S9), FAIL, or PASS'
        write(*,'(a)') 'Phase 4 (MPI_Comm_split + hipDevSync before) should be all PASS'
    end if

    ! ---------------------------------------------------------------
    ! Cleanup
    ! ---------------------------------------------------------------
    call MPI_Win_free(k_win, ierr_mpi)
    call MPI_Comm_free(shmem_comm_k, ierr_mpi)
    call MPI_Comm_free(ims_comm_x, ierr_mpi)
    call MPI_Comm_free(ims_comm_z, ierr_mpi)
    call MPI_Comm_free(ims_comm_xz, ierr_mpi)
    deallocate(staging, req, stat, maps_send_i, maps_recv_i)
    call MPI_Finalize(ierr_mpi)

end program vmpi_fabricdirect_commtest
