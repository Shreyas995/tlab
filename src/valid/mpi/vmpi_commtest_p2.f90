! vmpi_commtest_p2.f90
! Phase 2: comm = dup(ims_comm_x),  hipDevSync BEFORE IRECVs
! Expected: may fix the hang — if passes, the sync ORDER alone is the fix
! Output: fort.2XX where XX = rank (fort.200 to fort.247)
!
! Build:  hipcc -c hip_write_fence.hip -o hip_write_fence.o
!         ftn -DUSE_APU vmpi_commtest_p2.f90 hip_write_fence.o -lamdhip64 -o vmpi_commtest_p2.x

program vmpi_commtest_p2
    use mpi_f08
    use iso_c_binding
    implicit none

    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: NPRO_K = 8, NPRO_I = 6, CHUNK = 500
    integer, parameter :: FUNIT_BASE = 200   ! fort.200 + rank

#ifdef USE_APU
    interface
        function hipDeviceSynchronize() bind(C, name='hipDeviceSynchronize') result(r)
            use iso_c_binding; integer(c_int) :: r
        end function
    end interface
    integer(c_int) :: hip_err
#endif

    integer :: ims_rank, ims_nprocs, ims_pro_i, ims_pro_k
    integer :: coord(2), dims(2), ierr
    logical :: period(2), reorder, remain_dims(2)
    type(MPI_Comm) :: comm_xz, comm_x, comm_z, shmem_k, test_comm
    type(MPI_Win)  :: k_win
    type(c_ptr)    :: k_baseptr
    integer(MPI_ADDRESS_KIND) :: win_sz
    integer :: win_disp
    integer, allocatable :: maps_s(:), maps_r(:)
    real(dp), allocatable :: staging(:)
    type(MPI_Request), allocatable :: req(:)
    type(MPI_Status),  allocatable :: stat(:)
    integer :: ip, m, l, nr, ipr, ns, ips, nerr, i, funit, expected

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, ims_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, ims_nprocs, ierr)
    funit = FUNIT_BASE + ims_rank

    if (ims_nprocs /= NPRO_I * NPRO_K) then
        if (ims_rank == 0) write(*,*) 'ERROR: need 48 ranks, got', ims_nprocs
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    ! --- Cartesian topology (identical to tlab_mpi_procs.f90) ---
    dims(1) = NPRO_K; dims(2) = NPRO_I; period = .true.; reorder = .false.
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, reorder, comm_xz, ierr)
    call MPI_Cart_coords(comm_xz, ims_rank, 2, coord, ierr)
    ims_pro_k = coord(1); ims_pro_i = coord(2)

    remain_dims(1) = .false.; remain_dims(2) = .true.
    call MPI_Cart_sub(comm_xz, remain_dims, comm_x, ierr)

    remain_dims(1) = .true.; remain_dims(2) = .false.
    call MPI_Cart_sub(comm_xz, remain_dims, comm_z, ierr)

    ! --- Send/recv maps ---
    allocate(maps_s(NPRO_I), maps_r(NPRO_I))
    do ip = 0, NPRO_I - 1
        maps_s(ip+1) = ip
        maps_r(ip+1) = mod(NPRO_I - ip, NPRO_I)
    end do
    maps_s = cshift(maps_s,  ims_pro_i)
    maps_r = cshift(maps_r, -ims_pro_i)

    ! --- K-direction taint source ---
    call MPI_Comm_split_type(comm_z, MPI_COMM_TYPE_SHARED, ims_pro_k, &
                             MPI_INFO_NULL, shmem_k, ierr)
    call MPI_Win_allocate_shared( &
        int(CHUNK, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
        int(c_sizeof(1.0_dp)), MPI_INFO_NULL, shmem_k, k_baseptr, k_win, ierr)
    write(funit,'(a,i3,a,i2,a,i2,a)') '[INIT] rank', ims_rank, &
        ' k=', ims_pro_k, ' i=', ims_pro_i, ' K-shmem alloc done (taint active)'
    flush(funit)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! --- Test comm: dup(comm_x) ---
    call MPI_Comm_dup(comm_x, test_comm, ierr)
    write(funit,'(a)') '[COMM] dup(ims_comm_x) created'
    flush(funit)

    ! --- Staging buffer ---
    allocate(staging(2*CHUNK*NPRO_I), req(2*NPRO_I), stat(2*NPRO_I))
    do m = 0, NPRO_I-1
        staging(CHUNK*NPRO_I + CHUNK*m+1 : CHUNK*NPRO_I + CHUNK*(m+1)) = real(ims_rank, dp)
    end do
    staging(1:CHUNK*NPRO_I) = -999.0_dp

    ! --- IFR pattern: sync BEFORE IRECVs ---
#ifdef USE_APU
    hip_err = hipDeviceSynchronize()
    write(funit,'(a,i0)') '[Ssync] hipDevSync before IRECVs, err=', int(hip_err)
    flush(funit)
#endif

    l = 0
    do m = 1, NPRO_I
        nr = maps_r(m)+1; ipr = nr-1; l = l+1
        call MPI_IRECV(staging((nr-1)*CHUNK+1), CHUNK, MPI_DOUBLE_PRECISION, &
                       ipr, 77, test_comm, req(l), ierr)
    end do
    write(funit,'(a)') '[S4] IRECVs posted'
    flush(funit)

    do m = 1, NPRO_I
        ns = maps_s(m)+1; ips = ns-1; l = l+1
        call MPI_ISEND(staging(CHUNK*NPRO_I + (ns-1)*CHUNK+1), CHUNK, &
                       MPI_DOUBLE_PRECISION, ips, 77, test_comm, req(l), ierr)
    end do
    write(funit,'(a)') '[S8] ISENDs posted'
    flush(funit)

    call MPI_WAITALL(l, req, stat, ierr)
    write(funit,'(a)') '[S9] WAITALL done'
    flush(funit)

    ! --- Verify ---
    nerr = 0
    do m = 1, NPRO_I
        nr = maps_r(m)+1
        expected = ims_pro_k*NPRO_I + (nr-1)
        do i = 1, CHUNK
            if (abs(staging((nr-1)*CHUNK+i) - real(expected,dp)) > 0.5_dp) nerr = nerr+1
        end do
    end do
    if (nerr == 0) then
        write(funit,'(a)') '[RESULT] PASS'
    else
        write(funit,'(a,i0,a)') '[RESULT] FAIL  errors=', nerr
    end if
    flush(funit)

    call MPI_Comm_free(test_comm, ierr)
    call MPI_Win_free(k_win, ierr)
    call MPI_Finalize(ierr)
end program
