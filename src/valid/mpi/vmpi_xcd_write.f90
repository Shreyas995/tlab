! vmpi_xcd_write.f90
!
! Validates the session-2 FABRIC_DIRECT I-direction setup on Hunter (MI300A).
! Replicates the EXACT window-allocation order from tlab_mpi_transpose.f90 and
! tests two GPU-write approaches to diagnose what causes the signal-6 crash.
!
! What it tests
! -------------
! Comm/window order (mirrors session-2 production):
!   K-dup  -> I-window(ims_comm_x) -> I-dup -> K-split -> K-window
! This is the order that SHOULD give uniform VAs for ims_comm_x (no prior
! split-type on ims_comm_x when the I-window is allocated).
!
! VA uniformity check:
!   Every rank in the I-comm queries rank 0's VA. With uniform allocation,
!   rank0_VA + ims_pro_i * segsize must equal OWN segment VA for every rank.
!
! Phase 1 -- !$omp target write (production FABRIC_DIRECT exact replica):
!   Writes send_buf into every I-peer's recv slot via a single !$omp target
!   region over apu_async_all_i (contiguous span built from rank-0's VA).
!   If this crashes (signal 6): the shmem window is not GPU-accessible
!   without hipHostRegister, even with uniform VAs.
!   If this gives wrong data (no crash): cross-XCD cache incoherence.
!   If this passes: the session-2 production fix is complete.
!
! Phase 2 -- hipHostRegister + hip_write_with_fence (Gemini hypothesis fix):
!   Registers the full shared span with the HIP runtime, then uses the
!   hip_write_fence_kernel (__threadfence_system) for system-scope flush.
!   If Phase 1 failed but Phase 2 passes: integrate HIP fence into production.
!
! Build: make -f Makefile.xcd_write
! Run:   ./vmpi_xcd_write.x <npro_k> <npro_i> [chunk [phase]]
!        phase: 0=both (default), 1=Phase1 only, 2=Phase2 only
! Example: mpirun -np 48 ./vmpi_xcd_write.x 8 6 512
!          mpirun -np 48 ./vmpi_xcd_write.x 8 6 512 2   (skip P1, run P2 only)

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
    type(MPI_Comm)  :: mpi_comm_i, mpi_comm_k, shmem_comm_k
    type(MPI_Win)   :: win_i, win_k
    type(MPI_Group) :: dir_group_i, world_group

    integer :: ims_pro, ims_npro
    integer :: npro_i, npro_k
    integer :: ims_pro_i, ims_pro_k
    integer :: shmem_size_k
    integer :: ims_err, hip_reg_err

    ! -------------------------------------------------------------------
    ! Window pointers and VA check
    ! -------------------------------------------------------------------
    type(c_ptr) :: win_baseptr, own_seg_ptr, peer_i_0_cptr
    integer(MPI_ADDRESS_KIND) :: seg_size_abi, win_size_abi
    integer :: disp_unit_i

    integer(MPI_ADDRESS_KIND) :: rank0_VA, own_VA, expected_own_VA

    real(dp), pointer, contiguous :: recv_i(:) => null()  ! own I recv buffer
    real(dp), pointer, contiguous :: all_i(:)  => null()  ! full span all I-peers

    ! -------------------------------------------------------------------
    ! Test data
    ! -------------------------------------------------------------------
    real(dp), allocatable :: send_buf(:)
    integer :: chunk, seg_size_elems, flat_off
    integer :: phase_select   ! 0=both, 1=P1 only, 2=P2 only

    ! Verification
    integer, allocatable :: i_global_ranks(:), local_ranks(:)
    real(dp) :: expected, actual
    integer :: g_m, m, j
    integer :: errors_p1, errors_p2, total_errors, total_errors_global
    character(len=64) :: arg

    ! -------------------------------------------------------------------
    ! MPI init
    ! -------------------------------------------------------------------
    call MPI_Init(ims_err)
    call MPI_Comm_rank(MPI_COMM_WORLD, ims_pro, ims_err)
    call MPI_Comm_size(MPI_COMM_WORLD, ims_npro, ims_err)

    write(1000+ims_pro,*) 'MPI_Init done, rank=', ims_pro, ' npro=', ims_npro
    flush(1000+ims_pro)

    if (command_argument_count() < 2) then
        if (ims_pro == 0) write(*,*) &
            'Usage: vmpi_xcd_write.x <npro_k> <npro_i> [chunk [phase]]'
        call MPI_Finalize(ims_err); stop
    end if
    call get_command_argument(1, arg); read(arg,*) npro_k
    call get_command_argument(2, arg); read(arg,*) npro_i
    chunk = 512
    if (command_argument_count() >= 3) then
        call get_command_argument(3, arg); read(arg,*) chunk
    end if
    phase_select = 0
    if (command_argument_count() >= 4) then
        call get_command_argument(4, arg); read(arg,*) phase_select
    end if

    if (npro_i * npro_k /= ims_npro) then
        if (ims_pro == 0) write(*,'(a,3i6)') &
            'ERROR: npro_k*npro_i /= nproc:', npro_k, npro_i, ims_npro
        call MPI_Finalize(ims_err); stop
    end if
    write(1000+ims_pro,*) 'npro_k=', npro_k, ' npro_i=', npro_i, &
        ' chunk=', chunk, ' phase_select=', phase_select
    flush(1000+ims_pro)

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
    write(1000+ims_pro,*) 'pro_i=', ims_pro_i, ' pro_k=', ims_pro_k
    flush(1000+ims_pro)

    ! -------------------------------------------------------------------
    ! Session-2 production comm/window order:
    !   Step 1: K-dup      (before any window)
    !   Step 2: I-window   (on ims_comm_x -- NO prior split-type on this comm)
    !   Step 3: I-dup      (after I-window)
    !   Step 4: K-split    (MPI_Comm_split_type on ims_comm_z)
    !   Step 5: K-window   (on shmem_comm_k -- the potentially tainting alloc)
    ! -------------------------------------------------------------------

    ! Step 1: K-dup
    call MPI_Comm_dup(ims_comm_z, mpi_comm_k, ims_err)
    write(1000+ims_pro,*) 'Step1: K-dup done, err=', ims_err
    flush(1000+ims_pro)

    ! Step 2: I-window on full ims_comm_x
    ! Per-rank segment holds npro_i*chunk elements (slot per I-peer).
    seg_size_elems = npro_i * chunk
    win_size_abi = int(seg_size_elems, MPI_ADDRESS_KIND) &
                 * int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND)
    call MPI_Win_allocate_shared(win_size_abi, int(c_sizeof(1.0_dp)), MPI_INFO_NULL, &
                                  ims_comm_x, win_baseptr, win_i, ims_err)
    write(1000+ims_pro,*) 'Step2: I-window on ims_comm_x done, err=', ims_err
    flush(1000+ims_pro)

    ! Bind own recv_i via Win_shared_query (Bug-A fix: win_baseptr unreliable on Cray)
    call MPI_Win_shared_query(win_i, ims_pro_i, seg_size_abi, disp_unit_i, own_seg_ptr, ims_err)
    call c_f_pointer(own_seg_ptr, recv_i, [seg_size_elems])
    write(1000+ims_pro,*) 'I recv_i bound via Win_shared_query own rank'
    flush(1000+ims_pro)

    ! Step 3: I-dup (AFTER I-window, before K-split)
    call MPI_Comm_dup(ims_comm_x, mpi_comm_i, ims_err)
    write(1000+ims_pro,*) 'Step3: I-dup done, err=', ims_err
    flush(1000+ims_pro)

    ! Step 4: K-split (XCD-level shmem sub-comm of ims_comm_z)
    call MPI_Comm_split_type(ims_comm_z, MPI_COMM_TYPE_SHARED, ims_pro_k, &
                             MPI_INFO_NULL, shmem_comm_k, ims_err)
    call MPI_Comm_size(shmem_comm_k, shmem_size_k, ims_err)
    write(1000+ims_pro,*) 'Step4: K-split done, shmem_size_k=', shmem_size_k
    flush(1000+ims_pro)

    ! Step 5: K-window on shmem_comm_k (the potentially tainting allocation)
    win_size_abi = int(npro_k * chunk, MPI_ADDRESS_KIND) &
                 * int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND)
    call MPI_Win_allocate_shared(win_size_abi, int(c_sizeof(1.0_dp)), MPI_INFO_NULL, &
                                  shmem_comm_k, win_baseptr, win_k, ims_err)
    write(1000+ims_pro,*) 'Step5: K-window on shmem_comm_k done, err=', ims_err
    flush(1000+ims_pro)

    ! Global rank map for I-comm (for data verification)
    allocate(i_global_ranks(0:npro_i-1), local_ranks(0:npro_i-1))
    call MPI_Comm_group(ims_comm_x, dir_group_i, ims_err)
    call MPI_Comm_group(MPI_COMM_WORLD, world_group, ims_err)
    do m = 0, npro_i-1; local_ranks(m) = m; end do
    call MPI_Group_translate_ranks(dir_group_i, npro_i, local_ranks, &
                                   world_group, i_global_ranks, ims_err)

    allocate(send_buf(chunk))
    do j = 1, chunk
        send_buf(j) = dble(ims_pro) * 1000.0d0 + dble(j)
    end do
    write(1000+ims_pro,*) 'send_buf filled'
    flush(1000+ims_pro)

    ! -------------------------------------------------------------------
    ! VA uniformity check
    ! Query rank 0's segment VA from every rank in the I-comm.
    ! With uniform allocation: rank0_VA + ims_pro_i*segsize == own_VA.
    ! -------------------------------------------------------------------
    call MPI_Win_shared_query(win_i, 0, seg_size_abi, disp_unit_i, peer_i_0_cptr, ims_err)
    rank0_VA     = transfer(peer_i_0_cptr, rank0_VA)
    own_VA       = transfer(own_seg_ptr,   own_VA)
    expected_own_VA = rank0_VA &
        + int(ims_pro_i,   MPI_ADDRESS_KIND) &
        * int(seg_size_elems, MPI_ADDRESS_KIND) &
        * int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND)

    write(500+ims_pro,'(a,i4,3(a,i22),a,l1)') &
        '[VA_UNIFORM] PE', ims_pro, &
        '  rank0_VA=', rank0_VA, &
        '  own_VA=', own_VA, &
        '  expected_own=', expected_own_VA, &
        '  MATCH=', (own_VA == expected_own_VA)
    flush(500+ims_pro)
    call MPI_Barrier(mpi_comm_i, ims_err)

    ! Build all_i: contiguous pointer over all I-peers' segments.
    ! Valid only if VAs are uniform (peer_i_0_cptr is same value on all ranks).
    call c_f_pointer(peer_i_0_cptr, all_i, [seg_size_elems * npro_i])

    ! flat_off: byte-offset of OUR slot within any peer's segment (in elements)
    flat_off = ims_pro_i * chunk

    errors_p1 = 0; errors_p2 = 0

    ! ================================================================
    ! PHASE 1: !$omp target GPU write (production FABRIC_DIRECT replica)
    ! ================================================================
    if (phase_select == 0 .or. phase_select == 1) then
        recv_i = 0.0_dp
        call MPI_Barrier(MPI_COMM_WORLD, ims_err)
        write(1000+ims_pro,*) 'P1: barrier passed, starting !$omp target writes'
        flush(1000+ims_pro)

        ! Write OUR send_buf into every I-peer's slot at our pro_i offset.
        ! Mirrors production: apu_async_all_i(m*apu_async_size_i + flat_off + j)
        !$omp target teams distribute parallel do collapse(2)
        do m = 0, npro_i - 1
            do j = 1, chunk
                all_i(m * seg_size_elems + flat_off + j) = send_buf(j)
            end do
        end do
        !$omp end target teams distribute parallel do

        write(500+ims_pro,'(a,i4)') '[P1_OMP_WRITE] PE', ims_pro
        flush(500+ims_pro)
        call MPI_Barrier(mpi_comm_i, ims_err)
        write(500+ims_pro,'(a,i4)') '[P1_BARRIER_PASSED] PE', ims_pro
        flush(500+ims_pro)

        ! Verify: recv_i(m*chunk + j) == i_global_ranks(m)*1000 + j
        do m = 0, npro_i - 1
            g_m = i_global_ranks(m)
            do j = 1, chunk
                expected = dble(g_m) * 1000.0d0 + dble(j)
                actual   = recv_i(m * chunk + j)
                if (abs(actual - expected) > 1.0d-6) then
                    errors_p1 = errors_p1 + 1
                    if (errors_p1 <= 3) then
                        write(500+ims_pro,'(a,i4,a,i3,a,i5,a,g20.12,a,g20.12)') &
                            '[P1_ERR] PE', ims_pro, ' peer=', m, ' j=', j, &
                            '  exp=', expected, '  got=', actual
                        flush(500+ims_pro)
                    end if
                end if
            end do
        end do
        if (errors_p1 == 0) then
            write(500+ims_pro,'(a,i4)') '[P1_PASS] PE', ims_pro
        else
            write(500+ims_pro,'(a,i4,a,i6)') '[P1_FAIL] PE', ims_pro, '  errors=', errors_p1
        end if
        flush(500+ims_pro)
    end if

    ! ================================================================
    ! PHASE 2: hipHostRegister + hip_write_with_fence
    ! Tests the Gemini hypothesis: !$omp target may not issue a
    ! system-scope cache flush; __threadfence_system() in the HIP
    ! kernel guarantees all XCDs see the written data.
    ! ================================================================
    if (phase_select == 0 .or. phase_select == 2) then
        recv_i = 0.0_dp
        call MPI_Barrier(MPI_COMM_WORLD, ims_err)
        write(1000+ims_pro,*) 'P2: barrier passed, registering shmem with HIP...'
        flush(1000+ims_pro)

        ! Register the full contiguous I-window span with the HIP runtime.
        ! This enables GPU kernel access to the MPI shared-memory region.
        hip_reg_err = hipHostRegister(peer_i_0_cptr, &
            int(seg_size_elems * npro_i, c_size_t) &
            * int(c_sizeof(1.0_dp), c_size_t), 0_c_int)
        write(1000+ims_pro,*) 'hipHostRegister span, err=', hip_reg_err
        flush(1000+ims_pro)
        write(500+ims_pro,'(a,i4,a,i4)') '[P2_HIPREGISTER] PE', ims_pro, '  err=', hip_reg_err
        flush(500+ims_pro)

        ! Write OUR chunk into each I-peer's slot with system-scope fence.
        do m = 0, npro_i - 1
            call hip_write_with_fence(send_buf(1), &
                all_i(m * seg_size_elems + flat_off + 1), int(chunk, c_int))
        end do

        write(500+ims_pro,'(a,i4)') '[P2_HIP_WRITES_DONE] PE', ims_pro
        flush(500+ims_pro)
        call MPI_Barrier(mpi_comm_i, ims_err)
        write(500+ims_pro,'(a,i4)') '[P2_BARRIER_PASSED] PE', ims_pro
        flush(500+ims_pro)

        ! Verify
        do m = 0, npro_i - 1
            g_m = i_global_ranks(m)
            do j = 1, chunk
                expected = dble(g_m) * 1000.0d0 + dble(j)
                actual   = recv_i(m * chunk + j)
                if (abs(actual - expected) > 1.0d-6) then
                    errors_p2 = errors_p2 + 1
                    if (errors_p2 <= 3) then
                        write(500+ims_pro,'(a,i4,a,i3,a,i5,a,g20.12,a,g20.12)') &
                            '[P2_ERR] PE', ims_pro, ' peer=', m, ' j=', j, &
                            '  exp=', expected, '  got=', actual
                        flush(500+ims_pro)
                    end if
                end if
            end do
        end do
        if (errors_p2 == 0) then
            write(500+ims_pro,'(a,i4)') '[P2_PASS] PE', ims_pro
        else
            write(500+ims_pro,'(a,i4,a,i6)') '[P2_FAIL] PE', ims_pro, '  errors=', errors_p2
        end if
        flush(500+ims_pro)
    end if

    ! ================================================================
    ! Report
    ! ================================================================
    total_errors = errors_p1 + errors_p2
    call MPI_Reduce(total_errors, total_errors_global, 1, MPI_INTEGER, MPI_SUM, 0, &
                    MPI_COMM_WORLD, ims_err)

    write(500+ims_pro,'(a,i4,a,i6,a,i6)') &
        '[RESULT] PE', ims_pro, '  p1_errors=', errors_p1, '  p2_errors=', errors_p2
    flush(500+ims_pro)

    if (ims_pro == 0) then
        write(*,'(a)') '================================================'
        write(*,'(a,i4,a,i4,a,i5,a,i2)') &
            ' npro_k=', npro_k, '  npro_i=', npro_i, &
            '  chunk=', chunk, '  phase=', phase_select
        if (total_errors_global == 0) then
            write(*,'(a)') ' PASS: all tested phases correct'
        else
            write(*,'(a,i8,a)') ' FAIL: ', total_errors_global, ' total mismatches'
        end if
        write(*,'(a)') ''
        write(*,'(a)') ' Interpretation guide (check fort.500+rank for details):'
        write(*,'(a)') '  VA_UNIFORM MATCH=F -> ordering still gives non-uniform VA'
        write(*,'(a)') '  P1_OMP_WRITE in log but no P1_BARRIER -> crash in !$omp target'
        write(*,'(a)') '    -> shmem not GPU-accessible without hipHostRegister'
        write(*,'(a)') '    -> re-run with phase=2 to test HIP fence approach'
        write(*,'(a)') '  P1_FAIL + P2_PASS -> cache incoherence; HIP fence required'
        write(*,'(a)') '  P1_PASS -> session-2 !$omp target approach works as-is'
        write(*,'(a)') '================================================'
    end if

    ! -------------------------------------------------------------------
    ! Cleanup
    ! -------------------------------------------------------------------
    call MPI_Win_free(win_i, ims_err)
    call MPI_Win_free(win_k, ims_err)
    call MPI_Comm_free(shmem_comm_k, ims_err)
    call MPI_Comm_free(mpi_comm_i, ims_err)
    call MPI_Comm_free(mpi_comm_k, ims_err)
    nullify(recv_i, all_i)
    deallocate(send_buf, i_global_ranks, local_ranks)
    call MPI_Finalize(ims_err)

end program vmpi_xcd_write
