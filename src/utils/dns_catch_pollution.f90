#include "dns_error.h"

!########################################################################
!#
!# Cheap per-rank "polluted field" sentinel used to localize the random
!# one-step blow-up of the production run (CFL/dilatation -> ~1e20 in a
!# single iteration). Scans a(1:n) LOCALLY (no MPI) for values whose
!# magnitude exceeds a threshold or are NaN/Inf. If any are found, the
!# TRIPPING rank records the stage tag / iteration / substep / value to
!# tlab.err and aborts the whole job. On the happy path it does only one
!# local reduction (on the GPU under USE_APU) and returns -- no
!# communication and no file I/O, so it is safe to leave in production.
!#
!# Detection uses .not.(abs(a) <= thr) on purpose: it is .true. for huge
!# values AND for NaN/Inf (a plain abs(a) > thr would miss NaN, because
!# NaN > thr evaluates to .false.).
!#
!# DNS_CATCH_POLLUTION    uses the default physical-field threshold (1e6:
!#                        healthy |u|~O(10), |p|~O(100)).
!# DNS_CATCH_POLLUTION_HI takes an explicit threshold, for non-physical-
!#                        scale arrays (e.g. the Poisson forcing ~1e7 or
!#                        spectral intermediates ~1e10) where 1e6 would
!#                        false-trip; pass a high value (e.g. 1e15) that
!#                        still sits far below the 1e20-1e40 blow-up.
!# Both are plain external subroutines with all-required arguments, so no
!# explicit interface is needed at the call sites.
!#
!########################################################################
subroutine DNS_CATCH_POLLUTION(tag, a, n, isub)
    use TLab_Constants, only: wp, wi

    implicit none

    character(len=*), intent(in) :: tag         ! short stage/subroutine label
    integer(wi), intent(in) :: n                ! number of elements to scan
    real(wp), intent(in) :: a(n)                ! flattened field (hq, q, tmp1, ...)
    integer(wi), intent(in) :: isub             ! RK substep (-1 if not applicable)

    call DNS_CATCH_POLLUTION_HI(tag, a, n, isub, 1.0e6_wp)

    return

end subroutine DNS_CATCH_POLLUTION

!########################################################################
!########################################################################
subroutine DNS_CATCH_POLLUTION_HI(tag, a, n, isub, thr)
    use TLab_Constants, only: wp, wi, efile
    use TLab_Time, only: itime
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_pro
#endif

    implicit none

    character(len=*), intent(in) :: tag         ! short stage/subroutine label
    integer(wi), intent(in) :: n                ! number of elements to scan
    real(wp), intent(in) :: a(n)                ! flattened field (hq, q, tmp1, ...)
    integer(wi), intent(in) :: isub             ! RK substep (-1 if not applicable)
    real(wp), intent(in) :: thr                 ! magnitude threshold (huge OR NaN/Inf trips it)

#ifdef USE_APU
    !$omp requires unified_shared_memory
#endif

    ! -----------------------------------------------------------------------
    real(wp) vmax
    integer(wi) ij, nbad
    character(len=256) line
#ifndef USE_MPI
    integer, parameter :: ims_pro = 0
#endif

    ! #######################################################################
    vmax = 0.0_wp
    nbad = 0

#ifdef USE_APU
    !$omp target teams distribute parallel do private(ij) firstprivate(n, thr) &
    !$omp reduction(max:vmax) reduction(+:nbad)
#endif
    do ij = 1, n
        vmax = max(vmax, abs(a(ij)))
        if (.not. (abs(a(ij)) <= thr)) nbad = nbad + 1   ! huge OR NaN/Inf
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

    if (nbad > 0) then
        write (line, 1000) trim(adjustl(tag)), itime, isub, ims_pro, vmax, nbad
        call TLab_Write_ASCII(efile, line, .true.)  ! .true. => the tripping rank writes (any rank)
        call TLab_Stop(DNS_ERROR_POLLUTION)         ! writes code to tlab.err + MPI_Abort all ranks
    end if

    return

1000 format('[POLLUTION] tag=', a, ' it=', i7, ' sub=', i3, ' rank=', i5, ' max=', e13.6, ' nbad=', i12)

end subroutine DNS_CATCH_POLLUTION_HI

!########################################################################
!# PRINT-AND-CONTINUE probe: logs max|a| and the count of non-finite (NaN/Inf)
!# entries of a(1:n) to the per-rank flushed log (unit 500+rank / fort.5xx),
!# and ALWAYS returns -- it never aborts. Use it to flood the code with
!# write statements: run to the natural blow-up (DNS_CONTROL / NaN), then the
!# LAST line per rank whose max jumps to ~1e20 (or nbad>0) names the operation
!# that corrupted the field. No threshold => no false alarm.
!########################################################################
subroutine DNS_PRINT_MAXVAL(tag, a, n, isub)
    use TLab_Constants, only: wp, wi
    use TLab_Time, only: itime
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_pro
#endif

    implicit none

    character(len=*), intent(in) :: tag
    integer(wi), intent(in) :: n
    real(wp), intent(in) :: a(n)
    integer(wi), intent(in) :: isub

#ifdef USE_APU
    !$omp requires unified_shared_memory
#endif

    ! -----------------------------------------------------------------------
    real(wp) vmax
    integer(wi) ij, nbad, unit_num
#ifndef USE_MPI
    integer, parameter :: ims_pro = 0
#endif
#if defined(USE_APU) && defined(PROBE_SYNC_ONLY)
    ! HEISENBUG ISOLATION (build with -DPROBE_SYNC_ONLY): the probe does ONLY a
    ! GPU stream sync -- no array read, no log. Run the SAME restart three ways:
    !   (A) probe calls removed     -> crashes (the no-sentinel baseline)
    !   (B) full DNS_PRINT_MAXVAL   -> stable  (the masking sentinels)
    !   (C) this sync-only variant  -> ?
    ! C stable  => the masking agent is the device->host stream sync => the real
    !              bug is a MISSING hipDeviceSynchronize (a GPU->HBM coherency
    !              race), and the fix is to add that flush in production.
    ! C crashes => the masking agent is the full-array device READ (it faults/
    !              coheres the pages) => host->GPU page coherency or UB/OOB, NOT
    !              a GPU->HBM flush. Then hunt uninitialized/out-of-bounds memory.
    interface
        function hipDeviceSynchronize() bind(C, name='hipDeviceSynchronize') result(ierr)
            use iso_c_binding
            integer(c_int) :: ierr
        end function hipDeviceSynchronize
    end interface
    integer :: hip_sync_err
#endif

    ! #######################################################################
#ifndef DNS_DEBUG_PROBES
    return   ! DEBUG PROBES STRIPPED (production no-op: no GPU reduction, no flushed fort.5xx write).
             ! Build with -DDNS_DEBUG_PROBES to re-enable NWFR/NWBR/MAXVAL for validation.
#endif
#if defined(USE_APU) && defined(PROBE_SYNC_ONLY)
    hip_sync_err = hipDeviceSynchronize()
    return
#endif
    vmax = 0.0_wp
    nbad = 0

#ifdef USE_APU
    !$omp target teams distribute parallel do private(ij) firstprivate(n) &
    !$omp reduction(max:vmax) reduction(+:nbad)
#endif
    do ij = 1, n
        vmax = max(vmax, abs(a(ij)))
        if (.not. (abs(a(ij)) <= 1.0e30_wp)) nbad = nbad + 1   ! NaN/Inf/huge
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

    unit_num = 500 + ims_pro
    write (unit_num, 1100) trim(adjustl(tag)), itime, isub, vmax, nbad
    flush (unit_num)

    return

1100 format('[MAXVAL] tag=', a, ' it=', i7, ' sub=', i3, ' max=', e14.6, ' nbad=', i12)

end subroutine DNS_PRINT_MAXVAL

!########################################################################
!# CPU-ONLY twin of DNS_PRINT_MAXVAL. Scans a(1:n) on the HOST with NO
!# !$omp target -- it reads host-visible HBM directly, NOT through the GPU
!# L2 cache. Used to separate a writer-side fault from a reader-side stale
!# GPU read: place it on a shared-window recv buffer right after the closing
!# MPI_Win_fence, BEFORE any GPU read of that buffer. If this reads CLEAN
!# while the GPU read-back of the same buffer is garbage, the data in HBM is
!# correct and the GPU read-back returned stale cache (a reader-side bug).
!# Tag printed as [MAXVALC] so it greps apart from the GPU [MAXVAL] probe.
!########################################################################
subroutine DNS_PRINT_MAXVAL_CPU(tag, a, n, isub)
    use TLab_Constants, only: wp, wi
    use TLab_Time, only: itime
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_pro
#endif

    implicit none

    character(len=*), intent(in) :: tag
    integer(wi), intent(in) :: n
    real(wp), intent(in) :: a(n)
    integer(wi), intent(in) :: isub

    ! -----------------------------------------------------------------------
    real(wp) vmax
    integer(wi) ij, nbad, unit_num
#ifndef USE_MPI
    integer, parameter :: ims_pro = 0
#endif

    ! #######################################################################
#ifndef DNS_DEBUG_PROBES
    return   ! DEBUG PROBES STRIPPED (production no-op). Build -DDNS_DEBUG_PROBES to re-enable.
#endif
    vmax = 0.0_wp
    nbad = 0
    do ij = 1, n        ! pure host loop: reads HBM, never the GPU L2
        vmax = max(vmax, abs(a(ij)))
        if (.not. (abs(a(ij)) <= 1.0e30_wp)) nbad = nbad + 1
    end do

    unit_num = 500 + ims_pro
    write (unit_num, 1200) trim(adjustl(tag)), itime, isub, vmax, nbad
    flush (unit_num)

    return

1200 format('[MAXVALC] tag=', a, ' it=', i7, ' sub=', i3, ' max=', e14.6, ' nbad=', i12)

end subroutine DNS_PRINT_MAXVAL_CPU
