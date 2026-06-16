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
!########################################################################
subroutine DNS_CATCH_POLLUTION(tag, a, n, isub)
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

#ifdef USE_APU
    !$omp requires unified_shared_memory
#endif

    ! -----------------------------------------------------------------------
    real(wp), parameter :: thr = 1.0e6_wp       ! healthy |u|~O(10), |p|~O(100); blow-up ~1e20+
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
    !$omp target teams distribute parallel do private(ij) firstprivate(n) &
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

end subroutine DNS_CATCH_POLLUTION
