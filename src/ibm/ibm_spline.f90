#include "dns_const.h"
#include "dns_error.h"
#ifdef IBM_DEBUG
#ifdef USE_MPI

#endif
#endif
!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!# 2023/12/07 - Shreyas Deshpande
!#              Modified
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   cubic spline reconstruction in solid regions
!#
!#
!########################################################################
!# ARGUMENTS
!#
!#
!########################################################################
!# REQUIREMENTS
!#
!#
!########################################################################

subroutine IBM_SPLINE_XYZ(is, fld, fld_mod, g, isize_nob, isize_nob_be, nob, nob_b, nob_e, ibm_case)
    ! Steps C + D: geometry, LU(LHS), bisect indices, and writeback maps are precomputed in
    ! IBM_Spline_Cache (built once during IBM_INITIALIZE_GEOMETRY). At runtime this routine:
    !   1. picks the cache for the current direction (g%name)
    !   2. fills ya from fld using a precomputed source map
    !   3. calls CUBIC_SPLINE_PRECOMPUTED (skips LHS + TRIDFS + bisection)
    !   4. writes yb back to fld_mod via precomputed indices
    use TLab_Constants, only: efile, wp, wi
    use IBM_VARS, only: ibmscaljmin
    use IBM_Spline_Cache, only: ibm_spline_cache_t, cache_x, cache_y, cache_z, YA_SRC_SOLID
    use TLab_Memory, only: isize_field
    use FDM, only: fdm_dt
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Cubic_Splines, only: CUBIC_SPLINE_PRECOMPUTED
    use Tlab_Debug

    implicit none

    integer(wi), intent(in) :: is     ! scalar index; if 0, then velocity
    real(wp), dimension(isize_field), intent(in) :: fld
    real(wp), dimension(isize_field), intent(out) :: fld_mod
    type(fdm_dt), intent(in) :: g
    integer(wi), intent(in) :: isize_nob, isize_nob_be
    integer(wi), dimension(isize_nob), intent(in) :: nob
    integer(wi), dimension(isize_nob_be), intent(in) :: nob_b, nob_e
    integer(wi), dimension(isize_nob_be), intent(in) :: ibm_case

    type(ibm_spline_cache_t), pointer :: cache
    real(wp), allocatable :: ya_loc(:), yb_loc(:), wrk_loc(:)
    real(wp) :: ya_solid_val
    integer(wi) :: spl, k, ia, ib, idx
    integer(wi) :: max_ia, max_ib, wrk_size
    integer(wi), dimension(2) :: bc
    real(wp), dimension(2) :: bcval

    ! Aggregate debug counters — same set as the legacy reference so the diff comparison continues.
    integer(wi) :: dbg_total_ia, dbg_total_ib, dbg_total_iob
    real(wp)    :: dbg_sum_xa, dbg_sum_ya, dbg_sum_xb, dbg_sum_yb

    ! ================================================================== !

    ! Dispatch by direction
    select case (g%name)
    case ('x'); cache => cache_x
    case ('y'); cache => cache_y
    case ('z'); cache => cache_z
    case default
        call TLab_Write_ASCII(efile, 'IBM_SPLINE_XYZ: unknown direction '//trim(g%name))
        call TLab_Stop(DNS_ERROR_CUBIC_SPLINE)
        return
    end select

    !call TLab_Debug_Print_int('[I0_IBM_in] is=', is)
    !call TLab_Debug_Print_1D('[I0_IBM_in] sum(fld)=', fld)
    !call TLab_Debug_Print_int('[A0_skiplist] nactive=', cache%nactive)

    fld_mod = fld     ! never modify u,v,w,s directly !

    if (is /= 0) then
        ya_solid_val = ibmscaljmin(is)
    else
        ya_solid_val = 0.0_wp
    end if

    max_ia = cache%max_ia
    max_ib = cache%max_ib
    wrk_size = max(max_ia*11, 1)

    bc(:) = 2  ! fixed-first-derivative BC — matches what cache LU was factored for

    dbg_total_iob = 0
    dbg_total_ia  = 0; dbg_total_ib  = 0
    dbg_sum_xa    = 0.0_wp; dbg_sum_ya = 0.0_wp; dbg_sum_xb = 0.0_wp; dbg_sum_yb = 0.0_wp

    !$omp parallel &
    !$omp   private(ya_loc, yb_loc, wrk_loc) &
    !$omp   private(spl, k, ia, ib, idx, bcval)
    allocate (ya_loc(max(max_ia, 1)), yb_loc(max(max_ib, 1)), wrk_loc(wrk_size))
    !$omp do reduction(+:dbg_total_iob,dbg_total_ia,dbg_total_ib,dbg_sum_xa,dbg_sum_ya,dbg_sum_xb,dbg_sum_yb)
    do spl = 1, cache%n_splines
        ia = cache%spl_ia(spl)
        ib = cache%spl_ib(spl)

        ! Step C: fill ya from fld via the precomputed source map (velocity vs scalar)
        if (is == 0) then
            do k = 1, ia
                idx = cache%ya_src_idx_v(k, spl)
                if (idx > 0) then
                    ya_loc(k) = fld(idx)
                else
                    ya_loc(k) = ya_solid_val
                end if
            end do
        else
            do k = 1, ia
                idx = cache%ya_src_idx_s(k, spl)
                if (idx > 0) then
                    ya_loc(k) = fld(idx)
                else
                    ya_loc(k) = ya_solid_val
                end if
            end do
        end if

        ! Debug aggregates (geometry sums are static and could be precomputed; left inline for parity with reference path)
        dbg_total_iob = dbg_total_iob + 1
        dbg_total_ia  = dbg_total_ia  + ia
        dbg_total_ib  = dbg_total_ib  + ib
        dbg_sum_xa    = dbg_sum_xa    + sum(cache%spl_xa(1:ia, spl))
        dbg_sum_ya    = dbg_sum_ya    + sum(ya_loc(1:ia))
        dbg_sum_xb    = dbg_sum_xb    + sum(cache%spl_xb(1:ib, spl))

        if (.not. cache%spl_do_writeback(spl)) cycle  ! case 1 entries contribute to fingerprints but no spline math

        ! Step D: bcval still depends on current ya endpoints (RHS construction)
        bcval(1) = (ya_loc(2)  - ya_loc(1))     / cache%spl_dx(1,      spl)
        bcval(2) = (ya_loc(ia) - ya_loc(ia - 1)) / cache%spl_dx(ia - 1, spl)

        call CUBIC_SPLINE_PRECOMPUTED(bc, bcval, ia, ib, ya_loc(1:ia), &
                                      cache%spl_dx(1:ia - 1, spl), &
                                      cache%lu_aa(1:ia, spl), cache%lu_bb(1:ia, spl), cache%lu_cc(1:ia, spl), &
                                      cache%spl_bisect_idx(1:ib, spl), &
                                      cache%spl_xa(1:ia, spl), cache%spl_xb(1:ib, spl), &
                                      yb_loc(1:ib), wrk_loc)

        ! Force endpoint values to the physical BC (matches legacy override of yb_loc(1)/yb_loc(ib))
        yb_loc(1)  = ya_solid_val
        yb_loc(ib) = ya_solid_val

        ! Writeback via precomputed indices
        do k = 1, cache%spl_wb_count(spl)
            fld_mod(cache%spl_wb_fld_idx(k, spl)) = yb_loc(cache%spl_wb_yb_idx(k, spl))
        end do

        dbg_sum_yb = dbg_sum_yb + sum(yb_loc(1:ib))
    end do
    !$omp end do
    deallocate (ya_loc, yb_loc, wrk_loc)
    !$omp end parallel

    !call TLab_Debug_Print_int('[I1_IBM_agg] nactive=', cache%nactive)
    !call TLab_Debug_Print_int('[I1_IBM_agg] total_iob=', dbg_total_iob)
    !call TLab_Debug_Print_int('[I1_IBM_agg] total_ia=', dbg_total_ia)
    !call TLab_Debug_Print_int('[I1_IBM_agg] total_ib=', dbg_total_ib)
    !call TLab_Debug_Print_real('[I1_IBM_agg] sum_xa=', dbg_sum_xa)
    !call TLab_Debug_Print_real('[I1_IBM_agg] sum_ya=', dbg_sum_ya)
    !call TLab_Debug_Print_real('[I1_IBM_agg] sum_xb=', dbg_sum_xb)
    !call TLab_Debug_Print_real('[I1_IBM_agg] sum_yb=', dbg_sum_yb)
    !call TLab_Debug_Print_1D('[I2_IBM_out] sum(fld_mod)=', fld_mod)

    return
end subroutine IBM_SPLINE_XYZ

