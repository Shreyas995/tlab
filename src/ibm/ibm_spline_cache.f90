#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Step C + D cache for IBM_SPLINE_XYZ
!#
!# The IBM geometry is stationary, so per-spline geometry (xa, xb, dx),
!# the LHS tridiagonal LU factorization for CUBIC_SPLINE, and the bisect
!# index for each evaluation point can be computed once and reused.
!#
!# At runtime IBM_SPLINE_XYZ only:
!#   1. Fills ya from fld via ya_src_idx_v/s (Step C)
!#   2. Computes rhs, applies pre-factored TRIDSS, builds spline coeffs (Step D)
!#   3. Evaluates the cubic at the cached bisect indices
!#   4. Writes yb to fld_mod using the pre-baked writeback indices
!#
!# One cache per direction (x, y, z); IBM_SPLINE_XYZ dispatches by g%name.
!########################################################################

module IBM_Spline_Cache

    use TLab_Constants, only: efile, wp, wi
    use IBM_VARS, only: nflu, isize_wrk1d_ibm, nspl, ibmscaljmin
    use FDM, only: fdm_dt
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop

    implicit none

    private

    public :: ibm_spline_cache_t
    public :: IBM_Spline_Cache_Build
    public :: IBM_Spline_Cache_Build_All
    public :: cache_x, cache_y, cache_z

    ! Sentinel: ya source is "solid" (ya_solid = 0 for velocity, ibmscaljmin(is) for scalar).
    integer(wi), parameter, public :: YA_SRC_SOLID = -1

    type :: ibm_spline_cache_t
        logical :: built = .false.

        integer(wi) :: nactive = 0                       ! number of lines with nob/=0
        integer(wi) :: n_splines = 0                     ! total cached splines (iob entries)
        integer(wi) :: max_ia = 0                        ! max spline knot count over all splines
        integer(wi) :: max_ib = 0                        ! max evaluation-point count

        ! Per-spline metadata
        integer(wi), allocatable :: spl_ii(:)            ! (n_splines) line index ii
        integer(wi), allocatable :: spl_case(:)          ! (n_splines) case id 1..9
        integer(wi), allocatable :: spl_ia(:)            ! (n_splines) knot count
        integer(wi), allocatable :: spl_ib(:)            ! (n_splines) eval-point count
        integer(wi), allocatable :: spl_nob_b(:)         ! (n_splines) object start
        integer(wi), allocatable :: spl_nob_e(:)         ! (n_splines) object end
        logical, allocatable :: spl_do_writeback(:)      ! (n_splines) false for case 1 + iobs after case 1

        ! Geometry (Step C)
        real(wp), allocatable :: spl_xa(:, :)            ! (max_ia, n_splines)
        real(wp), allocatable :: spl_xb(:, :)            ! (max_ib, n_splines)
        real(wp), allocatable :: spl_dx(:, :)            ! (max_ia-1, n_splines)

        ! ya source map (Step C) — separate maps for velocity (is==0) vs scalar (is/=0)
        ! Positive value: fld index. YA_SRC_SOLID: use ya_solid for that scalar (constant per call).
        integer(wi), allocatable :: ya_src_idx_v(:, :)   ! (max_ia, n_splines)
        integer(wi), allocatable :: ya_src_idx_s(:, :)   ! (max_ia, n_splines)

        ! Pre-factored LHS tridiagonal for CUBIC_SPLINE (Step D)
        ! Built once per spline by calling CUBIC_SPLINE_LHS + TRIDFS at cache build.
        real(wp), allocatable :: lu_aa(:, :)             ! (max_ia, n_splines)
        real(wp), allocatable :: lu_bb(:, :)             ! (max_ia, n_splines)
        real(wp), allocatable :: lu_cc(:, :)             ! (max_ia, n_splines)

        ! Pre-computed CUBIC_SPLINE_BISECT result (Step D) — which xa interval each xb lies in
        integer(wi), allocatable :: spl_bisect_idx(:, :) ! (max_ib, n_splines)

        ! Writeback indices (Step D extension): for each spline, where to write yb into fld_mod
        integer(wi), allocatable :: spl_wb_count(:)              ! (n_splines)
        integer(wi), allocatable :: spl_wb_fld_idx(:, :)         ! (max_wb, n_splines)
        integer(wi), allocatable :: spl_wb_yb_idx(:, :)          ! (max_wb, n_splines)
        integer(wi) :: max_wb = 0
    end type

    type(ibm_spline_cache_t), save, target :: cache_x, cache_y, cache_z

contains

    !########################################################################
    subroutine IBM_Spline_Cache_Build_All()
        use IBM_VARS, only: isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, ibm_case_x
        use IBM_VARS, only: isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, ibm_case_y
        use IBM_VARS, only: isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, ibm_case_z
        use FDM, only: g

        call IBM_Spline_Cache_Build(cache_x, g(1), isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, ibm_case_x)
        call IBM_Spline_Cache_Build(cache_y, g(2), isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, ibm_case_y)
        call IBM_Spline_Cache_Build(cache_z, g(3), isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, ibm_case_z)
    end subroutine IBM_Spline_Cache_Build_All

    !########################################################################
    subroutine IBM_Spline_Cache_Build(cache, g, isize_nob, isize_nob_be, nob, nob_b, nob_e, ibm_case)
        use Cubic_Splines, only: CUBIC_SPLINE_LHS, CUBIC_SPLINE_BISECT

        type(ibm_spline_cache_t), intent(inout) :: cache
        type(fdm_dt), intent(in) :: g
        integer(wi), intent(in) :: isize_nob, isize_nob_be
        integer(wi), dimension(isize_nob), intent(in) :: nob
        integer(wi), dimension(isize_nob_be), intent(in) :: nob_b, nob_e, ibm_case

        integer(wi) :: ii, iob, ip, nlines, spl, k
        integer(wi) :: ia, ib, case_id
        integer(wi) :: ip_il, ip_ir, max_wb_local
        integer(wi) :: total_active, total_splines
        logical :: keep_iter
        real(wp), allocatable :: dx_tmp(:)
        real(wp), allocatable :: aa_tmp(:), bb_tmp(:), cc_tmp(:)
        integer(wi), dimension(2) :: bc

        nlines = isize_nob
        bc(:) = 2  ! fixed-first-derivative BC, matches what IBM_SPLINE_XYZ asks for

        ! First pass: count active lines and total splines, determine max_ia / max_ib / max_wb.
        total_active = 0
        total_splines = 0
        cache%max_ia = 0
        cache%max_ib = 0
        max_wb_local = 0
        do ii = 1, nlines
            if (nob(ii) == 0) cycle
            total_active = total_active + 1
            ip = 0
            keep_iter = .true.
            do iob = 1, nob(ii)
                if (.not. keep_iter) exit
                total_splines = total_splines + 1
                call cache_compute_sizes(ibm_case(ip + ii), nob_b(ip + ii), nob_e(ip + ii), g, &
                                        ia, ib, max_wb_local, cache%max_ia, cache%max_ib)
                if (ibm_case(ip + ii) == 1) keep_iter = .false.
                ip = ip + nlines
            end do
        end do

        cache%nactive = total_active
        cache%n_splines = total_splines
        cache%max_wb = max_wb_local

        ! Allocate
        allocate (cache%spl_ii(total_splines))
        allocate (cache%spl_case(total_splines))
        allocate (cache%spl_ia(total_splines))
        allocate (cache%spl_ib(total_splines))
        allocate (cache%spl_nob_b(total_splines))
        allocate (cache%spl_nob_e(total_splines))
        allocate (cache%spl_do_writeback(total_splines))
        allocate (cache%spl_xa(cache%max_ia, total_splines))
        allocate (cache%spl_xb(cache%max_ib, total_splines))
        allocate (cache%spl_dx(max(cache%max_ia - 1, 1), total_splines))
        allocate (cache%ya_src_idx_v(cache%max_ia, total_splines))
        allocate (cache%ya_src_idx_s(cache%max_ia, total_splines))
        allocate (cache%lu_aa(cache%max_ia, total_splines))
        allocate (cache%lu_bb(cache%max_ia, total_splines))
        allocate (cache%lu_cc(cache%max_ia, total_splines))
        allocate (cache%spl_bisect_idx(cache%max_ib, total_splines))
        allocate (cache%spl_wb_count(total_splines))
        allocate (cache%spl_wb_fld_idx(max(max_wb_local, 1), total_splines))
        allocate (cache%spl_wb_yb_idx(max(max_wb_local, 1), total_splines))

        allocate (dx_tmp(max(cache%max_ia - 1, 1)))
        allocate (aa_tmp(cache%max_ia), bb_tmp(cache%max_ia), cc_tmp(cache%max_ia))

        ! Second pass: populate cache
        spl = 0
        do ii = 1, nlines
            if (nob(ii) == 0) cycle
            ip = 0
            keep_iter = .true.
            do iob = 1, nob(ii)
                if (.not. keep_iter) exit
                spl = spl + 1
                case_id = ibm_case(ip + ii)
                ip_il = nob_b(ip + ii)
                ip_ir = nob_e(ip + ii)

                cache%spl_ii(spl) = ii
                cache%spl_case(spl) = case_id
                cache%spl_nob_b(spl) = ip_il
                cache%spl_nob_e(spl) = ip_ir
                cache%spl_do_writeback(spl) = (case_id /= 1)

                ! Build the per-case geometry, ya source map, and writeback indices.
                call cache_build_one(case_id, ip_il, ip_ir, ii, nlines, g, &
                                     cache%max_ia, cache%max_ib, max_wb_local, &
                                     cache%spl_xa(:, spl), cache%spl_xb(:, spl), &
                                     cache%ya_src_idx_v(:, spl), cache%ya_src_idx_s(:, spl), &
                                     cache%spl_wb_fld_idx(:, spl), cache%spl_wb_yb_idx(:, spl), &
                                     ia, ib, cache%spl_wb_count(spl))

                cache%spl_ia(spl) = ia
                cache%spl_ib(spl) = ib

                ! Step D: dx, LU(aa, bb, cc), bisect_idx
                do k = 1, ia - 1
                    dx_tmp(k) = cache%spl_xa(k + 1, spl) - cache%spl_xa(k, spl)
                end do
                cache%spl_dx(1:ia - 1, spl) = dx_tmp(1:ia - 1)

                call CUBIC_SPLINE_LHS(bc, ia, dx_tmp(1:ia - 1), aa_tmp(1:ia), bb_tmp(1:ia), cc_tmp(1:ia))
                ! Factorize in place via TRIDFS (fixed-first-deriv BC implies non-periodic path)
                call TRIDFS(ia, aa_tmp(1:ia), bb_tmp(1:ia), cc_tmp(1:ia))
                cache%lu_aa(1:ia, spl) = aa_tmp(1:ia)
                cache%lu_bb(1:ia, spl) = bb_tmp(1:ia)
                cache%lu_cc(1:ia, spl) = cc_tmp(1:ia)

                ! Bisect each xb(1:ib) into the xa(1:ia) intervals
                do k = 1, ib
                    call CUBIC_SPLINE_BISECT(ia, cache%spl_xa(1:ia, spl), cache%spl_xb(k, spl), &
                                             cache%spl_bisect_idx(k, spl))
                end do

                if (case_id == 1) keep_iter = .false.
                ip = ip + nlines
            end do
        end do

        deallocate (dx_tmp, aa_tmp, bb_tmp, cc_tmp)
        cache%built = .true.
    end subroutine IBM_Spline_Cache_Build

    !########################################################################
    ! First-pass sizing for one (case, ip_il, ip_ir).
    subroutine cache_compute_sizes(case_id, ip_il, ip_ir, g, ia, ib, max_wb, gmax_ia, gmax_ib)
        integer(wi), intent(in) :: case_id, ip_il, ip_ir
        type(fdm_dt), intent(in) :: g
        integer(wi), intent(out) :: ia, ib
        integer(wi), intent(inout) :: max_wb, gmax_ia, gmax_ib
        integer(wi) :: wb_count

        select case (case_id)
        case (1)
            ia = 2;             ib = ip_ir - ip_il + 1;        wb_count = 0
        case (2, 3, 4)
            ia = 2*nflu + 2;    ib = ip_ir - ip_il + 1;        wb_count = ib
        case (5, 6)
            ia = 2*nflu + 2;    ib = ip_ir - ip_il + 1;        wb_count = ib
        case (7)
            ia = 2*nflu + 2;    ib = (g%size - ip_il + 1) + ip_ir; wb_count = ib
        case (8)
            ia = 2*nflu + 2;    ib = 3;                        wb_count = ip_ir
        case (9)
            ia = (nflu + 1) + 2;ib = 3;                        wb_count = ip_ir
        case default
            ia = 0; ib = 0; wb_count = 0
        end select

        gmax_ia = max(gmax_ia, ia)
        gmax_ib = max(gmax_ib, ib)
        max_wb = max(max_wb, wb_count)
    end subroutine cache_compute_sizes

    !########################################################################
    ! Build xa, xb, ya source maps, and writeback indices for a single spline.
    subroutine cache_build_one(case_id, ip_il, ip_ir, ii, nlines, g, &
                                max_ia_buf, max_ib_buf, max_wb_buf, &
                                xa, xb, ya_v, ya_s, wb_fld, wb_yb, &
                                ia_out, ib_out, wb_count)
        integer(wi), intent(in) :: case_id, ip_il, ip_ir, ii, nlines
        type(fdm_dt), intent(in) :: g
        integer(wi), intent(in) :: max_ia_buf, max_ib_buf, max_wb_buf
        real(wp), intent(out) :: xa(max_ia_buf), xb(max_ib_buf)
        integer(wi), intent(out) :: ya_v(max_ia_buf), ya_s(max_ia_buf)
        integer(wi), intent(out) :: wb_fld(max_wb_buf), wb_yb(max_wb_buf)
        integer(wi), intent(out) :: ia_out, ib_out, wb_count

        integer(wi) :: ip_fl, iu_fl, iu_ir, ip_sol
        integer(wi) :: ia, ib, kflu, gap, l, n

        ia = 0; ib = 0; wb_count = 0
        xa(:) = 0.0_wp; xb(:) = 0.0_wp
        ya_v(:) = YA_SRC_SOLID; ya_s(:) = YA_SRC_SOLID

        select case (case_id)

        case (1) ! flagged "no spline" by IBM_SPLINE_XYZ; build the legacy dead work so debug sums match
            ia = ia + 1; xa(ia) = g%nodes(ip_il); ya_v(ia) = YA_SRC_SOLID; ya_s(ia) = YA_SRC_SOLID
            ia = ia + 1; xa(ia) = g%nodes(ip_ir); ya_v(ia) = YA_SRC_SOLID; ya_s(ia) = YA_SRC_SOLID
            do gap = ip_il, ip_ir
                ib = ib + 1; xb(ib) = g%nodes(gap)
            end do

        case (2) ! semi-immersed + periodic (left boundary)
            ip_fl = g%size - nflu
            iu_fl = ip_fl*nlines + ii
            iu_ir = (ip_ir - 1)*nlines + ii
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = -(g%scale - g%nodes(g%size - nflu + kflu))
                ya_v(ia) = iu_fl + (kflu - 1)*nlines
                ya_s(ia) = iu_fl + (kflu - 1)*nlines
            end do
            ia = ia + 1; xa(ia) = g%nodes(ip_il)
            ia = ia + 1; xa(ia) = g%nodes(ip_ir)
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(ip_ir + kflu)
                ya_v(ia) = iu_ir + kflu*nlines
                ya_s(ia) = iu_ir + kflu*nlines
            end do
            do gap = ip_il, ip_ir
                ib = ib + 1; xb(ib) = g%nodes(gap)
            end do
            call wb_default(ib, ip_il, nlines, ii, wb_fld, wb_yb, wb_count)

        case (3) ! semi-immersed + non-periodic (ground)
            ip_fl = ip_il - nflu
            iu_fl = (ip_fl - 1)*nlines + ii
            iu_ir = (ip_ir - 1)*nlines + ii
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = -g%nodes((nflu + 2) - kflu)
            end do
            ia = ia + 1; xa(ia) = g%nodes(ip_il)
            ia = ia + 1; xa(ia) = g%nodes(ip_ir)
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(ip_ir + kflu)
                ya_v(ia) = iu_ir + kflu*nlines
                ya_s(ia) = iu_ir + kflu*nlines
            end do
            do gap = ip_il, ip_ir
                ib = ib + 1; xb(ib) = g%nodes(gap)
            end do
            call wb_default(ib, ip_il, nlines, ii, wb_fld, wb_yb, wb_count)

        case (4) ! fully-immersed interior
            ip_fl = ip_il - nflu
            iu_fl = (ip_fl - 1)*nlines + ii
            iu_ir = (ip_ir - 1)*nlines + ii
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(ip_fl + (kflu - 1))
                ya_v(ia) = iu_fl + (kflu - 1)*nlines
                ya_s(ia) = iu_fl + (kflu - 1)*nlines
            end do
            ia = ia + 1; xa(ia) = g%nodes(ip_il)
            ia = ia + 1; xa(ia) = g%nodes(ip_ir)
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(ip_ir + kflu)
                ya_v(ia) = iu_ir + kflu*nlines
                ya_s(ia) = iu_ir + kflu*nlines
            end do
            do gap = ip_il, ip_ir
                ib = ib + 1; xb(ib) = g%nodes(gap)
            end do
            call wb_default(ib, ip_il, nlines, ii, wb_fld, wb_yb, wb_count)

        case (5) ! semi-immersed + periodic (right boundary)
            ip_fl = ip_il - nflu
            iu_fl = (ip_fl - 1)*nlines + ii
            iu_ir = ii
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(ip_fl + (kflu - 1))
                ya_v(ia) = iu_fl + (kflu - 1)*nlines
                ya_s(ia) = iu_fl + (kflu - 1)*nlines
            end do
            ia = ia + 1; xa(ia) = g%nodes(ip_il)
            ia = ia + 1; xa(ia) = g%nodes(ip_ir)
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(g%size) + g%nodes(kflu + 1)
                ya_v(ia) = iu_ir + (kflu - 1)*nlines
                ya_s(ia) = iu_ir + (kflu - 1)*nlines
            end do
            do gap = ip_il, ip_ir
                ib = ib + 1; xb(ib) = g%nodes(gap)
            end do
            call wb_default(ib, ip_il, nlines, ii, wb_fld, wb_yb, wb_count)

        case (6) ! semi-immersed + non-periodic (top)
            ip_fl = ip_il - nflu
            iu_fl = (ip_fl - 1)*nlines + ii
            iu_ir = (ip_ir - 1)*nlines + ii
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(ip_fl + (kflu - 1))
                ya_v(ia) = iu_fl + (kflu - 1)*nlines
                ya_s(ia) = iu_fl + (kflu - 1)*nlines
            end do
            ia = ia + 1; xa(ia) = g%nodes(ip_il)
            ia = ia + 1; xa(ia) = g%nodes(ip_ir)
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(g%size) + (g%nodes(g%size) - g%nodes(g%size - kflu))
            end do
            do gap = ip_il, ip_ir
                ib = ib + 1; xb(ib) = g%nodes(gap)
            end do
            call wb_default(ib, ip_il, nlines, ii, wb_fld, wb_yb, wb_count)

        case (7) ! fully-immersed wrapping periodic boundary
            ip_fl = ip_il - nflu
            iu_fl = (ip_fl - 1)*nlines + ii
            iu_ir = (ip_ir - 1)*nlines + ii
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(ip_fl + (kflu - 1))
                ya_v(ia) = iu_fl + (kflu - 1)*nlines
                ya_s(ia) = iu_fl + (kflu - 1)*nlines
            end do
            ia = ia + 1; xa(ia) = g%nodes(ip_il)
            ia = ia + 1; xa(ia) = g%nodes(ip_ir) + g%scale
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(ip_ir + kflu) + g%scale
                ya_v(ia) = iu_ir + kflu*nlines
                ya_s(ia) = iu_ir + kflu*nlines
            end do
            ip_sol = (g%size - ip_il + 1) + ip_ir
            n = 0
            do gap = 1, ip_sol
                ib = ib + 1
                if ((ip_il + gap - 1) <= g%size) then
                    xb(ib) = g%nodes(ip_il + gap - 1)
                else if ((ip_il + gap) >= g%size) then
                    xb(ib) = g%nodes(gap - ip_ir) + g%scale + (g%scale - g%nodes(g%size))
                else
                    call TLab_Write_ASCII(efile, 'IBM SPLINE_CACHE. Check case-7 gap vector.')
                    call TLab_Stop(DNS_ERROR_CUBIC_SPLINE)
                end if
            end do
            ! Writeback for case 7 (the wrap-aware loop in the legacy code)
            wb_count = 0
            n = 0
            do l = 1, ib
                if ((((ip_il - 1)*nlines + ii) + (l - 1)*nlines) <= (g%size*nlines)) then
                    n = n + 1
                    wb_count = wb_count + 1
                    wb_fld(wb_count) = ((ip_il - 1)*nlines + ii) + (l - 1)*nlines
                    wb_yb(wb_count) = l
                else
                    wb_count = wb_count + 1
                    wb_fld(wb_count) = ii + (l - n - 1)*nlines
                    wb_yb(wb_count) = l
                end if
            end do

        case (8) ! single solid point on non-periodic boundary (small stencil)
            ip_fl = ip_il - nflu
            iu_fl = (ip_fl - 1)*nlines + ii
            iu_ir = (ip_ir - 1)*nlines + ii
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = -g%nodes((nflu + 3) - kflu)
                ! For is==0 (velocity): use fld at mirrored position.
                ! For is/=0 (scalar): use ya_solid (ibmscaljmin).
                ya_v(ia) = ii + (ip_ir + nflu - kflu)*nlines
                ya_s(ia) = YA_SRC_SOLID
            end do
            ia = ia + 1
            xa(ia) = -g%nodes(nflu - 1)
            ya_v(ia) = (ip_il - 1) + ii
            ya_s(ia) = YA_SRC_SOLID
            ia = ia + 1; xa(ia) = g%nodes(ip_ir)
            do kflu = 1, nflu
                ia = ia + 1
                xa(ia) = g%nodes(ip_ir + kflu)
                ya_v(ia) = iu_ir + kflu*nlines
                ya_s(ia) = iu_ir + kflu*nlines
            end do
            xb(1) = -g%nodes(2); xb(2) = g%nodes(1); xb(3) = g%nodes(2)
            ib = 3
            do l = 1, ip_ir
                wb_count = wb_count + 1
                wb_fld(wb_count) = (l - 1)*nlines + ii
                wb_yb(wb_count) = l + 1
            end do

        case (9) ! single solid point on non-periodic boundary (larger stencil, no right half)
            ip_fl = ip_il - nflu
            iu_fl = (ip_fl - 1)*nlines + ii
            iu_ir = (ip_ir - 1)*nlines + ii
            do kflu = 1, nflu + 1
                ia = ia + 1
                xa(ia) = -g%nodes((nflu + 4) - kflu)
            end do
            ia = ia + 1; xa(ia) = -g%nodes(nflu)
            ia = ia + 1; xa(ia) = g%nodes(ip_ir)
            xb(1) = -g%nodes(3); xb(2) = -g%nodes(2); xb(3) = g%nodes(1)
            ib = 3
            do l = 1, ip_ir
                wb_count = wb_count + 1
                wb_fld(wb_count) = (l - 1)*nlines + ii
                wb_yb(wb_count) = l + 2
            end do

        end select

        ia_out = ia
        ib_out = ib
    end subroutine cache_build_one

    !########################################################################
    subroutine wb_default(ib, ip_il, nlines, ii, wb_fld, wb_yb, wb_count)
        integer(wi), intent(in) :: ib, ip_il, nlines, ii
        integer(wi), intent(out) :: wb_fld(:), wb_yb(:)
        integer(wi), intent(out) :: wb_count
        integer(wi) :: l, iu_il
        iu_il = (ip_il - 1)*nlines + ii
        wb_count = ib
        do l = 1, ib
            wb_fld(l) = iu_il + (l - 1)*nlines
            wb_yb(l) = l
        end do
    end subroutine wb_default

end module IBM_Spline_Cache
