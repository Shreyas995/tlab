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
    use TLab_Constants, only: efile, wp, wi
    ! xa/ya/xb/yb replaced by local allocatables — removes shared global state (Phase 1: thread-safety prep).
    use IBM_VARS, only: nspl, isize_wrk1d_ibm, ibmscaljmin
    use TLab_Memory, only: isize_field
    ! wrk1d replaced by local wrk_loc — same reason.
    use FDM, only: fdm_dt
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Cubic_Splines
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

    integer(wi) :: l, ii, ip, ia, ib, iob, iu_il, iu_ir, n, nlines
    integer(wi) :: iact, nactive
    logical :: splines
    integer(wi), dimension(2) :: bc
    real(wp), dimension(2) :: bcval
    real(wp) :: m1, m2
    ! Local scratch — own storage instead of module globals.
    ! CUBIC_SPLINE expects wrk(norg,11); passing 1D wrk_loc(nspl*11) is the same convention as wrk1d.
    real(wp), allocatable :: xa_loc(:), ya_loc(:), xb_loc(:), yb_loc(:), wrk_loc(:)
    ! Step A: skip-list of active ii indices (those with nob(ii)/=0).
    ! Built once at entry so the parallel loop only iterates over real work,
    ! avoiding stranded threads on empty lines and the per-iteration outer if.
    integer(wi), allocatable :: active_ii(:)
    ! Aggregate debug counters — accumulated across the ii loop, dumped at exit.
    ! These are loop-invariant fingerprints: if A+B refactor preserves logic, they match the reference run.
    integer(wi) :: dbg_nactive, dbg_total_ia, dbg_total_ib, dbg_total_iob
    real(wp)    :: dbg_sum_xa, dbg_sum_ya, dbg_sum_xb, dbg_sum_yb

    ! ================================================================== !
    ! cf. ibm_allocate.f90
    nlines = isize_nob

    call TLab_Debug_Print_int('[I0_IBM_in] is=', is)
    call TLab_Debug_Print_1D('[I0_IBM_in] sum(fld)=', fld)

    fld_mod = fld     ! never modify u,v,w,s directly !

    ! Build skip-list. Geometry-only: could be cached at init if profiling shows the build cost matters.
    allocate (active_ii(nlines))
    nactive = 0
    do ii = 1, nlines
        if (nob(ii) /= 0) then
            nactive = nactive + 1
            active_ii(nactive) = ii
        end if
    end do
    call TLab_Debug_Print_int('[A0_skiplist] nactive=', nactive)

    ! Reset aggregate debug counters
    dbg_nactive = 0; dbg_total_ia = 0; dbg_total_ib = 0; dbg_total_iob = 0
    dbg_sum_xa = 0.0_wp; dbg_sum_ya = 0.0_wp; dbg_sum_xb = 0.0_wp; dbg_sum_yb = 0.0_wp

    ! index ii (dummy index; for x,y,z: ii == jk,ik,ij)
    ! Each ii is independent: writes stride by nlines so different ii never alias.
    ! Scratch arrays are allocated once per thread (outside the ii loop) and reused across all iterations.
    ! Step A: loop over active_ii(1:nactive) instead of 1..nlines — removes the outer `if (nob/=0)` check
    ! and gives the scheduler dense work units, improving load balance on sparse-IBM lines.
    !$omp parallel &
    !$omp   private(xa_loc, ya_loc, xb_loc, yb_loc, wrk_loc) &
    !$omp   private(iact, ii, ip, iob, ia, ib, l, iu_il, iu_ir, n) &
    !$omp   private(splines, bc, bcval, m1, m2)
    allocate (xa_loc(nspl), ya_loc(nspl))
    allocate (xb_loc(isize_wrk1d_ibm), yb_loc(isize_wrk1d_ibm))
    allocate (wrk_loc(nspl*11))
    !$omp do reduction(+:dbg_nactive,dbg_total_ia,dbg_total_ib,dbg_total_iob,dbg_sum_xa,dbg_sum_ya,dbg_sum_xb,dbg_sum_yb)
    do iact = 1, nactive
        ii = active_ii(iact)
        dbg_nactive = dbg_nactive + 1
        splines = .true.
        ip = 0
        do iob = 1, nob(ii)    ! loop over immersed object(s)
                call IBM_SPLINE_VECTOR(is, ibm_case(ip + ii), fld, g, xa_loc, ya_loc, xb_loc, ia, ib, nob_b(ip + ii), nob_e(ip + ii), nlines, ii)
                if (ibm_case(ip + ii) == 1) splines = .false.
                dbg_total_iob = dbg_total_iob + 1
                dbg_total_ia  = dbg_total_ia  + ia
                dbg_total_ib  = dbg_total_ib  + ib
                dbg_sum_xa    = dbg_sum_xa    + sum(xa_loc(1:ia))
                dbg_sum_ya    = dbg_sum_ya    + sum(ya_loc(1:ia))
                dbg_sum_xb    = dbg_sum_xb    + sum(xb_loc(1:ib))
                ! ================================================================== !
                ! spline interpolation and fill gap in fld_ibm
                if (splines) then
                    ! generate splines (other possibility: natural boundary conditions)
                    bc(:) = 2 ! fixed first derivative at endpoints
                    m1 = (ya_loc(2) - ya_loc(1))/(xa_loc(2) - xa_loc(1)); bcval(1) = m1
                    m2 = (ya_loc(ia) - ya_loc(ia - 1))/(xa_loc(ia) - xa_loc(ia - 1)); bcval(2) = m2
!DIR$ INLINE CUBIC_SPLINE
                    call CUBIC_SPLINE(bc, bcval, ia, ib, xa_loc(1:ia), ya_loc(1:ia), xb_loc(1:ib), yb_loc(1:ib), wrk_loc)
                    ! force yb at interface to physical BCs again, to get exact boundary values here
                    if (is /= 0) then
                        yb_loc(1) = ibmscaljmin(is)
                        yb_loc(ib) = ibmscaljmin(is)
                    else
                        yb_loc(1) = 0.0_wp
                        yb_loc(ib) = 0.0_wp
                    end if
                    ! fld index of left interface
                    iu_il = (nob_b(ip + ii) - 1)*nlines + ii
                    iu_ir = (nob_e(ip + ii) - 1)*nlines + ii
                    n = 0
                    ! replace splines in solid gaps
                    if (((nob_e(ip + ii)) < (nob_b(ip + ii))) .and. (g%periodic .eqv. .true.)) then ! condition for case 7
                        do l = 1, ib
                            if ((iu_il + (l - 1)*nlines) <= (g%size*nlines)) then
                                n = n + 1
                                fld_mod(iu_il + (l - 1)*nlines) = yb_loc(l)
                            else if ((iu_il + (l - 1)*nlines) >= (g%size*nlines)) then
                                fld_mod(ii + (l - n - 1)*nlines) = yb_loc(l)
                            else
                                call TLab_Write_ASCII(efile, 'IBM SPLINE. Error in replacing spline in the solid.')
                                call TLab_Stop(DNS_ERROR_CUBIC_SPLINE)
                            end if
                        end do
                    else if (((nob_e(ip + ii)) == 1 + (nob_b(ip + ii))) .and. (g%periodic .eqv. .false.)) then ! condition for case 8
                        do l = 1, (nob_e(ip + ii))
                            fld_mod((l - 1)*nlines + ii) = yb_loc(l + 1)
                        end do
                    else if (((nob_e(ip + ii)) == (nob_b(ip + ii))) .and. (g%periodic .eqv. .false.)) then ! condition for case 9
                        do l = 1, (nob_e(ip + ii))
                            fld_mod((l - 1)*nlines + ii) = yb_loc(l + 2)
                        end do
                    else ! default execution
                        do l = 1, ib
                            fld_mod(iu_il + (l - 1)*nlines) = yb_loc(l)
                        end do
                    end if
                    dbg_sum_yb = dbg_sum_yb + sum(yb_loc(1:ib))
                end if
                ip = ip + nlines
            end do
    end do
    !$omp end do
    deallocate (xa_loc, ya_loc, xb_loc, yb_loc, wrk_loc)
    !$omp end parallel

    deallocate (active_ii)

    call TLab_Debug_Print_int('[I1_IBM_agg] nactive=', dbg_nactive)
    call TLab_Debug_Print_int('[I1_IBM_agg] total_iob=', dbg_total_iob)
    call TLab_Debug_Print_int('[I1_IBM_agg] total_ia=', dbg_total_ia)
    call TLab_Debug_Print_int('[I1_IBM_agg] total_ib=', dbg_total_ib)
    call TLab_Debug_Print_real('[I1_IBM_agg] sum_xa=', dbg_sum_xa)
    call TLab_Debug_Print_real('[I1_IBM_agg] sum_ya=', dbg_sum_ya)
    call TLab_Debug_Print_real('[I1_IBM_agg] sum_xb=', dbg_sum_xb)
    call TLab_Debug_Print_real('[I1_IBM_agg] sum_yb=', dbg_sum_yb)
    call TLab_Debug_Print_1D('[I2_IBM_out] sum(fld_mod)=', fld_mod)

    return
end subroutine IBM_SPLINE_XYZ

!########################################################################

subroutine IBM_SPLINE_VECTOR(is, case, fld, g, xa, ya, xb, ia, ib, ip_il, ip_ir, nlines, plane)

    ! Step B refactor: single per-case dispatch in place of the previous three successive select-case blocks
    ! (index setup, left/right half builds, gap build). Each case here builds its full xa/ya/xb in one
    ! contiguous, self-contained block — easier to read per-case, and removes two redundant branch dispatches.
    ! Behavior is bit-identical to the prior version (including the dead work done for case 1).

    use IBM_VARS, only: nflu, isize_wrk1d_ibm, nspl, ibmscaljmin
    use TLab_Memory, only: isize_field
    use FDM, only: fdm_dt
    use TLab_Constants, only: wp, wi, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop

    implicit none

    integer(wi), intent(in) :: is
    integer(wi), intent(in) :: case
    real(wp), dimension(isize_field), intent(in) :: fld
    type(fdm_dt), intent(in) :: g
    real(wp), dimension(nspl), intent(out) :: xa
    real(wp), dimension(nspl), intent(out) :: ya
    real(wp), dimension(isize_wrk1d_ibm), intent(out) :: xb
    integer(wi), intent(out) :: ia
    integer(wi), intent(out) :: ib
    integer(wi), intent(in) :: ip_il, ip_ir, nlines, plane

    integer(wi) :: kflu, gap, ip_sol
    integer(wi) :: ip_fl, iu_fl, iu_ir
    real(wp)    :: ya_solid  ! solid-region scalar value (ibmscaljmin or 0): same for both interfaces

    ! ================================================================== !
    ia = 0; ib = 0

    if (is /= 0) then
        ya_solid = ibmscaljmin(is)
    else
        ya_solid = 0.0_wp
    end if

    select case (case)

    case (1) ! flagged "no spline" by caller — caller sets splines=.false. and discards xa/ya/xb.
        ! Preserve the dead work the legacy 3-block dispatch produced for case 1, so debug fingerprints match.
        ip_fl = ip_il - nflu
        iu_fl = (ip_fl - 1)*nlines + plane
        iu_ir = (ip_ir - 1)*nlines + plane
        ia = ia + 1; xa(ia) = g%nodes(ip_il); ya(ia) = ya_solid     ! left interface
        ia = ia + 1; xa(ia) = g%nodes(ip_ir); ya(ia) = ya_solid     ! right interface
        do gap = ip_il, ip_ir
            ib = ib + 1; xb(ib) = g%nodes(gap)
        end do

    case (2) ! semi-immersed + periodic (left boundary)
        ip_fl = g%size - nflu
        iu_fl = ip_fl*nlines + plane
        iu_ir = (ip_ir - 1)*nlines + plane
        ! Left half: mirror nflu fluid points across left boundary via periodicity
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = -(g%scale - g%nodes(g%size - nflu + kflu))
            ya(ia) = fld(iu_fl + (kflu - 1)*nlines)
        end do
        ia = ia + 1; xa(ia) = g%nodes(ip_il); ya(ia) = ya_solid     ! left interface
        ia = ia + 1; xa(ia) = g%nodes(ip_ir); ya(ia) = ya_solid     ! right interface
        ! Right half: default fluid points
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_ir + kflu)
            ya(ia) = fld(iu_ir + kflu*nlines)
        end do
        do gap = ip_il, ip_ir
            ib = ib + 1; xb(ib) = g%nodes(gap)
        end do

    case (3) ! semi-immersed + non-periodic (ground)
        ip_fl = ip_il - nflu
        iu_fl = (ip_fl - 1)*nlines + plane
        iu_ir = (ip_ir - 1)*nlines + plane
        ! Left half: mirror nflu points on the ground + zeros at left interface
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = -g%nodes((nflu + 2) - kflu)
            ya(ia) = ya_solid
        end do
        ia = ia + 1; xa(ia) = g%nodes(ip_il); ya(ia) = ya_solid
        ia = ia + 1; xa(ia) = g%nodes(ip_ir); ya(ia) = ya_solid
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_ir + kflu)
            ya(ia) = fld(iu_ir + kflu*nlines)
        end do
        do gap = ip_il, ip_ir
            ib = ib + 1; xb(ib) = g%nodes(gap)
        end do

    case (4) ! fully-immersed interior, default fluid both sides
        ip_fl = ip_il - nflu
        iu_fl = (ip_fl - 1)*nlines + plane
        iu_ir = (ip_ir - 1)*nlines + plane
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_fl + (kflu - 1))
            ya(ia) = fld(iu_fl + (kflu - 1)*nlines)
        end do
        ia = ia + 1; xa(ia) = g%nodes(ip_il); ya(ia) = ya_solid
        ia = ia + 1; xa(ia) = g%nodes(ip_ir); ya(ia) = ya_solid
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_ir + kflu)
            ya(ia) = fld(iu_ir + kflu*nlines)
        end do
        do gap = ip_il, ip_ir
            ib = ib + 1; xb(ib) = g%nodes(gap)
        end do

    case (5) ! semi-immersed + periodic (right boundary)
        ip_fl = ip_il - nflu
        iu_fl = (ip_fl - 1)*nlines + plane
        iu_ir = plane
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_fl + (kflu - 1))
            ya(ia) = fld(iu_fl + (kflu - 1)*nlines)
        end do
        ia = ia + 1; xa(ia) = g%nodes(ip_il); ya(ia) = ya_solid
        ia = ia + 1; xa(ia) = g%nodes(ip_ir); ya(ia) = ya_solid
        ! Right half: mirror nflu points via periodicity
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(g%size) + g%nodes(kflu + 1)
            ya(ia) = fld(iu_ir + (kflu - 1)*nlines)
        end do
        do gap = ip_il, ip_ir
            ib = ib + 1; xb(ib) = g%nodes(gap)
        end do

    case (6) ! semi-immersed + non-periodic (top)
        ip_fl = ip_il - nflu
        iu_fl = (ip_fl - 1)*nlines + plane
        iu_ir = (ip_ir - 1)*nlines + plane
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_fl + (kflu - 1))
            ya(ia) = fld(iu_fl + (kflu - 1)*nlines)
        end do
        ia = ia + 1; xa(ia) = g%nodes(ip_il); ya(ia) = ya_solid
        ia = ia + 1; xa(ia) = g%nodes(ip_ir); ya(ia) = ya_solid
        ! Right half: mirror nflu solid points on top
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(g%size) + (g%nodes(g%size) - g%nodes(g%size - kflu))
            ya(ia) = ya_solid
        end do
        do gap = ip_il, ip_ir
            ib = ib + 1; xb(ib) = g%nodes(gap)
        end do

    case (7) ! fully-immersed object wrapping the periodic boundary
        ip_fl = ip_il - nflu
        iu_fl = (ip_fl - 1)*nlines + plane
        iu_ir = (ip_ir - 1)*nlines + plane
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_fl + (kflu - 1))
            ya(ia) = fld(iu_fl + (kflu - 1)*nlines)
        end do
        ia = ia + 1; xa(ia) = g%nodes(ip_il);             ya(ia) = ya_solid
        ia = ia + 1; xa(ia) = g%nodes(ip_ir) + g%scale;   ya(ia) = ya_solid    ! shifted across period
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_ir + kflu) + g%scale
            ya(ia) = fld(iu_ir + kflu*nlines)
        end do
        ! Gap wraps across the periodic boundary
        ip_sol = (g%size - ip_il + 1) + ip_ir
        do gap = 1, ip_sol
            ib = ib + 1
            if ((ip_il + gap - 1) <= g%size) then
                xb(ib) = g%nodes(ip_il + gap - 1)
            else if ((ip_il + gap) >= g%size) then
                xb(ib) = g%nodes(gap - ip_ir) + g%scale + (g%scale - g%nodes(g%size))
            else
                call TLab_Write_ASCII(efile, 'IBM SPLINE_VECTOR. Check gap vector.')
                call TLab_Stop(DNS_ERROR_CUBIC_SPLINE)
            end if
        end do

    case (8) ! single solid point on non-periodic boundary (small stencil)
        ip_fl = ip_il - nflu
        iu_fl = (ip_fl - 1)*nlines + plane
        iu_ir = (ip_ir - 1)*nlines + plane
        ! Left half: mirror nflu points; scalar uses fld mirrored across the boundary
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = -g%nodes((nflu + 3) - kflu)
            if (is /= 0) then
                ya(ia) = ibmscaljmin(is)
            else
                ya(ia) = fld(plane + (ip_ir + nflu - kflu)*nlines)
            end if
        end do
        ! Left interface (special: xa shifted, ya uses fld for is==0)
        ia = ia + 1
        xa(ia) = -g%nodes(nflu - 1)
        if (is /= 0) then
            ya(ia) = ibmscaljmin(is)
        else
            ya(ia) = fld((ip_il - 1) + plane)
        end if
        ! Right interface: default
        ia = ia + 1; xa(ia) = g%nodes(ip_ir); ya(ia) = ya_solid
        ! Right half: default
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_ir + kflu)
            ya(ia) = fld(iu_ir + kflu*nlines)
        end do
        ! Gap: explicit 3-point stencil
        xb(1) = -g%nodes(2); xb(2) = g%nodes(1); xb(3) = g%nodes(2)
        ib = 3

    case (9) ! single solid point on non-periodic boundary (larger stencil, no right half)
        ip_fl = ip_il - nflu
        iu_fl = (ip_fl - 1)*nlines + plane
        iu_ir = (ip_ir - 1)*nlines + plane
        ! Left half: nflu+1 mirrored solid points
        do kflu = 1, nflu + 1
            ia = ia + 1
            xa(ia) = -g%nodes((nflu + 4) - kflu)
            ya(ia) = ya_solid
        end do
        ! Left interface (special xa)
        ia = ia + 1; xa(ia) = -g%nodes(nflu); ya(ia) = ya_solid
        ! Right interface: default
        ia = ia + 1; xa(ia) = g%nodes(ip_ir); ya(ia) = ya_solid
        ! No right half for case 9
        ! Gap: explicit 3-point stencil
        xb(1) = -g%nodes(3); xb(2) = -g%nodes(2); xb(3) = g%nodes(1)
        ib = 3

    end select

    return
end subroutine IBM_SPLINE_VECTOR

!########################################################################
