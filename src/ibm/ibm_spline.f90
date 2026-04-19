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
    logical :: splines
    integer(wi), dimension(2) :: bc
    real(wp), dimension(2) :: bcval
    real(wp) :: m1, m2
    ! Local scratch — own storage instead of module globals.
    ! CUBIC_SPLINE expects wrk(norg,11); passing 1D wrk_loc(nspl*11) is the same convention as wrk1d.
    real(wp), allocatable :: xa_loc(:), ya_loc(:), xb_loc(:), yb_loc(:), wrk_loc(:)

    ! ================================================================== !
    ! cf. ibm_allocate.f90
    nlines = isize_nob

    fld_mod = fld     ! never modify u,v,w,s directly !

    ! index ii (dummy index; for x,y,z: ii == jk,ik,ij)
    ! Each ii is independent: writes stride by nlines so different ii never alias.
    !$omp parallel do &
    !$omp   private(xa_loc, ya_loc, xb_loc, yb_loc, wrk_loc) &
    !$omp   private(ip, iob, ia, ib, l, iu_il, iu_ir, n) &
    !$omp   private(splines, bc, bcval, m1, m2)
    do ii = 1, nlines          ! index of ii-plane, loop over plane and check for objects in each line
        if (nob(ii) /= 0) then ! if line contains immersed object(s) --yes-->  spline interpolation
            allocate (xa_loc(nspl), ya_loc(nspl))
            allocate (xb_loc(isize_wrk1d_ibm), yb_loc(isize_wrk1d_ibm))
            allocate (wrk_loc(nspl*11))
            splines = .true.
            ip = 0
            do iob = 1, nob(ii)    ! loop over immersed object(s)
                call IBM_SPLINE_VECTOR(is, ibm_case(ip + ii), fld, g, xa_loc, ya_loc, xb_loc, ia, ib, nob_b(ip + ii), nob_e(ip + ii), nlines, ii)
                if (ibm_case(ip + ii) == 1) splines = .false.
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
                end if
                ip = ip + nlines
            end do
            deallocate (xa_loc, ya_loc, xb_loc, yb_loc, wrk_loc)
        end if
    end do
    !$omp end parallel do

    return
end subroutine IBM_SPLINE_XYZ

!########################################################################

subroutine IBM_SPLINE_VECTOR(is, case, fld, g, xa, ya, xb, ia, ib, ip_il, ip_ir, nlines, plane)

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

    ! ================================================================== !
    ! indices
    ia = 0; ib = 0 ! remember current position in vectors
    select case (case)
    case (2) ! semi-immersed + periodic
        ip_fl = g%size - nflu                 ! k-axis index of most left fluid point
        iu_fl = ip_fl*nlines + plane  ! fld-index of most left fluid point
        iu_ir = (ip_ir - 1)*nlines + plane  ! fld-index of right interface  point
    case (5) ! semi-immersed + periodic
        ip_fl = ip_il - nflu
        iu_fl = (ip_fl - 1)*nlines + plane
        iu_ir = plane
    case default
        ip_fl = ip_il - nflu
        iu_fl = (ip_fl - 1)*nlines + plane
        iu_ir = (ip_ir - 1)*nlines + plane
    end select

    ! ================================================================== !
    ! vectors are built from left to right
    ! -----------------------------------------------------------------
    ! build left half of xa, ya
    select case (case)
    case (2) ! semi-immersed + periodic
        ! mirror nflu points
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = -(g%scale - g%nodes(g%size - nflu + kflu))
            ya(ia) = fld(iu_fl + (kflu - 1)*nlines)
        end do
    case (3) ! semi-immersed + non-periodic
        ! mirror nflu points on the ground + zeros at left interface
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = -g%nodes((nflu + 2) - kflu)
            if (is /= 0) then
                ya(ia) = ibmscaljmin(is)
            else
                ya(ia) = 0.0_wp
            end if
        end do
    case (4, 5, 6, 7)  ! Default case
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_fl + (kflu - 1))
            ya(ia) = fld(iu_fl + (kflu - 1)*nlines)
        end do
    case (8) ! Single solid point on non periodic boundary
        ! Construct 1 solid points and mirror 3 fluid points
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = -g%nodes((nflu + 3) - kflu)
            if (is /= 0) then
                ya(ia) = ibmscaljmin(is)
            else
                ya(ia) = 0.0_wp
                ya(ia) = fld(plane + (ip_ir + nflu - kflu)*nlines)
            end if
        end do
    case (9) ! Single solid point on non periodic boundary
        ! Construct 2 solid points and mirror 3 fluid points
        do kflu = 1, nflu + 1
            ia = ia + 1
            xa(ia) = -g%nodes((nflu + 4) - kflu)
            if (is /= 0) then
                ya(ia) = ibmscaljmin(is)
            else
                ya(ia) = 0.0_wp
            end if
        end do
    end select
    ! -----------------------------------------------------------------
    ! set interfaces (left and right)
    ia = ia + 1
    xa(ia) = g%nodes(ip_il)
    if (case == 9) xa(ia) = -g%nodes(nflu)
    if (case == 8) xa(ia) = -g%nodes(nflu - 1)
    if (is /= 0) then
        ya(ia) = ibmscaljmin(is)
    else
        ya(ia) = 0.0_wp
        if (case == 8) ya(ia) = fld((ip_il - 1) + plane)
    end if
    !
    ia = ia + 1
    xa(ia) = g%nodes(ip_ir)
    if (case == 7) xa(ia) = xa(ia) + g%scale
    if (is /= 0) then
        ya(ia) = ibmscaljmin(is)
    else
        ya(ia) = 0.0_wp
    end if
    ! -----------------------------------------------------------------
    ! build right half of xa, ya
    select case (case)
    case (2, 3, 4, 8)
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_ir + kflu)
            ya(ia) = fld(iu_ir + kflu*nlines)
        end do
    case (5) ! semi-immersed + periodic
        ! mirror nflu points
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(g%size) + g%nodes(kflu + 1)
            ya(ia) = fld(iu_ir + (kflu - 1)*nlines)
        end do
    case (6) ! semi-immersed + non-periodic
        ! mirror nflu points on the top + zeros at right interface
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(g%size) + (g%nodes(g%size) - g%nodes(g%size - kflu))
            if (is /= 0) then
                ya(ia) = ibmscaljmin(is)
            else
                ya(ia) = 0.0_wp
            end if
        end do
    case (7)
        do kflu = 1, nflu
            ia = ia + 1
            xa(ia) = g%nodes(ip_ir + kflu) + g%scale
            ya(ia) = fld(iu_ir + kflu*nlines)
        end do
    end select
    ! -----------------------------------------------------------------
    ! build gap vector where splines are evaluated (here: with interface points)
    select case (case)
    case (7)
        ip_sol = (g%size - ip_il + 1) + ip_ir
        do gap = 1, ip_sol
            ib = ib + 1
            if ((ip_il + gap - 1) <= g%size) then
                xb(ib) = g%nodes(ip_il + gap - 1)
            else if ((ip_il + gap) >= g%size) then
                xb(ib) = g%nodes(gap - ip_ir) + g%scale + (g%scale - g%nodes(g%size)) ! g%scale Warning!!
            else
                call TLab_Write_ASCII(efile, 'IBM SPLINE_VECTOR. Check gap vector.')
                call TLab_Stop(DNS_ERROR_CUBIC_SPLINE)
            end if
        end do
    case (8)
        xb(1) = -g%nodes(2)
        xb(2) = g%nodes(1)
        xb(3) = g%nodes(2)
        ib = 3
    case (9)
        xb(1) = -g%nodes(3)
        xb(2) = -g%nodes(2)
        xb(3) = g%nodes(1)
        ib = 3
    case default
        do gap = ip_il, ip_ir
            ib = ib + 1
            xb(ib) = g%nodes(gap)
        end do
    end select

    return
end subroutine IBM_SPLINE_VECTOR

!########################################################################
