#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!#
!# Assumes statistical homogeneity in xOz, so that the corresponding
!# partial derivative terms are assumed to be zero.
!#
!# In the incompressible case, the array p has been
!# pointed to dudz and the pressure field is stored there; do not
!# use array tmp3 until pressure block
!#
!# Reynolds and Favre averages
!#
!########################################################################

subroutine AVG_SCAL_XZ(is, q, s, s_local, dsdx, dsdy, dsdz, tmp1, tmp2, tmp3, mean2d)
    use TLab_Constants, only: MAX_AVG_TEMPORAL
    use TLab_Constants, only: efile, lfile, wp, wi
    use TLab_Memory, only: imax, jmax, kmax, inb_flow_array, inb_scal_array, inb_scal
    use TLab_Time, only: itime, rtime
    use TLab_Arrays, only: wrk1d
    use TLab_Pointers_3D, only: p_wrk3d, u, v, w, rho, vis
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use FDM, only: g
    use FDM, only: fdm_Int0
    use Thermodynamics, only: thermo_param, itransport, EQNS_TRANS_POWERLAW, EQNS_TRANS_SUTHERLAND
    use Thermodynamics, only: imixture, MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR
    use THERMO_ANELASTIC, only: Thermo_Anelastic_WEIGHT_INPLACE, Thermo_Anelastic_BUOYANCY, ribackground
    use THERMO_AIRWATER
    use NavierStokes
    use Gravity, only: buoyancy, bbackground, Gravity_Buoyancy, Gravity_Buoyancy_Source, EQNS_BOD_EXPLICIT
    use Rotation, only: coriolis
    use Radiation
    use Microphysics
    use FI_GRADIENT_EQN
    use OPR_Partial
    use IBM_VARS, only: imode_ibm, gamma_0, gamma_1, scal_bcs
    use Averages, only: AVG_IK_V

    implicit none

    integer, intent(IN) :: is
    real(wp), intent(IN) :: q(imax, jmax, kmax, inb_flow_array)
    real(wp), intent(IN) :: s(imax, jmax, kmax, inb_scal_array)
    real(wp), intent(IN) :: s_local(imax, jmax, kmax)
    real(wp), dimension(imax, jmax, kmax), intent(INOUT) :: dsdx, dsdy, dsdz, tmp1, tmp2, tmp3
    real(wp), intent(INOUT) :: mean2d(jmax, MAX_AVG_TEMPORAL)

    target q, tmp3

    ! -----------------------------------------------------------------------
    integer, parameter :: MAX_VARS_GROUPS = 10
    integer i ,j, k, bcs(2, 2), is_loc
    real(wp) diff, dummy, coefT, coefR, coefQ, c23

    integer ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nv, im

    character*32 name, groupname(MAX_VARS_GROUPS)
    character*250 line1, varname(MAX_VARS_GROUPS)

    ! Pointers to existing allocated space
    real(wp), dimension(:, :, :), pointer :: p_loc

    ! ###################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    ! Define pointers
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) then
        p_loc => q(:, :, :, 6)
    else
        p_loc => tmp3
    end if

    if (nse_diffusion == EQNS_NONE) then; diff = 0.0_wp
    else; diff = visc/schmidt(is)
    end if

    c23 = 2.0_wp/3.0_wp

    ! -----------------------------------------------------------------------
    ! Dependent variables
    ng = 1; ig(ng) = 1
#define rS(j)     mean2d(j,ig(1)  )
#define fS(j)     mean2d(j,ig(1)+1)
#define rS_y(j)   mean2d(j,ig(1)+2)
#define fS_y(j)   mean2d(j,ig(1)+3)
#define rQ(j)     mean2d(j,ig(1)+4)
#define fQ(j)     mean2d(j,ig(1)+5)
    sg(ng) = 6

    groupname(ng) = 'Mean'
    varname(ng) = 'rS fS rS_y fS_y rQ fQ'
    if (imode_ibm == 1) then
        varname(ng) = trim(adjustl(varname(ng)))//' eps_0 eps_1 Sbcs'
#define ep_0(j)   mean2d(j,ig(1)+6)
#define ep_1(j)   mean2d(j,ig(1)+7)
#define Sbcs(j)   mean2d(j,ig(1)+8)
        sg(ng) = sg(ng) + 3
    end if
    if (infraredProps%active(is)) then
        if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            varname(ng) = trim(adjustl(varname(ng)))//' rQrad rQradC'
        else
            varname(ng) = trim(adjustl(varname(ng)))//' rQrad rFrad'
        end if
        sg(ng) = sg(ng) + 2
    end if
    if (imixture == MIXT_TYPE_AIRWATER_LINEAR .or. imixture == MIXT_TYPE_AIRWATER) then
        varname(ng) = trim(adjustl(varname(ng)))//' rQeva'
        sg(ng) = sg(ng) + 1
    end if
    if (sedimentationProps%active(is)) then
        if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            varname(ng) = trim(adjustl(varname(ng)))//' rQtra rQtraC'
        else
            varname(ng) = trim(adjustl(varname(ng)))//' rQtra rFtra'
        end if
        sg(ng) = sg(ng) + 2
    end if

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsu(j)    mean2d(j,ig(2)  )
#define Rsv(j)    mean2d(j,ig(2)+1)
#define Rsw(j)    mean2d(j,ig(2)+2)
#define fS2(j)    mean2d(j,ig(2)+3)
#define fS3(j)    mean2d(j,ig(2)+4)
#define fS4(j)    mean2d(j,ig(2)+5)
#define rS2(j)    mean2d(j,ig(2)+6)
#define rS3(j)    mean2d(j,ig(2)+7)
#define rS4(j)    mean2d(j,ig(2)+8)
    sg(ng) = 9

    groupname(ng) = 'Fluctuations'
    varname(ng) = 'Rsu Rsv Rsw fS2 fS3 fS4 rS2 rS3 rS4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rss_t(j)  mean2d(j,ig(3)  )
#define Css(j)    mean2d(j,ig(3)+1)
#define Pss(j)    mean2d(j,ig(3)+2)
#define Ess(j)    mean2d(j,ig(3)+3)
#define Tssy1(j)  mean2d(j,ig(3)+4)
#define Tssy2(j)  mean2d(j,ig(3)+5)
#define Tssy_y(j) mean2d(j,ig(3)+6)
#define Dss(j)    mean2d(j,ig(3)+7)
#define Qss(j)    mean2d(j,ig(3)+8)
    sg(ng) = 9

    groupname(ng) = 'RssBudget'
    varname(ng) = 'Rss_t Css Pss Ess Tssy1 Tssy2 Tssy_y Dss Qss'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsu_t(j)  mean2d(j,ig(4)  )
#define Csu(j)    mean2d(j,ig(4)+1)
#define Psu(j)    mean2d(j,ig(4)+2)
#define Esu(j)    mean2d(j,ig(4)+3)
#define PIsu(j)   mean2d(j,ig(4)+4)
#define Tsuy1(j)  mean2d(j,ig(4)+5)
#define Tsuy2(j)  mean2d(j,ig(4)+6)
#define Tsuy_y(j) mean2d(j,ig(4)+7)
#define Dsu(j)    mean2d(j,ig(4)+8)
#define Gsu(j)    mean2d(j,ig(4)+9)
#define Bsu(j)    mean2d(j,ig(4)+10)
#define Fsu(j)    mean2d(j,ig(4)+11)
#define Qsu(j)    mean2d(j,ig(4)+12)
    sg(ng) = 13

    groupname(ng) = 'RsuBudget'
    varname(ng) = 'Rsu_t Csu Psu Esu PIsu Tsuy1 Tsuy2 Tsuy_y Dsu Gsu Bsu Fsu Qsu'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsv_t(j)  mean2d(j,ig(5)  )
#define Csv(j)    mean2d(j,ig(5)+1)
#define Psv(j)    mean2d(j,ig(5)+2)
#define Esv(j)    mean2d(j,ig(5)+3)
#define PIsv(j)   mean2d(j,ig(5)+4)
#define Tsvy1(j)  mean2d(j,ig(5)+5)
#define Tsvy2(j)  mean2d(j,ig(5)+6)
#define Tsvy3(j)  mean2d(j,ig(5)+7)
#define Tsvy_y(j) mean2d(j,ig(5)+8)
#define Dsv(j)    mean2d(j,ig(5)+9)
#define Gsv(j)    mean2d(j,ig(5)+10)
#define Bsv(j)    mean2d(j,ig(5)+11)
#define Fsv(j)    mean2d(j,ig(5)+12)
#define Qsv(j)    mean2d(j,ig(5)+13)
    sg(ng) = 14

    groupname(ng) = 'RsvBudget'
    varname(ng) = 'Rsv_t Csv Psv Esv PIsv Tsvy1 Tsvy2 Tsvy3 Tsvy_y Dsv Gsv Bsv Fsv Qsv'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rsw_t(j)  mean2d(j,ig(6)  )
#define Csw(j)    mean2d(j,ig(6)+1)
#define Psw(j)    mean2d(j,ig(6)+2)
#define Esw(j)    mean2d(j,ig(6)+3)
#define PIsw(j)   mean2d(j,ig(6)+4)
#define Tswy1(j)  mean2d(j,ig(6)+5)
#define Tswy2(j)  mean2d(j,ig(6)+6)
#define Tswy_y(j) mean2d(j,ig(6)+7)
#define Dsw(j)    mean2d(j,ig(6)+8)
#define Gsw(j)    mean2d(j,ig(6)+9)
#define Bsw(j)    mean2d(j,ig(6)+10)
#define Fsw(j)    mean2d(j,ig(6)+11)
#define Qsw(j)    mean2d(j,ig(6)+12)
    sg(ng) = 13

    groupname(ng) = 'RswBudget'
    varname(ng) = 'Rsw_t Csw Psw Esw PIsw Tswy1 Tswy2 Tswy_y Dsw Gsw Bsw Fsw Qsw'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define S_x2(j)  mean2d(j,ig(7)  )
#define S_y2(j)  mean2d(j,ig(7)+1)
#define S_z2(j)  mean2d(j,ig(7)+2)
#define S_x3(j)  mean2d(j,ig(7)+3)
#define S_y3(j)  mean2d(j,ig(7)+4)
#define S_z3(j)  mean2d(j,ig(7)+5)
#define S_x4(j)  mean2d(j,ig(7)+6)
#define S_y4(j)  mean2d(j,ig(7)+7)
#define S_z4(j)  mean2d(j,ig(7)+8)
    sg(ng) = 9

    groupname(ng) = 'DerivativeFluctuations'
    varname(ng) = 'S_x2 S_y2 S_z2 S_x3 S_y3 S_z3 S_x4 S_y4 S_z4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
    sg(ng) = 2*inb_scal_array

    groupname(ng) = 'CrossScalars'
    varname(ng) = ''
    do is_loc = 1, inb_scal_array
        write (name, *) is_loc
        varname(ng) = trim(adjustl(varname(ng)))//' Cs'//trim(adjustl(name))
        varname(ng) = trim(adjustl(varname(ng)))//' Css'//trim(adjustl(name))
    end do

    ! -----------------------------------------------------------------------
    ! Auxiliary variables depending on y and t; this last group is not written
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rR(j)       mean2d(j,ig(9))
#define rU(j)       mean2d(j,ig(9)+1)
#define rV(j)       mean2d(j,ig(9)+2)
#define rW(j)       mean2d(j,ig(9)+3)
#define fU(j)       mean2d(j,ig(9)+4)
#define fV(j)       mean2d(j,ig(9)+5)
#define fW(j)       mean2d(j,ig(9)+6)
#define fU_y(j)     mean2d(j,ig(9)+7)
#define fV_y(j)     mean2d(j,ig(9)+8)
#define fW_y(j)     mean2d(j,ig(9)+9)
#define rU_y(j)     mean2d(j,ig(9)+10)
#define rV_y(j)     mean2d(j,ig(9)+11)
#define rW_y(j)     mean2d(j,ig(9)+12)
#define Rvu(j)      mean2d(j,ig(9)+13)
#define Rvv(j)      mean2d(j,ig(9)+14)
#define Rvw(j)      mean2d(j,ig(9)+15)
#define Rss_y(j)    mean2d(j,ig(9)+16)
#define Rsu_y(j)    mean2d(j,ig(9)+17)
#define Rsv_y(j)    mean2d(j,ig(9)+18)
#define Rsw_y(j)    mean2d(j,ig(9)+19)
#define Fy(j)       mean2d(j,ig(9)+20)
#define Fy_y(j)     mean2d(j,ig(9)+21)
#define Tau_yy(j)   mean2d(j,ig(9)+22)
#define Tau_yy_y(j) mean2d(j,ig(9)+23)
#define Tau_yx(j)   mean2d(j,ig(9)+24)
#define Tau_yx_y(j) mean2d(j,ig(9)+25)
#define Tau_yz(j)   mean2d(j,ig(9)+26)
#define Tau_yz_y(j) mean2d(j,ig(9)+27)
#define rP(j)       mean2d(j,ig(9)+28)
#define aux(j)      mean2d(j,ig(9)+29)
    sg(ng) = 30

    ! -----------------------------------------------------------------------
    nv = ig(ng) + sg(ng) - 1
    if (MAX_AVG_TEMPORAL < nv) then
        call TLab_Write_ASCII(efile, 'AVG_SCAL_XZ. Not enough space in local arrays.')
        call TLab_Stop(DNS_ERROR_AVGTMP)
    end if
#ifndef USE_APU
    !$omp target teams distribute parallel do default(shared) private (i,j)
    do i = 1, nv
        do j = 1, jmax
            mean2d(j, i) = 0.0_wp
        end do
    end do
    !$omp end target teams distribute parallel do
    
    print *,'1', 'mean2d: ', sum(mean2d)
#else
    mean2d(:, 1:nv) = 0.0_wp
    print *,'1', 'mean2d: ', sum(mean2d)
#endif

    ng = ng - 1
    nv = ig(ng) + sg(ng) - 1 ! the last group is not written out

    ! #######################################################################
    write (line1, *) itime; line1 = 'Calculating scal statistics at It'//trim(adjustl(line1))//'...'
    call TLab_Write_ASCII(lfile, line1)

    ! #######################################################################
    ! Preliminary for IBM usage
    ! #######################################################################
    ! Asign gammas for conditional averages (c.f. Pope, p.170 [5.305])
    ! write out scalar boundary values applied in solids
    im = 5
    if (imode_ibm == 1) then
        ep_0(:) = gamma_0; ep_1(:) = gamma_1
        Sbcs(:) = scal_bcs(:, is)
        im = im + 3
    end if

    ! #######################################################################
    ! Preliminary data of velocity and density
    ! #######################################################################
    call AVG_IK_V(imax, jmax, kmax, u, rU(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, v, rV(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, w, rW(1), wrk1d)

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        rR(:) = 1.0_wp

        fU(:) = rU(:)
        fV(:) = rV(:)
        fW(:) = rW(:)

    else
        call AVG_IK_V(imax, jmax, kmax, rho, rR(1), wrk1d)

        dsdx = rho*u
        dsdy = rho*v
        dsdz = rho*w
        call AVG_IK_V(imax, jmax, kmax, dsdx, fU(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dsdy, fV(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dsdz, fW(1), wrk1d)
        fU(:) = fU(:)/rR(:)
        fV(:) = fV(:)/rR(:)
        fW(:) = fW(:)/rR(:)

    end if

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rU(1), rU_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rV(1), rV_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rW(1), rW_y(1))

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fU(1), fU_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fV(1), fV_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fW(1), fW_y(1))

    dsdx = v*u
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) dsdx = dsdx*rho
    dsdy = v*v
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) dsdy = dsdy*rho
    dsdz = v*w
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) dsdz = dsdz*rho
    call AVG_IK_V(imax, jmax, kmax, dsdx, Rvu(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdy, Rvv(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdz, Rvw(1), wrk1d)
    Rvu(:) = Rvu(:)/rR(:) - fV(:)*fU(:)
    Rvv(:) = Rvv(:)/rR(:) - fV(:)*fV(:)
    Rvw(:) = Rvw(:)/rR(:) - fV(:)*fW(:)

    ! #######################################################################
    ! Scalar
    ! #######################################################################
    call AVG_IK_V(imax, jmax, kmax, s_local, rS(1), wrk1d)

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        fS(:) = rS(:)
    else
        p_wrk3d = rho*s_local
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, fS(1), wrk1d)
        fS(:) = fS(:)/rR(:)
    end if

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rS(1), rS_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fS(1), fS_y(1))

    ! -----------------------------------------------------------------------
    ! Moments
#ifndef USE_APU
    !$omp target teams distribute parallel do collapse(3) default(shared) private(i,j,k)
    do j = 1, jmax !offload
        do i = 1, imax
            do k = 1, kmax
                p_wrk3d(i, j, k) = s_local(i, j, k) - rS(j)
            end do
        end do
    end do
    !$omp end target teams distribute parallel do
    print *,'2', 'p_wrk3d: ', sum(p_wrk3d)
#else
    do j = 1, jmax !offload
        p_wrk3d(:, j, :) = s_local(:, j, :) - rS(j)
    end do
    print *,'2', 'p_wrk3d: ', sum(p_wrk3d)
#endif
    tmp1 = p_wrk3d*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, tmp1, rS2(1), wrk1d)
    tmp1 = p_wrk3d*tmp1
    call AVG_IK_V(imax, jmax, kmax, tmp1, rS3(1), wrk1d)
    tmp1 = p_wrk3d*tmp1
    call AVG_IK_V(imax, jmax, kmax, tmp1, rS4(1), wrk1d)

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        fS2(:) = rS2(:)
        fS3(:) = rS3(:)
        fS4(:) = rS4(:)

    else
#ifndef USE_APU
        !$omp target teams distribute parallel do collapse(3) default(shared) private(i,j,k)
        do j = 1, jmax
            do i = 1, imax
                do k = 1, kmax
                    p_wrk3d(i, j, k) = s_local(i, j, k) - fS(j)
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
        print *,'3', 'p_wrk3d: ', sum(p_wrk3d)
#else
        do j = 1, jmax
            p_wrk3d(:, j, :) = s_local(:, j, :) - fS(j)
        end do
        print *,'3', 'p_wrk3d: ', sum(p_wrk3d)
#endif
        tmp1 = p_wrk3d*p_wrk3d*rho
        call AVG_IK_V(imax, jmax, kmax, tmp1, fS2(1), wrk1d)
        tmp1 = p_wrk3d*tmp1
        call AVG_IK_V(imax, jmax, kmax, tmp1, fS3(1), wrk1d)
        tmp1 = p_wrk3d*tmp1
        call AVG_IK_V(imax, jmax, kmax, tmp1, fS4(1), wrk1d)
        fS2(:) = fS2(:)/rR(:)
        fS3(:) = fS3(:)/rR(:)
        fS4(:) = fS4(:)/rR(:)

    end if

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fS2(1), Rss_y(1))

    ! -----------------------------------------------------------------------
    ! Cross terms
#ifndef USE_APU
        !$omp target teams distribute parallel do collapse(3) default(shared) private(i,j,k)
        do j = 1, jmax
            do i = 1, imax
                do k = 1, kmax
                    p_wrk3d(i, j, k) = s_local(i, j, k) - fS(j)
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
        print *,'4', 'p_wrk3d: ', sum(p_wrk3d)
#else
        do j = 1, jmax !offload
            p_wrk3d(:, j, :) = s_local(:, j, :) - fS(j)
        end do
        print *,'4', 'p_wrk3d: ', sum(p_wrk3d)
#endif
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) p_wrk3d = p_wrk3d*rho

#ifndef USE_APU
    !$omp target teams distribute parallel do collapse(3) default(shared) private(i,j,k)
    do j = 1, jmax !offload
        do i = 1, imax
            do k = 1, kmax
                dsdx(i, j, k) = p_wrk3d(i, j, k)*(u(i, j, k) - fU(j))
                dsdy(i, j, k) = p_wrk3d(i, j, k)*(v(i, j, k) - fV(j))
                dsdz(i, j, k) = p_wrk3d(i, j, k)*(w(i, j, k) - fW(j))
            end do
        end do
    end do
    !$omp end target teams distribute parallel do
    print *,'5', 'dsdx: ', sum(dsdx)
    print *,'5', 'dsdy: ', sum(dsdy)
    print *,'5', 'dsdz: ', sum(dsdz)
#else
    do j = 1, jmax !offload
        dsdx(:, j, :) = p_wrk3d(:, j, :)*(u(:, j, :) - fU(j))
        dsdy(:, j, :) = p_wrk3d(:, j, :)*(v(:, j, :) - fV(j))
        dsdz(:, j, :) = p_wrk3d(:, j, :)*(w(:, j, :) - fW(j))
    end do
    print *,'5', 'dsdx: ', sum(dsdx)
    print *,'5', 'dsdy: ', sum(dsdy)
    print *,'5', 'dsdz: ', sum(dsdz)
#endif
    call AVG_IK_V(imax, jmax, kmax, dsdx, Rsu(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdy, Rsv(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdz, Rsw(1), wrk1d)
    Rsu(:) = Rsu(:)/rR(:)
    Rsv(:) = Rsv(:)/rR(:)
    Rsw(:) = Rsw(:)/rR(:)

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rsu(1), Rsu_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rsv(1), Rsv_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rsw(1), Rsw_y(1))

    ! -----------------------------------------------------------------------
    ! turbulent transport terms
#ifndef USE_APU
    !$omp target teams distribute parallel do collapse(3) default(shared) private(i,j,k)
    do j = 1, jmax !offload
        do i = 1, imax
            do k = 1, kmax
                tmp1(i, j, k) = dsdy(i, j, k)*(s_local(i, j, k) - fS(j))
                dsdx(i, j, k) = dsdx(i, j, k)*(v(i, j, k) - fV(j))
                dsdy(i, j, k) = dsdy(i, j, k)*(v(i, j, k) - fV(j))
                dsdz(i, j, k) = dsdz(i, j, k)*(v(i, j, k) - fV(j))
            end do
        end do
    end do
    !$omp end target teams distribute parallel do
    print *,'6', 'tmp1: ', sum(tmp1)
    print *,'6', 'dsdx: ', sum(dsdx)
    print *,'6', 'dsdy: ', sum(dsdy)
    print *,'6', 'dsdz: ', sum(dsdz)
#else
    do j = 1, jmax !offload
        tmp1(:, j, :) = dsdy(:, j, :)*(s_local(:, j, :) - fS(j))
        dsdx(:, j, :) = dsdx(:, j, :)*(v(:, j, :) - fV(j))
        dsdy(:, j, :) = dsdy(:, j, :)*(v(:, j, :) - fV(j))
        dsdz(:, j, :) = dsdz(:, j, :)*(v(:, j, :) - fV(j))
    end do
    print *,'6', 'tmp1: ', sum(tmp1)
    print *,'6', 'dsdx: ', sum(dsdx)
    print *,'6', 'dsdy: ', sum(dsdy)
    print *,'6', 'dsdz: ', sum(dsdz)
#endif
    call AVG_IK_V(imax, jmax, kmax, tmp1, Tssy1(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdx, Tsuy1(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdy, Tsvy1(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdz, Tswy1(1), wrk1d)

    ! -----------------------------------------------------------------------
    ! Pressure terms in transport equations
    call AVG_IK_V(imax, jmax, kmax, p_loc, rP(1), wrk1d)

    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s_local, dsdx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s_local, dsdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s_local, dsdz)
#ifndef USE_APU
    !$omp target teams distribute parallel do collapse(3) default(shared) private(i,j,k)
    do j = 1, jmax
        do i = 1, imax
            do k = 1,kmax
                tmp1(i, j, k) = (p_loc(i, j, k) - rP(j))*(s_local(i, j, k) - fS(j))
                dsdx(i, j, k) = (p_loc(i, j, k) - rP(j))*dsdx(i, j, k)
                dsdy(i, j, k) = (p_loc(i, j, k) - rP(j))*(dsdy(i, j, k) - fS_y(j))
                dsdz(i, j, k) = (p_loc(i, j, k) - rP(j))*dsdz(i, j, k)
            end do
        end do
    end do
    !$omp end target teams distribute parallel do
    print *,'7', 'tmp1: ', sum(tmp1)
    print *,'7', 'dsdx: ', sum(dsdx)
    print *,'7', 'dsdy: ', sum(dsdy)
    print *,'7', 'dsdz: ', sum(dsdz)
#else
    do j = 1, jmax
        tmp1(:, j, :) = (p_loc(:, j, :) - rP(j))*(s_local(:, j, :) - fS(j))
        dsdx(:, j, :) = (p_loc(:, j, :) - rP(j))*dsdx(:, j, :)
        dsdy(:, j, :) = (p_loc(:, j, :) - rP(j))*(dsdy(:, j, :) - fS_y(j))
        dsdz(:, j, :) = (p_loc(:, j, :) - rP(j))*dsdz(:, j, :)
    end do
    print *,'7', 'tmp1: ', sum(tmp1)
    print *,'7', 'dsdx: ', sum(dsdx)
    print *,'7', 'dsdy: ', sum(dsdy)
    print *,'7', 'dsdz: ', sum(dsdz)
#endif
    call AVG_IK_V(imax, jmax, kmax, tmp1, Tsvy3(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdx, PIsu(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdy, PIsv(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdz, PIsw(1), wrk1d)

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rP(1), aux(1))
    Gsv(:) = (rS(:) - fS(:))*aux(:)

    ! #######################################################################
    ! Cross-scalar terms
    ! #######################################################################
    k = ig(8) - 1

    do is_loc = 1, inb_scal_array
        call AVG_IK_V(imax, jmax, kmax, s(:, :, :, is_loc), aux(1), wrk1d)
        do j = 1, jmax
            tmp1(:, j, :) = (s(:, j, :, is_loc) - aux(j))*(s_local(:, j, :) - fS(j))
            tmp2(:, j, :) = tmp1(:, j, :)*(s_local(:, j, :) - fS(j))
        end do
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, tmp1, mean2d(1, k), wrk1d)
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, tmp2, mean2d(1, k), wrk1d)
    end do

    ! #######################################################################
    ! Source terms
    ! #######################################################################
    if (infraredProps%active(is)) then       ! Radiation in tmp1 and dsdx
        call Radiation_Infrared_Y(infraredProps, imax, jmax, kmax, fdm_Int0, s, tmp1, tmp2, tmp3, dsdy, dsdx)
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_WEIGHT_INPLACE(imax, jmax, kmax, ribackground, tmp1)
        end if
    end if

    if (sedimentationProps%active(is)) then      ! Transport in tmp3 and dsdz
        call Microphysics_Sedimentation(sedimentationProps, imax, jmax, kmax, is, g(2), s, tmp3, dsdy, dsdz)
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_WEIGHT_INPLACE(imax, jmax, kmax, ribackground, tmp3)
        end if
    end if

    if (is > inb_scal) then     ! Diagnostic variables; I overwrite tmp1 and dsdx and recalculate them.
        if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            coefQ = 1.0_wp                        ! Coefficient in the evaporation term
            coefR = 0.0_wp                        ! Coefficient in the radiation term
            coefT = 0.0_wp                        ! Coefficient in the transport term
            if (is == inb_scal_array + 1) then ! Default values are for liquid; defining them for buoyancy
                coefQ = buoyancy%parameters(inb_scal_array)/froude
                coefR = buoyancy%parameters(inb_scal)/froude
                if (sedimentationProps%active(is)) then
                    do is_loc = 1, inb_scal
                        coefT = coefT + sedimentationProps%parameters(is_loc)/settling*buoyancy%parameters(is_loc)/froude
                    end do
                end if
            end if

            call THERMO_AIRWATER_LINEAR_SOURCE(imax*jmax*kmax, s, dsdx, dsdy, dsdz) ! calculate xi in dsdx
            call FI_GRADIENT(imax, jmax, kmax, dsdx, tmp2, tmp1)

            dummy = -diff*coefQ
            tmp2 = dsdz*tmp2*dummy         ! evaporation source

            if (sedimentationProps%active(is) .or. infraredProps%active(is)) then ! preparing correction terms into dsdz
                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), dsdx, tmp1)
                dsdz = dsdz*tmp1
            end if

            if (infraredProps%active(is)) then ! radiation source; needs dsdy
                ! only valid for IR_Bulk1D_Liquid, where tmp2, tmp3, dsdy are not used
                call Radiation_Infrared_Y(infraredProps, imax, jmax, kmax, fdm_Int0, s, tmp1, tmp2, tmp3, dsdy, dsdx)
                dummy = thermo_param(2)*coefQ
                tmp1 = tmp1*(coefR + dsdy*dummy)
                dsdx = dsdx*dsdz*dummy
            else
                tmp1 = 0.0_wp; dsdx = 0.0_wp
            end if

            if (sedimentationProps%active(is)) then ! transport source; needs dsdy
                dummy = coefQ
                tmp3 = tmp3*(coefT + dsdy*dummy)
                ! Correction term needs dsdz; the following call assumes incompressible mode
                call Microphysics_Sedimentation(sedimentationProps, imax, jmax, kmax, is, g(2), s, dsdy, dsdy, dsdy)
                dsdz = dsdy*dsdz*dummy
            else
                tmp3 = 0.0_wp; dsdz = 0.0_wp
            end if

        else
            if (buoyancy%type /= EQNS_EXPLICIT) then
                call FI_GRADIENT(imax, jmax, kmax, s, dsdx, dsdy)
                call Gravity_Buoyancy_Source(buoyancy, imax, jmax, kmax, s, dsdx, tmp1) ! dsdx contains gradient
                tmp1 = tmp1*diff/froude
            end if

        end if

    end if

    ! -----------------------------------------------------------------------
    ! Calculating averages
    k = ig(1) + im
    if (infraredProps%active(is)) then
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, tmp1, mean2d(1, k), wrk1d)
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, dsdx, mean2d(1, k), wrk1d) ! correction term or flux
    end if
    if (imixture == MIXT_TYPE_AIRWATER_LINEAR .or. imixture == MIXT_TYPE_AIRWATER) then           ! evaporation
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, tmp2, mean2d(1, k), wrk1d)
    end if
    if (sedimentationProps%active(is)) then
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, tmp3, mean2d(1, k), wrk1d)
        k = k + 1; call AVG_IK_V(imax, jmax, kmax, dsdz, mean2d(1, k), wrk1d) ! correction term or flux
    end if

    p_wrk3d = 0.0_wp ! total
    if (infraredProps%active(is)) then
        p_wrk3d = p_wrk3d + tmp1
    end if
    if (is > inb_scal) then
        if (imixture == MIXT_TYPE_AIRWATER_LINEAR .or. imixture == MIXT_TYPE_AIRWATER) then
            p_wrk3d = p_wrk3d + tmp2
        end if
    end if
    if (sedimentationProps%active(is)) then
        p_wrk3d = p_wrk3d + tmp3
    end if
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rQ(1), wrk1d)
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) p_wrk3d = p_wrk3d*rho
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, fQ(1), wrk1d)
    fQ(:) = fQ(:)/rR(:)

#ifndef USE_APU
    !$omp target teams distribute parallel do collapse(3) default(shared) private(i,j,k)
    do j = 1, jmax
        do i = 1, imax
            do k = 1, kmax
                tmp1(i, j, k) = (s_local(i, j, k) - fS(j))*p_wrk3d(i, j, k)
                dsdx(i, j, k) = (u(i, j, k) - fU(j))*p_wrk3d(i, j, k)
                dsdy(i, j, k) = (v(i, j, k) - fV(j))*p_wrk3d(i, j, k)
                dsdz(i, j, k) = (w(i, j, k) - fW(j))*p_wrk3d(i, j, k)
            end do
        end do
    end do
    !$omp end target teams distribute parallel do
    print *,'8', 'tmp1: ', sum(tmp1)
    print *,'8', 'dsdx: ', sum(dsdx)
    print *,'8', 'dsdy: ', sum(dsdy)
    print *,'8', 'dsdz: ', sum(dsdz)
#else
    do j = 1, jmax
        tmp1(:, j, :) = (s_local(:, j, :) - fS(j))*p_wrk3d(:, j, :)
        dsdx(:, j, :) = (u(:, j, :) - fU(j))*p_wrk3d(:, j, :)
        dsdy(:, j, :) = (v(:, j, :) - fV(j))*p_wrk3d(:, j, :)
        dsdz(:, j, :) = (w(:, j, :) - fW(j))*p_wrk3d(:, j, :)
    end do
    print *,'8', 'tmp1: ', sum(tmp1)
    print *,'8', 'dsdx: ', sum(dsdx)
    print *,'8', 'dsdy: ', sum(dsdy)
    print *,'8', 'dsdz: ', sum(dsdz)
#endif

    call AVG_IK_V(imax, jmax, kmax, tmp1, Qss(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdx, Qsu(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdy, Qsv(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dsdz, Qsw(1), wrk1d)
    Qss(:) = Qss(:)*2.0_wp/rR(:)
    Qsu(:) = Qsu(:)/rR(:)
    Qsv(:) = Qsv(:)/rR(:)
    Qsw(:) = Qsw(:)/rR(:)

    ! #######################################################################
    ! Derivatives
    ! #######################################################################
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s_local, dsdx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s_local, dsdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s_local, dsdz)

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdx*dsdx + dsdy*dsdy + dsdz*dsdz
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Ess(1), wrk1d)
    Ess(:) = Ess(:)*diff*2.0_wp

    ! -----------------------------------------------------------------------
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, tmp1)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, tmp2)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, tmp3)

    ! Transport term
    p_wrk3d = (tmp2*2.0_wp - tmp1 - tmp3)*c23*visc
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_yy(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = -(p_wrk3d(:, j, :) - Tau_yy(j))*(s_local(:, j, :) - fS(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tsvy2(1), wrk1d)

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_yy(1), Tau_yy_y(1))

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdx*((tmp1*2.0_wp - tmp2 - tmp3)*c23*visc + tmp1*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Esu(1), wrk1d)

    p_wrk3d = dsdy*((tmp2*2.0_wp - tmp1 - tmp3)*c23*visc + tmp2*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Esv(1), wrk1d)

    p_wrk3d = dsdz*((tmp3*2.0_wp - tmp1 - tmp2)*c23*visc + tmp3*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Esw(1), wrk1d)

    ! -----------------------------------------------------------------------
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, tmp2)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, tmp1)

    ! Transport term
    p_wrk3d = (tmp1 + tmp2)*visc
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_yx(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = -(p_wrk3d(:, j, :) - Tau_yx(j))*(s_local(:, j, :) - fS(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tsuy2(1), wrk1d)

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_yx(1), Tau_yx_y(1))

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdy*((tmp1 + tmp2)*visc + tmp1*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esu(:) = Esu(:) + aux(:)

    p_wrk3d = dsdx*((tmp1 + tmp2)*visc + tmp2*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esv(:) = Esv(:) + aux(:)

    ! -----------------------------------------------------------------------
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, tmp3)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, tmp2)

    ! Transport term
    p_wrk3d = (tmp3 + tmp2)*visc
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_yz(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = -(p_wrk3d(:, j, :) - Tau_yz(j))*(s_local(:, j, :) - fS(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tswy2(1), wrk1d)

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_yz(1), Tau_yz_y(1))

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdz*((tmp3 + tmp2)*visc + tmp2*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esv(:) = Esv(:) + aux(:)

    p_wrk3d = dsdy*((tmp3 + tmp2)*visc + tmp3*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esw(:) = Esw(:) + aux(:)

    ! -----------------------------------------------------------------------
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, tmp3)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, tmp1)

    ! Dissipation terms; mean terms substracted below
    p_wrk3d = dsdz*((tmp3 + tmp1)*visc + tmp1*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esu(:) = Esu(:) + aux(:)

    p_wrk3d = dsdx*((tmp3 + tmp1)*visc + tmp3*diff)
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, aux(1), wrk1d)
    Esw(:) = Esw(:) + aux(:)

    ! -----------------------------------------------------------------------
    ! Moments
#ifndef USE_APU
    !$omp target teams distribute parallel do collapse(3) default(shared) private(i,j,k)
    do j = 1, jmax !offload
        do i = 1, imax
            do k = 1, kmax
                p_wrk3d(i, j, k) = dsdy(i, j, k) - rS_y(j)
            end do
        end do
    end do
    !$omp end target teams distribute parallel do
    print *,'9', 'p_wrk3d: ', sum(p_wrk3d)
#else
    do j = 1, jmax
        p_wrk3d(:, j, :) = dsdy(:, j, :) - rS_y(j)
    end do
    print *,'9', 'p_wrk3d: ', sum(p_wrk3d)
#endif
    tmp1 = dsdx*dsdx
    tmp2 = p_wrk3d*p_wrk3d
    tmp3 = dsdz*dsdz
    call AVG_IK_V(imax, jmax, kmax, tmp1, S_x2(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp2, S_y2(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp3, S_z2(1), wrk1d)

    tmp1 = tmp1*dsdx
    tmp2 = tmp2*p_wrk3d
    tmp3 = tmp3*dsdz
    call AVG_IK_V(imax, jmax, kmax, tmp1, S_x3(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp2, S_y3(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp3, S_z3(1), wrk1d)

    tmp1 = tmp1*dsdx
    tmp2 = tmp2*p_wrk3d
    tmp3 = tmp3*dsdz
    call AVG_IK_V(imax, jmax, kmax, tmp1, S_x4(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp2, S_y4(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, tmp3, S_z4(1), wrk1d)

    ! -----------------------------------------------------------------------
    ! Molecular fluxes
    if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) dsdy = dsdy*vis
    call AVG_IK_V(imax, jmax, kmax, dsdy, Fy(1), wrk1d)

    ! Contribution to turbulent transport
#ifndef USE_APU
    !$omp target teams distribute parallel do collapse(3) default(shared) private(i,j,k)
    do j = 1, jmax
        do i = 1, imax
            do k = 1, kmax
                p_wrk3d(i, j, k) = (dsdy(i, j, k) - Fy(j))*(s_local(i, j, k) - fS(j))
                tmp1(i, j, k) = (dsdy(i, j, k) - Fy(j))*(u(i, j, k) - fU(j))
                tmp2(i, j, k) = (dsdy(i, j, k) - Fy(j))*(v(i, j, k) - fV(j))
                tmp3(i, j, k) = (dsdy(i, j, k) - Fy(j))*(w(i, j, k) - fW(j))
            end do
        end do
    end do
    !$omp end target teams distribute parallel do
    print *,'10', 'p_wrk3d: ', sum(p_wrk3d)
    print *,'10', 'tmp1: ', sum(tmp1)
    print *,'10', 'tmp2: ', sum(tmp2)
    print *,'10', 'tmp3: ', sum(tmp3)
#else
    do j = 1, jmax
        p_wrk3d(:, j, :) = (dsdy(:, j, :) - Fy(j))*(s_local(:, j, :) - fS(j))
        tmp1(:, j, :) = (dsdy(:, j, :) - Fy(j))*(u(:, j, :) - fU(j))
        tmp2(:, j, :) = (dsdy(:, j, :) - Fy(j))*(v(:, j, :) - fV(j))
        tmp3(:, j, :) = (dsdy(:, j, :) - Fy(j))*(w(:, j, :) - fW(j))
    end do
    print *,'10', 'p_wrk3d: ', sum(p_wrk3d)
    print *,'10', 'tmp1: ', sum(tmp1)
    print *,'10', 'tmp2: ', sum(tmp2)
    print *,'10', 'tmp3: ', sum(tmp3)
#endif
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tssy2(1), wrk1d)
    Tssy2(:) = -Tssy2(:)*diff*2.0_wp
    call AVG_IK_V(imax, jmax, kmax, tmp1, aux(1), wrk1d)
    Tsuy2(:) = Tsuy2(:) - aux(:)*diff
    call AVG_IK_V(imax, jmax, kmax, tmp2, aux(1), wrk1d)
    Tsvy2(:) = Tsvy2(:) - aux(:)*diff
    call AVG_IK_V(imax, jmax, kmax, tmp3, aux(1), wrk1d)
    Tswy2(:) = Tswy2(:) - aux(:)*diff

    Fy(:) = Fy(:)*diff
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Fy(1), Fy_y(1))

    ! Contribution to dissipation
    Ess(:) = (Ess(:) - Fy(:)*rS_y(:) - Fy(:)*rS_y(:))/rR(:)
    Esu(:) = (Esu(:) - Tau_yx(:)*rS_y(:) - Fy(:)*rU_y(:))/rR(:)
    Esv(:) = (Esv(:) - Tau_yy(:)*rS_y(:) - Fy(:)*rV_y(:))/rR(:)
    Esw(:) = (Esw(:) - Tau_yz(:)*rS_y(:) - Fy(:)*rW_y(:))/rR(:)

    ! #######################################################################
    ! Source terms in transport equations
    ! #######################################################################
    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        if (buoyancy%type == EQNS_BOD_EXPLICIT) then
            call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, p_wrk3d)
        else
            call Gravity_Buoyancy(buoyancy, imax, jmax, kmax, s, p_wrk3d, bbackground)
        end if
        dummy = 1.0_wp/froude
        p_wrk3d = p_wrk3d*dummy
    else
        p_wrk3d = rho*buoyancy%vector(2)
    end if
    do j = 1, jmax
        p_wrk3d(:, j, :) = (s_local(:, j, :) - fS(j))*p_wrk3d(:, j, :)
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Bsv(1), wrk1d)
    Bsv(:) = Bsv(:)/rR(:)

    ! #######################################################################
    ! Complete budget equations
    ! #######################################################################
    ! Transport terms
    aux(:) = Tssy1(:) + Tssy2(:)
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), aux(1), Tssy_y(1))
    aux(:) = Tsuy1(:) + Tsuy2(:)
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), aux(1), Tsuy_y(1))
    aux(:) = Tsvy1(:) + Tsvy2(:) + Tsvy3(:)
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), aux(1), Tsvy_y(1))
    aux(:) = Tswy1(:) + Tswy2(:)
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), aux(1), Tswy_y(1))

    ! Convective terms
    Css(:) = -fV(:)*Rss_y(:)
    Csu(:) = -fV(:)*Rsu_y(:)
    Csv(:) = -fV(:)*Rsv_y(:)
    Csw(:) = -fV(:)*Rsw_y(:)

    ! Production terms
    Pss(:) = -Rsv(:)*fS_y(:)*2.0_wp
    Psu(:) = -Rsv(:)*fU_y(:) - Rvu(:)*fS_y(:)
    Psv(:) = -Rsv(:)*fV_y(:) - Rvv(:)*fS_y(:)
    Psw(:) = -Rsv(:)*fW_y(:) - Rvw(:)*fS_y(:)

    ! Diffusion variable-density terms
    Dss(:) = (rS(:) - fS(:))*Fy_y(:)*2.0_wp
    Dsu(:) = (rS(:) - fS(:))*Tau_yx_y(:) + (rU(:) - fU(:))*Fy_y(:)
    Dsv(:) = (rS(:) - fS(:))*Tau_yy_y(:) + (rV(:) - fV(:))*Fy_y(:)
    Dsw(:) = (rS(:) - fS(:))*Tau_yz_y(:) + (rW(:) - fW(:))*Fy_y(:)

    ! Coriolis terms
    dummy = coriolis%vector(2)
    Fsu(:) = dummy*Rsw(:)
    Fsw(:) = -dummy*Rsu(:)

    ! Transient terms
    Rss_t(:) = Css(:) + Pss(:) - Ess(:) + Qss(:) + (Dss(:) - Tssy_y(:))/rR(:)
    Rsu_t(:) = Csu(:) + Psu(:) - Esu(:) + Bsu(:) - Fsu(:) + Qsu(:) + (PIsu(:) + Dsu(:) - Gsu(:) - Tsuy_y(:))/rR(:)
    Rsv_t(:) = Csv(:) + Psv(:) - Esv(:) + Bsv(:) - Fsv(:) + Qsv(:) + (PIsv(:) + Dsv(:) - Gsv(:) - Tsvy_y(:))/rR(:)
    Rsw_t(:) = Csw(:) + Psw(:) - Esw(:) + Bsw(:) - Fsw(:) + Qsw(:) + (PIsw(:) + Dsw(:) - Gsw(:) - Tswy_y(:))/rR(:)

    ! ###################################################################
    ! Output
    ! #######################################################################
    ! 11 t-dependent variables, for consistency with old format
    ! ng = ng +1
    ! groupname(ng) = ''
    ! varname(ng)   = 'dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy'
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''

    write (line1, *) is; line1 = 'avg'//trim(adjustl(line1))//'s'
    write (name, *) itime; name = trim(adjustl(line1))//trim(adjustl(name))
    call IO_WRITE_AVERAGES(name, itime, rtime, jmax, nv, ng, g(2)%nodes, varname, groupname, mean2d)

    return
end subroutine AVG_SCAL_XZ
