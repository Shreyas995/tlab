#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!#
!# Assumes statistical homogeneity in xOz, so that the corresponding
!# partial derivative terms are assumed to be zero.
!#
!# In the incompressible case, the array p has been
!# pointed to dudz and the pressure field is stored there; do not
!# use array dudz until pressure block
!#
!# Reynolds and Favre averages
!#
!########################################################################

subroutine AVG_FLOW_XZ(q, s, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, mean2d)
    use TLab_Constants, only: MAX_AVG_TEMPORAL
    use TLab_Constants, only: efile, lfile, wp, wi
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
#endif
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_flow_array, inb_scal_array, inb_scal
    use Tlab_Background, only: sbg, rbg
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Arrays, only: wrk1d
    use TLab_Pointers_3D, only: u, v, w, rho, T, e, rho, vis, p_wrk3d
    use FDM, only: g
    use Thermodynamics, only: imode_thermo, THERMO_TYPE_ANELASTIC, THERMO_TYPE_COMPRESSIBLE
    use Thermodynamics, only: itransport, EQNS_TRANS_POWERLAW
    use Thermodynamics, only: imixture, MIXT_TYPE_AIRWATER
    use Thermodynamics, only: PREF_1000, CRATIO_INV
    use Thermodynamics, only: Thermo_Psat_Polynomial
    use Thermo_Anelastic
    use THERMO_AIRWATER
    use THERMO_CALORIC
    use NavierStokes
    use Gravity, only: buoyancy, bbackground, Gravity_Buoyancy, Gravity_Buoyancy_Source, EQNS_BOD_NONE, EQNS_BOD_EXPLICIT
    use Rotation, only: coriolis
    use OPR_Partial
    use IBM_VARS, only: imode_ibm, gamma_0, gamma_1
    use Averages, only: AVG_IK_V

    implicit none

    real(wp), intent(IN) :: q(imax, jmax, kmax, inb_flow_array)
    real(wp), intent(IN) :: s(imax, jmax, kmax, inb_scal_array)
    real(wp), dimension(imax, jmax, kmax), intent(INOUT) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(wp), intent(INOUT) :: mean2d(jmax, MAX_AVG_TEMPORAL)

    target q, dudz

    ! -------------------------------------------------------------------
    integer, parameter :: MAX_VARS_GROUPS = 20
    integer i, j, k, bcs(2, 2)
    real(wp) dummy
    real(wp) c23

    integer ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nv

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
        p_loc => dudz
    end if

    c23 = 2.0_wp/3.0_wp

    ! Variable definition and memory management
    ! -----------------------------------------------------------------------
    ng = 1; ig(ng) = 1
#define rR(j)     mean2d(j,ig(1)  )
#define rU(j)     mean2d(j,ig(1)+1)
#define rV(j)     mean2d(j,ig(1)+2)
#define rW(j)     mean2d(j,ig(1)+3)
#define rP(j)     mean2d(j,ig(1)+4)
#define rT(j)     mean2d(j,ig(1)+5)
#define re(j)     mean2d(j,ig(1)+6)
#define rh(j)     mean2d(j,ig(1)+7)
#define rs(j)     mean2d(j,ig(1)+8)
#define rB(j)     mean2d(j,ig(1)+9)
#define fU(j)     mean2d(j,ig(1)+10)
#define fV(j)     mean2d(j,ig(1)+11)
#define fW(j)     mean2d(j,ig(1)+12)
#define fT(j)     mean2d(j,ig(1)+13)
#define fe(j)     mean2d(j,ig(1)+14)
#define fh(j)     mean2d(j,ig(1)+15)
#define fs(j)     mean2d(j,ig(1)+16)
    sg(ng) = 17

    groupname(ng) = 'Mean'
    varname(ng) = 'rR rU rV rW rP rT re rh rs rB fU fV fW fT fe fh fs'

    if (imode_ibm == 1) then
        varname(ng) = trim(adjustl(varname(ng)))//' eps_0 eps_1'
        sg(ng) = sg(ng) + 2
#define ep_0(j)   mean2d(j,ig(1)+17)
#define ep_1(j)   mean2d(j,ig(1)+18)
    end if

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Tke(j)    mean2d(j,ig(2)  )
#define Rxx(j)    mean2d(j,ig(2)+1)
#define Ryy(j)    mean2d(j,ig(2)+2)
#define Rzz(j)    mean2d(j,ig(2)+3)
#define Rxy(j)    mean2d(j,ig(2)+4)
#define Rxz(j)    mean2d(j,ig(2)+5)
#define Ryz(j)    mean2d(j,ig(2)+6)
#define rP2(j)    mean2d(j,ig(2)+7)
#define rR2(j)    mean2d(j,ig(2)+8)
#define rT2(j)    mean2d(j,ig(2)+9)
#define fT2(j)    mean2d(j,ig(2)+10)
#define re2(j)    mean2d(j,ig(2)+11)
#define fe2(j)    mean2d(j,ig(2)+12)
#define rh2(j)    mean2d(j,ig(2)+13)
#define fh2(j)    mean2d(j,ig(2)+14)
#define rs2(j)    mean2d(j,ig(2)+15)
#define fs2(j)    mean2d(j,ig(2)+16)
    sg(ng) = 17

    groupname(ng) = 'Fluctuations'
    varname(ng) = 'Tke Rxx Ryy Rzz Rxy Rxz Ryz rP2 rR2 rT2 fT2 re2 fe2 rh2 fh2 rs2 fs2'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define vortx(j)  mean2d(j,ig(3)  )
#define vorty(j)  mean2d(j,ig(3)+1)
#define vortz(j)  mean2d(j,ig(3)+2)
#define vortx2(j) mean2d(j,ig(3)+3)
#define vorty2(j) mean2d(j,ig(3)+4)
#define vortz2(j) mean2d(j,ig(3)+5)
    sg(ng) = 6

    groupname(ng) = 'Vorticity'
    varname(ng) = 'Wx Wy Wz Wx2 Wy2 Wz2'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rxx_t(j)  mean2d(j,ig(4)  )
#define Bxx(j)    mean2d(j,ig(4)+1)
#define Cxx(j)    mean2d(j,ig(4)+2)
#define Pxx(j)    mean2d(j,ig(4)+3)
#define Exx(j)    mean2d(j,ig(4)+4)
#define PIxx(j)   mean2d(j,ig(4)+5)
#define Fxx(j)    mean2d(j,ig(4)+6)
#define Txxy_y(j) mean2d(j,ig(4)+7)
#define Txxy(j)   mean2d(j,ig(4)+8)
#define Gxx(j)    mean2d(j,ig(4)+9)
#define Dxx(j)    mean2d(j,ig(4)+10)
    sg(ng) = 11

    groupname(ng) = 'RxxBudget'
    varname(ng) = 'Rxx_t Bxx Cxx Pxx Exx PIxx Fxx Txxy_y Txxy Gxx Dxx'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Ryy_t(j)  mean2d(j,ig(5)  )
#define Byy(j)    mean2d(j,ig(5)+1)
#define Cyy(j)    mean2d(j,ig(5)+2)
#define Pyy(j)    mean2d(j,ig(5)+3)
#define Eyy(j)    mean2d(j,ig(5)+4)
#define PIyy(j)   mean2d(j,ig(5)+5)
#define Fyy(j)    mean2d(j,ig(5)+6)
#define Tyyy_y(j) mean2d(j,ig(5)+7)
#define Tyyy(j)   mean2d(j,ig(5)+8)
#define Gyy(j)    mean2d(j,ig(5)+9)
#define Dyy(j)    mean2d(j,ig(5)+10)
    sg(ng) = 11

    groupname(ng) = 'RyyBudget'
    varname(ng) = 'Ryy_t Byy Cyy Pyy Eyy PIyy Fyy Tyyy_y Tyyy Gyy Dyy'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rzz_t(j)  mean2d(j,ig(6)  )
#define Bzz(j)    mean2d(j,ig(6)+1)
#define Czz(j)    mean2d(j,ig(6)+2)
#define Pzz(j)    mean2d(j,ig(6)+3)
#define Ezz(j)    mean2d(j,ig(6)+4)
#define PIzz(j)   mean2d(j,ig(6)+5)
#define Fzz(j)    mean2d(j,ig(6)+6)
#define Tzzy_y(j) mean2d(j,ig(6)+7)
#define Tzzy(j)   mean2d(j,ig(6)+8)
#define Gzz(j)    mean2d(j,ig(6)+9)
#define Dzz(j)    mean2d(j,ig(6)+10)
    sg(ng) = 11

    groupname(ng) = 'RzzBudget'
    varname(ng) = 'Rzz_t Bzz Czz Pzz Ezz PIzz Fzz Tzzy_y Tzzy Gzz Dzz'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rxy_t(j)  mean2d(j,ig(7)  )
#define Bxy(j)    mean2d(j,ig(7)+1)
#define Cxy(j)    mean2d(j,ig(7)+2)
#define Pxy(j)    mean2d(j,ig(7)+3)
#define Exy(j)    mean2d(j,ig(7)+4)
#define PIxy(j)   mean2d(j,ig(7)+5)
#define Fxy(j)    mean2d(j,ig(7)+6)
#define Txyy_y(j) mean2d(j,ig(7)+7)
#define Txyy(j)   mean2d(j,ig(7)+8)
#define Gxy(j)    mean2d(j,ig(7)+9)
#define Dxy(j)    mean2d(j,ig(7)+10)
    sg(ng) = 11

    groupname(ng) = 'RxyBudget'
    varname(ng) = 'Rxy_t Bxy Cxy Pxy Exy PIxy Fxy Txyy_y Txyy Gxy Dxy'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rxz_t(j)  mean2d(j,ig(8)  )
#define Bxz(j)    mean2d(j,ig(8)+1)
#define Cxz(j)    mean2d(j,ig(8)+2)
#define Pxz(j)    mean2d(j,ig(8)+3)
#define Exz(j)    mean2d(j,ig(8)+4)
#define PIxz(j)   mean2d(j,ig(8)+5)
#define Fxz(j)    mean2d(j,ig(8)+6)
#define Txzy_y(j) mean2d(j,ig(8)+7)
#define Txzy(j)   mean2d(j,ig(8)+8)
#define Gxz(j)    mean2d(j,ig(8)+9)
#define Dxz(j)    mean2d(j,ig(8)+10)
    sg(ng) = 11

    groupname(ng) = 'RxzBudget'
    varname(ng) = 'Rxz_t Bxz Cxz Pxz Exz PIxz Fxz Txzy_y Txzy Gxz Dxz'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Ryz_t(j)  mean2d(j,ig(9)  )
#define Byz(j)    mean2d(j,ig(9)+1)
#define Cyz(j)    mean2d(j,ig(9)+2)
#define Pyz(j)    mean2d(j,ig(9)+3)
#define Eyz(j)    mean2d(j,ig(9)+4)
#define PIyz(j)   mean2d(j,ig(9)+5)
#define Fyz(j)    mean2d(j,ig(9)+6)
#define Tyzy_y(j) mean2d(j,ig(9)+7)
#define Tyzy(j)   mean2d(j,ig(9)+8)
#define Gyz(j)    mean2d(j,ig(9)+9)
#define Dyz(j)    mean2d(j,ig(9)+10)
    sg(ng) = 11

    groupname(ng) = 'RyzBudget'
    varname(ng) = 'Ryz_t Byz Cyz Pyz Eyz PIyz Fyz Tyzy_y Tyzy Gyz Dyz'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Tke_t(j)  mean2d(j,ig(10)  )
#define Buo(j)    mean2d(j,ig(10)+1)
#define Con(j)    mean2d(j,ig(10)+2)
#define Prd(j)    mean2d(j,ig(10)+3)
#define Eps(j)    mean2d(j,ig(10)+4)
#define Pi(j)     mean2d(j,ig(10)+5)
#define Ty_y(j)   mean2d(j,ig(10)+6)
#define Ty1(j)    mean2d(j,ig(10)+7)
#define Ty2(j)    mean2d(j,ig(10)+8)
#define Ty3(j)    mean2d(j,ig(10)+9)
#define Ty1_y(j)  mean2d(j,ig(10)+10)
#define Ty2_y(j)  mean2d(j,ig(10)+11)
#define Ty3_y(j)  mean2d(j,ig(10)+12)
#define Gkin(j)   mean2d(j,ig(10)+13)
#define Dkin(j)   mean2d(j,ig(10)+14)
#define Phi(j)    mean2d(j,ig(10)+15)
#define ugradp(j) mean2d(j,ig(10)+16)
    sg(ng) = 17

    groupname(ng) = 'TkeBudget'
    varname(ng) = 'Tke_t Buo Con Prd Eps Pi Trp Trp1 Trp2 Trp3 Trp1_y Trp2_y Trp3_y G D Phi UgradP'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rU3(j)   mean2d(j,ig(11)  )
#define rU4(j)   mean2d(j,ig(11)+1)
#define rV3(j)   mean2d(j,ig(11)+2)
#define rV4(j)   mean2d(j,ig(11)+3)
#define rW3(j)   mean2d(j,ig(11)+4)
#define rW4(j)   mean2d(j,ig(11)+5)
    sg(ng) = 6

    groupname(ng) = 'HigherOrder'
    varname(ng) = 'rU3 rU4 rV3 rV4 rW3 rW4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define U_y1(j)  mean2d(j,ig(12)  )
#define V_y1(j)  mean2d(j,ig(12)+1)
#define W_y1(j)  mean2d(j,ig(12)+2)
#define U_ii2(j) mean2d(j,ig(12)+3)
#define U_x2(j)  mean2d(j,ig(12)+4)
#define U_y2(j)  mean2d(j,ig(12)+5)
#define U_z2(j)  mean2d(j,ig(12)+6)
#define V_x2(j)  mean2d(j,ig(12)+7)
#define V_y2(j)  mean2d(j,ig(12)+8)
#define V_z2(j)  mean2d(j,ig(12)+9)
#define W_x2(j)  mean2d(j,ig(12)+10)
#define W_y2(j)  mean2d(j,ig(12)+11)
#define W_z2(j)  mean2d(j,ig(12)+12)
#define U_x3(j)  mean2d(j,ig(12)+13)
#define U_y3(j)  mean2d(j,ig(12)+14)
#define U_z3(j)  mean2d(j,ig(12)+15)
#define V_x3(j)  mean2d(j,ig(12)+16)
#define V_y3(j)  mean2d(j,ig(12)+17)
#define V_z3(j)  mean2d(j,ig(12)+18)
#define W_x3(j)  mean2d(j,ig(12)+19)
#define W_y3(j)  mean2d(j,ig(12)+20)
#define W_z3(j)  mean2d(j,ig(12)+21)
#define U_x4(j)  mean2d(j,ig(12)+22)
#define U_y4(j)  mean2d(j,ig(12)+23)
#define U_z4(j)  mean2d(j,ig(12)+24)
#define V_x4(j)  mean2d(j,ig(12)+25)
#define V_y4(j)  mean2d(j,ig(12)+26)
#define V_z4(j)  mean2d(j,ig(12)+27)
#define W_x4(j)  mean2d(j,ig(12)+28)
#define W_y4(j)  mean2d(j,ig(12)+29)
#define W_z4(j)  mean2d(j,ig(12)+30)
    sg(ng) = 31

    groupname(ng) = 'DerivativeFluctuations'
    varname(ng) = 'U_y1 V_y1 W_y1 U_ii2 ' &
                  //'U_x2 U_y2 U_z2 V_x2 V_y2 V_z2 W_x2 W_y2 W_z2 ' &
                  //'U_x3 U_y3 U_z3 V_x3 V_y3 V_z3 W_x3 W_y3 W_z3 ' &
                  //'U_x4 U_y4 U_z4 V_x4 V_y4 V_z4 W_x4 W_y4 W_z4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rGamma(j)  mean2d(j,ig(13)  )
#define c2(j)      mean2d(j,ig(13)+1)
#define rho_ac(j)  mean2d(j,ig(13)+2)
#define rho_en(j)  mean2d(j,ig(13)+3)
#define T_ac(j)    mean2d(j,ig(13)+4)
#define T_en(j)    mean2d(j,ig(13)+5)
#define M_t(j)     mean2d(j,ig(13)+6)
#define rRP(j)     mean2d(j,ig(13)+7)
#define rRT(j)     mean2d(j,ig(13)+8)
    sg(ng) = 9

    groupname(ng) = 'Acoustics'
    varname(ng) = 'gamma C2 Rho_ac Rho_en T_ac T_en M_t rRP rRT'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rR2_flux_x(j) mean2d(j,ig(14)  )
#define rR2_flux_y(j) mean2d(j,ig(14)+1)
#define rR2_flux_z(j) mean2d(j,ig(14)+2)
#define rR2_dil1(j)   mean2d(j,ig(14)+3)
#define rR2_dil2(j)   mean2d(j,ig(14)+4)
#define rR2_trp(j)    mean2d(j,ig(14)+5)
#define rR2_prod(j)   mean2d(j,ig(14)+6)
#define rR2_conv(j)   mean2d(j,ig(14)+7)
    sg(ng) = 8

    groupname(ng) = 'RhoBudget'
    varname(ng) = 'RhoFluxX RhoFluxY RhoFluxZ RhoDil1 RhoDil2 RhoTrp RhoProd RhoConv'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Pot(j)       mean2d(j,ig(15)  )
#define rref(j)      mean2d(j,ig(15)+1)
#define tref(j)      mean2d(j,ig(15)+2)
#define bfreq_fr(j)  mean2d(j,ig(15)+3)
#define bfreq_eq(j)  mean2d(j,ig(15)+4)
#define lapse_fr(j)  mean2d(j,ig(15)+5)
#define lapse_eq(j)  mean2d(j,ig(15)+6)
#define potem_fr(j)  mean2d(j,ig(15)+7)
#define potem_eq(j)  mean2d(j,ig(15)+8)
#define psat(j)      mean2d(j,ig(15)+9)
#define pref(j)      mean2d(j,ig(15)+10)
#define relhum(j)    mean2d(j,ig(15)+11)
#define dewpoint(j)  mean2d(j,ig(15)+12)
#define lapse_dew(j) mean2d(j,ig(15)+13)
    sg(ng) = 14

    groupname(ng) = 'Stratification'
    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        varname(ng) = 'Pot rRref rTref BuoyFreq_fr BuoyFreq_eq LapseRate_fr LapseRate_eq ' &
                      //'PotTemp PotTemp_v SaturationPressure rPref RelativeHumidity Dewpoint LapseRate_dew'
    else
        varname(ng) = 'Pot rRref rTref BuoyFreq_fr BuoyFreq_eq LapseRate_fr LapseRate_eq ' &
                      //'PotTemp_fr PotTemp_eq SaturationPressure rPref RelativeHumidity Dewpoint LapseRate_dew'
    end if

    ! -----------------------------------------------------------------------
    ! Auxiliary variables depending on y and t; this last group is not written
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rUf(j)    mean2d(j,ig(16))
#define rVf(j)    mean2d(j,ig(16)+1)
#define rWf(j)    mean2d(j,ig(16)+2)

#define rU_y(j)   mean2d(j,ig(16)+3)
#define rV_y(j)   mean2d(j,ig(16)+4)
#define rW_y(j)   mean2d(j,ig(16)+5)
#define fU_y(j)   mean2d(j,ig(16)+6)
#define fV_y(j)   mean2d(j,ig(16)+7)
#define fW_y(j)   mean2d(j,ig(16)+8)
#define rP_y(j)   mean2d(j,ig(16)+9)
#define rR_y(j)   mean2d(j,ig(16)+10)
#define rT_y(j)   mean2d(j,ig(16)+11)
#define rB_y(j)   mean2d(j,ig(16)+12)

#define Rxx_y(j)  mean2d(j,ig(16)+13)
#define Ryy_y(j)  mean2d(j,ig(16)+14)
#define Rzz_y(j)  mean2d(j,ig(16)+15)
#define Rxy_y(j)  mean2d(j,ig(16)+16)
#define Rxz_y(j)  mean2d(j,ig(16)+17)
#define Ryz_y(j)  mean2d(j,ig(16)+18)
#define rR2_y(j)  mean2d(j,ig(16)+19)

#define Tau_yy(j)   mean2d(j,ig(16)+20)
#define Tau_xy(j)   mean2d(j,ig(16)+21)
#define Tau_yz(j)   mean2d(j,ig(16)+22)
#define Tau_xy_y(j) mean2d(j,ig(16)+23)
#define Tau_yy_y(j) mean2d(j,ig(16)+24)
#define Tau_yz_y(j) mean2d(j,ig(16)+25)
    sg(ng) = 26

    ! -----------------------------------------------------------------------
    nv = ig(ng) + sg(ng) - 1
    if (MAX_AVG_TEMPORAL < nv) then
        call TLab_Write_ASCII(efile, 'AVG_FLOW_XZ. Not enough space in local arrays.')
        call TLab_Stop(DNS_ERROR_AVGTMP)
    end if
    mean2d(:, 1:nv) = 0.0_wp

    ng = ng - 1
    nv = ig(ng) + sg(ng) - 1 ! the last group is not written out

    ! #######################################################################
    write (line1, *) itime; line1 = 'Calculating flow statistics at It'//trim(adjustl(line1))//'...'
    call TLab_Write_ASCII(lfile, line1)

    ! #######################################################################
    ! Preliminary for IBM usage
    ! #######################################################################
    ! Asign gammas for conditional averages (c.f. Pope, p.170 [5.305])
    if (imode_ibm == 1) then
        ep_0(:) = gamma_0; ep_1(:) = gamma_1
    end if

    ! ###################################################################
    ! Averages (do not overwrite dudz; it contains p for incompressible case)
    ! ###################################################################
#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 2')
#endif

    ! Velocity
    call AVG_IK_V(imax, jmax, kmax, u, rU(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, v, rV(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, w, rW(1), wrk1d)

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rU(1), rU_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rV(1), rV_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rW(1), rW_y(1))

    U_y1(:) = rU_y(:)
    V_y1(:) = rV_y(:)
    W_y1(:) = rW_y(:)

    ! Density and Favre avrages
    if (nse_eqns == DNS_EQNS_INCOMPRESSIBLE) then
        rR(:) = 1.0_wp ! I divide below by density; rbackground(:)

        fU(:) = rU(:); fV(:) = rV(:); fW(:) = rW(:)

    else if (nse_eqns == DNS_EQNS_ANELASTIC) then
        call Thermo_Anelastic_DENSITY(imax, jmax, kmax, s, dwdx, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dwdx, rR(1), wrk1d)

        fU(:) = rU(:); fV(:) = rV(:); fW(:) = rW(:)

    else
        call AVG_IK_V(imax, jmax, kmax, rho, rR(1), wrk1d)

        dwdx = rho*u
        dwdy = rho*v
        dwdz = rho*w
        call AVG_IK_V(imax, jmax, kmax, dwdx, fU(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dwdy, fV(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dwdz, fW(1), wrk1d)
        fU(:) = fU(:)/rR(:)
        fV(:) = fV(:)/rR(:)
        fW(:) = fW(:)/rR(:)

    end if

    rUf(:) = rU(:) - fU(:)
    rVf(:) = rV(:) - fV(:)
    rWf(:) = rW(:) - fW(:)

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rR(1), rR_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fU(1), fU_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fV(1), fV_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fW(1), fW_y(1))

    ! Pressure
    call AVG_IK_V(imax, jmax, kmax, p_loc, rP(1), wrk1d)
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rP(1), rP_y(1))

    ! #######################################################################
    ! Main covariances (do not overwrite dudz; it contains p for incompressible case)
    ! #######################################################################

    do j = 1, jmax
        dwdx(:, j, :) = u(:, j, :) - fU(j)
        dwdy(:, j, :) = v(:, j, :) - fV(j)
        dwdz(:, j, :) = w(:, j, :) - fW(j)
    end do

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        dvdx = dwdx*dwdx
        dvdy = dwdy*dwdy
        dvdz = dwdz*dwdz
    else
        dvdx = dwdx*dwdx*rho
        dvdy = dwdy*dwdy*rho
        dvdz = dwdz*dwdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, dvdx, Rxx(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdy, Ryy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdz, Rzz(1), wrk1d)
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) then
        Rxx(:) = Rxx(:)/rR(:)
        Ryy(:) = Ryy(:)/rR(:)
        Rzz(:) = Rzz(:)/rR(:)
    end if

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        dvdx = dwdx*dwdy
        dvdy = dwdx*dwdz
        dvdz = dwdy*dwdz
    else
        dvdx = dwdx*dwdy*rho
        dvdy = dwdx*dwdz*rho
        dvdz = dwdy*dwdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, dvdx, Rxy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdy, Rxz(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdz, Ryz(1), wrk1d)
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) then
        Rxy(:) = Rxy(:)/rR(:)
        Rxz(:) = Rxz(:)/rR(:)
        Ryz(:) = Ryz(:)/rR(:)
    end if

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rxx(1), Rxx_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Ryy(1), Ryy_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rzz(1), Rzz_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rxy(1), Rxy_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rxz(1), Rxz_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Ryz(1), Ryz_y(1))

    ! Density
    if (.not. (nse_eqns == DNS_EQNS_INCOMPRESSIBLE)) then
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_DENSITY(imax, jmax, kmax, s, dudx, p_wrk3d)
            do j = 1, jmax
                dudx(:, j, :) = dudx(:, j, :) - rR(j)
            end do
        else
            do j = 1, jmax
                dudx(:, j, :) = rho(:, j, :) - rR(j)
            end do
        end if
        dvdx = dudx*dudx
        call AVG_IK_V(imax, jmax, kmax, dvdx, rR2(1), wrk1d)

        call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rR2(1), rR2_y(1))
!$omp
        ! Density Fluctuations Budget

        do j = 1, jmax
            dvdx(:, j, :) = u(:, j, :) - rU(j)
            dvdy(:, j, :) = v(:, j, :) - rV(j)
            dvdz(:, j, :) = w(:, j, :) - rW(j)
        end do

        dvdx = dvdx*dudx
        dvdy = dvdy*dudx
        dvdz = dvdz*dudx
        call AVG_IK_V(imax, jmax, kmax, dvdx, rR2_flux_x(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdy, rR2_flux_y(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, rR2_flux_z(1), wrk1d)
        dvdy = dvdy*dudx
        call AVG_IK_V(imax, jmax, kmax, dvdy, rR2_trp(1), wrk1d)

    end if

    ! higher-order moments
    p_wrk3d = dwdx*dwdx*dwdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rU3(1), wrk1d)
    p_wrk3d = dwdx*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rU4(1), wrk1d)

    p_wrk3d = dwdy*dwdy*dwdy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rV3(1), wrk1d)
    p_wrk3d = dwdy*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rV4(1), wrk1d)

    p_wrk3d = dwdz*dwdz*dwdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rW3(1), wrk1d)
    p_wrk3d = dwdz*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rW4(1), wrk1d)

    ! Triple-velocity correlations
    dvdx = dwdx*dwdx*dwdy
    dvdy = dwdy*dwdy*dwdy
    dvdz = dwdz*dwdz*dwdy
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) then
        dvdx = dvdx*rho
        dvdy = dvdy*rho
        dvdz = dvdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, dvdx, Txxy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdy, Tyyy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdz, Tzzy(1), wrk1d)
    Ty1(:) = (Txxy(:) + Tyyy(:) + Tzzy(:))*0.5_wp

    dvdx = dwdx*dwdy*dwdy
    dvdy = dwdx*dwdy*dwdz
    dvdz = dwdy*dwdy*dwdz
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) then
        dvdx = dvdx*rho
        dvdy = dvdy*rho
        dvdz = dvdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, dvdx, Txyy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdy, Txzy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdz, Tyzy(1), wrk1d)

    ! Pressure
    do j = 1, jmax
        dvdz(:, j, :) = p_loc(:, j, :) - rP(j)
    end do
    p_wrk3d = dvdz*dvdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rP2(1), wrk1d)

    ! Pressure-velocity correlation in TKE transport terms
    dwdx = dwdx*dvdz
    dwdy = dwdy*dvdz
    dwdz = dwdz*dvdz
    call AVG_IK_V(imax, jmax, kmax, dwdx, wrk1d(1, 2), wrk1d)
    Txyy(:) = Txyy(:) + wrk1d(1:jmax, 2)
    call AVG_IK_V(imax, jmax, kmax, dwdy, Ty2(1), wrk1d)
    Tyyy(:) = Tyyy(:) + Ty2(:)*2.0_wp
    call AVG_IK_V(imax, jmax, kmax, dwdz, wrk1d(1, 2), wrk1d)
    Tyzy(:) = Tyzy(:) + wrk1d(1:jmax, 2)

    ! ###################################################################
    ! Pressure; array dudz containing p is used only up to this section
    !
    ! dudx = du/dx
    ! dudy = du/dy
    ! dudz = p
    ! dvdx = dv/dx
    ! dvdy = dv/dy
    ! dvdz = p_prime
    ! dwdx =       ; dp/dx
    ! dwdy =       ; dp/dy
    ! dwdz = dw/dz ; dp/dz
    ! ###################################################################
    ! Pressure convection term
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), p_loc, dwdx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), p_loc, dwdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), p_loc, dwdz)
    p_wrk3d = u*dwdx + v*dwdy + w*dwdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, ugradp(1), wrk1d)

    ! Pressure Strain Terms
    ! 9 derivatives are here recomputed; ok, this routine is not called that often
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, dudx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, dvdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, dwdz)
    dudx = dvdz*dudx ! dvdz contains the pressure fluctuation
    dvdy = dvdz*dvdy ! no need to substract rV_y
    dwdz = dvdz*dwdz
    call AVG_IK_V(imax, jmax, kmax, dudx, PIxx(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dvdy, PIyy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dwdz, PIzz(1), wrk1d)
    PIxx(:) = PIxx(:)*2.0_wp
    PIyy(:) = PIyy(:)*2.0_wp
    PIzz(:) = PIzz(:)*2.0_wp

    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, dudy)
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, dvdx)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, dwdz) !dudz not free
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, dwdx)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, dudx) !dvdz not free
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, dvdy)
    dudy = dvdz*(dudy + dvdx) ! no need to substract rU_y
    dwdz = dvdz*(dwdz + dwdx)
    dudx = dvdz*(dudx + dvdy) ! no need to substract rW_y
    call AVG_IK_V(imax, jmax, kmax, dudy, PIxy(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dwdz, PIxz(1), wrk1d)
    call AVG_IK_V(imax, jmax, kmax, dudx, PIyz(1), wrk1d)

    ! ###################################################################
    ! Thermodynamic variables
    !
    ! dudx = gamma
    ! dudy = dsdy
    ! dudz = dTdy
    ! dvdx = cp
    ! dvdy = drdy
    ! dvdz = psat
    ! dwdx = T
    ! dwdy = dpdy
    ! dwdz =
    ! ###################################################################
#define T_LOC(i,j,k)     dwdx(i,j,k)

    if (nse_eqns == DNS_EQNS_INCOMPRESSIBLE) then
        ! rT(:) = tbackground(:)

    else if (nse_eqns == DNS_EQNS_ANELASTIC) then
        call Thermo_Anelastic_TEMPERATURE(imax, jmax, kmax, s, T_LOC(1, 1, 1))
        call AVG_IK_V(imax, jmax, kmax, T_LOC(1, 1, 1), rT(1), wrk1d)

        do j = 1, jmax
            dvdx(:, j, :) = (T_LOC(:, j, :) - rT(j))**2
        end do
        call AVG_IK_V(imax, jmax, kmax, dvdx, rT2(1), wrk1d)

        call Thermo_Psat_Polynomial(imax*jmax*kmax, T_LOC(1, 1, 1), dvdz)
        call AVG_IK_V(imax, jmax, kmax, dvdz, psat(1), wrk1d)

        call Thermo_Anelastic_RELATIVEHUMIDITY(imax, jmax, kmax, s, dvdz, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, relhum(1), wrk1d)

        call Thermo_Anelastic_THETA(imax, jmax, kmax, s, dvdz, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, potem_fr(1), wrk1d)
        call Thermo_Anelastic_THETA_V(imax, jmax, kmax, s, dvdz, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, potem_eq(1), wrk1d)

        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), T_LOC(1, 1, 1), dudz)

        call Thermo_Anelastic_LAPSE_EQU(imax, jmax, kmax, s, dwdz, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, dwdz, lapse_eq(1), wrk1d)
        !frequency(:) = (lapse(:) + dTdy(:))/locT(:)
        p_wrk3d = (dwdz + dudz)/dwdx
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, bfreq_eq(1), wrk1d)
        bfreq_eq(:) = bfreq_eq(:)*buoyancy%vector(2)

        call Thermo_Anelastic_LAPSE_FR(imax, jmax, kmax, s, dwdz)
        call AVG_IK_V(imax, jmax, kmax, dwdz, lapse_fr(1), wrk1d)
        !frequency(:) = (lapse(:) + dTdy(:))/locT(:)
        p_wrk3d = (dwdz + dudz)/dwdx
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, bfreq_fr(1), wrk1d)
        bfreq_fr(:) = bfreq_fr(:)*buoyancy%vector(2)

        call Thermo_Anelastic_VAPOR_PRESSURE(imax, jmax, kmax, s, dudz)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), dudz, dudy)
        ! dwdz should contains lapse_fr, since lapse_dew = lapse_fr when saturated
        call Thermo_Anelastic_DEWPOINT(imax, jmax, kmax, s, dudy, p_wrk3d, dwdz)
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, dewpoint(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dwdz, lapse_dew(1), wrk1d)

#undef T_LOC

    else
#define GAMMA_LOC(i,j,k) dudx(i,j,k)
#define S_LOC(i,j,k)     dwdz(i,j,k)

        ! -------------------------------------------------------------------
        ! Main fields
        ! -------------------------------------------------------------------
        ! call THERMO_CALORIC_TEMPERATURE(imax*jmax*kmax, s, e, rho, T, p_wrk3d)
        call THERMO_GAMMA(imax*jmax*kmax, s, T, GAMMA_LOC(:, :, :))
        call THERMO_ENTROPY(imax*jmax*kmax, s, T, p_loc, S_LOC(1, 1, 1))

        call AVG_IK_V(imax, jmax, kmax, T, rT(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, e, re(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, S_LOC(1, 1, 1), rs(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, GAMMA_LOC(1, 1, 1), rGamma(1), wrk1d)

        ! Means
        dudy = rho*e
        dudz = e + CRATIO_INV*p_loc/rho             ! enthalpy
        dvdx = rho*dudz
        dvdy = rho*S_LOC(:, :, :)
        dvdz = rho*T
        p_wrk3d = GAMMA_LOC(:, :, :)*p_loc/rho      ! speed of sound
        call AVG_IK_V(imax, jmax, kmax, dudy, fe(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dudz, rh(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdx, fh(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdy, fs(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, fT(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, c2(1), wrk1d)
        fe(:) = fe(:)/rR(:)
        fh(:) = fh(:)/rR(:)
        fs(:) = fs(:)/rR(:)
        fT(:) = fT(:)/rR(:)

        call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rT(1), rT_y(1))

        ! Turbulent Mach number
        M_t(:) = sqrt((Rxx(:) + Ryy(:) + Rzz(:))/c2(:))

        ! Covariances
        do j = 1, jmax
            dudy(:, j, :) = (S_LOC(:, j, :) - rs(j))**2
            dudz(:, j, :) = rho(:, j, :)*(S_LOC(:, j, :) - fs(j))**2
            dvdx(:, j, :) = (T(:, j, :) - rT(j))**2
            dvdy(:, j, :) = rho(:, j, :)*(T(:, j, :) - fT(j))**2
            dvdz(:, j, :) = (rho(:, j, :) - rR(j))*(T(:, j, :) - fT(j))
            p_wrk3d(:, j, :) = (rho(:, j, :) - rR(j))*(p_loc(:, j, :) - rP(j))
        end do
        call AVG_IK_V(imax, jmax, kmax, dudy, rs2(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dudz, fs2(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdx, rT2(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdy, fT2(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, rRT(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rRP(1), wrk1d)
        fs2(:) = fs2(:)/rR(:)
        fT2(:) = fT2(:)/rR(:)

    do j = 1, jmax
        dudy(:, j, :) = (e(:, j, :) - re(j))**2
        dudz(:, j, :) = rho(:, j, :)*(e(:, j, :) - fe(j))**2
    end do
        call AVG_IK_V(imax, jmax, kmax, dudy, re2(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dudz, fe2(1), wrk1d)
        fe2(:) = fe2(:)/rR(:)

        p_wrk3d = e + CRATIO_INV*p_loc/rho
    do j = 1, jmax
        dudy(:, j, :) = (p_wrk3d(:, j, :) - rh(j))**2
        dudz(:, j, :) = rho(:, j, :)*(p_wrk3d(:, j, :) - fh(j))**2
    end do
        call AVG_IK_V(imax, jmax, kmax, dudy, rh2(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dudz, fh2(1), wrk1d)
        fh2(:) = fh2(:)/rR(:)

        ! Acoustic and entropic density and temperature fluctuations
    do j = 1, jmax
        dudy(:, j, :) = p_loc(:, j, :) - rP(j)                 ! pprime
        dudz(:, j, :) = dudy(:, j, :)/c2(j)               ! rho_ac
        dvdx(:, j, :) = rho(:, j, :) - rR(j) - dudz(:, j, :) ! rho_en = rprime - rho_ac
        dvdy(:, j, :) = (dudy(:, j, :)/rP(j) - dudz(:, j, :)/rR(j))*fT(j) ! T_ac
        dvdz(:, j, :) = T(:, j, :) - fT(j) - dvdy(:, j, :)               ! T_en = Tprime - T_ac
    end do
        dudz = dudz*dudz
        dvdx = dvdx*dvdx
        dvdy = dvdy*dvdy
        dvdz = dvdz*dvdz
        call AVG_IK_V(imax, jmax, kmax, dudz, rho_ac(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdx, rho_en(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdy, T_ac(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, T_en(1), wrk1d)

        ! -------------------------------------------------------------------
        ! Buoyancy frequency & saturation pressure
        ! -------------------------------------------------------------------
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), rho, dvdy)

        call Thermo_Psat_Polynomial(imax*jmax*kmax, T, dvdz)
        call THERMO_CP(imax*jmax*kmax, s, GAMMA_LOC(:, :, :), dvdx)

    do j = 1, jmax
        dudy(:, j, :) = dwdy(:, j, :)/p_loc(:, j, :)/GAMMA_LOC(:, j, :) - dvdy(:, j, :)/rho(:, j, :)
        dvdx(:, j, :) = 1.0_wp/dvdx(:, j, :)
        dvdy(:, j, :) = T(:, j, :)*((p_loc(:, j, :)/PREF_1000)**(1.0_wp/GAMMA_LOC(:, j, :) - 1.0_wp))
    end do
        call AVG_IK_V(imax, jmax, kmax, dudy, bfreq_fr(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdx, lapse_fr(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdy, potem_fr(1), wrk1d)
        call AVG_IK_V(imax, jmax, kmax, dvdz, psat(1), wrk1d)
        bfreq_fr(:) = -bfreq_fr(:)*buoyancy%vector(2)
        lapse_fr(:) = -lapse_fr(:)*buoyancy%vector(2)*CRATIO_INV

        if (imixture == MIXT_TYPE_AIRWATER) then
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), T, dudz)
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(1, 1, 1, 2), dudy)

            ! To be done
            ! call THERMO_AIRWATER_LAPSE_EQU(s, T, p_loc, dudz, dudy, dvdx, dvdy)!, dwdy, dwdz, p_wrk3d)
            ! call AVG_IK_V(imax, jmax, kmax, dvdy, bfreq_eq(1), wrk1d)
            ! call AVG_IK_V(imax, jmax, kmax, dvdx, lapse_eq(1), wrk1d)
            ! bfreq_eq(:) = -bfreq_eq(:)*buoyancy%vector(2)
            ! lapse_eq(:) = -lapse_eq(:)*buoyancy%vector(2)*CRATIO_INV

            call THERMO_AIRWATER_THETA_EQ(s, T, p_loc, dvdx, dvdy, dwdy)
            call AVG_IK_V(imax, jmax, kmax, dvdx, potem_eq(1), wrk1d)

        end if

#undef S_LOC
#undef GAMMA_LOC
    end if

    select case (imode_thermo)
    case (THERMO_TYPE_ANELASTIC)
        pref(:) = pbackground(:)
        tref(:) = tbackground(:)
        rref(:) = rbackground(:)
    case (THERMO_TYPE_COMPRESSIBLE)
        pref(:) = rP(:)
        tref(:) = rT(:)
        rref(:) = rR(:)
    end select

    ! ###################################################################
    ! Potential energy
    !
    ! dudx = buoyancy
    !
    ! ###################################################################
    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then

        if (buoyancy%type /= EQNS_BOD_NONE) then
            if (buoyancy%type == EQNS_BOD_EXPLICIT) then
                call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, dudx)
            else
                call Gravity_Buoyancy(buoyancy, imax, jmax, kmax, s, dudx, bbackground)
            end if

            call AVG_IK_V(imax, jmax, kmax, dudx, rB(1), wrk1d)
            do j = 1, jmax
                dvdx(:, j, :) = (u(:, j, :) - rU(j))*(dudx(:, j, :) - rB(j))
                dvdy(:, j, :) = (v(:, j, :) - rV(j))*(dudx(:, j, :) - rB(j))
                dvdz(:, j, :) = (w(:, j, :) - rW(j))*(dudx(:, j, :) - rB(j))
            end do
            call AVG_IK_V(imax, jmax, kmax, dvdx, Bxx(1), wrk1d)
            call AVG_IK_V(imax, jmax, kmax, dvdy, Byy(1), wrk1d)
            call AVG_IK_V(imax, jmax, kmax, dvdz, Bzz(1), wrk1d)
            Bxy(:) = Bxx(:)*buoyancy%vector(2) + Byy(:)*buoyancy%vector(1) ! buoyancy%vector includes the Froude
            Bxz(:) = Bxx(:)*buoyancy%vector(3) + Bzz(:)*buoyancy%vector(1)
            Byz(:) = Byy(:)*buoyancy%vector(3) + Bzz(:)*buoyancy%vector(2)

            Bxx(:) = 2.0_wp*Bxx(:)*buoyancy%vector(1)
            Byy(:) = 2.0_wp*Byy(:)*buoyancy%vector(2)
            Bzz(:) = 2.0_wp*Bzz(:)*buoyancy%vector(3)

            dummy = 1.0_wp/froude
            rB(:) = rB(:)*dummy

            !        pmod(:) =-rP_y(:) + SIGN(rB(:),buoyancy%vector(2))

            call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rB(1), rB_y(1))

        end if

    else ! Compressible case is not yet finished
        Bxx(:) = -rR(:)*rUf(:)*buoyancy%vector(1)
        Byy(:) = -rR(:)*rVf(:)*buoyancy%vector(2)
        Bzz(:) = -rR(:)*rWf(:)*buoyancy%vector(3)

        !     pmod(:) =-rP_y(:) +buoyancy%vector(2) *rR(:)

    end if

    ! ###################################################################
    ! # Array storage of velocity gradient tensor
    ! #
    ! # dudx = d U / d x
    ! # dudy = d U / d y
    ! # dudz = d U / d z
    ! # dvdx = d V / d x
    ! # dvdy = d V / d y
    ! # dvdz = d V / d z
    ! # dwdx = d W / d x
    ! # dwdy = d W / d y
    ! # dwdz = d W / d z
    ! ###################################################################
#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 3')
#endif

    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, dudx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, dudy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, dudz)

    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, dvdx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, dvdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, dvdz)

    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, dwdx)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, dwdy)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, dwdz)

    ! ###################################################################
    ! Vorticity
    ! ###################################################################
    p_wrk3d = dwdy - dvdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vortx(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = (p_wrk3d(:, j, :) - vortx(j))**2
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vortx2(1), wrk1d)

    p_wrk3d = dudz - dwdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vorty(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = (p_wrk3d(:, j, :) - vorty(j))**2
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vorty2(1), wrk1d)

    p_wrk3d = dvdx - dudy
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vortz(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = (p_wrk3d(:, j, :) - vortz(j))**2
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, vortz2(1), wrk1d)

    ! ###################################################################
    ! Derivatives Fluctuations
    ! ###################################################################
#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 11')
#endif

    ! -------------------------------------------------------------------
    ! Longitudinal terms
    p_wrk3d = dudx*dudx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_x2(1), wrk1d)
    p_wrk3d = p_wrk3d*dudx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_x3(1), wrk1d)
    p_wrk3d = p_wrk3d*dudx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_x4(1), wrk1d)

    do j = 1, jmax
        p_wrk3d(:, j, :) = (dvdy(:, j, :) - rV_y(j))*(dvdy(:, j, :) - rV_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_y2(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dvdy(:, j, :) - rV_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_y3(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dvdy(:, j, :) - rV_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_y4(1), wrk1d)

    p_wrk3d = dwdz*dwdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_z2(1), wrk1d)
    p_wrk3d = p_wrk3d*dwdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_z3(1), wrk1d)
    p_wrk3d = p_wrk3d*dwdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_z4(1), wrk1d)

    ! -------------------------------------------------------------------
    ! Lateral terms U
    do j = 1, jmax
        p_wrk3d(:, j, :) = (dudy(:, j, :) - rU_y(j))*(dudy(:, j, :) - rU_y(j))
    end do
    
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_y2(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dudy(:, j, :) - rU_y(j))
    end do
    
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_y3(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dudy(:, j, :) - rU_y(j))
    end do
    
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_y4(1), wrk1d)

    p_wrk3d = dudz*dudz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_z2(1), wrk1d)
    p_wrk3d = p_wrk3d*dudz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_z3(1), wrk1d)
    p_wrk3d = p_wrk3d*dudz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_z4(1), wrk1d)

    ! -------------------------------------------------------------------
    ! Lateral terms V
    p_wrk3d = dvdx*dvdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_x2(1), wrk1d)
    p_wrk3d = p_wrk3d*dvdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_x3(1), wrk1d)
    p_wrk3d = p_wrk3d*dvdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_x4(1), wrk1d)

    p_wrk3d = dvdz*dvdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_z2(1), wrk1d)
    p_wrk3d = p_wrk3d*dvdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_z3(1), wrk1d)
    p_wrk3d = p_wrk3d*dvdz
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, V_z4(1), wrk1d)

    ! -------------------------------------------------------------------
    ! Lateral terms W
    p_wrk3d = dwdx*dwdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_x2(1), wrk1d)
    p_wrk3d = p_wrk3d*dwdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_x3(1), wrk1d)
    p_wrk3d = p_wrk3d*dwdx
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_x4(1), wrk1d)

    do j = 1, jmax
        p_wrk3d(:, j, :) = (dwdy(:, j, :) - rW_y(j))*(dwdy(:, j, :) - rW_y(j))
    end do
    
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_y2(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dwdy(:, j, :) - rW_y(j))
    end do
    
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_y3(1), wrk1d)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dwdy(:, j, :) - rW_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, W_y4(1), wrk1d)

    ! -------------------------------------------------------------------
    ! Dilatation fluctuation
    p_wrk3d = dudx + dvdy + dwdz
    do j = 1, jmax
        p_wrk3d(:, j, :) = (p_wrk3d(:, j, :) - rV_y(j))*(p_wrk3d(:, j, :) - rV_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, U_ii2(1), wrk1d)

    ! ###################################################################
    ! Density Fluctuations Budget
    ! ###################################################################
    if (nse_eqns == DNS_EQNS_INTERNAL .or. nse_eqns == DNS_EQNS_TOTAL) then
        p_wrk3d = dudx + dvdy + dwdz
        do j = 1, jmax
            p_wrk3d(:, j, :) = (p_wrk3d(:, j, :) - rV_y(j))*(rho(:, j, :) - rR(j))
        end do
        
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rR2_dil1(1), wrk1d)

        do j = 1, jmax
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(rho(:, j, :) - rR(j))
        end do
        call AVG_IK_V(imax, jmax, kmax, p_wrk3d, rR2_dil2(1), wrk1d)

    end if

    ! ##################################################################
    ! Mean viscous dissipation rate
    ! ##################################################################
    p_wrk3d = dudx**2 + dvdy**2 + dwdz**2 + 0.5_wp*((dudy + dvdx)**2 + (dudz + dwdx)**2 + (dvdz + dwdy)**2) &
              - (dudx + dvdy + dwdz)**2/3.0_wp

    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Phi(1), wrk1d)
    Phi(:) = Phi(:)*visc*2.0_wp

    ! ###################################################################
    ! Dissipation Terms; final operation after viscous terms below
    ! ###################################################################
    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dudx*2.0_wp - p_wrk3d)*dudx + (dudy + dvdx)*dudy + (dudz + dwdx)*dudz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Exx(1), wrk1d)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dvdy*2.0_wp - p_wrk3d)*dvdy + (dudy + dvdx)*dvdx + (dvdz + dwdy)*dvdz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Eyy(1), wrk1d)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dwdz*2.0_wp - p_wrk3d)*dwdz + (dwdy + dvdz)*dwdy + (dwdx + dudz)*dwdx
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Ezz(1), wrk1d)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dudx*2.0_wp - p_wrk3d)*dvdx + (dudy + dvdx)*dvdy + (dudz + dwdx)*dvdz &
              + (dvdy*2.0_wp - p_wrk3d)*dudy + (dudy + dvdx)*dudx + (dvdz + dwdy)*dudz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Exy(1), wrk1d)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dudx*2.0_wp - p_wrk3d)*dwdx + (dudy + dvdx)*dwdy + (dudz + dwdx)*dwdz &
              + (dwdz*2.0_wp - p_wrk3d)*dudz + (dudz + dwdx)*dudx + (dvdz + dwdy)*dudy
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Exz(1), wrk1d)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dvdy*2.0_wp - p_wrk3d)*dwdy + (dudy + dvdx)*dwdx + (dvdz + dwdy)*dwdz &
              + (dwdz*2.0_wp - p_wrk3d)*dvdz + (dudz + dwdx)*dvdx + (dvdz + dwdy)*dvdy
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Eyz(1), wrk1d)

    ! ##################################################################
    ! Viscous shear-stress tensor
    ! ##################################################################
    p_wrk3d = dvdy*2.0_wp - dudx - dwdz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_yy(1), wrk1d)
    do j = 1, jmax ! fluctuation tau22'
        dvdy(:, j, :) = (p_wrk3d(:, j, :) - Tau_yy(j))*c23
    end do
    Tau_yy(:) = Tau_yy(:)*visc*c23

    p_wrk3d = dudy + dvdx
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_xy(1), wrk1d)
    do j = 1, jmax ! fluctuation tau12'
        dudy(:, j, :) = p_wrk3d(:, j, :) - Tau_xy(j)
    end do
    Tau_xy(:) = Tau_xy(:)*visc

    p_wrk3d = dvdz + dwdy
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, Tau_yz(1), wrk1d)
    do j = 1, jmax ! fluctuation tau23'
        dwdy(:, j, :) = p_wrk3d(:, j, :) - Tau_yz(j)
    end do
    Tau_yz(:) = Tau_yz(:)*visc

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_xy(1), Tau_xy_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_yy(1), Tau_yy_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_yz(1), Tau_yz_y(1))

    ! -------------------------------------------------------------------
    ! Contribution to turbulent transport terms
    do j = 1, jmax
        p_wrk3d(:, j, :) = dudy(:, j, :)*(u(:, j, :) - fU(j)) ! -2*u'*tau12'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, wrk1d(1, 2), wrk1d)
    Txxy(:) = Txxy(:) - wrk1d(1:jmax, 2)*visc*2.0_wp
    Ty3(:) = Ty3(:) - wrk1d(1:jmax, 2)*visc

    do j = 1, jmax
        p_wrk3d(:, j, :) = dvdy(:, j, :)*(v(:, j, :) - fV(j)) ! -2*v'*tau22'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, wrk1d(1, 2), wrk1d)
    Tyyy(:) = Tyyy(:) - wrk1d(1:jmax, 2)*visc*2.0_wp
    Ty3(:) = Ty3(:) - wrk1d(1:jmax, 2)*visc

    do j = 1, jmax
        p_wrk3d(:, j, :) = dwdy(:, j, :)*(w(:, j, :) - fW(j)) ! -2*w'*tau23'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, wrk1d(1, 2), wrk1d)
    Tzzy(:) = Tzzy(:) - wrk1d(1:jmax, 2)*visc*2.0_wp
    Ty3(:) = Ty3(:) - wrk1d(1:jmax, 2)*visc

    do j = 1, jmax
        p_wrk3d(:, j, :) = dvdy(:, j, :)*(u(:, j, :) - fU(j)) + dudy(:, j, :)*(v(:, j, :) - fV(j))! -u'*tau22' -v'*tau12'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, wrk1d(1, 2), wrk1d)
    Txyy(:) = Txyy(:) - wrk1d(1:jmax, 2)*visc
    do j = 1, jmax
        p_wrk3d(:, j, :) = dwdy(:, j, :)*(u(:, j, :) - fU(j)) + dudy(:, j, :)*(w(:, j, :) - fW(j))! -u'*tau23' -w'*tau12'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, wrk1d(1, 2), wrk1d)
    Txzy(:) = Txzy(:) - wrk1d(1:jmax, 2)*visc

    do j = 1, jmax
        p_wrk3d(:, j, :) = dwdy(:, j, :)*(v(:, j, :) - fV(j)) + dvdy(:, j, :)*(w(:, j, :) - fW(j))! -v'*tau23' -w'*tau22'
    end do
    call AVG_IK_V(imax, jmax, kmax, p_wrk3d, wrk1d(1, 2), wrk1d)
    Tyzy(:) = Tyzy(:) - wrk1d(1:jmax, 2)*visc

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Txxy(1), Txxy_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tyyy(1), Tyyy_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tzzy(1), Tzzy_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Txyy(1), Txyy_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Txzy(1), Txzy_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tyzy(1), Tyzy_y(1))

    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Ty1(1), Ty1_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Ty2(1), Ty2_y(1))
    call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Ty3(1), Ty3_y(1))

    ! -------------------------------------------------------------------
    ! Contribution to dissipation
    Exx(:) = (Exx(:)*visc - Tau_xy(:)*rU_y(:))*2.0_wp
    Eyy(:) = (Eyy(:)*visc - Tau_yy(:)*rV_y(:))*2.0_wp
    Ezz(:) = (Ezz(:)*visc - Tau_yz(:)*rW_y(:))*2.0_wp
    Exy(:) = Exy(:)*visc - Tau_xy(:)*rV_y(:) - Tau_yy(:)*rU_y(:)
    Exz(:) = Exz(:)*visc - Tau_xy(:)*rW_y(:) - Tau_yz(:)*rU_y(:)
    Eyz(:) = Eyz(:)*visc - Tau_yy(:)*rW_y(:) - Tau_yz(:)*rV_y(:)

    ! ###################################################################
    ! Complete budget equations
    ! ###################################################################
    ! Density fluctuations budget equation
    if (nse_eqns == DNS_EQNS_INTERNAL .or. nse_eqns == DNS_EQNS_TOTAL) then
        rR2_prod(:) = -2.0_wp*(rR2_flux_y(:)*rR_y(:) + rR2(:)*rV_y(:))
        rR2_conv(:) = -rV(:)*rR2_y(:)
        rR2_dil1(:) = 2.0_wp*rR(:)*rR2_dil1(:)
    end if

    ! Rij Convective Terms
    Cxx(:) = -fV(:)*Rxx_y(:)
    Cyy(:) = -fV(:)*Ryy_y(:)
    Czz(:) = -fV(:)*Rzz_y(:)
    Cxy(:) = -fV(:)*Rxy_y(:)
    Cxz(:) = -fV(:)*Rxz_y(:)
    Cyz(:) = -fV(:)*Ryz_y(:)

    ! Rij Production Terms
    Pxx(:) = -2.0_wp*Rxy(:)*fU_y(:)
    Pyy(:) = -2.0_wp*Ryy(:)*fV_y(:)
    Pzz(:) = -2.0_wp*Ryz(:)*fW_y(:)
    Pxy(:) = -(Rxy(:)*fV_y(:) + Ryy(:)*fU_y(:))
    Pxz(:) = -(Rxy(:)*fW_y(:) + Ryz(:)*fU_y(:))
    Pyz(:) = -(Ryy(:)*fW_y(:) + Ryz(:)*fV_y(:))

    ! Rij Pressure Variable-Density Terms
    Gxx(:) = 0.0_wp
    Gyy(:) = 2.0_wp*rVf(:)*rP_y(:)
    Gzz(:) = 0.0_wp
    Gxy(:) = rUf(:)*rP_y(:)
    Gxz(:) = 0.0_wp
    Gyz(:) = rWf(:)*rP_y(:)

    ! Rij Viscous Variable-Density Terms
    Dxx(:) = 2.0_wp*rUf(:)*Tau_xy_y(:)
    Dyy(:) = 2.0_wp*rVf(:)*Tau_yy_y(:)
    Dzz(:) = 2.0_wp*rWf(:)*Tau_yz_y(:)
    Dxy(:) = rUf(:)*Tau_yy_y(:) + rVf(:)*Tau_xy_y(:)
    Dxz(:) = rUf(:)*Tau_yz_y(:) + rWf(:)*Tau_xy_y(:)
    Dyz(:) = rVf(:)*Tau_yz_y(:) + rWf(:)*Tau_yy_y(:)

    ! Rij Coriolis Terms
    if (coriolis%active(1) .and. coriolis%active(3)) then ! contribution from angular velocity Oy
        dummy = coriolis%vector(2)
        Fxx(:) = dummy*2.0_wp*Rxz(:)
        Fyy(:) = 0.0_wp
        Fzz(:) = -dummy*2.0_wp*Rxz(:)
        Fxy(:) = dummy*Ryz(:)
        Fxz(:) = dummy*(Rzz(:) - Rxx(:))
        Fyz(:) = -dummy*Rxy(:)
    end if

    ! Rij Buoyancy Terms; Calculated in Section Potential Energy

    ! Rij Transient terms
    Rxx_t(:) = -Fxx(:) + Bxx(:) + Cxx(:) + Pxx(:) - Exx(:) + (PIxx(:) - Txxy_y(:) - Gxx(:) + Dxx(:))/rR(:)
    Ryy_t(:) = -Fyy(:) + Byy(:) + Cyy(:) + Pyy(:) - Eyy(:) + (PIyy(:) - Tyyy_y(:) - Gyy(:) + Dyy(:))/rR(:)
    Rzz_t(:) = -Fzz(:) + Bzz(:) + Czz(:) + Pzz(:) - Ezz(:) + (PIzz(:) - Tzzy_y(:) - Gzz(:) + Dzz(:))/rR(:)
    Rxy_t(:) = -Fxy(:) + Bxy(:) + Cxy(:) + Pxy(:) - Exy(:) + (PIxy(:) - Txyy_y(:) - Gxy(:) + Dxy(:))/rR(:)
    Rxz_t(:) = -Fxz(:) + Bxz(:) + Cxz(:) + Pxz(:) - Exz(:) + (PIxz(:) - Txzy_y(:) - Gxz(:) + Dxz(:))/rR(:)
    Ryz_t(:) = -Fyz(:) + Byz(:) + Cyz(:) + Pyz(:) - Eyz(:) + (PIyz(:) - Tyzy_y(:) - Gyz(:) + Dyz(:))/rR(:)

    ! Kinetic energy equation
    Tke(:) = 0.5_wp*(Rxx(:) + Ryy(:) + Rzz(:))

    Buo(:) = 0.5_wp*(Bxx(:) + Byy(:) + Bzz(:))
    Con(:) = 0.5_wp*(Cxx(:) + Cyy(:) + Czz(:))
    Prd(:) = 0.5_wp*(Pxx(:) + Pyy(:) + Pzz(:))
    Pi(:) = 0.5_wp*(PIxx(:) + PIyy(:) + PIzz(:))
    Eps(:) = 0.5_wp*(Exx(:) + Eyy(:) + Ezz(:))
    Ty_y(:) = 0.5_wp*(Txxy_y(:) + Tyyy_y(:) + Tzzy_y(:))
    Gkin(:) = 0.5_wp*(Gxx(:) + Gyy(:) + Gzz(:))
    Dkin(:) = 0.5_wp*(Dxx(:) + Dyy(:) + Dzz(:))

    Tke_t(:) = Buo(:) + Con(:) + Prd(:) - Eps(:) + (-Ty_y(:) + Pi(:) - Gkin(:) + Dkin(:))/rR(:)

    ! Potential energy equation
    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == nse_eqns)) then
        Pot(:) = -rB(:)*(g(2)%nodes(:) - sbg(inb_scal)%ymean)

    else
        Pot(:) = -rR(:)*(g(2)%nodes(:) - rbg%ymean)*buoyancy%vector(2)

    end if

    ! ###################################################################
    ! Output
    ! ###################################################################
    ! 14 t-dependent variables, for consistency with old format
    ! ng = ng +1
    ! groupname(ng) = ''
    ! varname(ng)   = 'dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy'
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''

    write (name, *) itime; name = 'avg'//trim(adjustl(name))
    call IO_WRITE_AVERAGES(name, itime, rtime, jmax, nv, ng, g(2)%nodes, varname, groupname, mean2d)

    return
end subroutine AVG_FLOW_XZ
