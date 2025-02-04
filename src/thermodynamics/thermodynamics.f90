#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# Dimensional and nondimensional formulations. In the latter case:
!# - Compressible formulations use p nondimensionalized by reference dynamic pressure rho_0 U_0^2
!# - Incompressible and anelastic formulations use p nondimensionalized by reference thermodynamic pressure p_0
!#
!# Mixture:
!# inb_scal          # of prognostic scalars (transported during the simulation and saved)
!# inb_scal_array    # of scalars in array s (prognositc+diagnostic, normally = inb_scal)
!# NSP               # of species in the mixture (NSP>=inb_scal)
!#
!# General multispecies formulation retains NSP-1 species in scalar array s, Y_i=s_i, the last one is obtained by sum Y_i=1
!# There might be additional scalars at the end of array s, e.g., conserved scalars or diagnostic variables
!# Multispecies admits the case Y_i=f_i(s_j), s_j could be conserved scalar, e.g. using total water or mixture fraction
!#
!# Saturation pressure implies that reference R_0 in non-dimensionalization is such that reference pressure is 1 bar.
!#
!########################################################################
module Thermodynamics
    use TLab_Constants, only: wp, wi, efile, lfile, MAX_PROF
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    integer, public :: imode_thermo
    integer, parameter, public :: THERMO_TYPE_NONE = 0
    integer, parameter, public :: THERMO_TYPE_COMPRESSIBLE = 1
    integer, parameter, public :: THERMO_TYPE_ANELASTIC = 2
    integer, parameter, public :: THERMO_TYPE_LINEAR = 3

    ! Thermodynamic description of the mixture
    integer(wi), public :: imixture
    character*128 :: chemkin_file                       ! File with thermodynamic data, if used

    integer, parameter, public :: MAX_NSP = 10          ! Maximum number of components (species) in a mixture
    integer(wi), public :: NSP = 0                      ! Number of components (species) in a mixture
    character(len=32), public :: THERMO_SPNAME(MAX_NSP) = ''

    ! Thermal data
    real(wp), public :: THERMO_R(MAX_NSP)               ! Gas constants

    ! Caloric data
    integer(wi), parameter :: MAX_NCP = 7               ! Caloric data; cp(T), formation enthalpy, formation entropy
    integer(wi), public :: NCP                          ! Number of terms in polynomial fit to cp
    real(wp), public :: THERMO_AI(MAX_NCP, 2, MAX_NSP)  ! Polynomial coefficients. Last to terms are formation enthalpy and formation entropy
    !                                                   The second index indicates each of the 2 temperature intervals considered in the fit
    real(wp), public :: THERMO_TLIM(3, MAX_NSP)         ! Temperature limits of the two temperature intervals in the polynomial fits.

    ! Saturation vapor pressures
    integer(wi), parameter :: MAX_NPSAT = 10            ! Polynomial fit to saturation pressure
    integer(wi), public :: NPSAT
    real(wp), public :: THERMO_PSAT(MAX_NPSAT), NEWTONRAPHSON_ERROR

    ! Compressible formulation, different combinations of parameters \gamma0 and mach to save calculations
    real(wp), public :: gama0                           ! Specific heat ratio, Cp0/Cv0 = Cp0/(Cp0-R0)
    !                                                     For imixture=NONE, I only need gama0 and it is set in tlab.ini
    !                                                     Otherwise, I need the thermodynamic data that is given in thermo_initialize, and gama0 is derived.
    real(wp), public :: RRATIO                          ! 1/(gama0 mach^2) = R0/(U0^2/T0)
    real(wp), public :: CRATIO_INV                      ! (gamma0-1)*mach^2 = (U0^2/T0)/Cp0
    real(wp), public :: RRATIO_INV                      ! gama0 mach^2 = (U0^2/T0)/R0, inverse of RRATIO to save computational time in some routines
    ! Anelastic and incompressible formulation
    real(wp), public :: GRATIO                          ! (gama0-1)/gama0 = R0/Cp0
    real(wp), public :: scaleheightinv                  ! Inverse of pressure scale height. Equivalent to 1/(Fr*RRATIO) in compressible formulation
    ! nondimensional parameters; to be moved to navierstokes
    real(wp), public :: mach                                ! compressibility


    ! Nondimensional formulation
    logical, public :: nondimensional = .true.          ! consider nondimensional formulation
    !                                                   A dimensional formulation can be imposed by setting RRATIO=CRATIO_INV=1, or GRATIO=1

    real(wp), public :: thermo_param(MAX_PROF)          ! Additional data
    real(wp), public :: dsmooth                         ! Smoothing factor for derivative discontinuity in inifinitely fast chemistry and saturation adjustment

    ! Derived parameters, for clarity in airwater formulation
    real(wp), public :: Rv, Rd, Rdv, Cd, Cl, Cdv, Cvl, Cdl, Lv0, Ld, Lv, Ldv, Lvl, Ldl, rd_ov_rv, rd_ov_cd, PREF_1000

    ! Transport phenomena
    integer, public :: itransport                       ! variable viscosity

    public :: Thermodynamics_Initialize_Parameters
    public :: Thermo_Psat_Polynomial, Thermo_dPsat_Polynomial

contains
    !########################################################################
    !########################################################################
    subroutine Thermodynamics_Initialize_Parameters(inifile)
        use TLab_Memory, only: inb_scal, inb_scal_array
        ! use THERMO_ANELASTIC, only: scaleheightinv, GRATIO

        character(len=*), intent(in), optional :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=512) sRes
        integer(wi) idummy

        real(wp) WGHT(MAX_NSP)                              ! Molar masses
        real(wp) TREF_LOC, HREF_LOC(MAX_NSP), SREF_LOC(MAX_NSP)
        integer(wi) icp, is, im
        real(wp) WRK1D_LOC(MAX_NPSAT)
        integer(wi) ipsat, i, j
        real(wp) tmp1, tmp2
        character*46 str
        logical :: molar_data = .true.

        real(wp), parameter :: RGAS = 8314_wp               ! J /kg /K, universal gas constant,
        ! Reference values for the nondimensionalization; in principle, these could be changed, if desired.
        real(wp), parameter :: TREF = 298.0_wp              ! K, reference temperature T_0
        real(wp), parameter :: PREF = 1e5_wp                ! Pa, reference pressure p_0
        integer, parameter :: ISPREF = 2                    ! Species 2 is taken as reference species for CPREF and RREF
        real(wp) CPREF, RREF                                ! Reference cp, referece gas constant R_0
        !                                                     Reference density results from rho_0=p_0/(T_0R_0)
        ! Reference to calculate potential temperatures. Global variable.
        PREF_1000 = 1e5_wp                                  ! 1000 hPa,

        !########################################################################
        if (present(inifile)) then
            bakfile = trim(adjustl(inifile))//'.bak'
            block = 'Thermodynamics'

            call TLab_Write_ASCII(bakfile, '#')
            call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
            call TLab_Write_ASCII(bakfile, '#Type=<compressible/anelastic>')
            call TLab_Write_ASCII(bakfile, '#Gama=<value>')
            call TLab_Write_ASCII(bakfile, '#Mixture=<value>')
            call TLab_Write_ASCII(bakfile, '#Transport=<constant/powerlaw/sutherland>')
            call TLab_Write_ASCII(bakfile, '#Parameters=<value>')
            call TLab_Write_ASCII(bakfile, '#Nondimensional=<yes,no>')

            call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
            if (trim(adjustl(sRes)) == 'none') then; imode_thermo = THERMO_TYPE_NONE
            else if (trim(adjustl(sRes)) == 'compressible') then; imode_thermo = THERMO_TYPE_COMPRESSIBLE
            else if (trim(adjustl(sRes)) == 'anelastic') then; imode_thermo = THERMO_TYPE_ANELASTIC
            else if (trim(adjustl(sRes)) == 'linear') then; imode_thermo = THERMO_TYPE_LINEAR
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Error in Thermodynamics.Type.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

            call ScanFile_Real(bakfile, inifile, block, 'HeatCapacityRatio', '1.4', gama0)      ! needed in compressible formulation
            call ScanFile_Real(bakfile, inifile, block, 'ScaleHeight', '0.0', scaleheightinv)   ! needed in anelastic formulation
            if (scaleheightinv > 0.0_wp) scaleheightinv = 1.0_wp/scaleheightinv

            call ScanFile_Char(bakfile, inifile, block, 'Mixture', 'None', sRes)
            if (trim(adjustl(sRes)) == 'none') &
                call ScanFile_Char(bakfile, inifile, 'Main', 'Mixture', 'None', sRes)                   ! backwards compatibility, to be removed
            if (trim(adjustl(sRes)) == 'none') then; imixture = MIXT_TYPE_NONE
            else if (trim(adjustl(sRes)) == 'air') then; imixture = MIXT_TYPE_AIR
            else if (trim(adjustl(sRes)) == 'airvapor') then; imixture = MIXT_TYPE_AIRVAPOR
            else if (trim(adjustl(sRes)) == 'airwater') then; imixture = MIXT_TYPE_AIRWATER
            else if (trim(adjustl(sRes)) == 'airwaterlinear') then; imixture = MIXT_TYPE_AIRWATER_LINEAR
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Error in Thermodynamics.Mixture.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

            call ScanFile_Char(bakfile, inifile, block, 'Transport', 'None', sRes)
            if (trim(adjustl(sRes)) == 'sutherland') then; itransport = EQNS_TRANS_SUTHERLAND; 
            elseif (trim(adjustl(sRes)) == 'powerlaw') then; itransport = EQNS_TRANS_POWERLAW; 
            else; itransport = EQNS_NONE; end if

            if (imixture /= EQNS_NONE) then
                thermo_param(:) = 0.0_wp
                call ScanFile_Char(bakfile, inifile, block, 'Parameters', '1.0', sRes)
                idummy = MAX_PROF
                call LIST_REAL(sRes, idummy, thermo_param)

            end if

            if (imixture == MIXT_TYPE_AIRWATER) then
                call ScanFile_Real(bakfile, inifile, block, 'SmoothFactor', '0.1', dsmooth)
            end if

            call ScanFile_Char(bakfile, inifile, block, 'Nondimensional', 'yes', sRes)
            if (trim(adjustl(sRes)) == 'yes') then; nondimensional = .true.
            else if (trim(adjustl(sRes)) == 'no') then; nondimensional = .false.
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Error in Thermodynamics.Nondimensional')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

        end if

        if (imode_thermo == THERMO_TYPE_NONE) return

        ! ###################################################################
        ! Species tags
        ! Thermal equation, molar masses in kg/kmol
        ! ###################################################################
        WGHT(:) = 1.0_wp            ! We devide by WGTH below even when mxiture is none

        select case (imixture)
            ! -------------------------------------------------------------------
            ! Burke-Schuman case: transport just mixture fraction, and then equilibrium
            ! 4 species + Nitrogen + Conserved Scalar
            ! -------------------------------------------------------------------
        case (MIXT_TYPE_BS, MIXT_TYPE_QUASIBS)
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'CH4'; WGHT(NSP) = 16.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'O2 '; WGHT(NSP) = 32.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'H2O'; WGHT(NSP) = 18.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'CO2'; WGHT(NSP) = 44.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'N2 '; WGHT(NSP) = 28.0_wp

            ! -------------------------------------------------------------------
            ! Peters Mechanism for Methane
            ! 7 species + Nitrogen + Conserved Scalar
            ! -------------------------------------------------------------------
        case (MIXT_TYPE_PETERS1991, MIXT_TYPE_PETERS1988)
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'CH4'; WGHT(NSP) = 16.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'O2 '; WGHT(NSP) = 32.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'H2O'; WGHT(NSP) = 18.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'CO2'; WGHT(NSP) = 44.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'CO '; WGHT(NSP) = 28.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'H2 '; WGHT(NSP) = 2.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'H  '; WGHT(NSP) = 1.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'N2 '; WGHT(NSP) = 28.0_wp

            ! -------------------------------------------------------------------
            ! Unimolecular decomposition flame
            ! 1 reactant + 1 product + Conserved Scalar
            ! -------------------------------------------------------------------
        case (MIXT_TYPE_UNIDECOMP)
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Reactant'; WGHT(NSP) = 32.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Product'; WGHT(NSP) = 32.0_wp

            ! -------------------------------------------------------------------
            ! Unimolecular decomposition flame
            ! 2 reactant + 1 product + Conserved Scalar
            ! -------------------------------------------------------------------
        case (MIXT_TYPE_ONESTEP)
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Reactant'; WGHT(NSP) = 32.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Oxidizer'; WGHT(NSP) = 32.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Product'; WGHT(NSP) = 32.0_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Inert'; WGHT(NSP) = 32.0_wp

            ! -------------------------------------------------------------------
            ! Air; data from Iribarne and Godson, 1981
            ! -------------------------------------------------------------------
        case (MIXT_TYPE_AIR)
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Air'; WGHT(NSP) = 28.9644_wp
            WGHT(2) = 28.9644_wp  ! needed for the nondimensionalization, based on species 2; to be fixed

            ! -------------------------------------------------------------------
            ! Water vapor and air; data from Iribarne and Godson, 1981
            ! -------------------------------------------------------------------
        case (MIXT_TYPE_AIRVAPOR)
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'H2Ov'; WGHT(NSP) = 18.015_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Air'; WGHT(NSP) = 28.9644_wp

            ! -------------------------------------------------------------------
            ! Water vapor, air and liquid water; data from Iribarne and Godson, 1981
            ! Compressible:   Transport    q_t, and q_l from equilibrium; add space for q_l
            ! Incompressible: Transport h, q_t, and q_l from equilibrium; add space for q_l
            ! -------------------------------------------------------------------
        case (MIXT_TYPE_AIRWATER)
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'H2Ov'; WGHT(NSP) = 18.015_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Air'; WGHT(NSP) = 28.9644_wp
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'H2Ol'; WGHT(NSP) = 18.015_wp

            inb_scal_array = inb_scal + 1   ! assumes phase equilibrium so that q_l is not prognostic

            ! -------------------------------------------------------------------
            ! Linearized thermodynamics for stratocumulus case; mainly for the tags
            ! -------------------------------------------------------------------
        case (MIXT_TYPE_AIRWATER_LINEAR)
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Chi'   ! scaled total water (mixing fraction)
            NSP = NSP + 1; THERMO_SPNAME(NSP) = 'Psi'   ! enthalpy deviations from the pure mixing case (described by mixing fraction only)
            do is = 3, inb_scal
                NSP = NSP + 1; write (THERMO_SPNAME(NSP), *) is; THERMO_SPNAME(NSP) = 'Scalar'//trim(adjustl(THERMO_SPNAME(NSP)))
            end do

            inb_scal_array = inb_scal + 1
            NSP = inb_scal_array; THERMO_SPNAME(NSP) = 'Liquid' ! normalized concentration of liquid.

        end select

        ! ###################################################################
        ! Caloric equations
        !
        ! General formulation is CHEMKIN format: 7-coefficient NASA
        ! polynomials (see Burcat&Ruscic):
        !
        ! THERMO_AI(i,im,k) = a_i of species k at
        !         im=1 high temperature
        !         im=2 low   temperature
        !
        ! C_{p,i} = \sum_1^5 a_i T^{i-1}
        ! h_i     = \sum_1^5 a_i T^i/i + a_6
        ! s_{T,i} = a_1 ln(T) + \sum_2^5 a_i T^{i-1}/(i-1) + a_7
        !
        ! i.e., dh_i = C_{p,i}dT and ds_{T,i}=C_{p,i}dT/T, where a_6 is
        ! related to the formation enthalpy, and a_7 to the formation entropy.
        ! HREF_LOC and SREF_LOC (at TREF_LOC) are used to fix last Cpi 6-7 coefficients.
        !
        ! Note that Burcat&Ruscic give values devided by R^0
        !
        ! NCP gives the number of coefficients a_i.
        ! The simplified situation NCP=1 assumes also only one range, which then gives C_p constant.
        ! This is used to expedite the calculation of T in the energy formulation.
        !
        ! The pressure contribution to the entropy still needs to be added
        !
        ! ###################################################################
        THERMO_AI(:, :, :) = 0.0_wp               ! Initialize to zero
        HREF_LOC(:) = 0.0_wp
        SREF_LOC(:) = 0.0_wp

        THERMO_TLIM(1, :) = 200.0_wp    ! Default T limits for the 2 intervals of polynomial fits
        THERMO_TLIM(2, :) = 5000.0_wp   ! These intervals are currently not used
        THERMO_TLIM(3, :) = 5000.0_wp

        NCP = 1                         ! Default order of heat capacity polynomial

        select case (imixture)
        case (MIXT_TYPE_BS, MIXT_TYPE_QUASIBS)      ! Methane
            TREF_LOC = 298.0_wp         ! K, auxiliar value to define caloric data

            ! Enthalpy of Formation in Jules/Kmol
            HREF_LOC(1) = -74.0e+6_wp
            HREF_LOC(2) = 0.0_wp
            HREF_LOC(3) = -241.82e+6_wp
            HREF_LOC(4) = -393.51e+6_wp
            HREF_LOC(5) = 0.0_wp

            ! Entropy of Formation in Jules/(Kelvin Kmol)
            SREF_LOC(1) = 186.37e+3_wp
            SREF_LOC(2) = 205.15e+3_wp
            SREF_LOC(3) = 188.83e+3_wp
            SREF_LOC(4) = 213.78e+3_wp
            SREF_LOC(5) = 191.61e+3_wp

            ! Heat capacity polynomial. Values taken from fit with T-TREF_LOC instead T
            NCP = 2
            THERMO_AI(1, :, 1) = 35.70e+3_wp - 42.4833_wp*TREF_LOC; THERMO_AI(2, :, 1) = 42.4833_wp
            THERMO_AI(1, :, 2) = 28.96e+3_wp - 6.21666_wp*TREF_LOC; THERMO_AI(2, :, 2) = 6.21666_wp
            THERMO_AI(1, :, 3) = 32.76e+3_wp - 11.9570_wp*TREF_LOC; THERMO_AI(2, :, 3) = 11.9570_wp
            THERMO_AI(1, :, 4) = 37.22e+3_wp - 17.6500_wp*TREF_LOC; THERMO_AI(2, :, 4) = 17.6500_wp
            THERMO_AI(1, :, 5) = 28.88e+3_wp - 4.70833_wp*TREF_LOC; THERMO_AI(2, :, 5) = 4.70833_wp

            ! -------------------------------------------------------------------
        case (MIXT_TYPE_UNIDECOMP)      ! Simple reaction, R -> P
            TREF_LOC = 298.0_wp         ! K, auxiliar value to define caloric data

            ! Enthalpy of Formation in Jules/Kmol
            HREF_LOC(1) = 0.0_wp
            HREF_LOC(2) = -86.71502e+6_wp

            ! Entropy of Formation in Jules/(Kelvin Kmol)
            SREF_LOC(1) = 205.15e+3_wp
            SREF_LOC(2) = 205.15e+3_wp

            ! Heat capacity polynomial
            THERMO_AI(1, :, 1) = 29.099e+3_wp
            THERMO_AI(1, :, 2) = 29.099e+3_wp   ! = 59.099e+3_wp

            ! -------------------------------------------------------------------
        case (MIXT_TYPE_ONESTEP)        ! Simple reaction, R + O -> P
            TREF_LOC = 298.0_wp         ! K, auxiliar value to define caloric data

            ! Enthalpy of Formation in Jules/Kmol
            HREF_LOC(1) = 0.0_wp
            HREF_LOC(2) = 0.0_wp
            HREF_LOC(3) = -86.71502e+6_wp
            HREF_LOC(4) = 0.0_wp

            ! Entropy of Formation in Jules/(Kelvin Kmol)
            SREF_LOC(1) = 205.15e+3_wp
            SREF_LOC(2) = 205.15e+3_wp
            SREF_LOC(3) = 205.15e+3_wp
            SREF_LOC(4) = 205.15e+3_wp

            ! Heat capacity polynomial
            THERMO_AI(1, :, 1) = 29.099e+3_wp
            THERMO_AI(1, :, 2) = 29.099e+3_wp
            THERMO_AI(1, :, 3) = 29.099e+3_wp
            THERMO_AI(1, :, 4) = 29.099e+3_wp

            ! -------------------------------------------------------------------
        case (MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR)
            TREF_LOC = 273.15_wp

            molar_data = .false. ! this block gives data already in spefic

            ! Enthalpy of Formation in J /kg
            HREF_LOC(1) = 1870.0_wp*TREF_LOC                 ! values s.t. THERMO_AI(6,im,1:2) = 0, i.e., liquid-water enthalpy
            HREF_LOC(2) = 1007.0_wp*TREF_LOC
            HREF_LOC(3) = 1870.0_wp*TREF_LOC - 2501600_wp    ! latent heat of vaporization at 273.15 K is 2501.6 kJ/kg

            ! Entropy of Formation in J /kg /K
            SREF_LOC(1) = 0.0_wp
            SREF_LOC(2) = 0.0_wp
            SREF_LOC(3) = -2501600_wp/TREF_LOC

            ! Heat capacity polynomial in J /kg /K; using only the low temperature range
            THERMO_AI(1, :, 1) = 1870.0_wp     ! water vapor
            THERMO_AI(1, :, 2) = 1007.0_wp     ! dry air
            THERMO_AI(1, :, 3) = 4217.6_wp     ! liquid water

            ! -------------------------------------------------------------------
        case (MIXT_TYPE_CHEMKIN) ! Load thermodynamic data from chemkin file
            call THERMO_READ_CHEMKIN(chemkin_file)
            THERMO_AI = THERMO_AI*RGAS
            NCP = 5             ! Fifth-order polynomial

        end select

        ! -------------------------------------------------------------------
        ! Final calculations

        ! 6th and 7th coefficients of heat capacity polynomial are calculated from reference enthalpy
        do is = 1, NSP
            THERMO_AI(6, :, is) = HREF_LOC(is) - THERMO_AI(1, :, is)*TREF_LOC - THERMO_AI(2, :, is)*TREF_LOC*TREF_LOC*0.5_wp
            THERMO_AI(7, :, is) = SREF_LOC(is) - THERMO_AI(2, :, is)*TREF_LOC
        end do

        if (molar_data) then ! Change heat capacities from molar to mass specific, i.e., J /K /kg
            do is = 1, NSP
                THERMO_AI(:, :, is) = THERMO_AI(:, :, is)/WGHT(is)
            end do
        end if

        ! ###################################################################
        ! Phase change. Polynomial fit to saturation vapor pressure
        !
        ! p_sat(T) = \sum_1^9 a_i T^{i-1}
        !
        ! ###################################################################
        THERMO_PSAT(:) = 0.0_wp ! Initialize to zero
        NPSAT = 0

        select case (imixture)

        case (MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR)

            ! Flatau et al., J. Applied Meteorol., 1507-1513, 1992
            NPSAT = 9
            WRK1D_LOC(1) = 0.611213476e+3_wp    ! Pa
            WRK1D_LOC(2) = 0.444007856e+2_wp
            WRK1D_LOC(3) = 0.143064234e+1_wp
            WRK1D_LOC(4) = 0.264461437e-1_wp
            WRK1D_LOC(5) = 0.305930558e-3_wp
            WRK1D_LOC(6) = 0.196237241e-5_wp
            WRK1D_LOC(7) = 0.892344772e-8_wp
            WRK1D_LOC(8) = -0.373208410e-10_wp
            WRK1D_LOC(9) = 0.209339997e-13_wp

            ! going from powers of (T-T_ref) to T, with T_ref = 273.15 K
            TREF_LOC = 273.15_wp
            do ipsat = 1, NPSAT
                THERMO_PSAT(ipsat) = 0.0_wp
                do i = ipsat, NPSAT
                    tmp1 = 1.0_wp
                    do j = i - 1, i - ipsat + 1, -1
                        tmp1 = tmp1*real(j, wp)
                    end do
                    THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat) + &
                                         WRK1D_LOC(i)*TREF_LOC**(i - 1)*tmp1*(-1)**(i - ipsat)
                end do
                tmp2 = 1.0_wp
                do j = ipsat - 1, 1, -1
                    tmp2 = tmp2*real(j, wp)
                end do
                THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)/tmp2/TREF_LOC**(ipsat - 1)
            end do

        case default

        end select

        ! ###################################################################
        ! Final calculations
        ! ###################################################################
        THERMO_R(:) = RGAS/WGHT(:)              ! Specific gas constants, J /kg /K
        RREF = THERMO_R(ISPREF)                 ! Reference value R_0

        CPREF = 0.0_wp                          ! Reference heat capacity Cp0, J /kg /K
        do icp = NCP, 1, -1
            CPREF = CPREF*TREF + THERMO_AI(icp, 2, ISPREF)
        end do

        if (imixture /= MIXT_TYPE_NONE) then    ! Reference heat capacity ratio; othewise, gama0 is read in tlab.ini
            gama0 = CPREF/(CPREF - RREF)
        end if

        ! -------------------------------------------------------------------
        ! Nondimensionalization
        !- Compressible formulations use p nondimensionalized by reference dynamic pressure rho_0 U_0^2
        !- Incompressible and anelastic formulations use p nondimensionalized by reference thermodynamic pressure p_0
        if (nondimensional) then
            ! Thermal equation of state
            THERMO_R(:) = THERMO_R(:)/RREF              ! normalized gas constants (Inverse of molar masses)

            ! Caloric equations of state
            do is = 1, NSP
                do im = 1, 2
                    THERMO_AI(6, im, is) = THERMO_AI(6, im, is)/(CPREF*TREF)    ! Formation enthalpy
                    THERMO_AI(7, im, is) = THERMO_AI(7, im, is)/CPREF           ! Formation entropy
                    do icp = 1, NCP
                        THERMO_AI(icp, im, is) = THERMO_AI(icp, im, is)*(TREF**(icp - 1))/CPREF
                    end do
                end do
            end do
            THERMO_TLIM = THERMO_TLIM/TREF              ! Temperature limis for polynomial fits to cp

            ! Pressures
            PREF_1000 = PREF_1000/PREF                  ! 1000 hPa, reference to calculate potential temperatures

            THERMO_PSAT(:) = THERMO_PSAT(:)/PREF        ! Saturation vapor pressure
            do ipsat = 1, NPSAT
                THERMO_PSAT(ipsat) = THERMO_PSAT(ipsat)*(TREF**(ipsat - 1))
            end do

        else
            if (imixture == MIXT_TYPE_NONE) then
                call TLab_Write_ASCII(efile, __FILE__//'. Single species formulation must be nondimensional.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if

        ! Parameters in the evolution equations
        RRATIO = 1.0_wp                                 ! Compressible formulation uses RRATIO, CRATIO_INV, but RRATIO also used below
        CRATIO_INV = 1.0_wp
        GRATIO = 1.0_wp                                 ! Anelastic formulation uses GRATIO, but GRATIO also used below
        if (nondimensional) then
            if (imode_thermo == THERMO_TYPE_COMPRESSIBLE) then
                RRATIO = 1.0_wp/(gama0*mach*mach)       ! (R_0T_0)/U_0^2 = p_0/(rho_0U_0^2), a scaled reference pressure
                CRATIO_INV = (gama0 - 1.0_wp)*mach*mach
            end if
            GRATIO = (gama0 - 1.0_wp)/gama0             ! R_0/C_{p,0}
        end if

        ! Derived parameters to save operations
        PREF_1000 = PREF_1000*RRATIO                    ! scaling by dynamic reference pressure U0^2/T0 for compressible mode
        THERMO_PSAT(:) = THERMO_PSAT(:)*RRATIO          ! this assumes that RRATIO is 1 in anelastic, incompressible mode
        THERMO_R(:) = THERMO_R(:)*RRATIO
        RRATIO_INV = 1.0_wp/RRATIO

        ! -------------------------------------------------------------------
        ! Definitions for clarity in the thermodynamics code
        select case (imixture)
        case (MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR)
            Rv = THERMO_R(1)
            Rd = THERMO_R(2)
            Rdv = THERMO_R(1) - THERMO_R(2)

            Cd = THERMO_AI(1, 1, 2)
            Cl = THERMO_AI(1, 1, 3)
            Cdv = THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 2)
            Cvl = THERMO_AI(1, 1, 3) - THERMO_AI(1, 1, 1)
            Cdl = THERMO_AI(1, 1, 3) - THERMO_AI(1, 1, 2)
            Lv0 = -THERMO_AI(6, 1, 3)
            Ld = THERMO_AI(6, 1, 2)
            Lv = THERMO_AI(6, 1, 1)
            Ldv = THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 2)
            Lvl = THERMO_AI(6, 1, 3) - THERMO_AI(6, 1, 1)
            Ldl = THERMO_AI(6, 1, 3) - THERMO_AI(6, 1, 2)
            rd_ov_rv = Rd/Rv
            rd_ov_cd = Rd/Cd*GRATIO

        end select

        ! ###################################################################
        ! Output
        ! ###################################################################
        if (imixture /= MIXT_TYPE_NONE) then        ! othewise, gama0 is read in tlab.ini
            call TLab_Write_ASCII(lfile, 'Thermodynamic properties have been initialized.')
            do is = 1, NSP
                write (str, *) is; str = 'Setting Species'//trim(adjustl(str))//'='//trim(adjustl(THERMO_SPNAME(is)))
                call TLab_Write_ASCII(lfile, str)
            end do
            write (str, 1010) 'Setting RREF = ', RREF
            call TLab_Write_ASCII(lfile, str)
            write (str, 1010) 'Setting TREF = ', TREF
            call TLab_Write_ASCII(lfile, str)
            write (str, 1010) 'Setting PREF = ', PREF
            call TLab_Write_ASCII(lfile, str)
            if (NPSAT > 0) then
                write (str, 1010) 'Setting RHOREF = ', PREF/(RREF*TREF)
                call TLab_Write_ASCII(lfile, str)
            end if
            write (str, 1020) 'Setting CPREF = ', CPREF
            call TLab_Write_ASCII(lfile, str)
            write (str, 1020) 'Setting Gama0 = ', gama0
            call TLab_Write_ASCII(lfile, str)
        end if

        return

1010    format(A14, 1x, G_FORMAT_R)
1020    format(A15, 1x, G_FORMAT_R)

    end subroutine Thermodynamics_Initialize_Parameters

    !########################################################################
    !########################################################################
    subroutine THERMO_READ_CHEMKIN(name)    ! Read chemkin mixtures ! To be checked
        character(len=*), intent(in) :: name

        ! -----------------------------------------------------------------------
        real(wp) T1, T2, T3
        integer(wi) il, i, is
        logical frun
        character*15 token
        character*225 wline
        character*80 line, line1, line2, line3
        integer(wi) THERMO_FLAG(MAX_NSP)

        integer, parameter :: i23 = 23

        ! #######################################################################
        ! Initialize thermodynamic data structure
        do is = 1, NSP
            THERMO_FLAG(is) = 0
        end do

        ! Read Thermodynamic file
        open (i23, file=name, status='old')

        rewind (i23)

        ! Read Header
        read (i23, *) line
        call TLab_Write_ASCII(lfile, line)

        if (trim(adjustl(line)) /= 'THERMO') then
            call TLab_Write_ASCII(efile, 'THERMO_READ_CHEMKIN. Thermodynamic file format error')
            call TLab_Stop(DNS_ERROR_THERMOFORMAT)
        end if

        ! Read Temperature ranges
        read (i23, *) T1, T2, T3
        write (wline, *) T1, T2, T3
        call TLab_Write_ASCII(lfile, wline(1:80))

        ! Remove comments
        frun = .true.
        do while (frun)
            read (i23, '(A80)', end=50) line
            if (line(1:1) /= '!') frun = .false.
        end do

        frun = .true.
        do while (frun)
            !    Check for end of file
            read (line, '(A15)', end=50) token
            if (trim(adjustl(token)) == 'END') then
                frun = .false.
                goto 50
            end if

            !    Read all relevant information
            read (i23, '(A80)', end=50) line1
            read (i23, '(A80)', end=50) line2
            read (i23, '(A80)', end=50) line3

            !    Process lines
            do is = 1, NSP
                if (trim(adjustl(token)) == THERMO_SPNAME(is)) then
                    call TLab_Write_ASCII(lfile, line)
                    call TLab_Write_ASCII(lfile, line1)
                    call TLab_Write_ASCII(lfile, line2)
                    call TLab_Write_ASCII(lfile, line3)

                    !          Required species found, process information
                    !          Get limit temperatures
                    do i = 1, 225
                        wline(i:i) = ' '
                    end do
                    wline = line(46:75)
                    read (wline, *) (THERMO_TLIM(i, is), i=1, 3)

                    !          Concatenate lines so read is simpler
                    wline = line1(1:75)//line2(1:75)//line3(1:75)

                    do i = 1, 14
                        il = (i - 1)*15 + 1
                        read (wline(il:il + 14), *) THERMO_AI(i, 1, is)
                    end do

                    THERMO_FLAG(is) = 1
                end if
            end do

            !    Read next line
            read (i23, '(A80)', end=50) line

        end do

50      close (i23)

        do is = 1, NSP
            if (THERMO_FLAG(is) == 0) then
                call TLab_Write_ASCII(efile, 'THERMO_READ_CHEMKIN. Not all thermodynamic data contained in thermo file')
                call TLab_Stop(DNS_ERROR_THERMOCONT)
            end if
        end do

        return
    end subroutine THERMO_READ_CHEMKIN

    ! ###################################################################
    ! ###################################################################
    subroutine Thermo_Psat_Polynomial(ijmax, T, p)
        integer(wi) ijmax
        real(wp) T(ijmax)
        real(wp) p(ijmax)

        ! -------------------------------------------------------------------
        integer(wi) i, ipsat

        ! ###################################################################
        if (NPSAT > 0) then
            do i = 1, ijmax
                p(i) = THERMO_PSAT(NPSAT)
                do ipsat = NPSAT - 1, 1, -1
                    p(i) = p(i)*T(i) + THERMO_PSAT(ipsat)
                end do
            end do
        else
            p(:) = 0.0_wp
        end if

        return
    end subroutine Thermo_Psat_Polynomial

    ! ###################################################################
    ! ###################################################################
    subroutine Thermo_dPsat_Polynomial(ijmax, T, dp)
        integer(wi) ijmax
        real(wp) T(ijmax)
        real(wp) dp(ijmax)

        ! -------------------------------------------------------------------
        integer(wi) i, ipsat

        ! ###################################################################
        if (NPSAT > 0) then
            do i = 1, ijmax
                dp(i) = 0.0_wp
                do ipsat = NPSAT - 1, 1, -1
                    dp(i) = dp(i)*T(i) + THERMO_PSAT(ipsat + 1)*real(ipsat, wp)
                end do
            end do
        else
            dp(:) = 0.0_wp
        end if

        return
    end subroutine Thermo_dPsat_Polynomial

end module Thermodynamics
