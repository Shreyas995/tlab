#include "dns_const.h"

! Implementation of airwater cases
! Airwater considers the special case where qt=qv+ql and ql instead of qv and ql are used as prognostic

! we use explicit-shape arrays in arguments for the routines to be callable with different array ranks

module THERMO_AIRWATER
    use TLAB_VARS, only: inb_scal
    use TLAB_CONSTANTS, only: wp, wi, small_wp, big_wp
    use THERMO_VARS, only: imixture, WGHT_INV
    use THERMO_VARS, only: gama0, CRATIO_INV, NCP, THERMO_AI, THERMO_TLIM
    use THERMO_VARS, only: MRATIO, RRATIO, GRATIO
    use THERMO_VARS, only: THERMO_PSAT, NPSAT
    use THERMO_VARS, only: Cd, Cdv, Cvl, Cdl, Ld, Lv, Ldv, Lvl, Ldl, Rd, Rdv, Rv, rd_ov_rv
    use THERMO_VARS, only: dsmooth, NEWTONRAPHSON_ERROR
    use THERMO_VARS, only: thermo_param
    implicit none
    private

    integer(wi) ij, is
    real(wp) RMEAN, CPMEAN, T_LOC, psat, qsat, dsmooth_loc, dummy

    public :: THERMO_AIRWATER_PT
    public :: THERMO_AIRWATER_RE
    public :: THERMO_AIRWATER_LINEAR
    public :: THERMO_AIRWATER_LINEAR_SOURCE

contains
    !########################################################################
    !# Calculate liquid content from p, T and water content.
    !########################################################################
    subroutine THERMO_AIRWATER_PT(ijmax, s, p, T)
        integer(wi) ijmax
        real(wp), intent(in) :: p(ijmax), T(ijmax)
        real(wp), intent(out) :: s(ijmax, *)

        real(wp) dqldqt
        
        ! ###################################################################
        call THERMO_POLYNOMIAL_PSAT(ijmax, T, s(1, 2))
        do ij = 1, ijmax
            ! this is really the vapor content
            qsat = 1.0_wp/(MRATIO*p(ij)/s(ij, 2) - 1.0_wp)*rd_ov_rv*(1 - s(ij, 1))
            if (qsat >= s(ij, 1)) then
                s(ij, 2) = 0.0_wp
            else
                s(ij, 2) = s(ij, 1) - qsat
            end if
            if (dsmooth > 0.0_wp) then
                qsat = qsat/(1 - s(ij, 1))
                dqldqt = 1.0_wp + qsat
                ! this is the real qsat
                qsat = qsat/(1.0_wp + qsat)
                dsmooth_loc = dsmooth*qsat
                s(ij, 2) = dsmooth_loc*dqldqt &
                           *log(exp((s(ij, 1) - qsat)/dsmooth_loc) + 1.0_wp)
            end if
        end do

        return
    end subroutine THERMO_AIRWATER_PT

!########################################################################
!# Calculate T and liquid content from rho, e and water content using caloric equation of state.
!########################################################################
    subroutine THERMO_AIRWATER_RE(ijmax, s, e, rho, T, dqldqt)
        integer(wi), intent(in) :: ijmax
        real(wp), intent(in) :: e(ijmax), rho(ijmax)
        real(wp), intent(out) :: T(ijmax)
        real(wp), intent(inout) :: s(ijmax, *), dqldqt(ijmax)

        ! -------------------------------------------------------------------
        integer(wi) i, inr, nrmax, ipsat
        real(wp) HEAT_CAPACITY_LD, HEAT_CAPACITY_LV, HEAT_CAPACITY_VD
        real(wp) B_LOC(10), FUN, DER, B_LOC_CONST_2, B_LOC_CONST_3
        real(wp) alpha, dummy1, dummy2

        ! ###################################################################
        NEWTONRAPHSON_ERROR = 0.0_wp
        ! maximum number of iterations in Newton-Raphson
        nrmax = 3

        ! reference case q_l = 0
        do ij = 1, ijmax
            T(ij) = e(ij) - (Ld + s(ij, 1)*Ldv)
            CPMEAN = Cd + s(ij, 1)*Cdv
            RMEAN = Rd + s(ij, 1)*Rdv
            T(ij) = (e(ij) - (Ld + s(ij, 1)*Ldv))/(CPMEAN - RMEAN*CRATIO_INV)
        end do

        ! -------------------------------------------------------------------
        ! calculate saturation specific humidity q_s(\rho,e)
        ! -------------------------------------------------------------------
        if (dsmooth <= 0.0_wp) then
            ! calculate saturation specific humidity, in array s(1,2).
            ! THERMO_POLYNOMIAL_PSAT is duplicated here to avoid array calls
            do i = 1, ijmax
                psat = 0.0_wp
                do ipsat = NPSAT, 1, -1
                    psat = psat*T(i) + THERMO_PSAT(ipsat)
                end do
                s(i, 2) = psat/(rho(i)*T(i)*WGHT_INV(1))
            end do

        else
            ! initialize homogeneous data
            HEAT_CAPACITY_VD = Rdv*CRATIO_INV - Cdv
            do i = 1, 9
                B_LOC(i) = -THERMO_PSAT(i)*Ldv
            end do
            B_LOC(10) = 0.0_wp
            do i = 2, 10
                B_LOC(i) = B_LOC(i) + THERMO_PSAT(i - 1)*HEAT_CAPACITY_VD
            end do
            B_LOC_CONST_2 = B_LOC(2)
            B_LOC_CONST_3 = B_LOC(3)

            ! loop on all points
            do i = 1, ijmax
                B_LOC(2) = B_LOC_CONST_2 + rho(i)*WGHT_INV(1)*(e(i) - Lv)
                B_LOC(3) = B_LOC_CONST_3 - rho(i)*WGHT_INV(1)*(Cd - GRATIO*WGHT_INV(2))
                ! Newton-Raphson
                t_loc = T(i)
                do inr = 1, nrmax
                    FUN = B_LOC(10)
                    DER = 0.0_wp
                    do is = 9, 1, -1
                        FUN = FUN*t_loc + B_LOC(is)
                        DER = DER*t_loc + B_LOC(is + 1)*real(is, wp)
                    end do
                    t_loc = t_loc - FUN/DER
                end do
                NEWTONRAPHSON_ERROR = max(NEWTONRAPHSON_ERROR, abs(FUN/DER)/t_loc)
                ! calculate saturation specific humidity, in array s(1,2).
                ! THERMO_POLYNOMIAL_PSAT is duplicated here to avoid routine calls
                psat = 0.0_wp
                do ipsat = NPSAT, 1, -1
                    psat = psat*t_loc + THERMO_PSAT(ipsat)
                end do
                s(i, 2) = psat/(rho(i)*t_loc*WGHT_INV(1))

                ! calculate dqldqt
                qsat = s(i, 2)
                alpha = -Lvl - (Cvl + GRATIO*WGHT_INV(1))*t_loc
                alpha = alpha/(t_loc*GRATIO*WGHT_INV(1)) - 1.0_wp

                dummy1 = -Ldl - (Cdl + GRATIO*WGHT_INV(2))*t_loc
                dummy1 = dummy1*qsat*alpha

                dummy2 = -Ldv + t_loc*GRATIO*WGHT_INV(1)*alpha*(alpha + 1.0_wp)
                dummy2 = dummy2*qsat + e(i) - Lv

                dqldqt(i) = 1.0_wp - dummy1/dummy2
            end do

        end if

        ! -------------------------------------------------------------------
        ! calculate final T and q_l
        ! -------------------------------------------------------------------
        ! initialize homogeneous data
        HEAT_CAPACITY_LV = Cvl + GRATIO*WGHT_INV(1)
        HEAT_CAPACITY_LD = Cdl + GRATIO*WGHT_INV(2)
        HEAT_CAPACITY_VD = HEAT_CAPACITY_LD - HEAT_CAPACITY_LV
        do i = 1, 9
            B_LOC(i) = THERMO_PSAT(i)*Lvl
        end do
        B_LOC(10) = 0.0_wp
        do i = 2, 10
            B_LOC(i) = B_LOC(i) + THERMO_PSAT(i - 1)*HEAT_CAPACITY_LV
        end do
        B_LOC_CONST_2 = B_LOC(2)
        B_LOC_CONST_3 = B_LOC(3)

        ! loop on all points
        do i = 1, ijmax
            qsat = s(i, 2)

            if (qsat >= s(i, 1)) then
                s(i, 2) = 0.0_wp
                if (dsmooth > 0.0_wp) then
                    dsmooth_loc = dsmooth*qsat
                    s(i, 2) = dsmooth_loc*dqldqt(i) &
                              *log(exp((s(i, 1) - qsat)/dsmooth_loc) + 1.0_wp)
                    ! change T consistently
                    T(i) = (e(i) - s(i, 1)*Ldv - Lv - s(i, 2)*Lvl)/ &
                           (s(i, 1)*HEAT_CAPACITY_VD + Cd - GRATIO*WGHT_INV(2) + s(i, 2)*HEAT_CAPACITY_LV)
                end if

                ! if q_s < q_t, then we have to repeat calculation of T
            else
                B_LOC(2) = B_LOC_CONST_2 + rho(i)*WGHT_INV(1)* &
                           (e(i) - s(i, 1)*(THERMO_AI(6, 1, 3) - THERMO_AI(6, 1, 2)) - THERMO_AI(6, 1, 2))
                B_LOC(3) = B_LOC_CONST_3 - rho(i)*WGHT_INV(1)* &
                           (s(i, 1)*HEAT_CAPACITY_LD + THERMO_AI(1, 1, 2) - GRATIO*WGHT_INV(2))
                ! IF ( dsmooth .GT. 0.0_wp ) THEN
                ! dsmooth_loc = dsmooth*qsat
                ! alpha =-( dsmooth_loc*LOG(EXP((s(i,1)-qsat)/dsmooth_loc)+1.0_wp)
                ! $                 -(s(i,1)-qsat) )*dqldqt(i)
                ! B_LOC(2) = B_LOC(2) + rho(i)*WGHT_INV(1)*alpha*Lvl
                ! B_LOC(3) = B_LOC(3) + rho(i)*WGHT_INV(1)*alpha*HEAT_CAPACITY_LV
                ! ENDIF

                ! Newton-Raphson
                do inr = 1, nrmax
                    FUN = B_LOC(10)
                    DER = 0.0_wp
                    do is = 9, 1, -1
                        FUN = FUN*T(i) + B_LOC(is)
                        DER = DER*T(i) + B_LOC(is + 1)*real(is, wp)
                    end do
                    T(i) = T(i) - FUN/DER
                end do
                NEWTONRAPHSON_ERROR = max(NEWTONRAPHSON_ERROR, abs(FUN/DER)/T(i))

                ! calculate saturation specific humidity, in array s(1,2).
                ! THERMO_POLYNOMIAL_PSAT is duplicated here to avoid routine calls
                psat = 0.0_wp
                do ipsat = NPSAT, 1, -1
                    psat = psat*T(i) + THERMO_PSAT(ipsat)
                end do
                s(i, 2) = psat/(rho(i)*T(i)*WGHT_INV(1))

                ! liquid content
                s(i, 2) = s(i, 1) - s(i, 2)
                if (dsmooth > 0.0_wp) then
                    dsmooth_loc = dsmooth*qsat
                    s(i, 2) = dsmooth_loc*dqldqt(i) &
                              *log(exp((s(i, 1) - qsat)/dsmooth_loc) + 1.0_wp) &
                              + s(i, 2) - dqldqt(i)*(s(i, 1) - qsat)
                    ! change T consistently
                    T(i) = (e(i) - s(i, 1)*Ldv - Lv - s(i, 2)*Lvl)/ &
                           (s(i, 1)*HEAT_CAPACITY_VD + Cd - GRATIO*WGHT_INV(2) &
                            + s(i, 2)*HEAT_CAPACITY_LV)
                end if
            end if

        end do

        return
    end subroutine THERMO_AIRWATER_RE

    !########################################################################
    !# Calculating the equilibrium liquid from the mixture fraction and the
    !# enthalpy deviation according to the linearized equilibrium thermodynamics
    !########################################################################
    subroutine THERMO_AIRWATER_LINEAR(ijmax, s, l)
        integer(wi), intent(IN) :: ijmax
        real(wp), dimension(ijmax, inb_scal), intent(IN) :: s     ! chi, psi
        real(wp), dimension(ijmax), intent(OUT) :: l     ! normalized liquid

        ! -------------------------------------------------------------------
        real(wp) dummy2

        ! ###################################################################
        ! Calculating \xi
        if (inb_scal == 1) then
            l = 1.0_wp + thermo_param(1)*s(:, 1)
        else
            l = 1.0_wp + thermo_param(1)*s(:, 1) + thermo_param(2)*s(:, 2)
        end if

        ! Calculating liquid; \xi is overwritten in this routine
        if (abs(thermo_param(inb_scal + 1)) < small_wp) then
            do ij = 1, ijmax
                l(ij) = max(l(ij), 0.0_wp)
            end do

        else
            dummy = thermo_param(inb_scal + 1)
            dummy2 = 1.0_wp/dummy
            !     l = dummy *LOG( EXP(dummy2 *l) +1.0_wp )
            do ij = 1, ijmax
                l(ij) = dummy*log(exp(dummy2*l(ij)) + 1.0_wp)
            end do

        end if

        return
    end subroutine THERMO_AIRWATER_LINEAR

    !########################################################################
    !########################################################################
    subroutine THERMO_AIRWATER_LINEAR_SOURCE(ijmax, s, xi, der1, der2)
        integer(wi), intent(IN) :: ijmax
        real(wp), intent(IN) :: s(ijmax, inb_scal)        ! chi, psi
        real(wp), intent(OUT) :: xi(ijmax), der1(ijmax), der2(ijmax)

        ! ###################################################################
        if (inb_scal == 1) then
            xi(:) = 1.0_wp + thermo_param(1)*s(:, 1)
        else
            xi(:) = 1.0_wp + thermo_param(1)*s(:, 1) + thermo_param(2)*s(:, 2)
        end if

        if (abs(thermo_param(inb_scal + 1)) < small_wp) then
            der2(:) = big_wp

            do ij = 1, ijmax
                if (xi(ij) <= 0.0_wp) then
                    der1(ij) = 0.0_wp
                else
                    der1(ij) = 1.0_wp
                end if
            end do

        else
            ! Formulation in terms of the exponential might lead to NAN because of too large numbers!
            !     dummy=-1.0_wp/thermo_param(inb_scal+1)
            !     der1 = 1.0_wp / ( 1.0_wp + EXP(dummy *xi) )

            !     der2 = (der1-1.0_wp) *der1 *dummy

            dummy = 0.5_wp/thermo_param(inb_scal + 1)
            der1(:) = 0.5_wp*(1.0_wp + tanh(dummy*xi(:)))

            dummy = 1.0_wp/thermo_param(inb_scal + 1)
            der2(:) = (1.0_wp - der1(:))*der1(:)*dummy

        end if

        return
    end subroutine THERMO_AIRWATER_LINEAR_SOURCE

end module THERMO_AIRWATER
