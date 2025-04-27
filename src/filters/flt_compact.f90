#include "dns_const.h"

! Compact filters following Lele, JCP, 1992
module Filters_Compact
    use TLab_Constants, only: wp, wi
    use FLT_Base
    implicit none
    private

    ! constants for the cutoff filter
    real(wp), parameter :: C4_ALPHA = 0.6522474d0
    real(wp), parameter :: C4_BETA = 0.1702929d0
    real(wp), parameter :: C4_A = 0.9891856d0
    real(wp), parameter :: C4_BD2 = 0.66059d0    ! 1.321180 /2
    real(wp), parameter :: C4_CD2 = 0.1666774d0  ! 0.3333548 /2
    real(wp), parameter :: C4_DD2 = 0.679925d-3  ! 0.001359850 /2

    integer(wi) i

    public :: FLT_C4_LHS, FLT_C4_RHS_COEFFS, FLT_C4_RHS
    public :: FLT_C4_CUTOFF_LHS, FLT_C4_CUTOFF_RHS, FLT_C4P_CUTOFF_LHS, FLT_C4P_CUTOFF_RHS

contains
!########################################################################
!#
!# 4th-Order Compact Filter from Lele, J. Comp. Phys., 1992, eqn C.2.4
!#
!# uf_i + alpha*(uf_i-1 + uf_i+1) = a*u_i + b*(u_i-1 + u_i+1) + c*(u_i-2 + u_i+2)
!#
!########################################################################
    subroutine FLT_C4_LHS(imax, bcsimin, bcsimax, alpha, a, b, c)
        integer(wi), intent(in) :: imax, bcsimin, bcsimax
        real(wp), intent(in) :: alpha
        real(wp), intent(out) :: a(imax), b(imax), c(imax)

        a(:) = alpha
        b(:) = 1.0_wp
        c(:) = alpha

        if (bcsimin == DNS_FILTER_BCS_ZERO) then
            c(1) = 0.0_wp
        end if

        if (bcsimax == DNS_FILTER_BCS_ZERO) then
            a(imax) = 0.0_wp
        end if

        return
    end subroutine FLT_C4_LHS

    subroutine FLT_C4_RHS_COEFFS(imax, alpha, periodic, dx, cxi)
        integer(wi), intent(in) :: imax
        real(wp), intent(in) :: alpha
        logical, intent(in) :: periodic
        real(wp), dimension(imax), intent(IN) :: dx
        real(wp), dimension(imax, 5), intent(out) :: cxi

        ! -----------------------------------------------------------------------
        real(wp) ac_loc

        ! #######################################################################
        ! Calculate constants for RHS, interior point
        ! #######################################################################
        ac_loc = (5.0_wp + 6.0_wp*alpha)/8.0_wp

        i = 1
        cxi(i, 1) = (ac_loc - 1.0_wp)*dx(i)*dx(i + 1)*(dx(i + 1) + dx(i + 2))/ &
                    (dx(i)*(dx(i) + dx(i) + dx(i + 1))* &
                     (dx(i) + dx(i) + dx(i + 1) + dx(i + 2)))
        cxi(i, 2) = alpha &
                    + (1.0_wp - ac_loc)*(dx(i) + dx(1))**2*dx(i + 1)*(dx(i + 1) + dx(i + 2))/ &
                    (dx(i)*dx(i)*(dx(i) + dx(i + 1))*(dx(i) + dx(i) + dx(i + 1) + dx(i + 2))) &
                    + (ac_loc - 1.0_wp)*(dx(i) + dx(i))*dx(i + 1)*(dx(i + 1) + dx(i + 2))**2/ &
                    (dx(i)*(dx(i) + dx(i + 1))*(dx(i) + dx(i + 1) + dx(i + 2))* &
                     (dx(i) + dx(i) + dx(i + 1) + dx(i + 2)))
        cxi(i, 3) = ac_loc
        cxi(i, 4) = alpha + (1.0_wp - ac_loc)*dx(i) &
                    *(dx(i) + dx(1))*(dx(i + 1) + dx(i + 2))/ &
                    (dx(i + 2)*(dx(i) + dx(i + 1))*(dx(i) + dx(1) + dx(i + 1)))
        cxi(i, 5) = (ac_loc - 1.0_wp)*dx(i)*(dx(i) + dx(i))*dx(i + 1)/ &
                    (dx(i + 2)*(dx(i) + dx(i + 1) + dx(i + 2))*(dx(i) + dx(i) + &
                                                                dx(i + 1) + dx(i + 2)))

        do i = 2, imax - 2
            cxi(i, 1) = (ac_loc - 1.0_wp)*dx(i)*dx(i + 1)*(dx(i + 1) + dx(i + 2))/ &
                        (dx(i - 1)*(dx(i) + dx(i - 1) + dx(i + 1))* &
                         (dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 2)))
            cxi(i, 2) = alpha + (1.0_wp - ac_loc) &
                        *(dx(i) + dx(i - 1))**2*dx(i + 1)*(dx(i + 1) + dx(i + 2))/ &
                        (dx(i)*dx(i - 1)*(dx(i) + dx(i + 1))*(dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 2))) &
                        + (ac_loc - 1.0_wp)*(dx(i) + dx(i - 1))*dx(i + 1)*(dx(i + 1) + dx(i + 2))**2/ &
                        (dx(i)*(dx(i) + dx(i + 1))*(dx(i) + dx(i + 1) + dx(i + 2))* &
                         (dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 2)))
            cxi(i, 3) = ac_loc
            cxi(i, 4) = alpha + (1.0_wp - ac_loc) &
                        *dx(i)*(dx(i) + dx(i - 1))*(dx(i + 1) + dx(i + 2))/ &
                        (dx(i + 2)*(dx(i) + dx(i + 1))*(dx(i) + dx(i - 1) + dx(i + 1)))
            cxi(i, 5) = (ac_loc - 1.0_wp)*dx(i)*(dx(i) + dx(i - 1))*dx(i + 1)/ &
                        (dx(i + 2)*(dx(i) + dx(i + 1) + dx(i + 2))* &
                         (dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 2)))
        end do

        i = imax - 1
        cxi(i, 1) = (ac_loc - 1.0_wp)*dx(i)*dx(i + 1)*(dx(i + 1) + dx(i + 1))/ &
                    (dx(i - 1)*(dx(i) + dx(i - 1) + dx(i + 1))* &
                     (dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 1)))
        cxi(i, 2) = alpha &
                    + (1.0_wp - ac_loc)*(dx(i) + dx(i - 1))**2*dx(i + 1)*(dx(i + 1) + dx(i + 1))/ &
                    (dx(i)*dx(i - 1)*(dx(i) + dx(i + 1))*(dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 1))) &
                    + (ac_loc - 1.0_wp)*(dx(i) + dx(i - 1))*dx(i + 1)*(dx(i + 1) + dx(i + 1))**2/ &
                    (dx(i)*(dx(i) + dx(i + 1))*(dx(i) + dx(i + 1) + dx(i + 1))* &
                     (dx(i) + dx(i - 1) + dx(i + 1) + dx(i + 1)))
        cxi(i, 3) = ac_loc
        cxi(i, 4) = alpha + (1.0_wp - ac_loc) &
                    *dx(i)*(dx(i) + dx(i - 1))*(dx(i + 1) + dx(i + 1))/ &
                    (dx(i + 1)*(dx(i) + dx(i + 1))*(dx(i) + dx(i - 1) + dx(i + 1)))
        cxi(i, 5) = (ac_loc - 1.0_wp)*dx(i)*(dx(i) + dx(i - 1))*dx(i + 1)/ &
                    (dx(i + 1)*(dx(i) + dx(i + 1) + dx(i + 1))*(dx(i) + dx(i - 1) + &
                                                                dx(i + 1) + dx(i + 1)))

        i = imax
        cxi(i, 1) = (ac_loc - 1.0_wp)*dx(i)*dx(i)*(dx(i) + dx(i))/ &
                    (dx(i - 1)*(dx(i) + dx(i - 1) + dx(i))* &
                     (dx(i) + dx(i - 1) + dx(i) + dx(i)))
        cxi(i, 2) = alpha &
                    + (1.0_wp - ac_loc)*(dx(i) + dx(i - 1))**2*dx(i)*(dx(i) + dx(i))/ &
                    (dx(i)*dx(i - 1)*(dx(i) + dx(i))*(dx(i) + dx(i - 1) + dx(i) + dx(i))) &
                    + (ac_loc - 1.0_wp)*(dx(i) + dx(i - 1))*dx(i)*(dx(i) + dx(i))**2/ &
                    (dx(i)*(dx(i) + dx(i))*(dx(i) + dx(i) + dx(i))* &
                     (dx(i) + dx(i - 1) + dx(i) + dx(i)))
        cxi(i, 3) = ac_loc
        cxi(i, 4) = alpha + (1.0_wp - ac_loc)*dx(i)*(dx(i) + dx(i - 1))*(dx(i) + dx(i))/ &
                    (dx(i)*(dx(i) + dx(i))*(dx(i) + dx(i - 1) + dx(i)))
        cxi(i, 5) = (ac_loc - 1.0_wp)*dx(i)*(dx(i) + dx(i - 1))*dx(i)/ &
                    (dx(i)*(dx(i) + dx(i) + dx(i))*(dx(i) + dx(i - 1) + &
                                                    dx(i) + dx(i)))

        ! #######################################################################
        ! biased formulation at the BCs
        ! #######################################################################
        if (.not. periodic) then

            ! -----------------------------------------------------------------------
            i = 1
            ac_loc = (15.0_wp + alpha)/16.0_wp

            cxi(i, 1) = ac_loc
            cxi(i, 2) = alpha + (1.0_wp - ac_loc) &
                        *(dx(2) + dx(3))*(dx(2) + dx(3) + dx(4))* &
                        (dx(2) + dx(3) + dx(4) + dx(5))/(dx(3)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5)))
            cxi(i, 3) = (ac_loc - 1.0_wp) &
                        *dx(2)*(dx(2) + dx(3) + dx(4))*(dx(2) + dx(3) + dx(4) + dx(5))/ &
                        (dx(3)*dx(4)*(dx(4) + dx(5)))
            cxi(i, 4) = (1.0_wp - ac_loc)*dx(2)*(dx(2) + dx(3))*(dx(2) + dx(3) + dx(4) + dx(5))/ &
                        (dx(4)*dx(5)*(dx(3) + dx(4)))
            cxi(i, 5) = (ac_loc - 1.0_wp)*dx(2)*(dx(2) + dx(3))*(dx(2) + dx(3) + dx(4))/ &
                        (dx(5)*(dx(4) + dx(5))*(dx(3) + dx(4) + dx(5)))

            ! -----------------------------------------------------------------------
            i = 2
            ac_loc = (3.0_wp + 2.0_wp*alpha)/4.0_wp

            cxi(i, 1) = alpha + (1.0_wp &
                                 - ac_loc)*dx(3)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/ &
                        ((dx(2) + dx(3))*(dx(2) + dx(3) + dx(4))*dx(5)) &
                        + (-1.0_wp + ac_loc) &
                        *dx(3)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/ &
                        ((dx(2) + dx(3))*dx(5)*(dx(2) + dx(3) + dx(4) + dx(5)))
            cxi(i, 2) = ac_loc
            cxi(i, 3) = alpha + dx(2)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/ &
                        ((dx(2) + dx(3))*dx(4)*(dx(4) + dx(5))) &
                        - ac_loc*dx(2)*(dx(3) + dx(4))*(dx(3) + dx(4) + dx(5))/ &
                        ((dx(2) + dx(3))*dx(4)*(dx(4) + dx(5)))
            cxi(i, 4) = (-1.0_wp + ac_loc)*dx(2)*dx(3)*(dx(3) + dx(4) + dx(5))/ &
                        (dx(4)*(dx(2) + dx(3) + dx(4))*dx(5))
            cxi(i, 5) = (1.0_wp - ac_loc)*dx(2)*dx(3)*(dx(3) + dx(4))/ &
                        (dx(5)*(dx(4) + dx(5))*(dx(2) + dx(3) + dx(4) + dx(5)))

            ! -----------------------------------------------------------------------
            i = imax - 1
            ac_loc = (3.0_wp + 2.0_wp*alpha)/4.0_wp

            cxi(i, 5) = alpha + (1.0_wp - ac_loc)* &
                        dx(i)*(dx(i) + dx(i - 1))*(dx(i) + dx(i - 1) + dx(i - 2))/ &
                        ((dx(i) + dx(i + 1))*(dx(i) + dx(i - 1) + dx(i + 1))* &
                         (dx(i) + dx(i - 1) + dx(i - 2) + dx(i + 1)))
            cxi(i, 4) = ac_loc
            cxi(i, 3) = alpha + (1.0_wp - ac_loc)* &
                        (dx(i) + dx(i - 1))*(dx(i) + dx(i - 1) + dx(i - 2))*dx(i + 1)/ &
                        (dx(i - 1)*dx(i - 2)*(dx(i) + dx(i + 1))) &
                        + (-1.0_wp + ac_loc)*(dx(i) + dx(i - 1))* &
                        (dx(i) + dx(i - 1) + dx(i - 2))*dx(i + 1)/ &
                        (dx(i - 2)*(dx(i - 1) + dx(i - 2))*(dx(i) + dx(i + 1)))
            cxi(i, 2) = (-1.0_wp + ac_loc) &
                        *dx(i)*(dx(i) + dx(i - 1) + dx(i - 2))*dx(i + 1)/ &
                        (dx(i - 1)*dx(i - 2)*(dx(i) + dx(i - 1) + dx(i + 1)))
            cxi(i, 1) = (1.0_wp - ac_loc)*dx(i)*(dx(i) + dx(i - 1))*dx(i + 1)/ &
                        (dx(i - 2)*(dx(i - 1) + dx(i - 2))*(dx(i) + dx(i - 1) + dx(i - 2) + &
                                                            dx(i + 1)))

            ! -----------------------------------------------------------------------
            i = imax
            ac_loc = (15.0_wp + alpha)/16.0_wp

            cxi(i, 5) = ac_loc
            cxi(i, 4) = alpha + (1.0_wp - ac_loc) &
                        *(dx(i) + dx(i - 1))*(dx(i) + dx(i - 1) + dx(i - 2))* &
                        (dx(i) + dx(i - 1) + dx(i - 2) + dx(i - 3))/ &
                        (dx(i - 1)*(dx(i - 1) + dx(i - 2))*(dx(i - 1) + dx(i - 2) + dx(i - 3)))
            cxi(i, 3) = (ac_loc - 1.0_wp)*dx(i)*(dx(i) + dx(i - 1) + dx(i - 2))* &
                        (dx(i) + dx(i - 1) + dx(i - 2) + dx(i - 3))/ &
                        (dx(i - 1)*dx(i - 2)*(dx(i - 2) + dx(i - 3)))
            cxi(i, 2) = (1.0_wp - ac_loc)*dx(i) &
                        *(dx(i) + dx(i - 1))*(dx(i) + dx(i - 1) + dx(i - 2) + dx(i - 3))/ &
                        (dx(i - 2)*dx(i - 3)*(dx(i - 1) + dx(i - 2)))
            cxi(i, 1) = (ac_loc - 1.0_wp)*dx(i) &
                        *(dx(i) + dx(i - 1))*(dx(i) + dx(i - 1) + dx(i - 2))/ &
                        (dx(i - 3)*(dx(i - 2) + dx(i - 3))*(dx(i - 1) + dx(i - 2) + dx(i - 3)))

        end if

        return
    end subroutine FLT_C4_RHS_COEFFS

    subroutine FLT_C4_RHS(imax, jkmax, periodic, bcsimin, bcsimax, cxi, u, rhs)
        logical, intent(in) :: periodic
        integer(wi), intent(IN) :: imax, jkmax, bcsimin, bcsimax
        real(wp), dimension(jkmax, imax), intent(IN) :: u
        real(wp), dimension(jkmax, imax), intent(OUT) :: rhs
        real(wp), dimension(imax, 5), intent(IN) :: cxi

! -----------------------------------------------------------------------

! #######################################################################
! Set up right hand side and DO forward/backward substitution
! #######################################################################
        if (periodic) then
            rhs(:, 1) = cxi(1, 1)*u(:, imax - 1) + cxi(1, 2)*u(:, imax) + cxi(1, 3)*u(:, 1) + &
                       cxi(1, 4)*u(:, 2) + cxi(1, 5)*u(:, 3)

            rhs(:, 2) = cxi(2, 1)*u(:, imax) + cxi(2, 2)*u(:, 1) + cxi(2, 3)*u(:, 2) + &
                       cxi(2, 4)*u(:, 3) + cxi(2, 5)*u(:, 4)

            rhs(:, imax) = cxi(imax, 1)*u(:, imax - 2) + cxi(imax, 2)*u(:, imax - 1) + &
                          cxi(imax, 3)*u(:, imax) + cxi(imax, 4)*u(:, 1) + cxi(imax, 5)*u(:, 2)

            rhs(:, imax - 1) = cxi(imax - 1, 1)*u(:, imax - 3) + cxi(imax - 1, 2)*u(:, imax - 2) + &
                              cxi(imax - 1, 3)*u(:, imax - 1) + cxi(imax - 1, 4)*u(:, imax) + cxi(imax - 1, 5)*u(:, 1)

        else ! biased
            rhs(:, 1) = cxi(1, 1)*u(:, 1) + &
                       cxi(1, 2)*u(:, 2) + &
                       cxi(1, 3)*u(:, 3) + &
                       cxi(1, 4)*u(:, 4) + &
                       cxi(1, 5)*u(:, 5)
            rhs(:, 2) = cxi(2, 1)*u(:, 1) + &
                       cxi(2, 2)*u(:, 2) + &
                       cxi(2, 3)*u(:, 3) + &
                       cxi(2, 4)*u(:, 4) + &
                       cxi(2, 5)*u(:, 5)

            rhs(:, imax - 1) = cxi(imax - 1, 5)*u(:, imax) + &
                              cxi(imax - 1, 4)*u(:, imax - 1) + &
                              cxi(imax - 1, 3)*u(:, imax - 2) + &
                              cxi(imax - 1, 2)*u(:, imax - 3) + &
                              cxi(imax - 1, 1)*u(:, imax - 4)
            rhs(:, imax) = cxi(imax, 5)*u(:, imax) + &
                          cxi(imax, 4)*u(:, imax - 1) + &
                          cxi(imax, 3)*u(:, imax - 2) + &
                          cxi(imax, 2)*u(:, imax - 3) + &
                          cxi(imax, 1)*u(:, imax - 4)

            if (bcsimin == DNS_FILTER_BCS_ZERO) then ! No filter at i=1
                rhs(:, 1) = u(:, 1)
            end if
            if (bcsimax == DNS_FILTER_BCS_ZERO) then ! No filter at i=1
                rhs(:, imax) = u(:, imax)
            end if

        end if

        do i = 3, imax - 2
            rhs(:, i) = cxi(i, 1)*u(:, i - 2) + cxi(i, 2)*u(:, i - 1) + cxi(i, 3)*u(:, i) + &
                       cxi(i, 4)*u(:, i + 1) + cxi(i, 5)*u(:, i + 2)
        end do

        return
    end subroutine FLT_C4_RHS

!########################################################################
!#
!# 4th-Order Compact Filter from Lele, J. Comp. Phys., 1992, eqn C.2.10b
!#
!# uf_i + alpha*(uf_i-1 + uf_i+1) = a*u_i + b*(u_i-1 + u_i+1) + c*(u_i-2 + u_i+2)
!#
!########################################################################

    subroutine FLT_C4P_CUTOFF_LHS(imax, a, b, c, d, e)
        integer(wi), intent(in) :: imax
        real(wp), intent(out) :: a(imax), b(imax), c(imax), d(imax), e(imax)

        a(:) = C4_BETA
        b(:) = C4_ALPHA
        c(:) = 1.0_wp
        d(:) = C4_ALPHA
        e(:) = C4_BETA

        return
    end subroutine FLT_C4P_CUTOFF_LHS

    subroutine FLT_C4_CUTOFF_LHS(imax, a, b, c, d, e)
        integer(wi), intent(in) :: imax
        real(wp), intent(out) :: a(imax), b(imax), c(imax), d(imax), e(imax)

        a(:) = C4_BETA
        b(:) = C4_ALPHA
        c(:) = 1.0_wp
        d(:) = C4_ALPHA
        e(:) = C4_BETA

        a(1:3) = 0.0_wp; b(1:3) = 0.0_wp; d(1:3) = 0.0_wp; e(1:3) = 0.0_wp
        a(imax - 2:imax) = 0.0_wp; b(imax - 2:imax) = 0.0_wp; d(imax - 2:imax) = 0.0_wp; e(imax - 2:imax) = 0.0_wp

        return
    end subroutine FLT_C4_CUTOFF_LHS

    subroutine FLT_C4P_CUTOFF_RHS(imax, jkmax, u, rhs)
        integer(wi), intent(in) :: imax, jkmax
        real(wp), intent(in) :: u(jkmax, imax)
        real(wp), intent(out) :: rhs(jkmax, imax)

        integer i
        integer im3, im2, im1, ip1, ip2, ip3

        do i = 1, imax
            im3 = mod(i + imax - 4, imax) + 1   ! i-3
            im2 = mod(i + imax - 3, imax) + 1   ! i-2
            im1 = mod(i + imax - 2, imax) + 1   ! i-1
            ip1 = mod(i, imax) + 1              ! i+1
            ip2 = mod(i + 1, imax) + 1          ! i+2
            ip3 = mod(i + 2, imax) + 1          ! i+3
            rhs(:, i) = C4_BD2*(u(:, ip1) + u(:, im1)) + &
                        C4_CD2*(u(:, ip2) + u(:, im2)) + &
                        C4_DD2*(u(:, ip3) + u(:, im3)) + &
                        C4_A*u(:, i)
        end do

        return
    end subroutine FLT_C4P_CUTOFF_RHS

    subroutine FLT_C4_CUTOFF_RHS(imax, jkmax, u, rhs)
        integer(wi), intent(in) :: imax, jkmax
        real(wp), intent(in) :: u(jkmax, imax)
        real(wp), intent(out) :: rhs(jkmax, imax)

        ! rhs(:, 1:3) = u(:, 1:3)
        ! rhs(:, imax - 2:imax) = u(:, imax - 2:imax)

        rhs(:, 1) = (15.0_wp*u(:, 1) + 4.0_wp*u(:, 2) - 6.0_wp*u(:, 3) + 4.0_wp*u(:, 4) - u(:, 5))/16.0_wp
        rhs(:, 2) = (12.0_wp*u(:, 2) + u(:, 1) + 6.0_wp*u(:, 3) - 4.0_wp*u(:, 4) + u(:, 5))/16.0_wp
        rhs(:, 3) = (10.0_wp*u(:, 3) - u(:, 1) + 4.0_wp*u(:, 2) + 4.0_wp*u(:, 4) - u(:, 5))/16.0_wp

       rhs(:, imax - 2) = (10.0_wp*u(:, imax - 2) - u(:, imax) + 4.0_wp*u(:, imax - 1) + 4.0_wp*u(:, imax - 3) - u(:, imax - 4))/16.0_wp
       rhs(:, imax - 1) = (12.0_wp*u(:, imax - 1) + u(:, imax) + 6.0_wp*u(:, imax - 2) - 4.0_wp*u(:, imax - 3) + u(:, imax - 4))/16.0_wp
     rhs(:, imax) = (15.0_wp*u(:, imax) + 4.0_wp*u(:, imax - 1) - 6.0_wp*u(:, imax - 2) + 4.0_wp*u(:, imax - 3) - u(:, imax - 4))/16.0_wp

        do i = 4, imax - 3
            rhs(:, i) = C4_BD2*(u(:, i + 1) + u(:, i - 1)) + &
                        C4_CD2*(u(:, i + 2) + u(:, i - 2)) + &
                        C4_DD2*(u(:, i + 3) + u(:, i - 3)) + &
                        C4_A*u(:, i)
        end do

        return
    end subroutine FLT_C4_CUTOFF_RHS

end module Filters_Compact
