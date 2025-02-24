#include "dns_const.h"

! # Solvers for linear ordinary differential equations with constant coefficients
module OPR_ODES
    use TLab_Constants, only: wp, wi, BCS_MIN, BCS_MAX, BCS_BOTH
    use FDM, only: fdm_dt, FDM_COM6_JACOBIAN, FDM_COM6_DIRECT, FDM_COM6_JACOBIAN_PENTA
    use FDM_Integral, only: fdm_integral_dt
    use FDM_MatMul
    implicit none
    private

    integer(wi) i
    real(wp) dummy, lambda
    real(wp), dimension(:), pointer :: a, b, c, d, e, g, h, ep, em

    integer, parameter :: i1 = 1, i2 = 2

    ! First-order ODEs
    public :: OPR_Integral1             ! Using derived type grid_dp to pass information, typically initialized for lambda = 0
    ! public :: OPR_ODE1                ! Passing lu and rhs directly, which could be initialized for lambda /= 0

    ! Second-order ODEs
    ! public :: OPR_ODE2
    public :: OPR_ODE2_1_SINGULAR_DD    ! Reducing problem to a system of 2 first-order equations
    public :: OPR_ODE2_1_SINGULAR_DN
    public :: OPR_ODE2_1_SINGULAR_DN_New
    public :: OPR_ODE2_1_SINGULAR_ND
    public :: OPR_ODE2_1_SINGULAR_ND_New
    public :: OPR_ODE2_1_SINGULAR_NN
    public :: OPR_ODE2_1_SINGULAR_NN_New
    public :: OPR_ODE2_1_REGULAR_DD     ! Dirichlet/Dirichlet boundary conditions
    public :: OPR_ODE2_1_REGULAR_NN     ! Neumann/Neumann boundary conditions

contains
!########################################################################
!#
!#     u'_i + \lamba u_i = f_i  N-1 eqns
!#     u_1 or u_N given         1   eqn
!#     Au' = Bu                 N   eqns
!#
!# See FDM_Int1_Initialize. Similar to OPR_PARTIAL1
!# I wonder if this one and OPR_PARTIAL1 should be in FDM module.
!########################################################################
    subroutine OPR_Integral1(nlines, fdmi, f, result, wrk2d, du_boundary)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi
        real(wp), intent(in) :: f(nlines, size(fdmi%lhs, 1))
        real(wp), intent(inout) :: result(nlines, size(fdmi%lhs, 1))   ! contains bcs
        real(wp), intent(inout) :: wrk2d(nlines, 2)
        real(wp), intent(out), optional :: du_boundary(nlines)

        integer(wi) :: nx
        integer(wi) :: idl, ndl, idr, ndr, ic

        ! ###################################################################
        nx = size(f, 2)

        ndl = size(fdmi%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(fdmi%rhs, 2)
        idr = ndr/2 + 1

        select case (fdmi%bc)
        case (BCS_MIN)
            result(:, nx) = f(:, nx)
        case (BCS_MAX)
            result(:, 1) = f(:, 1)
        end select

        select case (size(fdmi%rhs, 2))
        case (3)
            call MatMul_3d(nx, nlines, fdmi%rhs(:, 1), fdmi%rhs(:, 3), f, result, &
                           BCS_BOTH, rhs_b=fdmi%rhs_b(1:3, 0:3), rhs_t=fdmi%rhs_t(0:2, 1:4), bcs_b=wrk2d(:, 1), bcs_t=wrk2d(:, 2))
        case (5)
            call MatMul_5d(nx, nlines, fdmi%rhs(:, 1), fdmi%rhs(:, 2), fdmi%rhs(:, 4), fdmi%rhs(:, 5), f, result, &
                           BCS_BOTH, rhs_b=fdmi%rhs_b(1:4, 0:5), rhs_t=fdmi%rhs_t(0:3, 1:6), bcs_b=wrk2d(:, 1), bcs_t=wrk2d(:, 2))
        end select

        select case (ndl)
        case (3)
            call TRIDSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), result(:, 2:))
        case (5)
            call PENTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), result(:, 2:))
        case (7)
            call HEPTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), fdmi%lhs(2:, 6), fdmi%lhs(2:, 7), result(:, 2:))
        end select

        if (any([BCS_MAX] == fdmi%bc)) then
            result(:, 1) = wrk2d(:, 1)
            do ic = 1, idl - 1
                result(:, 1) = result(:, 1) + fdmi%lhs(1, idl + ic)*result(:, 1 + ic)
            end do
            result(:, 1) = result(:, 1) + fdmi%lhs(1, 1)*result(:, 1 + ic)

            if (present(du_boundary)) then      ! calculate u'n
            end if

        end if

        if (any([BCS_MIN] == fdmi%bc)) then
            result(:, nx) = wrk2d(:, 2)
            do ic = 1, idl - 1
                result(:, nx) = result(:, nx) + fdmi%lhs(nx, idl - ic)*result(:, nx - ic)
            end do
            result(:, nx) = result(:, nx) + fdmi%lhs(nx, ndl)*result(:, nx - ic)

            if (present(du_boundary)) then      ! calculate u'1
                du_boundary(:) = fdmi%lhs(1, idl)*result(:, 1)
                do ic = 1, idl - 1
                    du_boundary(:) = du_boundary(:) + fdmi%lhs(1, idl + ic)*result(:, 1 + ic)
                end do
                do ic = 1, idr - 1
                    du_boundary(:) = du_boundary(:) + fdmi%rhs(1, idr + ic)*f(:, 1 + ic)
                end do
            end if

        end if

        return
    end subroutine OPR_Integral1

!########################################################################
!#
!#     u'_i + \lamba u_i = f_i  N-1 eqns
!#     u_1 or u_N given         1   eqn
!#     Au' = Bu                 N   eqns
!#
!#  Same as before, but passing information explicitly instead of using derived type fdm_dt
!########################################################################
    ! subroutine OPR_ODE1(nlines, lu, rhs, rhs_b, rhs_t, f, result, ibc)
    !     return
    ! end subroutine OPR_ODE1

!########################################################################
!#
!#     u''_i + \lamba u_i = f_i             N-2 eqns
!#     u_1, u'_1, u_N, or u'_N given        2   eqn
!#     Au'' = Bu                            N   eqns
!#
!########################################################################
    ! subroutine OPR_ODE2(nlines, lu, rhs, rhs_b, rhs_t, f, result, ibc)
    !     return
    ! end subroutine OPR_ODE2

!########################################################################
!#
!#     (u')'_i - \lambda^2 u_i = f_i        N - 2 equtions (singular case when \lambda = 0)
!#     u_1, u'_1, u_N, or u'_N given        2   eqn
!#     Au' = Bu                             N   eqns
!#     A(u')' = Bu'                         N   eqns
!#
!########################################################################
!########################################################################
!#     (u')' = f                            singular cases
!Dirichlet/Neumann boundary conditions
!########################################################################
    subroutine OPR_ODE2_1_SINGULAR_DN(imode_fdm, nx, nlines, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi), intent(in) :: imode_fdm, nx, nlines
        real(wp), intent(in) :: dx(nx)
        real(wp), intent(out) :: u(nlines, nx)          ! solution
        real(wp), intent(inout) :: f(nlines, nx)        ! forcing
        real(wp), intent(inout) :: tmp1(nlines, nx)     ! first derivative
        real(wp), intent(inout) :: bcs(nlines, 2)
        real(wp), intent(inout), target :: wrk1d(nx, 9)

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 8)
        h => wrk1d(:, 9)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v in v' = f , v_nx given
! -----------------------------------------------------------------------
        f(:, 1) = 0.0_wp
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(nx, i2, a, b, c, d, e)
            call INT_C1N6_RHS(nx, nlines, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = a(1)*dx(1) ! for v^1
            call PENTADFS(nx - 1, a, b, c, d, e)

!   obtain v^0, array tmp1
            call PENTADSS(nx - 1, nlines, a, b, c, d, e, tmp1)
            tmp1(:, nx) = 0.0_wp
            do i = 1, nx
                tmp1(:, i) = tmp1(:, i) + bcs(:, 2) ! add v_N to free array bcs(:,2)
            end do

!   obtain v^1, array wrk1d(:,6)
            call PENTADSS(nx - 1, i1, a, b, c, d, e, wrk1d(1, 6))
            wrk1d(nx, 6) = 0.0_wp
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(nx, i2, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(nx, nlines, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = b(1)*dx(1) ! for v^1
            call HEPTADFS(nx - 1, a, b, c, d, e, g, h)

!   obtain v^0, array tmp1
            call HEPTADSS(nx - 1, nlines, a, b, c, d, e, g, h, tmp1)
            tmp1(:, nx) = 0.0_wp
            do i = 1, nx
                tmp1(:, i) = tmp1(:, i) + bcs(:, 2) ! add v_N to free array bcs(:,2)
            end do

!   obtain v^1, array wrk1d(:,6)
            call HEPTADSS(nx - 1, i1, a, b, c, d, e, g, h, wrk1d(1, 6))
            wrk1d(nx, 6) = 0.0_wp
        end if

! -----------------------------------------------------------------------
! solve for u in u' = v, u_1 given
! -----------------------------------------------------------------------
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(nx, i1, a, b, c, d, e)
            call INT_C1N6_RHS(nx, nlines, i1, dx, tmp1, u)
            call PENTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:))

!   obtain u^0
            call PENTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), u(1, 2))
            bcs(:, 2) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 2) = (bcs(:, 2) + c(1)*u(:, 1) + d(1)*u(:, 2) + e(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6_RHS(nx, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call PENTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + c(1)*wrk1d(1, 7) + d(1)*wrk1d(2, 7) + e(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(nx, i1, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(nx, nlines, i1, dx, tmp1, u)
            call HEPTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))

!   obtain u^0
            call HEPTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), u(1, 2))
            bcs(:, 2) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 2) = (bcs(:, 2) + d(1)*u(:, 1) + e(1)*u(:, 2) + g(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6M_RHS(nx, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call HEPTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + d(1)*wrk1d(1, 7) + e(1)*wrk1d(2, 7) + g(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        end if

! Constraint
        dummy = 1.0_wp/(dummy - wrk1d(1, 6))
        bcs(:, 2) = (tmp1(:, 1) - bcs(:, 2))*dummy

! Result
        do i = 1, nx
            u(:, i) = u(:, i) + bcs(:, 2)*wrk1d(i, 7) + bcs(:, 1)
            tmp1(:, i) = tmp1(:, i) + bcs(:, 2)*wrk1d(i, 6)
        end do

        return
    end subroutine OPR_ODE2_1_SINGULAR_DN

    subroutine OPR_ODE2_1_SINGULAR_DN_NEW(nlines, fdmi, u, f, bcs, tmp1, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(in) :: f(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(inout) :: u(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(in) :: bcs(nlines, 2)
        real(wp), intent(inout) :: tmp1(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(inout) :: wrk2d(nlines, 2)

        integer(wi) nx

! #######################################################################
        nx = size(fdmi(1)%lhs, 1)

! -----------------------------------------------------------------------
! solve for v in v' = f , v_n given
! -----------------------------------------------------------------------
        tmp1(:, nx) = bcs(:, 2)
        call OPR_Integral1(nlines, fdmi(BCS_MAX), f, tmp1, wrk2d)

! -----------------------------------------------------------------------
! solve for u in u' = v, u_1 given
! -----------------------------------------------------------------------
        u(:, 1) = bcs(:, 1)
        call OPR_Integral1(nlines, fdmi(BCS_MIN), tmp1, u, wrk2d)

        return
    end subroutine OPR_ODE2_1_SINGULAR_DN_NEW

!########################################################################
!Neumann/Dirichlet boundary conditions
!########################################################################
    subroutine OPR_ODE2_1_SINGULAR_ND(imode_fdm, nx, nlines, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, nx, nlines
        real(wp), dimension(nx) :: dx
        real(wp), dimension(nx, 9), target :: wrk1d
        real(wp), dimension(nlines, nx) :: u, f, tmp1
        real(wp), dimension(nlines, 2) :: bcs

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 8)
        h => wrk1d(:, 9)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v in v' = f , v_1 given
! -----------------------------------------------------------------------
        f(:, nx) = 0.0_wp
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(nx, i1, a, b, c, d, e)
            call INT_C1N6_RHS(nx, nlines, i1, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(nx, 6) = dx(nx); wrk1d(nx - 1, 6) = e(nx)*dx(nx) ! for v^1
            call PENTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:))

!   obtain v^0, array tmp1
            call PENTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), tmp1(1, 2))
            tmp1(:, 1) = 0.0_wp
            do i = 1, nx
                tmp1(:, i) = tmp1(:, i) + bcs(:, 1) ! add v_1 to free array bcs(:,1)
            end do

!   obtain v^1, array wrk1d(:,6)
            call PENTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 6))
            wrk1d(1, 6) = 0.0_wp
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(nx, i1, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(nx, nlines, i1, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(nx, 6) = dx(nx); wrk1d(nx - 1, 6) = g(nx)*dx(nx) ! for v^1
            call HEPTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))

!   obtain v^0, array tmp1
            call HEPTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), tmp1(1, 2))
            tmp1(:, 1) = 0.0_wp
            do i = 1, nx
                tmp1(:, i) = tmp1(:, i) + bcs(:, 1) ! add v_1 to free array bcs(:,1)
            end do

!   obtain v^1, array wrk1d(:,6)
            call HEPTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 6))
            wrk1d(1, 6) = 0.0_wp
        end if

! -----------------------------------------------------------------------
! solve for u in u' = v, u_N given
! -----------------------------------------------------------------------
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(nx, i2, a, b, c, d, e)
            call INT_C1N6_RHS(nx, nlines, i2, dx, tmp1, u)
            call PENTADFS(nx - 1, a, b, c, d, e)

!   obtain u^0
            call PENTADSS(nx - 1, nlines, a, b, c, d, e, u)
            bcs(:, 1) = u(:, nx); u(:, nx) = 0.0_wp
            bcs(:, 1) = (bcs(:, 1) + a(nx)*u(:, nx - 2) + b(nx)*u(:, nx - 1) + c(nx)*u(:, nx))/dx(nx) !u^(0)'_nx

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6_RHS(nx, i1, i2, dx, wrk1d(1, 6), wrk1d(1, 7))
            call PENTADSS(nx - 1, i1, a, b, c, d, e, wrk1d(1, 7))
            dummy = wrk1d(nx, 7); wrk1d(nx, 7) = 0.0_wp
            dummy = (dummy + a(nx)*wrk1d(nx - 2, 7) + b(nx)*wrk1d(nx - 1, 7) + c(nx)*wrk1d(nx, 7))/dx(nx) ! u^(1)'_nx
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(nx, i2, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(nx, nlines, i2, dx, tmp1, u)
            call HEPTADFS(nx - 1, a, b, c, d, e, g, h)

!   obtain u^0
            call HEPTADSS(nx - 1, nlines, a, b, c, d, e, g, h, u)
            bcs(:, 1) = u(:, nx); u(:, nx) = 0.0_wp
            bcs(:, 1) = (bcs(:, 1) + b(nx)*u(:, nx - 2) + c(nx)*u(:, nx - 1) + d(nx)*u(:, nx))/dx(nx) !u^(0)'_nx

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6M_RHS(nx, i1, i2, dx, wrk1d(1, 6), wrk1d(1, 7))
            call HEPTADSS(nx - 1, i1, a, b, c, d, e, g, h, wrk1d(1, 7))
            dummy = wrk1d(nx, 7); wrk1d(nx, 7) = 0.0_wp
            dummy = (dummy + b(nx)*wrk1d(nx - 2, 7) + c(nx)*wrk1d(nx - 1, 7) + d(nx)*wrk1d(nx, 7))/dx(nx) ! u^(1)'_nx
        end if

! Constraint
        dummy = 1.0_wp/(dummy - wrk1d(nx, 6))
        bcs(:, 1) = (tmp1(:, nx) - bcs(:, 1))*dummy

! Result
        do i = 1, nx
            u(:, i) = u(:, i) + bcs(:, 1)*wrk1d(i, 7) + bcs(:, 2)
            tmp1(:, i) = tmp1(:, i) + bcs(:, 1)*wrk1d(i, 6)
        end do

        return
    end subroutine OPR_ODE2_1_SINGULAR_ND

    subroutine OPR_ODE2_1_SINGULAR_ND_NEW(nlines, fdmi, u, f, bcs, tmp1, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(in) :: f(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(inout) :: u(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(in) :: bcs(nlines, 2)
        real(wp), intent(inout) :: tmp1(nlines, size(fdmi(1)%lhs, 1))
        ! real(wp), intent(inout) :: wrk1d(size(fdmi(1)%lhs, 2))
        real(wp), intent(inout) :: wrk2d(nlines, 2)

        integer(wi) nx

! #######################################################################
        nx = size(fdmi(1)%lhs, 1)

! #define f1(:) wrk1d(:,1)
! #define v1(:) wrk1d(:,2)
! #define u1(:) wrk1d(:,3)

! -----------------------------------------------------------------------
! solve for v in v' = f , v_1 given
! -----------------------------------------------------------------------
        ! f(:, nx) = 0.0_wp
        tmp1(:, 1) = bcs(:, 1)
        call OPR_Integral1(nlines, fdmi(BCS_MIN), f, tmp1, wrk2d)
!   solve for v1
        ! f1(:) = 0.0_wp; f1(nx) = 1.0_wp
        ! v1(1) = 0.0_wp
        ! call OPR_Integral1(nlines, fdmi(BCS_MIN), f1, v1, wrk2d)

! -----------------------------------------------------------------------
! solve for u in u' = v, u_n given
! -----------------------------------------------------------------------
        u(:, nx) = bcs(:, 2)
        call OPR_Integral1(nlines, fdmi(BCS_MAX), tmp1, u, wrk2d)
!   solve for u1
        ! u1(:, nx) = 0.0_wp
        ! call OPR_Integral1(nlines, fdmi(BCS_MAX), v1, u1, wrk2d)

! ! Constraint
!         dummy = 1.0_wp/(du1(nx) - v1(nx))
!         bcs(:, 1) = (tmp1(:, nx) - du(:, 1))*dummy

! ! Result
!         do i = 1, nx
!             u(:, i) = u(:, i) + bcs(:, 1)*u1(i)
!             tmp1(:, i) = tmp1(:, i) + bcs(:, 1)*v1(i)
!         end do

! #undef f1
! #undef v1
! #undef u1

        return
    end subroutine OPR_ODE2_1_SINGULAR_ND_NEW

!########################################################################
!Dirichlet/Dirichlet boundary conditions
!########################################################################
    subroutine OPR_ODE2_1_SINGULAR_DD(imode_fdm, nx, nlines, x, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, nx, nlines
        real(wp), dimension(nx) :: dx, x
        real(wp), dimension(nx, 9), target :: wrk1d
        real(wp), dimension(nlines, nx) :: u, f, tmp1
        real(wp), dimension(nlines, 3) :: bcs

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 8)
        h => wrk1d(:, 9)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v = u' in (u')' = f , u'_nx given
! -----------------------------------------------------------------------
        f(:, 1) = 0.0_wp
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(nx, i2, a, b, c, d, e)
            call INT_C1N6_RHS(nx, nlines, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = a(1)*dx(1) ! for v^1
            call PENTADFS(nx - 1, a, b, c, d, e)

!   obtain v^0, array tmp1
            call PENTADSS(nx - 1, nlines, a, b, c, d, e, tmp1)
            tmp1(:, nx) = 0.0_wp

!   obtain v^1, array wrk1d(:,6)
            call PENTADSS(nx - 1, i1, a, b, c, d, e, wrk1d(1, 6))
            wrk1d(nx, 6) = 0.0_wp
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(nx, i2, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(nx, nlines, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = b(1)*dx(1) ! for v^1
            call HEPTADFS(nx - 1, a, b, c, d, e, g, h)

!     obtain v^0, array tmp1
            call HEPTADSS(nx - 1, nlines, a, b, c, d, e, g, h, tmp1)
            tmp1(:, nx) = 0.0_wp

!     obtain v^1, array wrk1d(:,6)
            call HEPTADSS(nx - 1, i1, a, b, c, d, e, g, h, wrk1d(1, 6))
            wrk1d(nx, 6) = 0.0_wp
        end if
! -----------------------------------------------------------------------
! solve for u in u' v f, u_1 given
! -----------------------------------------------------------------------
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(nx, i1, a, b, c, d, e)
            call INT_C1N6_RHS(nx, nlines, i1, dx, tmp1, u)
            call PENTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:))

!   obtain u^0
            call PENTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), u(1, 2))
            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + c(1)*u(:, 1) + d(1)*u(:, 2) + e(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6_RHS(nx, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call PENTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + c(1)*wrk1d(1, 7) + d(1)*wrk1d(2, 7) + e(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(nx, i1, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(nx, nlines, i1, dx, tmp1, u)
            call HEPTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))

!   obtain u^0
            call HEPTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), u(1, 2))
            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + d(1)*u(:, 1) + e(1)*u(:, 2) + g(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^1, array wrk1d(:,7)
            call INT_C1N6M_RHS(nx, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call HEPTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + c(1)*wrk1d(1, 7) + e(1)*wrk1d(2, 7) + g(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        end if

! Constraint
        dummy = 1.0_wp/(dummy - wrk1d(1, 6))
        bcs(:, 3) = (tmp1(:, 1) - bcs(:, 3))*dummy

! BCs
        dummy = 1.0_wp/(x(nx) - x(1))
        bcs(:, 2) = (bcs(:, 2) - u(:, nx) - bcs(:, 3)*wrk1d(nx, 7) - bcs(:, 1))*dummy

! Result
        do i = 1, nx
            u(:, i) = u(:, i) + bcs(:, 3)*wrk1d(i, 7) + bcs(:, 2)*(x(i) - x(1)) + bcs(:, 1)
            tmp1(:, i) = tmp1(:, i) + bcs(:, 3)*wrk1d(i, 6) + bcs(:, 2)
        end do

        return
    end subroutine OPR_ODE2_1_SINGULAR_DD

    subroutine OPR_ODE2_1_SINGULAR_DD_NEW(nlines, fdmi, u, f, bcs, tmp1, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(in) :: f(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(inout) :: u(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(in) :: bcs(nlines, 3)
        real(wp), intent(inout) :: tmp1(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(inout) :: wrk2d(nlines, 2)

        integer(wi) nx

! #######################################################################
        nx = size(fdmi(1)%lhs, 1)

! -----------------------------------------------------------------------
! solve for v in v' = f , v_max = 0 given
! -----------------------------------------------------------------------
        tmp1(:, nx) = 0.0_wp
        call OPR_Integral1(nlines, fdmi(BCS_MAX), f, tmp1, wrk2d)

! -----------------------------------------------------------------------
! solve for u in u' = v, u_1 given
! -----------------------------------------------------------------------
        u(:, 1) = 0.0_wp
        call OPR_Integral1(nlines, fdmi(BCS_MIN), tmp1, u, wrk2d)

        return
    end subroutine OPR_ODE2_1_SINGULAR_DD_NEW

    !########################################################################
!Neumann/Neumann boundary conditions; must be compatible!
!########################################################################
    subroutine OPR_ODE2_1_SINGULAR_NN(imode_fdm, nx, nlines, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, nx, nlines
        real(wp), dimension(nx) :: dx
        real(wp), dimension(nx, 7), target :: wrk1d
        real(wp), dimension(nlines, nx) :: u, f, tmp1
        real(wp), dimension(nlines, 2) :: bcs

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 6)
        h => wrk1d(:, 7)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v in v' = f , v_1 given
! -----------------------------------------------------------------------
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS(nx, i1, a, b, c, d, e)
            call INT_C1N6_RHS(nx, nlines, i1, dx, f, tmp1)
            call PENTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:))
            call PENTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), tmp1(1, 2))
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS(nx, i1, a, b, c, d, e, g, h)
            call INT_C1N6M_RHS(nx, nlines, i1, dx, f, tmp1)
            call HEPTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))
            call HEPTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), tmp1(1, 2))
        end if
        tmp1(:, 1) = 0.0_wp
        do i = 1, nx
            tmp1(:, i) = tmp1(:, i) + bcs(:, 1) ! this step assumes compatible problem
        end do
! -----------------------------------------------------------------------
! solve for u in u' = v, u_1 given
! -----------------------------------------------------------------------
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
!   same l.h.s. as before
            call INT_C1N6_RHS(nx, nlines, i1, dx, tmp1, u)
!   same l.h.s. as before
            call PENTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), u(1, 2))
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
!   same l.h.s. as before
            call INT_C1N6M_RHS(nx, nlines, i1, dx, tmp1, u)
!   same l.h.s. as before
            call HEPTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), u(1, 2))
        end if
        u(:, 1) = 0.0_wp ! this integration constant is free and set to zero

        return
    end subroutine OPR_ODE2_1_SINGULAR_NN

    subroutine OPR_ODE2_1_SINGULAR_NN_NEW(nlines, fdmi, u, f, bcs, tmp1, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(in) :: f(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(inout) :: u(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(in) :: bcs(nlines, 2)
        real(wp), intent(inout) :: tmp1(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(inout) :: wrk2d(nlines, 2)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v in v' = f , v_1 given
! -----------------------------------------------------------------------
        tmp1(:, 1) = bcs(:, 1)   ! this step assumes compatible problem
        call OPR_Integral1(nlines, fdmi(BCS_MIN), f, tmp1, wrk2d)

! -----------------------------------------------------------------------
! solve for u in u' = v, u_1 given
! -----------------------------------------------------------------------
        u(:, 1) = 0.0_wp        ! this integration constant is free and set to zero
        call OPR_Integral1(nlines, fdmi(BCS_MIN), tmp1, u, wrk2d)

        return
    end subroutine OPR_ODE2_1_SINGULAR_NN_NEW

!########################################################################
!Neumann/Neumann boundary conditions
!########################################################################
    subroutine OPR_ODE2_1_REGULAR_NN(imode_fdm, nx, nlines, cst, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, nx, nlines
        real(wp) cst
        real(wp), dimension(nx) :: dx
        real(wp), dimension(nx, 12), target :: wrk1d
        real(wp), dimension(nlines, nx) :: u, f, tmp1
        real(wp), dimension(nlines, 3) :: bcs

! -----------------------------------------------------------------------
        real(wp) g_1, g_2, a_1, b_1, deti

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 11)
        h => wrk1d(:, 12)

        ep => wrk1d(:, 9)
        em => wrk1d(:, 10)

        lambda = sqrt(cst)

! #######################################################################
! -----------------------------------------------------------------------
! 1st step; solve for v^(0) and v^(1)
! -----------------------------------------------------------------------
        dummy = -lambda
        f(:, 1) = 0.0_wp
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS_E(nx, i2, dx, dummy, a, b, c, d, e, ep)
            call INT_C1N6_RHS(nx, nlines, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = a(1)*dx(1) ! for v^1
            call PENTADFS(nx - 1, a, b, c, d, e)

!   obtain e^(+), array ep
            call PENTADSS(nx - 1, i1, a, b, c, d, e, ep)

!   obtain v^(0), array tmp1
            call PENTADSS(nx - 1, nlines, a, b, c, d, e, tmp1)
            tmp1(:, nx) = 0.0_wp

!   obtain v^(1), array wrk1d(:,6)
            call PENTADSS(nx - 1, i1, a, b, c, d, e, wrk1d(1, 6))
            wrk1d(nx, 6) = 0.0_wp
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS_E(nx, i2, dx, dummy, a, b, c, d, e, g, h, ep)
            call INT_C1N6M_RHS(nx, nlines, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = b(1)*dx(1) ! for v^1
            call HEPTADFS(nx - 1, a, b, c, d, e, g, h)

!   obtain e^(+), array ep
            call HEPTADSS(nx - 1, i1, a, b, c, d, e, g, h, ep)

!   obtain v^(0), array tmp1
            call HEPTADSS(nx - 1, nlines, a, b, c, d, e, g, h, tmp1)
            tmp1(:, nx) = 0.0_wp

!   obtain v^(1), array wrk1d(:,6)
            call HEPTADSS(nx - 1, i1, a, b, c, d, e, g, h, wrk1d(1, 6))
            wrk1d(nx, 6) = 0.0_wp
        end if
! -----------------------------------------------------------------------
! 2nd step; solve for u^(0) and u^(1) and u^(2:)
! -----------------------------------------------------------------------
        dummy = lambda
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS_E(nx, i1, dx, dummy, a, b, c, d, e, em)
            call INT_C1N6_RHS(nx, nlines, i1, dx, tmp1, u)
            call PENTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:))

!   obtain e^(m), array em
            call PENTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), em(2:))
            g_1 = (c(1)*em(1) + d(1)*em(2) + e(1)*em(3))/dx(1)/lambda + 1.0_wp ! e^(-)'_1/\lambda + 1

!   obtain u^(2), array wrk1d(:,8)
            call INT_C1N6_RHS(nx, i1, i1, dx, ep, wrk1d(1, 8))
            call PENTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 8))
            g_2 = wrk1d(1, 8); wrk1d(1, 8) = 0.0_wp
            g_2 = (g_2 + c(1)*wrk1d(1, 8) + d(1)*wrk1d(2, 8) + e(1)*wrk1d(3, 8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1

!   obtain u^(0), array u
            call PENTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), u(1, 2))

!   BCs; intermediate step to save memory space
            dummy = 1.0_wp - wrk1d(nx, 8)*lambda
            bcs(:, 3) = tmp1(:, 1) - bcs(:, 1)
            bcs(:, 2) = u(:, nx)*lambda + bcs(:, 2)
            bcs(:, 1) = bcs(:, 2) + bcs(:, 3)*em(nx)! a_0 *det
            bcs(:, 2) = bcs(:, 2)*ep(1) + bcs(:, 3)*dummy   ! b_0 *det

            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + c(1)*u(:, 1) + d(1)*u(:, 2) + e(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^(1), array wrk1d(:,7)
            call INT_C1N6_RHS(nx, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call PENTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + c(1)*wrk1d(1, 7) + d(1)*wrk1d(2, 7) + e(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS_E(nx, i1, dx, dummy, a, b, c, d, e, g, h, em)
            call INT_C1N6M_RHS(nx, nlines, i1, dx, tmp1, u)
            call HEPTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))

!   obtain e^(m), array em
            call HEPTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), em(2:))
            g_1 = (d(1)*em(1) + e(1)*em(2) + g(1)*em(3))/dx(1)/lambda + 1.0_wp ! e^(-)'_1/\lambda + 1

!   obtain u^(2), array wrk1d(:,8)
            call INT_C1N6M_RHS(nx, i1, i1, dx, ep, wrk1d(1, 8))
            call HEPTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 8))
            g_2 = wrk1d(1, 8); wrk1d(1, 8) = 0.0_wp
            g_2 = (g_2 + d(1)*wrk1d(1, 8) + e(1)*wrk1d(2, 8) + g(1)*wrk1d(3, 8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1

!   obtain u^(0), array u
            call HEPTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), u(1, 2))

!   BCs; intermediate step to save memory space
            dummy = 1.0_wp - wrk1d(nx, 8)*lambda
            bcs(:, 3) = tmp1(:, 1) - bcs(:, 1)
            bcs(:, 2) = u(:, nx)*lambda + bcs(:, 2)
            bcs(:, 1) = bcs(:, 2) + bcs(:, 3)*em(nx)! a_0 *det
            bcs(:, 2) = bcs(:, 2)*ep(1) + bcs(:, 3)*dummy   ! b_0 *det

            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + d(1)*u(:, 1) + e(1)*u(:, 2) + g(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^(1), array wrk1d(:,7)
            call INT_C1N6M_RHS(nx, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call HEPTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + d(1)*wrk1d(1, 7) + e(1)*wrk1d(2, 7) + g(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        end if

! BCs; final step
        deti = 1.0_wp/(1.0_wp - ep(1)*em(nx) - wrk1d(nx, 8)*lambda) ! inverse of determinant
        bcs(:, 1) = bcs(:, 1)*deti ! a_0
        bcs(:, 2) = bcs(:, 2)*deti ! b_0
        a_1 = (lambda*wrk1d(nx, 7) + wrk1d(1, 6)*em(nx))*deti
        b_1 = (lambda*wrk1d(nx, 7)*ep(1) + wrk1d(1, 6)*(1.0_wp - wrk1d(nx, 8)*lambda))*deti

! Constraint
        dummy = 1.0_wp/(dummy + b_1*g_1 - wrk1d(1, 6) + a_1*g_2)
        bcs(:, 3) = (tmp1(:, 1) - bcs(:, 3) - bcs(:, 2)*g_1 - bcs(:, 1)*g_2)*dummy

        dummy = 1.0_wp/lambda
        bcs(:, 1) = bcs(:, 1) + bcs(:, 3)*a_1
        bcs(:, 2) = (bcs(:, 2) + bcs(:, 3)*b_1)*dummy

! Result
        do i = 1, nx
            u(:, i) = u(:, i) + bcs(:, 3)*wrk1d(i, 7) + bcs(:, 1)*wrk1d(i, 8) + bcs(:, 2)*em(i)
            tmp1(:, i) = tmp1(:, i) + bcs(:, 3)*wrk1d(i, 6) + bcs(:, 1)*ep(i) - lambda*u(:, i)
        end do

        return
    end subroutine OPR_ODE2_1_REGULAR_NN

!########################################################################
!Dirichlet/Dirichlet boundary conditions
!########################################################################
    subroutine OPR_ODE2_1_REGULAR_DD(imode_fdm, nx, nlines, cst, dx, u, f, bcs, tmp1, wrk1d)
        integer(wi) imode_fdm, nx, nlines
        real(wp) cst
        real(wp), dimension(nx) :: dx
        real(wp), dimension(nx, 12), target :: wrk1d
        real(wp), dimension(nlines, nx) :: u, f, tmp1
        real(wp), dimension(nlines, 3) :: bcs

! -----------------------------------------------------------------------
        real(wp) g_1, g_2, a_1, deti

! #######################################################################
        a => wrk1d(:, 1)
        b => wrk1d(:, 2)
        c => wrk1d(:, 3)
        d => wrk1d(:, 4)
        e => wrk1d(:, 5)
! additional diagonals
        g => wrk1d(:, 11)
        h => wrk1d(:, 12)

        ep => wrk1d(:, 9)
        em => wrk1d(:, 10)

        lambda = sqrt(cst)

! #######################################################################
! -----------------------------------------------------------------------
! 1st step; solve for v^(0) and v^(1)
! -----------------------------------------------------------------------
        dummy = -lambda
        f(:, 1) = 0.0_wp
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS_E(nx, i2, dx, dummy, a, b, c, d, e, ep)
            call INT_C1N6_RHS(nx, nlines, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = a(1)*dx(1) ! for v^1
            call PENTADFS(nx - 1, a, b, c, d, e)

!   obtain e^(+), array ep
            call PENTADSS(nx - 1, i1, a, b, c, d, e, ep)

!   obtain v^(0), array tmp1
            call PENTADSS(nx - 1, nlines, a, b, c, d, e, tmp1)
            tmp1(:, nx) = 0.0_wp

!   obtain v^(1), array wrk1d(:,6)
            call PENTADSS(nx - 1, i1, a, b, c, d, e, wrk1d(1, 6))
            wrk1d(nx, 6) = 0.0_wp
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS_E(nx, i2, dx, dummy, a, b, c, d, e, g, h, ep)
            call INT_C1N6M_RHS(nx, nlines, i2, dx, f, tmp1)
            wrk1d(:, 6) = 0.0_wp; wrk1d(1, 6) = dx(1); wrk1d(2, 6) = b(1)*dx(1) ! for v^1
            call HEPTADFS(nx - 1, a, b, c, d, e, g, h)

!   obtain e^(+), array ep
            call HEPTADSS(nx - 1, i1, a, b, c, d, e, g, h, ep)

!   obtain v^(0), array tmp1
            call HEPTADSS(nx - 1, nlines, a, b, c, d, e, g, h, tmp1)
            tmp1(:, nx) = 0.0_wp

!   obtain v^(1), array wrk1d(:,6)
            call HEPTADSS(nx - 1, i1, a, b, c, d, e, g, h, wrk1d(1, 6))
            wrk1d(nx, 6) = 0.0_wp
        end if

! -----------------------------------------------------------------------
! 2nd step; solve for u^(0) and u^(1) and u^(2)
! -----------------------------------------------------------------------
        dummy = lambda
        if (imode_fdm == FDM_COM6_JACOBIAN .or. imode_fdm == FDM_COM6_DIRECT) then
            call INT_C1N6_LHS_E(nx, i1, dx, dummy, a, b, c, d, e, em)
            call INT_C1N6_RHS(nx, nlines, i1, dx, tmp1, u)
            call PENTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:))

!   obtain e^(m), array em
            call PENTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), em(2:))
            g_1 = (c(1)*em(1) + d(1)*em(2) + e(1)*em(3))/dx(1)/lambda + 1.0_wp ! e^(-)'_1/\lambda + 1

!   obtain u^(2), array wrk1d(:,8)
            call INT_C1N6_RHS(nx, i1, i1, dx, ep, wrk1d(1, 8))
            call PENTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 8))
            g_2 = wrk1d(1, 8); wrk1d(1, 8) = 0.0_wp
            g_2 = (g_2 + c(1)*wrk1d(1, 8) + d(1)*wrk1d(2, 8) + e(1)*wrk1d(3, 8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1

!   obtain u^(0), array u
            call PENTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), u(1, 2))
            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + c(1)*u(:, 1) + d(1)*u(:, 2) + e(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^(1), array wrk1d(:,7)
            call INT_C1N6_RHS(nx, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call PENTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + c(1)*wrk1d(1, 7) + d(1)*wrk1d(2, 7) + e(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        elseif (imode_fdm == FDM_COM6_JACOBIAN_PENTA) then
            call INT_C1N6M_LHS_E(nx, i1, dx, dummy, a, b, c, d, e, g, h, em)
            call INT_C1N6M_RHS(nx, nlines, i1, dx, tmp1, u)
            call HEPTADFS(nx - 1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:))

!   obtain e^(m), array em
            call HEPTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), em(2:))
            g_1 = (d(1)*em(1) + e(1)*em(2) + g(1)*em(3))/dx(1)/lambda + 1.0_wp ! e^(-)'_1/\lambda + 1

!   obtain u^(2), array wrk1d(:,8)
            call INT_C1N6M_RHS(nx, i1, i1, dx, ep, wrk1d(1, 8))
            call HEPTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 8))
            g_2 = wrk1d(1, 8); wrk1d(1, 8) = 0.0_wp
            g_2 = (g_2 + d(1)*wrk1d(1, 8) + e(1)*wrk1d(2, 8) + g(1)*wrk1d(3, 8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1

!   obtain u^(0), array u
            call HEPTADSS(nx - 1, nlines, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), u(1, 2))
            bcs(:, 3) = u(:, 1); u(:, 1) = 0.0_wp
            bcs(:, 3) = (bcs(:, 3) + d(1)*u(:, 1) + e(1)*u(:, 2) + g(1)*u(:, 3))/dx(1) !u^(0)'_1

!   obtain u^(1), array wrk1d(:,7)
            call INT_C1N6M_RHS(nx, i1, i1, dx, wrk1d(1, 6), wrk1d(1, 7))
            call HEPTADSS(nx - 1, i1, a(2:), b(2:), c(2:), d(2:), e(2:), g(2:), h(2:), wrk1d(2, 7))
            dummy = wrk1d(1, 7); wrk1d(1, 7) = 0.0_wp
            dummy = (dummy + d(1)*wrk1d(1, 7) + e(1)*wrk1d(2, 7) + g(1)*wrk1d(3, 7))/dx(1) ! u^(1)'_1
        end if

! BCs; final step
        deti = 1.0_wp/wrk1d(nx, 8) ! inverse of determinant
        bcs(:, 2) = (bcs(:, 2) - u(:, nx) - bcs(:, 1)*em(nx))*deti ! a_0
! bcs(:,1) = bcs(:,1)                                       ! b_0/lambda
        a_1 = -wrk1d(nx, 7)*deti
! b_1 = 0.0_wp

! Constraint
        g_1 = g_1*lambda
        dummy = 1.0_wp/(dummy - wrk1d(1, 6) + a_1*g_2)
        bcs(:, 3) = (tmp1(:, 1) - bcs(:, 3) - bcs(:, 1)*g_1 - bcs(:, 2)*g_2)*dummy

        bcs(:, 2) = bcs(:, 2) + bcs(:, 3)*a_1 ! a = v_N
! bcs(:,1) = bcs(:,1)                ! b = u_1

! Result
        do i = 1, nx
            u(:, i) = u(:, i) + bcs(:, 3)*wrk1d(i, 7) + bcs(:, 2)*wrk1d(i, 8) + bcs(:, 1)*em(i)
            tmp1(:, i) = tmp1(:, i) + bcs(:, 3)*wrk1d(i, 6) + bcs(:, 2)*ep(i) - lambda*u(:, i)
        end do

        return
    end subroutine OPR_ODE2_1_REGULAR_DD

end module OPR_ODES
