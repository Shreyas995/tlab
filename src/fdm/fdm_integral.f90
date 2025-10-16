#include "dns_error.h"

!########################################################################
! Boundary-value problems and integrals based on the compact schemes.
! Should we move this to OPR_ODE in operators?
!########################################################################

module FDM_Integral
    use TLab_Constants, only: wp, wi, efile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: imax, jmax, kmax, isize_line
    use FDM_MatMul
    use FDM_Derivative, only: fdm_derivative_dt
    use FDM_Base
    use Tlab_Type, only: fdm_integral_dt, fdm_integral_dt2
    use LinearDss

    implicit none
    private

    ! type, public :: fdm_integral_dt
    !     sequence
    !     integer mode_fdm                            ! original finite-difference method; only informative
    !     real(wp) :: lambda                          ! constant of the equation
    !     integer :: bc                               ! type of boundary condition, [ BCS_MIN, BCS_MAX ]
    !     real(wp) :: rhs_b(1:5, 0:7), rhs_t(0:4, 8)  ! # of diagonals is 7, # rows is 7/2+1
    !     real(wp), allocatable :: lhs(:, :)          ! Often overwritten to LU decomposition.
    !     real(wp), allocatable :: rhs(:, :)
    ! end type fdm_integral_dt
    ! This type is used in elliptic operators for different eigenvalues. This can lead to fragmented memory.
    ! One could use pointers instead of allocatable for lhs and rhs, and point the pointers to the
    ! corresponding memory space.

    public FDM_Int1_Initialize                      ! Prepare to solve u' +\lambda u = f
    public FDM_Int1_CreateSystem
    public FDM_Int1_Solve
    public FDM_Int1_Initialize_APU
    public FDM_Int1_CreateSystem_APU
    ! public FDM_Int1_Solve_APU

    public FDM_Int2_Initialize                      ! Prepare to solve (u')' - \lamba^2 u = f
    ! public FDM_Int2_CreateSystem
    public FDM_Int2_Solve
    public FDM_Int2_CreateSystem_APU
    public FDM_Int2_Initialize_APU
    public FDM_Int2_Solve_APU

contains
    !########################################################################
    !#
    !#     u'_i + \lamba u_i = f_i  N-1 eqns
    !#     u_1 or u_N given         1   eqn
    !#     Au' = Bu                 N   eqns
    !#
    !# starting from generic diagonal matrices A (lhs) and B (rhs).
    !#
    !# The system of N-1 eqns:
    !#
    !#                    (B + \lambda A)u = Af
    !#
    !# is established in this routine (see notes).
    !#
    !# System normalized s.t. RHS diagonals are O(1), LHS diagonals O(h^2)
    !# System normalized s.t. 1. upper-diagonal in B is 1 (except at boundaries)
    !#
    !########################################################################
    subroutine FDM_Int1_Initialize(x, g, lambda, ibc, fdmi)
        real(wp), intent(in) :: x(:)                    ! node positions
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda                  ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) nx, nd

        !########################################################################
        call FDM_Int1_CreateSystem(x, g, lambda, ibc, fdmi)

        ! LU decomposition
        nx = size(fdmi%lhs, 1)              ! # of grid points
        nd = size(fdmi%lhs, 2)              ! # of diagonals

        select case (nd)
        case (3)
            call TRIDFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3))
        case (5)
            call PENTADFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5))
        case (7)
            call HEPTADFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), fdmi%lhs(2:, 6), fdmi%lhs(2:, 7))
        end select

        return
    end subroutine FDM_Int1_Initialize

    subroutine FDM_Int1_Initialize_APU(k, i, x, g, lambda, ibc, fdmi)
        integer, intent(in) :: k, i
        real(wp), intent(in) :: x(:)                    ! node positions
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda                  ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt2), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) nx, nd

        !########################################################################
        call FDM_Int1_CreateSystem_APU(k, i, x, g, lambda, ibc, fdmi)

        ! LU decomposition
        nx = size(fdmi%lhs, 1)              ! # of grid points
        nd = size(fdmi%lhs, 4)              ! # of diagonals

        select case (nd)
        case (3)
            call TRIDFS(nx - 2, fdmi%lhs(2:nx-1, k, i, 1), fdmi%lhs(2:nx-1, k, i, 2), fdmi%lhs(2:nx-1, k, i, 3))
        case (5)
            call PENTADFS(nx - 2, fdmi%lhs(2:nx-1, k, i, 1), fdmi%lhs(2:nx-1, k, i, 2), fdmi%lhs(2:nx-1, k, i, 3), &
                          fdmi%lhs(2:nx-1, k, i, 4), fdmi%lhs(2:nx-1, k, i, 5))
        case (7)
            call HEPTADFS(nx - 2, fdmi%lhs(2:nx-1, k, i, 1), fdmi%lhs(2:nx-1, k, i, 2), fdmi%lhs(2:nx-1, k, i, 3), &
                          fdmi%lhs(2:nx-1, k, i, 4), fdmi%lhs(2:nx-1, k, i, 5), fdmi%lhs(2:nx-1, k, i, 6), fdmi%lhs(2:nx-1, k, i, 7))
        end select

        return
    end subroutine FDM_Int1_Initialize_APU

    !########################################################################
    !#####################################################################

    subroutine FDM_Int1_CreateSystem_APU(ik, ii, x, g, lambda, ibc, fdmi)
        integer, intent(in) :: ik, ii
        real(wp), intent(in) :: x(:)                    ! node positions
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda                  ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt2), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) idl, ndl, idr, ndr, ir, nx, i
        real(wp) dummy, rhsr_b(5, 0:7), rhsr_t(0:4, 8)

        ! ###################################################################
        ndl = g%nb_diag(1)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = g%nb_diag(2)
        idr = ndr/2 + 1             ! center diagonal in rhs
        nx = g%size                 ! # grid points

        ! check sizes
        if (abs(idl - idr) > 1) then
            call TLab_Write_ASCII(efile, __FILE__//'. lhs and rhs cannot differ by more than 2 diagonals.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        fdmi%mode_fdm = g%mode_fdm
        fdmi%lambda(ik, ii) = lambda
        fdmi%bc(ik, ii) = ibc

        ! if (allocated(fdmi%lhs)) deallocate (fdmi%lhs)
        ! if (allocated(fdmi%rhs)) deallocate (fdmi%rhs)
        ! allocate (fdmi%lhs(nx, ndr))
        ! allocate (fdmi%rhs(nx, ndl))

        ! -------------------------------------------------------------------
        ! new rhs diagonals (array A), independent of lambda
        fdmi%rhs(:, ik, ii, 1:ndl) = g%lhs(:, 1:ndl)

        call FDM_Bcs_Reduce_APU(nx, ndl, ndr, fdmi%bc(ik, ii), fdmi%rhs(1:nx, ik, ii, 1:ndl), g%rhs(:,1:ndr), rhsr_b, rhsr_t)

        fdmi%rhs_b(1:5, ik, ii, 0:7) = 0.0_wp
        fdmi%rhs_t(0:4, ik, ii, 1:8) = 0.0_wp
        select case (fdmi%bc(ik, ii))
        case (BCS_MIN)
            fdmi%rhs_b(1:idl + 1, ik, ii, 1:ndl) = fdmi%rhs(1:idl + 1, ik, ii, 1:ndl)
            do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
                fdmi%rhs_b(1 + ir, ik, ii, idl - ir) = -rhsr_b(1 + ir, idr - ir)
            end do

        case (BCS_MAX)
            fdmi%rhs_t(0:idl, ik, ii, 1:ndl) = fdmi%rhs(nx - idl:nx, ik, ii, 1:ndl)
            do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
                fdmi%rhs_t(idl - ir, ik, ii, idl + ir) = -rhsr_t(idr - ir, idr + ir)
            end do

        end select

        ! -------------------------------------------------------------------
        ! new lhs diagonals (array C = B + h \lambda A), dependent on lambda
        fdmi%lhs(:, ik, ii, :) = g%rhs(:, 1:ndr)

        fdmi%lhs(:, ik, ii, idr) = fdmi%lhs(:, ik, ii, idr) + lambda*g%lhs(:, idl)                ! center diagonal
        do i = 1, idl - 1                                                       ! off-diagonals
            fdmi%lhs(1 + i:nx, ik, ii, idr - i) = fdmi%lhs(1 + i:nx, ik, ii, idr - i) + lambda*g%lhs(1 + i:nx, idl - i)
            fdmi%lhs(1:nx - i, ik, ii, idr + i) = fdmi%lhs(1:nx - i, ik, ii, idr + i) + lambda*g%lhs(1:nx - i, idl + i)
        end do

        select case (fdmi%bc(ik, ii))
        case (BCS_MIN)
            fdmi%lhs(1:idr, ik, ii, 1:ndr) = rhsr_b(1:idr, 1:ndr)
            fdmi%lhs(1, ik, ii, idr + 1:idr + idl - 1) = fdmi%lhs(1, ik, ii, idr + 1:idr + idl - 1) - lambda*fdmi%rhs_b(1, ik, ii, idl + 1:ndl)
            ! fdmi%lhs(2:idr, 1:ndr) = rhsr_b(2:idr, 1:ndr)
            do ir = 1, idr - 1
                fdmi%lhs(1 + ir, ik, ii, idr - idl + 1:idr + idl - 1) = fdmi%lhs(1 + ir, ik, ii, idr - idl + 1:idr + idl - 1) + lambda*fdmi%rhs_b(1 + ir, ik, ii, 1:ndl)
            end do
        case (BCS_MAX)
            fdmi%lhs(nx - idr + 1:nx, ik, ii, 1:ndr) = rhsr_t(1:idr, 1:ndr)
            fdmi%lhs(nx, ik, ii, idr - idl + 1:idr - 1) = fdmi%lhs(nx, ik, ii, idr - idl + 1:idr - 1) - lambda*fdmi%rhs_t(idl, ik, ii, 1:idl - 1)
            ! fdmi%lhs(nx - idr + 1:nx - 1, 1:ndr) = rhsr_t(1:idr - 1, 1:ndr)
            do ir = 1, idr - 1
                fdmi%lhs(nx - ir, ik, ii, idr - idl + 1:idr + idl - 1) = fdmi%lhs(nx - ir, ik, ii, idr - idl + 1:idr + idl - 1) + lambda*fdmi%rhs_t(idl - ir, ik, ii, 1:ndl)
            end do
        end select

        ! -------------------------------------------------------------------
        ! normalization such that new central diagonal in rhs is 1
        do ir = 1, max(idr, idl + 1)
            dummy = 1.0_wp/fdmi%rhs(ir, ik, ii, idl)
            fdmi%rhs_b(ir, ik, ii, 0:ndl) = fdmi%rhs_b(ir, ik, ii, 0:ndl)*dummy

            dummy = 1.0_wp/fdmi%rhs(nx - ir + 1, ik, ii, idl)
            fdmi%rhs_t(idl - ir + 1, ik, ii, 1:ndl + 1) = fdmi%rhs_t(idl - ir + 1, ik, ii, 1:ndl + 1)*dummy

            dummy = 1.0_wp/fdmi%rhs(ir, ik, ii, idl)
            fdmi%rhs(ir, ik, ii, 1:ndl) = fdmi%rhs(ir, ik, ii, 1:ndl)*dummy
            fdmi%lhs(ir, ik, ii, 1:ndr) = fdmi%lhs(ir, ik, ii, 1:ndr)*dummy

            dummy = 1.0_wp/fdmi%rhs(nx - ir + 1, ik, ii, idl)
            fdmi%rhs(nx - ir + 1, ik, ii, 1:ndl) = fdmi%rhs(nx - ir + 1, ik, ii, 1:ndl)*dummy
            fdmi%lhs(nx - ir + 1, ik, ii, 1:ndr) = fdmi%lhs(nx - ir + 1, ik, ii, 1:ndr)*dummy

        end do

        ! interior points: normalization such that 1. upper-diagonal is 1
        do ir = max(idr, idl + 1) + 1, nx - max(idr, idl + 1)
            dummy = 1.0_wp/fdmi%rhs(ir, ik, ii, idl + 1)

            fdmi%rhs(ir, ik, ii, 1:ndl) = fdmi%rhs(ir, ik, ii, 1:ndl)*dummy
            fdmi%lhs(ir, ik, ii, 1:ndr) = fdmi%lhs(ir, ik, ii, 1:ndr)*dummy

        end do

        ! -------------------------------------------------------------------
        ! reducing system in the opposite end to account for the case of extended stencils
        ! to move it up, you need to recalculate the expression for p_1 and p_n because they assume division by a_11
        select case (fdmi%bc(ik, ii))
        case (BCS_MIN)
            call FDM_Bcs_Reduce_APU(ik, ndl, ndr, BCS_MAX, fdmi%lhs(1:nx, ik, ii, 1:ndl), fdmi%rhs(1:nx, ik, ii, 1:ndr), rhs_t=fdmi%rhs_t(0:4, ik, ii, 1:8))
        case (BCS_MAX)
            call FDM_Bcs_Reduce_APU(ik, ndl, ndr, BCS_MIN, fdmi%lhs(1:nx, ik, ii, 1:ndl), fdmi%rhs(1:nx, ik, ii, 1:ndr), rhs_b=fdmi%rhs_b(1:5, ik, ii, 0:7))
        end select

        return
    end subroutine FDM_Int1_CreateSystem_APU

    !########################################################################
    !########################################################################
    subroutine FDM_Int1_CreateSystem(x, g, lambda, ibc, fdmi)
        real(wp), intent(in) :: x(:)                    ! node positions
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda                  ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) i
        integer(wi) idl, ndl, idr, ndr, ir, nx
        real(wp) dummy, rhsr_b(5, 0:7), rhsr_t(0:4, 8)

        ! ###################################################################
        ndl = g%nb_diag(1)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = g%nb_diag(2)
        idr = ndr/2 + 1             ! center diagonal in rhs
        nx = g%size                 ! # grid points

        ! check sizes
        if (abs(idl - idr) > 1) then
            call TLab_Write_ASCII(efile, __FILE__//'. lhs and rhs cannot differ by more than 2 diagonals.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        fdmi%mode_fdm = g%mode_fdm
        fdmi%lambda = lambda
        fdmi%bc = ibc

        if (allocated(fdmi%lhs)) deallocate (fdmi%lhs)
        if (allocated(fdmi%rhs)) deallocate (fdmi%rhs)
        allocate (fdmi%lhs(nx, ndr))
        allocate (fdmi%rhs(nx, ndl))

        ! -------------------------------------------------------------------
        ! new rhs diagonals (array A), independent of lambda
        fdmi%rhs(:, :) = g%lhs(:, 1:ndl)

        call FDM_Bcs_Reduce(fdmi%bc, fdmi%rhs, g%rhs(:,1:ndr), rhsr_b, rhsr_t)

        fdmi%rhs_b = 0.0_wp
        fdmi%rhs_t = 0.0_wp
        select case (fdmi%bc)
        case (BCS_MIN)
            fdmi%rhs_b(1:idl + 1, 1:ndl) = fdmi%rhs(1:idl + 1, 1:ndl)
            do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
                fdmi%rhs_b(1 + ir, idl - ir) = -rhsr_b(1 + ir, idr - ir)
            end do

        case (BCS_MAX)
            fdmi%rhs_t(0:idl, 1:ndl) = fdmi%rhs(nx - idl:nx, 1:ndl)
            do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
                fdmi%rhs_t(idl - ir, idl + ir) = -rhsr_t(idr - ir, idr + ir)
            end do

        end select

        ! -------------------------------------------------------------------
        ! new lhs diagonals (array C = B + h \lambda A), dependent on lambda
        fdmi%lhs(:, :) = g%rhs(:, 1:ndr)

        fdmi%lhs(:, idr) = fdmi%lhs(:, idr) + lambda*g%lhs(:, idl)                ! center diagonal
        do i = 1, idl - 1                                                       ! off-diagonals
            fdmi%lhs(1 + i:nx, idr - i) = fdmi%lhs(1 + i:nx, idr - i) + lambda*g%lhs(1 + i:nx, idl - i)
            fdmi%lhs(1:nx - i, idr + i) = fdmi%lhs(1:nx - i, idr + i) + lambda*g%lhs(1:nx - i, idl + i)
        end do

        select case (fdmi%bc)
        case (BCS_MIN)
            fdmi%lhs(1:idr, 1:ndr) = rhsr_b(1:idr, 1:ndr)
            fdmi%lhs(1, idr + 1:idr + idl - 1) = fdmi%lhs(1, idr + 1:idr + idl - 1) - lambda*fdmi%rhs_b(1, idl + 1:ndl)
            ! fdmi%lhs(2:idr, 1:ndr) = rhsr_b(2:idr, 1:ndr)
            do ir = 1, idr - 1
                fdmi%lhs(1 + ir, idr - idl + 1:idr + idl - 1) = fdmi%lhs(1 + ir, idr - idl + 1:idr + idl - 1) + lambda*fdmi%rhs_b(1 + ir, 1:ndl)
            end do
        case (BCS_MAX)
            fdmi%lhs(nx - idr + 1:nx, 1:ndr) = rhsr_t(1:idr, 1:ndr)
            fdmi%lhs(nx, idr - idl + 1:idr - 1) = fdmi%lhs(nx, idr - idl + 1:idr - 1) - lambda*fdmi%rhs_t(idl, 1:idl - 1)
            ! fdmi%lhs(nx - idr + 1:nx - 1, 1:ndr) = rhsr_t(1:idr - 1, 1:ndr)
            do ir = 1, idr - 1
                fdmi%lhs(nx - ir, idr - idl + 1:idr + idl - 1) = fdmi%lhs(nx - ir, idr - idl + 1:idr + idl - 1) + lambda*fdmi%rhs_t(idl - ir, 1:ndl)
            end do
        end select

        ! -------------------------------------------------------------------
        ! normalization such that new central diagonal in rhs is 1
        do ir = 1, max(idr, idl + 1)
            dummy = 1.0_wp/fdmi%rhs(ir, idl)
            fdmi%rhs_b(ir, 0:ndl) = fdmi%rhs_b(ir, 0:ndl)*dummy

            dummy = 1.0_wp/fdmi%rhs(nx - ir + 1, idl)
            fdmi%rhs_t(idl - ir + 1, 1:ndl + 1) = fdmi%rhs_t(idl - ir + 1, 1:ndl + 1)*dummy

            dummy = 1.0_wp/fdmi%rhs(ir, idl)
            fdmi%rhs(ir, 1:ndl) = fdmi%rhs(ir, 1:ndl)*dummy
            fdmi%lhs(ir, 1:ndr) = fdmi%lhs(ir, 1:ndr)*dummy

            dummy = 1.0_wp/fdmi%rhs(nx - ir + 1, idl)
            fdmi%rhs(nx - ir + 1, 1:ndl) = fdmi%rhs(nx - ir + 1, 1:ndl)*dummy
            fdmi%lhs(nx - ir + 1, 1:ndr) = fdmi%lhs(nx - ir + 1, 1:ndr)*dummy

        end do

        ! interior points: normalization such that 1. upper-diagonal is 1
        do ir = max(idr, idl + 1) + 1, nx - max(idr, idl + 1)
            dummy = 1.0_wp/fdmi%rhs(ir, idl + 1)

            fdmi%rhs(ir, 1:ndl) = fdmi%rhs(ir, 1:ndl)*dummy
            fdmi%lhs(ir, 1:ndr) = fdmi%lhs(ir, 1:ndr)*dummy

        end do

        ! -------------------------------------------------------------------
        ! reducing system in the opposite end to account for the case of extended stencils
        ! to move it up, you need to recalculate the expression for p_1 and p_n because they assume division by a_11
        select case (fdmi%bc)
        case (BCS_MIN)
            call FDM_Bcs_Reduce(BCS_MAX, fdmi%lhs, fdmi%rhs, rhs_t=fdmi%rhs_t)
        case (BCS_MAX)
            call FDM_Bcs_Reduce(BCS_MIN, fdmi%lhs, fdmi%rhs, rhs_b=fdmi%rhs_b)
        end select

        return
    end subroutine FDM_Int1_CreateSystem

    !########################################################################
    !########################################################################
    ! Allow to pass separate rhs because this part does not depend on lambda
    subroutine FDM_Int1_Solve(nlines, fdmi, rhsi, f, result, wrk2d, du_boundary)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi
        real(wp), intent(in) :: rhsi(:, :)
        real(wp), intent(in) :: f(nlines, size(fdmi%lhs, 1))
        real(wp), intent(inout) :: result(nlines, size(fdmi%lhs, 1))   ! contains bcs
        real(wp), intent(inout) :: wrk2d(nlines, 2)
        real(wp), intent(out), optional :: du_boundary(nlines)

        ! -------------------------------------------------------------------
        integer(wi) :: nx
        integer(wi) :: idl, ndl, idr, ndr, ic

        ! ###################################################################
        nx = size(fdmi%lhs, 1)

        ndl = size(fdmi%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(rhsi, 2)
        idr = ndr/2 + 1

        select case (fdmi%bc)
        case (BCS_MIN)
            result(:, nx) = f(:, nx)
        case (BCS_MAX)
            result(:, 1) = f(:, 1)
        end select

        select case (ndr)
        case (3)
            call MatMul_3d(rhsi(:, 1:3), f, result, &
                           BCS_BOTH, rhs_b=fdmi%rhs_b(1:3, 0:3), rhs_t=fdmi%rhs_t(0:2, 1:4), bcs_b=wrk2d(:, 1), bcs_t=wrk2d(:, 2))
        case (5)
            call MatMul_5d(rhsi(:, 1:5), f, result, &
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
            result(:, 1) = wrk2d(:, 1) !*fdmi%lhs(1, idl)
            do ic = 1, idl - 1
                result(:, 1) = result(:, 1) + fdmi%lhs(1, idl + ic)*result(:, 1 + ic)
            end do
            result(:, 1) = result(:, 1) + fdmi%lhs(1, 1)*result(:, 1 + ic)
            ! result(:, 1) = (result(:, 1) + fdmi%lhs(1, 1)*result(:, 1 + ic))*fdmi%lhs(1, idl)

            if (present(du_boundary)) then      ! calculate u'n
                du_boundary(:) = fdmi%lhs(nx, idl)*result(:, nx)
                do ic = 1, idl - 1
                    du_boundary(:) = du_boundary(:) + fdmi%lhs(nx, idl - ic)*result(:, nx - ic)
                end do
                ic = idl                        ! longer stencil at the boundary
                du_boundary(:) = du_boundary(:) + fdmi%lhs(nx, ndl)*result(:, nx - ic)

                do ic = 1, idr - 1
                    du_boundary(:) = du_boundary(:) + rhsi(nx, idr - ic)*f(:, nx - ic)
                end do

            end if

        end if

        if (any([BCS_MIN] == fdmi%bc)) then
            result(:, nx) = wrk2d(:, 2) !*fdmi%lhs(nx, idl)
            do ic = 1, idl - 1
                result(:, nx) = result(:, nx) + fdmi%lhs(nx, idl - ic)*result(:, nx - ic)
            end do
            result(:, nx) = result(:, nx) + fdmi%lhs(nx, ndl)*result(:, nx - ic)
            ! result(:, nx) = (result(:, nx) + fdmi%lhs(nx, ndl)*result(:, nx - ic))*fdmi%lhs(nx, idl)

            if (present(du_boundary)) then      ! calculate u'1
                du_boundary(:) = fdmi%lhs(1, idl)*result(:, 1)
                do ic = 1, idl - 1
                    du_boundary(:) = du_boundary(:) + fdmi%lhs(1, idl + ic)*result(:, 1 + ic)
                end do
                ic = idl                        ! longer stencil at the boundary
                du_boundary(:) = du_boundary(:) + fdmi%lhs(1, 1)*result(:, 1 + ic)

                do ic = 1, idr - 1
                    du_boundary(:) = du_boundary(:) + rhsi(1, idr + ic)*f(:, 1 + ic)
                end do

            end if

        end if

        return
    end subroutine FDM_Int1_Solve

    !########################################################################

    ! subroutine FDM_Int1_Solve_APU(nlines, ilines, klines, fdmi, rhsi, f, result, wrk2d, du_boundary)
    !     integer(wi) nlines, ilines, klines
    !     type(fdm_integral_dt2), intent(in) :: fdmi
    !     real(wp), intent(in) :: rhsi(:, :)
    !     real(wp), intent(in) :: f(1:2*size(fdmi%lhs, 1), 1:klines, 1:ilines)
    !     real(wp), intent(inout) :: result(1:nlines, 1:size(fdmi%lhs, 1), 1:klines, 1:ilines)   ! contains bcs
    !     real(wp), intent(inout) :: wrk2d(nlines, 2, klines, ilines)
    !     real(wp), intent(out), optional :: du_boundary(nlines)

    !     ! -------------------------------------------------------------------
    !     integer(wi) :: nx, k, i
    !     integer(wi) :: idl, ndl, idr, ndr, ic

    !     ! ###################################################################
    !     nx = size(fdmi%lhs, 1)

    !     ndl = size(fdmi%lhs, 4)
    !     idl = ndl/2 + 1
    !     ndr = size(rhsi, 2)
    !     idr = ndr/2 + 1
        
    !     do i = 1, ilines
    !         do k = 1, klines
    !             select case (fdmi%bc(k, i))
    !             case (BCS_MIN)
    !                 result(1, nx, k, i) = f(1, k, i)
    !                 result(2, nx, k, i) = f(nx+1, k, i)
    !             case (BCS_MAX)
    !                 result(1, 1, k, i) = f(1, k, i)
    !                 result(2, 1, k, i) = f(nx+1, k, i)
    !             end select
    !         end do
    !     end do

    !     select case (ndr)
    !     case (3)
    !         call MatMul_3d_APU(nlines, klines, ilines, nx, fdmi, rhsi(:, 1:3), f(1:2*size(fdmi%lhs, 1), 1:klines, 1:ilines), &
    !         result(1:nlines, 1:nx, 1:klines, 1:ilines), BCS_BOTH, bcs_b=wrk2d(:, 1, :, :), bcs_t=wrk2d(:, 2, :, :))
    !     case (5)
    !         call MatMul_5d_APU(nlines, ilines, klines, nx, fdmi, rhsi(:, 1:5), f(1:2*size(fdmi%lhs, 1), 1:klines, 1:ilines), &
    !         result(1:nlines, 1:nx, 1:klines, 1:ilines), BCS_BOTH, bcs_b=wrk2d(:, 1, :, :), bcs_t=wrk2d(:, 2, :, :))
    !     end select

    !     select case (ndl)
    !     case (3)
    !         call TRIDSS_APU(nlines, nx, klines, ilines, fdmi, result(1:nlines, 1:nx, 1:klines, 1:ilines))
    !     case (5)
    !         call PENTADSS_APU(nlines, nx, klines, ilines, fdmi, result(1:nlines, 1:nx, 1:klines, 1:ilines))
    !     case (7)
    !         call HEPTADSS_APU(nlines, nx, klines, ilines, fdmi, result(1:nlines, 1:nx, 1:klines, 1:ilines))
    !     end select

    !     do i = 1, ilines
    !         do k = 1, klines

    !             if (any([BCS_MAX] == fdmi%bc(k,i))) then
    !                 result(:, 1, k, i) = wrk2d(:, 1, k, i) !*fdmi%lhs(1, idl)
    !                 do ic = 1, idl - 1
    !                     result(:, 1, k, i) = result(:, 1, k, i) + fdmi%lhs(1, k, i, idl + ic)*result(:, 1 + ic, k, i)
    !                 end do
    !                 result(:, 1, k, i) = result(:, 1, k, i) + fdmi%lhs(1, k, i, 1)*result(:, 1 + ic, k, i)
    !                 ! result(:, 1) = (result(:, 1) + fdmi%lhs(1, 1)*result(:, 1 + ic))*fdmi%lhs(1, idl)

    !                 if (present(du_boundary)) then      ! calculate u'n
    !                     du_boundary(:) = fdmi%lhs(nx, k, i, idl)*result(:, nx, k ,i)
    !                     do ic = 1, idl - 1
    !                         du_boundary(:) = du_boundary(:) + fdmi%lhs(nx, idl - ic, k, i)*result(:, nx - ic, k, i)
    !                     end do
    !                     ic = idl                        ! longer stencil at the boundary
    !                     du_boundary(:) = du_boundary(:) + fdmi%lhs(nx, k, i, ndl)*result(:, nx - ic, k, i)

    !                     do ic = 1, idr - 1
    !                         du_boundary(:) = du_boundary(:) + rhsi(nx, idr - ic)*f(:, nx - ic, k, i)
    !                     end do

    !                 end if

    !             end if

    !             if (any([BCS_MIN] == fdmi%bc)) thens
    !                 result(:, nx, k, i) = wrk2d(:, 2, k, i) !*fdmi%lhs(nx, idl)
    !                 do ic = 1, idl - 1
    !                     result(:, nx, k, i) = result(:, nx, k, i) + fdmi%lhs(nx, k, i, idl - ic)*result(:, nx - ic, k, i)
    !                 end do
    !                 result(:, nx, k, i) = result(:, nx, k, i) + fdmi%lhs(nx, k, i, ndl)*result(:, nx - ic, k, i)
    !                 ! result(:, nx) = (result(:, nx) + fdmi%lhs(nx, ndl)*result(:, nx - ic))*fdmi%lhs(nx, idl)

    !                 if (present(du_boundary)) then      ! calculate u'1
    !                     du_boundary(:) = fdmi%lhs(1, k, i, idl)*result(:, 1, k, i)
    !                     do ic = 1, idl - 1
    !                         du_boundary(:) = du_boundary(:) + fdmi%lhs(1, k, i, idl + ic)*result(:, 1 + ic, k, i)
    !                     end do
    !                     ic = idl                        ! longer stencil at the boundary
    !                     du_boundary(:) = du_boundary(:) + fdmi%lhs(1, k, i, 1)*result(:, 1 + ic, k, i)

    !                     do ic = 1, idr - 1
    !                         du_boundary(:) = du_boundary(:) + rhsi(1, idr + ic)*f(:, 1 + ic, k, i)
    !                     end do

    !                 end if

    !             end if
    !         end do
    !     end do

    !     return
    ! end subroutine FDM_Int1_Solve_APU

    !########################################################################
    !#
    !#     u''_i - \lamba^2 u_i = f_i  N-2 eqns
    !#     u_1 and u_N given           2   eqns
    !#     Au'' = Bu                   N   eqns
    !#
    !# starting from the matrices A (lhs, tridiagonal) and B (rhs, pentadiagonal, with unitary central diagonal).
    !#
    !# The system of N-2 eqns:
    !#
    !#                    (B - \lambda^2 A)u = Af = g
    !#
    !# is established in this routine, giving diagonals a-e and g (see notes).
    !#
    !# System normalized s.t. RHS diagonals are O(1), LHS diagonals O(h^2)
    !# System normalized s.t. 1. upper-diagonal in B is 1 (except at boundaries)
    !#
    !########################################################################
    subroutine FDM_Int2_Initialize(x, g, lambda2, ibc, fdmi)
        real(wp), intent(in) :: x(:)                    ! node positions
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda2                 ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) nx, nd

        !########################################################################
        call FDM_Int2_CreateSystem(x, g, lambda2, ibc, fdmi)

        ! LU decomposition
        nx = size(fdmi%lhs, 1)              ! # of grid points
        nd = size(fdmi%lhs, 2)              ! # of diagonals

        select case (nd)
        case (3)
            call TRIDFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3))
        case (5)
            ! We rely on this routines not changing a(2:3), b(2), e(ny-2:ny-1), d(ny-1)
            call PENTADFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5))
        case (7)
            call HEPTADFS(nx - 2, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), &
                          fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), fdmi%lhs(2:, 6), fdmi%lhs(2:, 7))
        end select

        return
    end subroutine FDM_Int2_Initialize

    !########################################################################
    !########################################################################
    ! Follows FDM_Int1_CreateSystem very much (see notes)
    subroutine FDM_Int2_CreateSystem(x, g, lambda2, ibc, fdmi)
        real(wp), intent(in) :: x(:)                    ! node positions
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda2                 ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt), intent(inout) :: fdmi    ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) idl, ndl, idr, ndr, ir, nx, i
        real(wp) dummy, rhsr_b(5, 0:7), rhsr_t(0:4, 8)
        real(wp) coef(5)

        ! ###################################################################
        ndl = g%nb_diag(1)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = g%nb_diag(2)
        idr = ndr/2 + 1             ! center diagonal in rhs
        nx = g%size                 ! # grid points

        ! check sizes
        if (abs(idl - idr) > 1) then
            call TLab_Write_ASCII(efile, __FILE__//'. lhs and rhs cannot differ by more than 2 diagonals.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        fdmi%mode_fdm = g%mode_fdm
        fdmi%lambda = lambda2
        fdmi%bc = ibc

        if (allocated(fdmi%lhs)) deallocate (fdmi%lhs)
        if (allocated(fdmi%rhs)) deallocate (fdmi%rhs)
        allocate (fdmi%lhs(nx, ndr))
        allocate (fdmi%rhs(nx, ndl))

        ! -------------------------------------------------------------------
        ! new rhs diagonals (array A22R), independent of lambda
        fdmi%rhs(:, :) = g%lhs(:, 1:ndl)

        call FDM_Bcs_Reduce(BCS_BOTH, fdmi%rhs, g%rhs(:, 1:ndr), rhsr_b, rhsr_t)

        fdmi%rhs_b = 0.0_wp
        fdmi%rhs_t = 0.0_wp

        fdmi%rhs_b(1:idl + 1, 1:ndl) = fdmi%rhs(1:idl + 1, 1:ndl)
        do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
            fdmi%rhs_b(1 + ir, idl - ir) = -rhsr_b(1 + ir, idr - ir)
        end do

        fdmi%rhs_t(0:idl, 1:ndl) = fdmi%rhs(nx - idl:nx, 1:ndl)
        do ir = 1, idr - 1              ! change sign in b^R_{2n} for nonzero bc
            fdmi%rhs_t(idl - ir, idl + ir) = -rhsr_t(idr - ir, idr + ir)
        end do

        ! -------------------------------------------------------------------
        ! new lhs diagonals (array C22R); remember rhs center diagonal is not saved because it was 1
        fdmi%lhs(:, :) = g%rhs(:, 1:ndr)

        fdmi%lhs(:, idr) = fdmi%lhs(:, idr) - lambda2*g%lhs(:, idl)               ! center diagonal
        do i = 1, idl - 1                                                       ! off-diagonals
            fdmi%lhs(1 + i:nx, idr - i) = fdmi%lhs(1 + i:nx, idr - i) - lambda2*g%lhs(1 + i:nx, idl - i)
            fdmi%lhs(1:nx - i, idr + i) = fdmi%lhs(1:nx - i, idr + i) - lambda2*g%lhs(1:nx - i, idl + i)
        end do

        ! fdmi%lhs(1:idr, 1:ndr) = rhsr_b(1:idr, 1:ndr)
        ! fdmi%lhs(1, idr + 1:idr + idl - 1) = fdmi%lhs(1, idr + 1:idr + idl - 1) + lambda2*fdmi%rhs_b(1, idl + 1:ndl)
        fdmi%lhs(2:idr, 1:ndr) = rhsr_b(2:idr, 1:ndr)
        do ir = 1, idr - 1
            fdmi%lhs(1 + ir, idr - idl + 1:idr + idl - 1) = fdmi%lhs(1 + ir, idr - idl + 1:idr + idl - 1) - lambda2*fdmi%rhs_b(1 + ir, 1:ndl)
        end do

        ! fdmi%lhs(nx - idr + 1:nx, 1:ndr) = rhsr_t(1:idr, 1:ndr)
        ! fdmi%lhs(nx, idr - idl + 1:idr - 1) = fdmi%lhs(nx, idr - idl + 1:idr - 1) + lambda2*fdmi%rhs_t(idl, 1:idl - 1)
        fdmi%lhs(nx - idr + 1:nx - 1, 1:ndr) = rhsr_t(1:idr - 1, 1:ndr)
        do ir = 1, idr - 1
            fdmi%lhs(nx - ir, idr - idl + 1:idr + idl - 1) = fdmi%lhs(nx - ir, idr - idl + 1:idr + idl - 1) - lambda2*fdmi%rhs_t(idl - ir, 1:ndl)
        end do

        ! -------------------------------------------------------------------
        ! Corrections to the BCS_DD to account for Neumann using third-order fdm for derivative at the boundary
        if (any([BCS_ND, BCS_NN] == fdmi%bc)) then
            ! Coefficients in FDM p'_1 = b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
            coef(:) = 0.0_wp
            ! coef(1:3) = coef_e1n2_biased(x, 1)                  ! second-order
            ! coef(1:4) = coef_e1n3_biased(x, 1)                  ! third-order
            coef(1:5) = coef_c1n4_biased(x, 1)                  ! fourth-order

            ! Solve for p_1 (see notes)
            fdmi%lhs(1, :) = 0.0_wp
            fdmi%lhs(1, 1:3) = -coef(2:4)/coef(1)               ! vector d_2
            fdmi%rhs_b(1, :) = 0.0_wp
            fdmi%rhs_b(1, idl) = 1.0_wp/coef(1)                 ! coefficient d_1
            fdmi%rhs_b(1, idl + 1) = -coef(5)/coef(1)           ! vector e_2, only 1 component

            ! Construct vector d + lambda^2h^2 e, e only 1 component
            fdmi%lhs(1, 1) = fdmi%lhs(1, 1) + lambda2*fdmi%rhs_b(1, idl + 1)

            ! Derived coefficients; contribution from -b^R_{21} (see notes)
            ! fdmi%lhs(2, 3:5) = fdmi%lhs(2, 3:5) - fdmi%rhs_b(1 + 1, idl - 1)*fdmi%lhs(1, 1:3)           ! in reduced C matrix
            ! fdmi%lhs(3, 2:4) = fdmi%lhs(3, 2:4) - fdmi%rhs_b(1 + 2, idl - 2)*fdmi%lhs(1, 1:3)

            ! fdmi%rhs_b(2, 2) = fdmi%rhs_b(2, 2) + fdmi%rhs_b(1 + 1, idl - 1)*fdmi%rhs_b(1, idl + 1)     ! in reduced A matrix
            ! fdmi%rhs_b(3, 1) = fdmi%rhs_b(3, 1) + fdmi%rhs_b(1 + 2, idl - 2)*fdmi%rhs_b(1, idl + 1)

            do ir = 1, idr - 1
                fdmi%lhs(1 + ir, idr - ir + 1:idr - ir + 1 + 2) = fdmi%lhs(1 + ir, idr - ir + 1:idr - ir + 1 + 2) &
                                                                  - fdmi%rhs_b(1 + ir, idl - ir)*fdmi%lhs(1, 1:3)       ! in reduced C matrix

                fdmi%rhs_b(1 + ir, idl - ir + 1) = fdmi%rhs_b(1 + ir, idl - ir + 1) &
                                                   + fdmi%rhs_b(1 + ir, idl - ir)*fdmi%rhs_b(1, idl + 1)                ! in reduced A matrix

                fdmi%rhs_b(1 + ir, idl - ir) = fdmi%rhs_b(1 + ir, idl - ir)*fdmi%rhs_b(1, idl)                          ! d_1 b^R_{21}
            end do

        end if

        if (any([BCS_DN, BCS_NN] == fdmi%bc)) then
            ! Coefficients in FDM p'_n = b_1 p_n + b_2 p_{n-1} + b_3 p_{n-2} +...
            coef(:) = 0.0_wp
            ! coef(1:3) = coef_e1n2_biased(x, nx, backwards=.true.)
            ! coef(1:4) = coef_e1n3_biased(x, nx, backwards=.true.)
            coef(1:5) = coef_c1n4_biased(x, nx, backwards=.true.)

            ! Solve for p_n (see notes)
            fdmi%lhs(nx, :) = 0.0_wp
            fdmi%lhs(nx, ndr - 2:ndr) = -coef([4, 3, 2])/coef(1)  ! vector d_n-1
            fdmi%rhs_t(idl, :) = 0.0_wp
            fdmi%rhs_t(idl, idl) = 1.0_wp/coef(1)                 ! coefficient d_n
            fdmi%rhs_t(idl, idl - 1) = -coef(5)/coef(1)           ! vector e_n-1, only 1 component

            ! Construct vector d + lambda^2h^2 e, e only 1 component
            fdmi%lhs(nx, ndr) = fdmi%lhs(nx, ndr) + lambda2*fdmi%rhs_t(idl, idl - 1)

            ! Derived coefficients; contribution from -b^R_{2n} (see notes)
            ! fdmi%lhs(nx - 1, 1:3) = fdmi%lhs(nx - 1, 1:3) - fdmi%rhs_t(idl - 1, idl + 1)*fdmi%lhs(nx, ndr - 2:ndr)              ! in reduced C matrix
            ! fdmi%lhs(nx - 2, 2:4) = fdmi%lhs(nx - 2, 2:4) - fdmi%rhs_t(idl - 2, idl + 2)*fdmi%lhs(nx, ndr - 2:ndr)

            ! fdmi%rhs_t(idl - 1, idl + 0) = fdmi%rhs_t(idl - 1, idl + 0) + fdmi%rhs_t(idl - 1, idl + 1)*fdmi%rhs_t(idl, idl - 1) ! in reduced A matrix
            ! fdmi%rhs_t(idl - 2, idl + 1) = fdmi%rhs_t(idl - 2, idl + 1) + fdmi%rhs_t(idl - 2, idl + 2)*fdmi%rhs_t(idl, idl - 1)

            do ir = 1, idr - 1
                fdmi%lhs(nx - ir, ir - 1 + 1:ir - 1 + 3) = fdmi%lhs(nx - ir, ir - 1 + 1:ir - 1 + 3) &
                                                           - fdmi%rhs_t(idl - ir, idl + ir)*fdmi%lhs(nx, ndr - 2:ndr)              ! in reduced C matrix

                fdmi%rhs_t(idl - ir, idl + ir - 1) = fdmi%rhs_t(idl - ir, idl + ir - 1) &
                                                     + fdmi%rhs_t(idl - ir, idl + ir)*fdmi%rhs_t(idl, idl - 1)      ! in reduced A matrix

                fdmi%rhs_t(idl - ir, idl + ir) = fdmi%rhs_t(idl - ir, idl + ir)*fdmi%rhs_t(idl, idl)                ! d_n b^R_{2n}
            end do

        end if

        ! -------------------------------------------------------------------
        ! normalization such that new central diagonal in rhs is 1
        do ir = 2, max(idr, idl + 1)
            dummy = 1.0_wp/fdmi%rhs(ir, idl)
            fdmi%rhs_b(ir, 0:ndl) = fdmi%rhs_b(ir, 0:ndl)*dummy

            dummy = 1.0_wp/fdmi%rhs(nx - ir + 1, idl)
            fdmi%rhs_t(idl - ir + 1, 1:ndl + 1) = fdmi%rhs_t(idl - ir + 1, 1:ndl + 1)*dummy

            dummy = 1.0_wp/fdmi%rhs(ir, idl)
            fdmi%rhs(ir, 1:ndl) = fdmi%rhs(ir, 1:ndl)*dummy
            fdmi%lhs(ir, 1:ndr) = fdmi%lhs(ir, 1:ndr)*dummy

            dummy = 1.0_wp/fdmi%rhs(nx - ir + 1, idl)
            fdmi%rhs(nx - ir + 1, 1:ndl) = fdmi%rhs(nx - ir + 1, 1:ndl)*dummy
            fdmi%lhs(nx - ir + 1, 1:ndr) = fdmi%lhs(nx - ir + 1, 1:ndr)*dummy

        end do

        ! interior points: normalization such that 1. upper-diagonal is 1
        do ir = max(idr, idl + 1) + 1, nx - max(idr, idl + 1)
            dummy = 1.0_wp/fdmi%rhs(ir, idl + 1)

            fdmi%rhs(ir, 1:ndl) = fdmi%rhs(ir, 1:ndl)*dummy
            fdmi%lhs(ir, 1:ndr) = fdmi%lhs(ir, 1:ndr)*dummy

        end do

        return
    contains
        !########################################################################
        ! 1. derivatie of interpolation polynomial between equations (15) and (16)
        !    p'_1= b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
        !
        ! Notation in Shukla and Zhong (2005), JCP, 204, 404–429 for the interpolation:
        !
        !       +                    I_n: set of points where the function and derivatives are given
        !   +---+---+---+---...
        !   +       +   +            I_m: set of points where only the function is given.
        !########################################################################
        function coef_c1n4_biased(x, i, backwards) result(coef)
            real(wp), intent(in) :: x(:)
            integer(wi), intent(in) :: i
            logical, intent(in), optional :: backwards
            real(wp) coef(5)

            real(wp) a2, b1, b2, b3, b4
            real(wp) dx1, dx3, dx4
            real(wp) D
            integer(wi) set_m(3), i1, i2, i3, i4

            i1 = i
            if (present(backwards)) then
                ! same as fowards, but changing the signs of the increments w.r.t. i
                ! To understand it, e.g., define a new variable k = -j, where k is the
                ! discrete variable moving around i
                i2 = i - 1
                i3 = i - 2
                i4 = i - 3
            else
                i2 = i + 1
                i3 = i + 2
                i4 = i + 3
            end if
            dx1 = x(i2) - x(i1)
            dx3 = x(i2) - x(i3)
            dx4 = x(i2) - x(i4)
            set_m = [i1, i3, i4]

            ! -------------------------------------------------------------------
            a2 = 0.5_wp*(Pi(x, i1, set_m) - dx1*Pi_p(x, i1, set_m))/Pi_p(x, i2, set_m)

            b2 = Pi_p(x, i1, set_m)*(2.0_wp*Pi_p(x, i2, set_m) + dx1*Pi_pp_3(x, i2, set_m)) &
                 - Pi(x, i1, set_m)*Pi_pp_3(x, i2, set_m)
            b2 = 0.5_wp*b2/Pi(x, i2, set_m)/Pi_p(x, i2, set_m)

            ! -------------------------------------------------------------------
            D = Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m)
            b1 = Lag(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + 2*dx1*Lag_p(x, i2, i1, set_m)) &
                 - dx1*Lag_p(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m))
            b1 = -b1/dx1/D

            D = Lag(x, i2, i3, set_m) + dx3*Lag_p(x, i2, i3, set_m)
            b3 = Lag(x, i1, i3, set_m)*(Lag(x, i2, i3, set_m) + 2*dx1*Lag_p(x, i2, i3, set_m)) &
                 - dx1*Lag_p(x, i1, i3, set_m)*(Lag(x, i2, i3, set_m) + dx1*Lag_p(x, i2, i3, set_m))
            b3 = -b3/dx3/D

            D = Lag(x, i2, i4, set_m) + dx4*Lag_p(x, i2, i4, set_m)
            b4 = Lag(x, i1, i4, set_m)*(Lag(x, i2, i4, set_m) + 2*dx1*Lag_p(x, i2, i4, set_m)) &
                 - dx1*Lag_p(x, i1, i4, set_m)*(Lag(x, i2, i4, set_m) + dx1*Lag_p(x, i2, i4, set_m))
            b4 = -b4/dx4/D

            coef = [b1, b2, b3, b4, a2]

            ! if uniform, we should have ( -29/6 54/6 -27/6 2/6 )/h and 3h
            ! print*, [b1, b2, b3, b4]*(x(2)-x(1))
            ! print*, a2/(x(2)-x(1))

            return
        end function

    end subroutine FDM_Int2_CreateSystem

    subroutine FDM_Int2_Initialize_APU(k, i, x, g, lambda2, ibc, fdmi_int2)
        integer, intent(in) :: k, i
        real(wp), intent(in) :: x(:)                        ! node positions
        type(fdm_derivative_dt), intent(in) :: g            ! derivative plan to be inverted
        real(wp), intent(in) :: lambda2                     ! system constant
        integer, intent(in) :: ibc                          ! type of boundary condition
        type(fdm_integral_dt2), intent(inout) :: fdmi_int2  ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) nx, nd

        !########################################################################
        call FDM_Int2_CreateSystem_APU(k, i, x, g, lambda2, ibc, fdmi_int2)

        ! LU decomposition
        nx = size(fdmi_int2%lhs, 1)              ! # of grid points
        nd = size(fdmi_int2%lhs, 4)              ! # of diagonals

        select case (nd)
        case (3)
            call TRIDFS(nx - 2, fdmi_int2%lhs(2:nx-1, k, i, 1), fdmi_int2%lhs(2:nx-1, k, i, 2), fdmi_int2%lhs(2:nx-1, k, i, 3))
        case (5)
            ! We rely on this routines not changing a(2:3), b(2), e(ny-2:ny-1), d(ny-1)
            call PENTADFS(nx - 2, fdmi_int2%lhs(2:nx-1, k, i, 1), fdmi_int2%lhs(2:nx-1, k, i, 2), fdmi_int2%lhs(2:nx-1, k, i, 3), &
                          fdmi_int2%lhs(2:nx-1, k, i, 4), fdmi_int2%lhs(2:nx-1, k, i, 5))
        case (7)
            call HEPTADFS(nx - 2, fdmi_int2%lhs(2:nx-1, k, i, 1), fdmi_int2%lhs(2:nx-1, k, i, 2), fdmi_int2%lhs(2:nx-1, k, i, 3), &
                          fdmi_int2%lhs(2:nx-1, k, i, 4), fdmi_int2%lhs(2:nx-1, k, i, 5), fdmi_int2%lhs(2:nx-1, k, i, 6), fdmi_int2%lhs(2:nx-1, k, i, 7))
        end select

        return
    end subroutine FDM_Int2_Initialize_APU
    !########################################################################
    !########################################################################

    subroutine FDM_Int2_CreateSystem_APU(ik, ii, x, g, lambda2, ibc, fdmi_int2)
        integer, intent(in) :: ik, ii
        real(wp), intent(in) :: x(:)                    ! node positions
        type(fdm_derivative_dt), intent(in) :: g        ! derivative plan to be inverted
        real(wp), intent(in) :: lambda2                 ! system constant
        integer, intent(in) :: ibc                      ! type of boundary condition
        type(fdm_integral_dt2), intent(inout) :: fdmi_int2   ! int_plan to be created; inout because otherwise allocatable arrays are deallocated

        ! -------------------------------------------------------------------
        integer(wi) idl, ndl, idr, ndr, ir, nx, i
        real(wp) dummy, rhsr_b(5, 0:7), rhsr_t(0:4, 8)
        real(wp) coef(5)

        ! ###################################################################
        ndl = g%nb_diag(1)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = g%nb_diag(2)
        idr = ndr/2 + 1             ! center diagonal in rhs
        nx = g%size                 ! # grid points

        ! check sizes
        if (abs(idl - idr) > 1) then
            call TLab_Write_ASCII(efile, __FILE__//'. lhs and rhs cannot differ by more than 2 diagonals.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        fdmi_int2%mode_fdm = g%mode_fdm
        fdmi_int2%lambda(ik,ii) = lambda2
        fdmi_int2%bc(ik,ii) = ibc

        ! -------------------------------------------------------------------
        ! new rhs diagonals (array A22R), independent of lambda
        fdmi_int2%rhs(:, ik, ii, 1:ndl) = g%lhs(:, 1:ndl)

        call FDM_Bcs_Reduce_APU(nx, ndl, ndr, BCS_BOTH, fdmi_int2%rhs(1:nx, ik, ii, 1:ndl), g%rhs(:, 1:ndr), rhsr_b, rhsr_t)

        fdmi_int2%rhs_b = 0.0_wp
        fdmi_int2%rhs_t = 0.0_wp

        fdmi_int2%rhs_b(1:idl + 1, ik, ii, 1:ndl) = fdmi_int2%rhs(1:idl + 1, ik , ii, 1:ndl)
        do ir = 1, idr - 1              ! change sign in b^R_{21} for nonzero bc
            fdmi_int2%rhs_b(1 + ir, ik, ii, idl - ir) = -rhsr_b(1 + ir, idr - ir)
        end do

        fdmi_int2%rhs_t(0:idl, ik, ii, 1:ndl) = fdmi_int2%rhs(nx - idl:nx, ik, ii, 1:ndl)
        do ir = 1, idr - 1              ! change sign in b^R_{2n} for nonzero bc
            fdmi_int2%rhs_t(idl - ir, ik, ii, idl + ir) = -rhsr_t(idr - ir, idr + ir)
        end do

        ! -------------------------------------------------------------------
        ! new lhs diagonals (array C22R); remember rhs center diagonal is not saved because it was 1
        fdmi_int2%lhs(:, ik, ii, :) = g%rhs(:, 1:ndr)

        fdmi_int2%lhs(:, ik, ii, idr) = fdmi_int2%lhs(:, ik, ii, idr) - lambda2*g%lhs(:, idl)               ! center diagonal
        do i = 1, idl - 1                                                       ! off-diagonals
            fdmi_int2%lhs(1 + i:nx, ik, ii, idr - i) = fdmi_int2%lhs(1 + i:nx, ik, ii, idr - i) - lambda2*g%lhs(1 + i:nx, idl - i)
            fdmi_int2%lhs(1:nx - i, ik, ii, idr + i) = fdmi_int2%lhs(1:nx - i, ik, ii, idr + i) - lambda2*g%lhs(1:nx - i, idl + i)
        end do

        ! fdmi%lhs(1:idr, 1:ndr) = rhsr_b(1:idr, 1:ndr)
        ! fdmi%lhs(1, idr + 1:idr + idl - 1) = fdmi%lhs(1, idr + 1:idr + idl - 1) + lambda2*fdmi%rhs_b(1, idl + 1:ndl)
        fdmi_int2%lhs(2:idr, ik, ii, 1:ndr) = rhsr_b(2:idr, 1:ndr)
        do ir = 1, idr - 1
            fdmi_int2%lhs(1 + ir, ik, ii, idr - idl + 1:idr + idl - 1) = fdmi_int2%lhs(1 + ir, ik, ii, idr - idl + 1:idr + idl - 1) - lambda2*fdmi_int2%rhs_b(1 + ir, ik, ii, 1:ndl)
        end do

        ! fdmi%lhs(nx - idr + 1:nx, 1:ndr) = rhsr_t(1:idr, 1:ndr)
        ! fdmi%lhs(nx, idr - idl + 1:idr - 1) = fdmi%lhs(nx, idr - idl + 1:idr - 1) + lambda2*fdmi%rhs_t(idl, 1:idl - 1)
        fdmi_int2%lhs(nx - idr + 1:nx - 1, ik, ii, 1:ndr) = rhsr_t(1:idr - 1, 1:ndr)
        do ir = 1, idr - 1
            fdmi_int2%lhs(nx - ir, ik, ii, idr - idl + 1:idr + idl - 1) = fdmi_int2%lhs(nx - ir, ik, ii, idr - idl + 1:idr + idl - 1) - lambda2*fdmi_int2%rhs_t(idl - ir, ik, ii, 1:ndl)
        end do

        ! -------------------------------------------------------------------
        ! Corrections to the BCS_DD to account for Neumann using third-order fdm for derivative at the boundary
        if (any([BCS_ND, BCS_NN] == fdmi_int2%bc(ik, ii))) then
            ! Coefficients in FDM p'_1 = b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
            coef(:) = 0.0_wp
            ! coef(1:3) = coef_e1n2_biased(x, 1)                  ! second-order
            ! coef(1:4) = coef_e1n3_biased(x, 1)                  ! third-order
            coef(1:5) = coef_c1n4_biased(x, 1)                  ! fourth-order

            ! Solve for p_1 (see notes)
            fdmi_int2%lhs(1, ik, ii, :) = 0.0_wp
            fdmi_int2%lhs(1, ik, ii, 1:3) = -coef(2:4)/coef(1)               ! vector d_2
            fdmi_int2%rhs_b(1, ik, ii, :) = 0.0_wp
            fdmi_int2%rhs_b(1, ik, ii, idl) = 1.0_wp/coef(1)                 ! coefficient d_1
            fdmi_int2%rhs_b(1, ik, ii, idl + 1) = -coef(5)/coef(1)           ! vector e_2, only 1 component

            ! Construct vector d + lambda^2h^2 e, e only 1 component
            fdmi_int2%lhs(1, ik, ii, 1) = fdmi_int2%lhs(1, ik, ii, 1) + lambda2*fdmi_int2%rhs_b(1, ik, ii, idl + 1)

            ! Derived coefficients; contribution from -b^R_{21} (see notes)
            ! fdmi_int2%lhs(2, 3:5) = fdmi_int2%lhs(2, 3:5) - fdmi_int2%rhs_b(1 + 1, idl - 1)*fdmi_int2%lhs(1, 1:3)           ! in reduced C matrix
            ! fdmi_int2%lhs(3, 2:4) = fdmi_int2%lhs(3, 2:4) - fdmi_int2%rhs_b(1 + 2, idl - 2)*fdmi_int2%lhs(1, 1:3)

            ! fdmi_int2%rhs_b(2, 2) = fdmi_int2%rhs_b(2, 2) + fdmi_int2%rhs_b(1 + 1, idl - 1)*fdmi_int2%rhs_b(1, idl + 1)     ! in reduced A matrix
            ! fdmi_int2%rhs_b(3, 1) = fdmi_int2%rhs_b(3, 1) + fdmi_int2%rhs_b(1 + 2, idl - 2)*fdmi_int2%rhs_b(1, idl + 1)

            do ir = 1, idr - 1
                fdmi_int2%lhs(1 + ir, ik, ii, idr - ir + 1:idr - ir + 1 + 2) = fdmi_int2%lhs(1 + ir, ik, ii, idr - ir + 1:idr - ir + 1 + 2) &
                                                   - fdmi_int2%rhs_b(1 + ir, ik, ii, idl - ir)*fdmi_int2%lhs(1, ik, ii, 1:3)       ! in reduced C matrix

                fdmi_int2%rhs_b(1 + ir, ik, ii, idl - ir + 1) = fdmi_int2%rhs_b(1 + ir, ik, ii, idl - ir + 1) &
                                                   + fdmi_int2%rhs_b(1 + ir, ik, ii, idl - ir)*fdmi_int2%rhs_b(1, ik, ii, idl + 1)                ! in reduced A matrix

                fdmi_int2%rhs_b(1 + ir, ik, ii, idl - ir) = fdmi_int2%rhs_b(1 + ir, ik, ii, idl - ir)*fdmi_int2%rhs_b(1, ik, ii, idl)                          ! d_1 b^R_{21}
            end do

        end if

        if (any([BCS_DN, BCS_NN] == fdmi_int2%bc(ik, ii))) then
            ! Coefficients in FDM p'_n = b_1 p_n + b_2 p_{n-1} + b_3 p_{n-2} +...
            coef(:) = 0.0_wp
            ! coef(1:3) = coef_e1n2_biased(x, nx, backwards=.true.)
            ! coef(1:4) = coef_e1n3_biased(x, nx, backwards=.true.)
            coef(1:5) = coef_c1n4_biased(x, nx, backwards=.true.)

            ! Solve for p_n (see notes)
            fdmi_int2%lhs(nx, ik ,ii, :) = 0.0_wp
            fdmi_int2%lhs(nx, ik, ii, ndr - 2:ndr) = -coef([4, 3, 2])/coef(1)  ! vector d_n-1
            fdmi_int2%rhs_t(idl, ik, ii, :) = 0.0_wp
            fdmi_int2%rhs_t(idl, ik, ii, idl) = 1.0_wp/coef(1)                 ! coefficient d_n
            fdmi_int2%rhs_t(idl, ik, ii, idl - 1) = -coef(5)/coef(1)           ! vector e_n-1, only 1 component

            ! Construct vector d + lambda^2h^2 e, e only 1 component
            fdmi_int2%lhs(nx, ik, ii, ndr) = fdmi_int2%lhs(nx, ik, ii, ndr) + lambda2*fdmi_int2%rhs_t(idl, ik, ii, idl - 1)

            ! Derived coefficients; contribution from -b^R_{2n} (see notes)
            ! fdmi_int2%lhs(nx - 1, 1:3) = fdmi_int2%lhs(nx - 1, 1:3) - fdmi_int2%rhs_t(idl - 1, idl + 1)*fdmi_int2%lhs(nx, ndr - 2:ndr)              ! in reduced C matrix
            ! fdmi_int2%lhs(nx - 2, 2:4) = fdmi_int2%lhs(nx - 2, 2:4) - fdmi_int2%rhs_t(idl - 2, idl + 2)*fdmi_int2%lhs(nx, ndr - 2:ndr)

            ! fdmi_int2%rhs_t(idl - 1, idl + 0) = fdmi_int2%rhs_t(idl - 1, idl + 0) + fdmi_int2%rhs_t(idl - 1, idl + 1)*fdmi_int2%rhs_t(idl, idl - 1) ! in reduced A matrix
            ! fdmi_int2%rhs_t(idl - 2, idl + 1) = fdmi_int2%rhs_t(idl - 2, idl + 1) + fdmi_int2%rhs_t(idl - 2, idl + 2)*fdmi_int2%rhs_t(idl, idl - 1)

            do ir = 1, idr - 1
                fdmi_int2%lhs(nx - ir, ik, ii, ir - 1 + 1:ir - 1 + 3) = fdmi_int2%lhs(nx - ir, ik, ii, ir - 1 + 1:ir - 1 + 3) &
                                                     - fdmi_int2%rhs_t(idl - ir, ik, ii, idl + ir)*fdmi_int2%lhs(nx, ik, ii, ndr - 2:ndr)              ! in reduced C matrix

                fdmi_int2%rhs_t(idl - ir, ik, ii, idl + ir - 1) = fdmi_int2%rhs_t(idl - ir, ik, ii, idl + ir - 1) &
                                                     + fdmi_int2%rhs_t(idl - ir, ik, ii, idl + ir)*fdmi_int2%rhs_t(idl, ik, ii, idl - 1)      ! in reduced A matrix

                fdmi_int2%rhs_t(idl - ir, ik, ii, idl + ir) = fdmi_int2%rhs_t(idl - ir, ik, ii, idl + ir)*fdmi_int2%rhs_t(idl, ik, ii, idl)                ! d_n b^R_{2n}
            end do

        end if

        ! -------------------------------------------------------------------
        ! normalization such that new central diagonal in rhs is 1
        do ir = 2, max(idr, idl + 1)
            dummy = 1.0_wp/fdmi_int2%rhs(ir, ik, ii, idl)
            fdmi_int2%rhs_b(ir, ik, ii, 0:ndl) = fdmi_int2%rhs_b(ir, ik, ii, 0:ndl)*dummy

            dummy = 1.0_wp/fdmi_int2%rhs(nx - ir + 1, ik, ii, idl)
            fdmi_int2%rhs_t(idl - ir + 1, ik, ii, 1:ndl + 1) = fdmi_int2%rhs_t(idl - ir + 1, ik, ii, 1:ndl + 1)*dummy

            dummy = 1.0_wp/fdmi_int2%rhs(ir, ik, ii, idl)
            fdmi_int2%rhs(ir, ik, ii, 1:ndl) = fdmi_int2%rhs(ir, ik, ii, 1:ndl)*dummy
            fdmi_int2%lhs(ir, ik, ii, 1:ndr) = fdmi_int2%lhs(ir, ik, ii, 1:ndr)*dummy

            dummy = 1.0_wp/fdmi_int2%rhs(nx - ir + 1, ik, ii, idl)
            fdmi_int2%rhs(nx - ir + 1, ik, ii, 1:ndl) = fdmi_int2%rhs(nx - ir + 1, ik, ii, 1:ndl)*dummy
            fdmi_int2%lhs(nx - ir + 1, ik, ii, 1:ndr) = fdmi_int2%lhs(nx - ir + 1, ik, ii, 1:ndr)*dummy

        end do

        ! interior points: normalization such that 1. upper-diagonal is 1
        do ir = max(idr, idl + 1) + 1, nx - max(idr, idl + 1)
            dummy = 1.0_wp/fdmi_int2%rhs(ir, ik, ii, idl + 1)

            fdmi_int2%rhs(ir, ik, ii, 1:ndl) = fdmi_int2%rhs(ir, ik, ii, 1:ndl)*dummy
            fdmi_int2%lhs(ir, ik, ii, 1:ndr) = fdmi_int2%lhs(ir, ik, ii, 1:ndr)*dummy

        end do

        return
    contains
        !########################################################################
        ! 1. derivatie of interpolation polynomial between equations (15) and (16)
        !    p'_1= b_1 p_1 + b_2 p_2 + b_3 p_3 + b_4 p_4 + a_2 p''_2
        !
        ! Notation in Shukla and Zhong (2005), JCP, 204, 404–429 for the interpolation:
        !
        !       +                    I_n: set of points where the function and derivatives are given
        !   +---+---+---+---...
        !   +       +   +            I_m: set of points where only the function is given.
        !########################################################################
        function coef_c1n4_biased(x, i, backwards) result(coef)
            real(wp), intent(in) :: x(:)
            integer(wi), intent(in) :: i
            logical, intent(in), optional :: backwards
            real(wp) coef(5)

            real(wp) a2, b1, b2, b3, b4
            real(wp) dx1, dx3, dx4
            real(wp) D
            integer(wi) set_m(3), i1, i2, i3, i4

            i1 = i
            if (present(backwards)) then
                ! same as fowards, but changing the signs of the increments w.r.t. i
                ! To understand it, e.g., define a new variable k = -j, where k is the
                ! discrete variable moving around i
                i2 = i - 1
                i3 = i - 2
                i4 = i - 3
            else
                i2 = i + 1
                i3 = i + 2
                i4 = i + 3
            end if
            dx1 = x(i2) - x(i1)
            dx3 = x(i2) - x(i3)
            dx4 = x(i2) - x(i4)
            set_m = [i1, i3, i4]

            ! -------------------------------------------------------------------
            a2 = 0.5_wp*(Pi(x, i1, set_m) - dx1*Pi_p(x, i1, set_m))/Pi_p(x, i2, set_m)

            b2 = Pi_p(x, i1, set_m)*(2.0_wp*Pi_p(x, i2, set_m) + dx1*Pi_pp_3(x, i2, set_m)) &
                 - Pi(x, i1, set_m)*Pi_pp_3(x, i2, set_m)
            b2 = 0.5_wp*b2/Pi(x, i2, set_m)/Pi_p(x, i2, set_m)

            ! -------------------------------------------------------------------
            D = Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m)
            b1 = Lag(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + 2*dx1*Lag_p(x, i2, i1, set_m)) &
                 - dx1*Lag_p(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m))
            b1 = -b1/dx1/D

            D = Lag(x, i2, i3, set_m) + dx3*Lag_p(x, i2, i3, set_m)
            b3 = Lag(x, i1, i3, set_m)*(Lag(x, i2, i3, set_m) + 2*dx1*Lag_p(x, i2, i3, set_m)) &
                 - dx1*Lag_p(x, i1, i3, set_m)*(Lag(x, i2, i3, set_m) + dx1*Lag_p(x, i2, i3, set_m))
            b3 = -b3/dx3/D

            D = Lag(x, i2, i4, set_m) + dx4*Lag_p(x, i2, i4, set_m)
            b4 = Lag(x, i1, i4, set_m)*(Lag(x, i2, i4, set_m) + 2*dx1*Lag_p(x, i2, i4, set_m)) &
                 - dx1*Lag_p(x, i1, i4, set_m)*(Lag(x, i2, i4, set_m) + dx1*Lag_p(x, i2, i4, set_m))
            b4 = -b4/dx4/D

            coef = [b1, b2, b3, b4, a2]

            ! if uniform, we should have ( -29/6 54/6 -27/6 2/6 )/h and 3h
            ! print*, [b1, b2, b3, b4]*(x(2)-x(1))
            ! print*, a2/(x(2)-x(1))

            return
        end function

    end subroutine FDM_Int2_CreateSystem_APU

    !########################################################################
    !########################################################################
    ! Allow to pass separate rhs because this part does not depend on lambda
    subroutine FDM_Int2_Solve(nlines, fdmi, rhsi, f, result, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi
        real(wp), intent(in) :: rhsi(:, :)
        real(wp), intent(in) :: f(nlines, size(fdmi%lhs, 1))
        real(wp), intent(inout) :: result(nlines, size(fdmi%lhs, 1))   ! contains bcs
        real(wp), intent(inout) :: wrk2d(nlines, 2)

        ! -------------------------------------------------------------------
        integer(wi) :: nx, ndl, ndr

        ! ###################################################################
        nx = size(fdmi%lhs, 1)
        ndl = size(fdmi%lhs, 2)
        ndr = size(rhsi, 2)

        select case (ndr)
        case (3)
            call MatMul_3d(rhsi(:, 1:3), f, result, &
                           BCS_BOTH, rhs_b=fdmi%rhs_b(1:3, 0:3), rhs_t=fdmi%rhs_t(0:2, 1:4), bcs_b=wrk2d(:, 1), bcs_t=wrk2d(:, 2))
        case (5)
            call MatMul_5d(rhsi(:, 1:5), f, result, &
                           BCS_BOTH, rhs_b=fdmi%rhs_b(1:4, 0:5), rhs_t=fdmi%rhs_t(0:3, 1:6), bcs_b=wrk2d(:, 1), bcs_t=wrk2d(:, 2))

        end select

        ! Solve pentadiagonal linear system
        select case (ndl)
        case (3)
            call TRIDSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), result(:, 2:))
        case (5)
            call PENTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), result(:, 2:))
        case (7)
            call HEPTADSS(nx - 2, nlines, fdmi%lhs(2:, 1), fdmi%lhs(2:, 2), fdmi%lhs(2:, 3), fdmi%lhs(2:, 4), fdmi%lhs(2:, 5), fdmi%lhs(2:, 6), fdmi%lhs(2:, 7), result(:, 2:))
        end select

        !   Corrections to the BCS_DD to account for Neumann
        if (any([BCS_ND, BCS_NN] == fdmi%bc)) then
            result(:, 1) = wrk2d(:, 1) &
                           + fdmi%lhs(1, 1)*result(:, 2) + fdmi%lhs(1, 2)*result(:, 3) + fdmi%lhs(1, 3)*result(:, 4)
        end if

        if (any([BCS_DN, BCS_NN] == fdmi%bc)) then
            result(:, nx) = wrk2d(:, 2) &
                            + fdmi%lhs(nx, ndl)*result(:, nx - 1) + fdmi%lhs(nx, ndl - 1)*result(:, nx - 2) + fdmi%lhs(nx, ndl - 2)*result(:, nx - 3)
        end if

        return
    end subroutine FDM_Int2_Solve

    subroutine FDM_Int2_Solve_APU(nlines, ilines, klines, fdmi_int2, rhsi, f, result, wrk2d)
        use TLab_Time, only: fdm_solve2_time
        integer(wi) nlines, ilines, klines
        type(fdm_integral_dt2), intent(in) :: fdmi_int2
        real(wp), intent(in) :: rhsi(:, :)
        real(wp), intent(in) :: f(1:2*size(fdmi_int2%lhs, 1), 1:klines, 1:ilines)
        real(wp), intent(inout) :: result(1:nlines, 1:size(fdmi_int2%lhs, 1), 1:klines, 1:ilines)   ! contains bcs
        real(wp), intent(inout) :: wrk2d(nlines, 2, klines, ilines)

        ! -------------------------------------------------------------------
        integer(wi) :: nx, ndl, ndr, i, j, k, bcs
        integer(wi) :: clock_0, clock_1, clock_cycle

        ! ###################################################################
        nx = size(fdmi_int2%lhs, 1)
        ndl = size(fdmi_int2%lhs, 4)
        ndr = size(rhsi, 2)

        select case (ndr)
        case (3)
            call MatMul_3d_APU(nlines, klines, ilines, nx, fdmi_int2, rhsi(:, 1:3), f(1:2*size(fdmi_int2%lhs, 1), 1:klines, 1:ilines), &
            result(1:nlines, 1:nx, 1:klines, 1:ilines), BCS_BOTH, bcs_b=wrk2d(:, 1, :, :), bcs_t=wrk2d(:, 2, :, :))
        case (5)
            call MatMul_5d_APU(nlines, ilines, klines, nx, fdmi_int2, rhsi(:, 1:5), f(1:2*size(fdmi_int2%lhs, 1), 1:klines, 1:ilines), &
            result(1:nlines, 1:nx, 1:klines, 1:ilines), BCS_BOTH, bcs_b=wrk2d(:, 1, :, :), bcs_t=wrk2d(:, 2, :, :))
        end select

        ! Solve pentadiagonal linear system
        select case (ndl)
        case (3)
            call TRIDSS_APU(nlines, nx, klines, ilines, fdmi_int2, result(1:nlines, 1:nx, 1:klines, 1:ilines))
        case (5)
            call PENTADSS_APU(nlines, nx, klines, ilines, fdmi_int2, result(1:nlines, 1:nx, 1:klines, 1:ilines))  !%lhs(2:, 1), fdmi_int2%lhs(2:, 2), fdmi_int2%lhs(2:, 3), fdmi_int2%lhs(2:, 4), fdmi_int2%lhs(2:, 5), result(:, 2:))
        case (7)
            call HEPTADSS_APU(nlines, nx, klines, ilines, fdmi_int2, result(1:nlines, 1:nx, 1:klines, 1:ilines))
        end select

        CALL SYSTEM_CLOCK(clock_0,clock_cycle) 
        !$omp target teams distribute parallel do collapse(2) &
        !$omp private(i,j,k,bcs) &
        !$omp shared(ilines,klines,nlines,fdmi_int2%bc,fdmi_int2%lhs,BCS_ND,BCS_NN,BCS_DN,result,wrk2d,nx,ndl)
        do i = 1, ilines
            do k = 1, klines
                bcs = fdmi_int2%bc(k,i)
                !   Corrections to the BCS_DD to account for Neumann
                if ((BCS_ND == bcs) .or. (BCS_NN == bcs) ) then
                    do j = 1, nlines
                        result(j, 1, k, i) = wrk2d(j, 1, k, i) &
                                + fdmi_int2%lhs(1, k, i, 1)*result(j, 2, k, i) + fdmi_int2%lhs(1, k, i, 2)*result(j, 3, k, i) + fdmi_int2%lhs(1, k, i, 3)*result(j, 4, k, i)
                    end do
                end if

                if ((BCS_DN == bcs) .or. (BCS_NN == bcs)) then
                    do j = 1, nlines
                        result(j, nx, k, i) = wrk2d(j, 2, k, i) &
                                    + fdmi_int2%lhs(nx, k, i, ndl)*result(j, nx - 1, k, i) + fdmi_int2%lhs(nx, k, i, ndl - 1)*result(j, nx - 2, k, i) &
                                    + fdmi_int2%lhs(nx, k, i, ndl - 2)*result(j, nx - 3, k, i)
                    end do
                end if
            end do
        end do
        !$omp end target teams distribute parallel do
        CALL SYSTEM_CLOCK(clock_1,clock_cycle)
        fdm_solve2_time = fdm_solve2_time + real(clock_1 - clock_0)/real(clock_cycle)
        return
    end subroutine FDM_Int2_Solve_APU
    
end module FDM_Integral
