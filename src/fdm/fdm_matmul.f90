!########################################################################
! Matrix multiplication of n-diagonal matrix with a vector with special boundary conditions
! The boundary conditions can extend over n/2 + 2 points
! The 1. upper-diagonal in interior points is equal to 1
! This allows use to handle systems A y = B x in which A amd B differ by up to 2 diagonals (see notes)
!########################################################################
module FDM_MatMul
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use Tlab_Type, only : fdm_integral_dt
    implicit none
    private

    ! generic cases
    public MatMul_3d            ! Calculate f = B u, assuming B is tridiagonal
    public MatMul_3d_test
    public MatMul_3d_add        ! Calculate f = f + B u, assuming B is tridiagonal
    ! special cases where coefficients are constant in the interior points
    public MatMul_3d_antisym    ! Calculate f = B u, assuming B is tridiagonal, antisymmetric
    public MatMul_3d_sym        ! Calculate f = B u, assuming B is tridiagonal, symmetric

    ! generic cases
    public MatMul_5d            ! Calculate f = B u, assuming B is pentadiagonal
    public MatMul_5d_add        ! Calculate f = f + B u, assuming B is pentadiagonal
    ! special cases where coefficients are constant in the interior points
    public MatMul_5d_antisym    ! Calculate f = B u, assuming B is pentadiagonal, antisymmetric
    public MatMul_5d_sym        ! Calculate f = B u, assuming B is pentadiagonal, symmetric

    ! generic cases
    ! tbd when needed
    ! special cases where coefficients are constant in the interior points
    public MatMul_7d_antisym    ! Calculate f = B u, assuming B is heptadiagonal, antisymmetric
    public MatMul_7d_sym        ! Calculate f = B u, assuming B is heptadiagonal, symmetric

#define r0_b(j) rhs_b(j,0)
#define r1_b(j) rhs_b(j,1)
#define r2_b(j) rhs_b(j,2)
#define r3_b(j) rhs_b(j,3)
#define r4_b(j) rhs_b(j,4)
#define r5_b(j) rhs_b(j,5)
#define r6_b(j) rhs_b(j,6)
#define r7_b(j) rhs_b(j,7)

#define r1_t(j) rhs_t(j,1)
#define r2_t(j) rhs_t(j,2)
#define r3_t(j) rhs_t(j,3)
#define r4_t(j) rhs_t(j,4)
#define r5_t(j) rhs_t(j,5)
#define r6_t(j) rhs_t(j,6)
#define r7_t(j) rhs_t(j,7)

! to be implemented and checked
#define r1_i(j) rhs(j,1)
#define r2_i(j) rhs(j,2)
#define r3_i(j) rhs(j,3)
#define r4_i(j) rhs(j,4)
#define r5_i(j) rhs(j,5)
#define r6_i(j) rhs(j,6)
#define r7_i(j) rhs(j,7)

#define r0b(i,j) fdmi(i)%rhs_b(j,0)
#define r1b(i,j) fdmi(i)%rhs_b(j,1)
#define r2b(i,j) fdmi(i)%rhs_b(j,2)
#define r3b(i,j) fdmi(i)%rhs_b(j,3)
#define r4b(i,j) fdmi(i)%rhs_b(j,4)
#define r5b(i,j) fdmi(i)%rhs_b(j,5)
#define r6b(i,j) fdmi(i)%rhs_b(j,6)
#define r7b(i,j) fdmi(i)%rhs_b(j,7)

#define r1t(i,j) fdmi(i)%rhs_t(j,1)
#define r2t(i,j) fdmi(i)%rhs_t(j,2)
#define r3t(i,j) fdmi(i)%rhs_t(j,3)
#define r4t(i,j) fdmi(i)%rhs_t(j,4)
#define r5t(i,j) fdmi(i)%rhs_t(j,5)
#define r6t(i,j) fdmi(i)%rhs_t(j,6)
#define r7t(i,j) fdmi(i)%rhs_t(j,7)

contains
    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is tri-diagonal and 1. upper-diagonal in interior points is equal to 1
    ! Special boundary conditions restricted to 3 points:
    ! r_11 r_12 r_13
    !      r_21 r_22 r_23
    !      r_30 r_31 r_32 r_33
    !                r_41  1.  r_43         <- interior points start here
    !                     ...  ...  ...
    subroutine MatMul_3d(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        use TLab_Time, only: mat3d_time, t_compute
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u
        integer, intent(in), optional :: ibc
        real(wp), intent(in), optional :: rhs_b(1:3, 0:3), rhs_t(0:2, 1:4)  ! Special bcs at bottom and top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx, len
        integer ibc_loc
        ! APU offloading 
#ifdef USE_APU
        integer(wi) l
#endif
        integer clock_0, clock_1, clock_2, clock_cycle, clock_3

        ! -----------------------------------------------------------------------
        ! Profiling
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_0,clock_cycle) 

        ! #######################################################################
        nx = size(rhs, 1)
        len = size(f,1)
        ! print *, nd = size(rhs, 2)    ! # diagonals, should be 3
        ! size(u,2) and size(f,2) should be nx
        ! size(u,1) and size(f,1) and size(bcs_b) and size(bcs_t) should be the same (number of equations to solve)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        ! -------------------------------------------------------------------
        ! Boundary; the first 3/2+1+1=3 rows might be different
        if (any([BCS_MIN, BCS_BOTH] == ibc_loc)) then
            if (present(bcs_b)) bcs_b(:) = f(:, 1)*r2_b(1) + u(:, 2)*r3_b(1) + u(:, 3)*r1_b(1) ! r1(1) contains extended stencil
            ! f(1) contains the boundary condition
            f(:, 2) = f(:, 1)*r1_b(2) + u(:, 2)*r2_b(2) + u(:, 3)*r3_b(2)
            f(:, 3) = f(:, 1)*r0_b(3) + u(:, 2)*r1_b(3) + u(:, 3)*r2_b(3) + u(:, 4)*r3_b(3)
        else
            f(:, 1) = u(:, 1)*r2_i(1) + u(:, 2)*r3_i(1) + u(:, 3)*r1_i(1)   ! r1(1) contains extended stencil
            f(:, 2) = u(:, 1)*r1_i(2) + u(:, 2)*r2_i(2) + u(:, 3)*r3_i(2)
            f(:, 3) = u(:, 2)*r1_i(3) + u(:, 3)*r2_i(3) + u(:, 4)*r3_i(3)
        end if
#ifndef USE_APU
        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 4, nx - 3
            f(:, n) = u(:, n - 1)*r1_i(n) + u(:, n)*r2_i(n) + u(:, n + 1)
        end do
#else
        ! -----------------------------------------------------------------------
        ! With APU ACCELERATION 
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_2) 
        !$omp target teams distribute parallel do collapse(2) default(none) &
        !$omp private(l,n) & 
        !$omp shared(len,u,f,rhs,nx)
            do n = 4, nx - 3
                do l = 1, len
                    f(l, n) = u(l, n - 1)*rhs(n,1) + u(l, n)*rhs(n,2) + u(l, n + 1)
                end do
            end do
        !$omp end target teams distribute parallel do
        call SYSTEM_CLOCK(clock_3) 
#endif
        ! -------------------------------------------------------------------
        ! Boundary; the last 3/2+1+1=3 rows might be different
        if (any([BCS_MAX, BCS_BOTH] == ibc_loc)) then
            ! f(nx) contains the boundary condition
            f(:, nx - 2) = u(:, nx - 3)*r1_t(0) + u(:, nx - 2)*r2_t(0) + u(:, nx - 1)*r3_t(0) + f(:, nx)*r4_t(0)
            f(:, nx - 1) = u(:, nx - 2)*r1_t(1) + u(:, nx - 1)*r2_t(1) + f(:, nx)*r3_t(1)
            if (present(bcs_t)) bcs_t(:) = u(:, nx - 2)*r3_t(2) + u(:, nx - 1)*r1_t(2) + f(:, nx)*r2_t(2) ! r3(nx) contains extended stencil
        else
            f(:, nx - 2) = u(:, nx - 3)*r1_i(nx - 2) + u(:, nx - 2)*r2_i(nx - 2) + u(:, nx - 1)*r3_i(nx - 2)
            f(:, nx - 1) = u(:, nx - 2)*r1_i(nx - 1) + u(:, nx - 1)*r2_i(nx - 1) + u(:, nx)*r3_i(nx - 1)
            f(:, nx) = u(:, nx - 2)*r3_i(nx) + u(:, nx - 1)*r1_i(nx) + u(:, nx)*r2_i(nx) ! r3(nx) contains extended stencil
        end if

        ! -----------------------------------------------------------------------
        ! Profiling
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_1)
        mat3d_time = mat3d_time + real(clock_1 - clock_0)/ real(clock_cycle) 
        t_compute  = t_compute + real(clock_3 - clock_2) / real(clock_cycle)
        return
    end subroutine MatMul_3d

    subroutine MatMul_3d_test(nlines, ilines, rhs, u, f, ibc, fdmi, bcs_b, bcs_t)
        use TLab_Time, only: mat3d_time, t_compute
        integer(wi) nlines, ilines
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :, :)                                 ! vector f = B u
        integer, intent(in), optional :: ibc
        type(fdm_integral_dt), intent(in) :: fdmi(:) !rhs_b(1:3, 0:3), rhs_t(0:2, 1:4)  ! Special bcs at bottom and top
        real(wp), intent(out), optional :: bcs_b(:,:), bcs_t(:,:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx, len, i
        integer(wi) pa, pb, pc, pd, pe, pf
        integer(wi) lp0, lp1, lp2, lp3, lp4, lp5, lp6, lp7
        integer ibc_loc
        ! APU offloading 

        integer clock_0, clock_1, clock_2, clock_cycle, clock_3

        ! -----------------------------------------------------------------------
        ! Profiling
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_0,clock_cycle) 

        ! #######################################################################
        nx = size(rhs, 1)
        len = size(f,1)

        lp0 = 2*nx; lp1 = 2*nx - 1; lp2 = 2*nx - 2; lp3 = 2*nx - 3 
        lp5 = 2*nx - 5; lp4 = 2*nx - 4; lp7 = 2*nx - 7; lp6 = 2*nx - 6 

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        do i = 1, ilines
            ! -------------------------------------------------------------------
            ! Boundary; the first 3/2+1+1=3 rows might be different
            if (any([BCS_MIN, BCS_BOTH] == ibc_loc)) then
                if (present(bcs_b)) bcs_b(:,i) = f(:, 1, i)*r2b(i, 1) + u(3:4, i)*r3b(i, 1) + u(5:6, i)*r1b(i, 1) ! r1(1) contains extended stencil
                ! f(1) contains the boundary condition
                f(:, 2, i) = f(:, 1, i)*r1b(i, 2) + u(3:4, i)*r2b(i, 2) + u(5:6, i)*r3b(i, 2)
                f(:, 3, i) = f(:, 1, i)*r0b(i, 3) + u(3:4, i)*r1b(i, 3) + u(5:6, i)*r2b(i, 3) + u(7:8, i)*r3b(i, 3)
            else
                f(:, 1, i) = u(1:2, i)*r2_i(1) + u(3:4, i)*r3_i(1) + u(5:6, i)*r1_i(1)   ! r1(1) contains extended stencil
                f(:, 2, i) = u(1:2, i)*r1_i(2) + u(3:4, i)*r2_i(2) + u(5:6, i)*r3_i(2)
                f(:, 3, i) = u(3:4, i)*r1_i(3) + u(5:6, i)*r2_i(3) + u(7:8, i)*r3_i(3)
            end if
            ! -------------------------------------------------------------------
            ! Interior points; accelerate
                do n = 4, nx - 3
                    pa = 2*n - 3; pb = 2*n - 2; pc = 2*n - 1; pd = 2*n; pe = 2*n + 1; pf = 2*n + 2
                    f(:, n, i) = u(pa:pb, i)*r1_i(n) + u(pc:pd, i)*r2_i(n) + u(pe:pf, i)
                end do
            ! -------------------------------------------------------------------
            ! Boundary; the last 3/2+1+1=3 rows might be different
            if (any([BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(nx) contains the boundary condition
                f(:, nx - 2, i) = u(lp7:lp6, i)*r1t(i, 0) + u(lp5:lp4, i)*r2t(i, 0) + u(lp3:lp2, i)*r3t(i, 0) + f(:, nx, i)*r4t(i, 0)
                f(:, nx - 1, i) = u(lp5:lp4, i)*r1t(i, 1) + u(lp3:lp2, i)*r2t(i, 1) + f(:, nx, i)*r3t(i, 1)
                if (present(bcs_t)) bcs_t(:,i) = u(lp5:lp4, i)*r3t(i,2) + u(lp3:lp2, i)*r1t(i,2) + f(:, nx, i)*r2t(i, 2) ! r3(nx) contains extended stencil
            else
                f(:, nx - 2, i) = u(lp7:lp6, i)*r1_i(nx - 2) + u(lp5:lp4, i)*r2_i(nx - 2) + u(lp3:lp2, i)*r3_i(nx - 2)
                f(:, nx - 1, i) = u(lp5:lp4, i)*r1_i(nx - 1) + u(lp3:lp2, i)*r2_i(nx - 1) + u(lp1:lp0, i)*r3_i(nx - 1)
                f(:, nx, i) = u(lp5:lp4, i)*r3_i(nx) + u(lp3:lp2, i)*r1_i(nx) + u(lp1:lp0, i)*r2_i(nx) ! r3(nx) contains extended stencil
            end if
        end do
        ! -----------------------------------------------------------------------
        ! Profiling
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_1)
        mat3d_time = mat3d_time + real(clock_1 - clock_0)/ real(clock_cycle) 
        t_compute  = t_compute + real(clock_3 - clock_2) / real(clock_cycle)
        return
    end subroutine MatMul_3d_test

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = f + B u, assuming B is tri-diagonal
    subroutine MatMul_3d_add(nx, nlines, r1, r2, r3, u, f)
        use Tlab_Time, only: mat3dadd_time
        integer(wi), intent(in) :: nx, nlines                       ! nlines linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx)              ! RHS diagonals
        real(wp), intent(in) :: u(nlines, nx)                       ! function u
        real(wp), intent(inout) :: f(nlines, nx)                    ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! APU offloading 
#ifdef USE_APU
        integer(wi) l
#endif
        integer clock_0, clock_1, clock_cycle

        ! -----------------------------------------------------------------------
        ! Profiling
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_0,clock_cycle) 
        
#ifndef USE_APU

        ! -------------------------------------------------------------------
        ! Boundary
        n = 1
        f(:, n) = f(:, n) + u(:, n)*r2(n) + u(:, n + 1)*r3(n) + u(:, n + 2)*r1(n)   ! r1(1) contains extended stencil

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = f(:, n) + u(:, n - 1)*r1(n) + u(:, n)*r2(n) + u(:, n + 1)*r3(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        n = nx
        f(:, n) = f(:, n) + u(:, n - 2)*r3(n) + u(:, n - 1)*r1(n) + u(:, n)*r2(n)   ! r3(nx) contains extended stencil

        ! -----------------------------------------------------------------------
        ! With APU ACCELERATION 
        ! -----------------------------------------------------------------------
#else
        ! -------------------------------------------------------------------
        ! Boundary
        n = 1
        !$omp target teams distribute parallel do default(none) &
        !$omp private(l) &
        !$omp shared(nlines,f,u,n,r1,r2,r3)
        do l = 1, nlines
            f(l, n) = f(l, n) + u(l, n)*r2(n) + u(l, n + 1)*r3(n) + u(l, n + 2)*r1(n)   ! r1(1) contains extended stencil
        end do
        !$omp end target teams distribute parallel do
        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        !$omp target teams distribute parallel do collapse(2) default(none) &
        !$omp private(n,l) &
        !$omp shared(nlines,f,u,nx,r1,r2,r3)
        do n = 2, nx - 1
            do l = 1, nlines
                f(l, n) = f(l, n) + u(l, n - 1)*r1(n) + u(l, n)*r2(n) + u(l, n + 1)*r3(n)
            end do
        end do
        !$omp end target teams distribute parallel do

        ! -------------------------------------------------------------------
        ! Boundary
        n = nx
        !$omp target teams distribute parallel do default(none) &
        !$omp private(l) &
        !$omp shared(nlines,f,u,n,r1,r2,r3)
        do l = 1, nlines
            f(l, n) = f(l, n) + u(l, n - 2)*r3(n) + u(l, n - 1)*r1(n) + u(l, n)*r2(n)   ! r3(nx) contains extended stencil
        end do
        !$omp end target teams distribute parallel do

#endif

        ! -----------------------------------------------------------------------
        ! Profiling
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_1)
        mat3dadd_time = mat3dadd_time + real(clock_1 - clock_0)/ clock_cycle 

        return
    end subroutine MatMul_3d_add

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_3d_antisym(nx, nlines, r1, r2, r3, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        integer(wi), intent(in) :: nx, nlines                       ! nlines linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx)              ! RHS diagonals
        real(wp), intent(in) :: u(nlines, nx)                       ! function u
        real(wp), intent(inout) :: f(nlines, nx)                    ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc
        real(wp), intent(in), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(out), optional :: bcs_b(nlines), bcs_t(nlines)

        ! -------------------------------------------------------------------
        integer(wi) n
        integer ibc_loc

        ! -------------------------------------------------------------------
        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx)
            f(:, 2) = u(:, 3) - u(:, 1)

        else
            if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc_loc)) then
                ! f(1) contains the boundary condition
                if (present(bcs_b)) bcs_b(:) = f(:, 1)*r2_b(1) + u(:, 2)*r3_b(1) + u(:, 3)*r1_b(1) ! r1(1) contains extended stencil
                f(:, 2) = f(:, 1)*r1_b(2) + u(:, 2)*r2_b(2) + u(:, 3)*r3_b(2)

            else
                f(:, 1) = u(:, 1)*r2(1) + u(:, 2)*r3(1) + u(:, 3)*r1(1) ! r1(1) contains extended stencil
                f(:, 2) = u(:, 1)*r1(2) + u(:, 2)*r2(2) + u(:, 3)*r3(2)

            end if

        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 3, nx - 2
            f(:, n) = u(:, n + 1) - u(:, n - 1)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2)
            f(:, nx) = u(:, 1) - u(:, nx - 1)

        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(nx) contains the boundary condition
                f(:, nx - 1) = u(:, nx - 2)*r1_t(1) + u(:, nx - 1)*r2_t(1) + f(:, nx)*r3_t(1)
                if (present(bcs_t)) bcs_t(:) = u(:, nx - 2)*r3_t(2) + u(:, nx - 1)*r1_t(2) + f(:, nx)*r2_t(2) ! r3(nx) contains extended stencil

            else
                f(:, nx - 1) = u(:, nx - 2)*r1(nx - 1) + u(:, nx - 1)*r2(nx - 1) + u(:, nx)*r3(nx - 1)
                f(:, nx) = u(:, nx - 2)*r3(nx) + u(:, nx - 1)*r1(nx) + u(:, nx)*r2(nx)  ! r3(nx) contains extended stencil

            end if

        end if

        return
    end subroutine MatMul_3d_antisym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_3d_sym(nx, nlines, r1, r2, r3, u, f, periodic, ibc)
        integer(wi), intent(in) :: nx, nlines                   ! nlines linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx)          ! RHS diagonals
        real(wp), intent(in) :: u(nlines, nx)                   ! function u
        real(wp), intent(out) :: f(nlines, nx)                  ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n
        integer ibc_loc

        ! -------------------------------------------------------------------
        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) + u(:, nx) + u(:, 1)*r2(1)

        else
            f(:, 1) = u(:, 1)*r2(1) + u(:, 2)*r3(1) + u(:, 3)*r1(1) ! r1(1) contains extended stencil

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = u(:, n + 1) + u(:, n - 1) + u(:, n)*r2(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, nx) = u(:, 1) + u(:, nx - 1) + u(:, nx)*r2(nx)

        else
            f(:, nx) = u(:, nx - 2)*r3(nx) + u(:, nx - 1)*r1(nx) + u(:, nx)*r2(nx) ! r3(nx) contains extended stencil

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_3d_sym

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is penta-diagonal and 1. upper-diagonal in interior points is equal to 1
    subroutine MatMul_5d(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u
        integer, intent(in), optional :: ibc
        real(wp), intent(in), optional :: rhs_b(1:4, 0:5), rhs_t(0:3, 1:6)  ! Special bcs at bottom and top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx

        integer ibc_loc

        ! #######################################################################
        nx = size(rhs, 1)
        ! print *, nd = size(rhs, 2)    ! # diagonals, should be 5
        ! size(u,2) and size(f,2) should be nx
        ! size(u,1) and size(f,1) and size(bcs_b) and size(bcs_t) should be the same (number of equations to solve)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        ! -------------------------------------------------------------------
        ! Boundary; the first 5/2+1+1=4 rows might be different
        if (any([BCS_MIN, BCS_BOTH] == ibc_loc)) then
            ! f(1) contains the boundary condition
            if (present(bcs_b)) bcs_b(:) = f(:, 1)*r3_b(1) + u(:, 2)*r4_b(1) + u(:, 3)*r5_b(1) + u(:, 4)*r1_b(1) ! r1(1) contains extended stencil
            f(:, 2) = f(:, 1)*r2_b(2) + u(:, 2)*r3_b(2) + u(:, 3)*r4_b(2) + u(:, 4)*r5_b(2)
            f(:, 3) = f(:, 1)*r1_b(3) + u(:, 2)*r2_b(3) + u(:, 3)*r3_b(3) + u(:, 4)*r4_b(3) + u(:, 5)*r5_b(3)
            f(:, 4) = f(:, 1)*r0_b(4) + u(:, 2)*r1_b(4) + u(:, 3)*r2_b(4) + u(:, 4)*r3_b(4) + u(:, 5)*r4_b(4) + u(:, 6)*r5_b(4)
        else
            f(:, 1) = u(:, 1)*r3_i(1) + u(:, 2)*r4_i(1) + u(:, 3)*r5_i(1) + u(:, 4)*r1_i(1)   ! r1(1) contains extended stencil
            f(:, 2) = u(:, 1)*r2_i(2) + u(:, 2)*r3_i(2) + u(:, 3)*r4_i(2) + u(:, 4)*r5_i(2)
            f(:, 3) = u(:, 1)*r1_i(3) + u(:, 2)*r2_i(3) + u(:, 3)*r3_i(3) + u(:, 4)*r4_i(3) + u(:, 5)*r5_i(3)
            f(:, 4) = u(:, 2)*r1_i(4) + u(:, 3)*r2_i(4) + u(:, 4)*r3_i(4) + u(:, 5)*r4_i(4) + u(:, 6)*r5_i(4)
        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 5, nx - 4
            f(:, n) = u(:, n - 2)*r1_i(n) + u(:, n - 1)*r2_i(n) + u(:, n)*r3_i(n) + u(:, n + 1) + u(:, n + 2)*r5_i(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary; the last 5/2+1+1=4 rows might be different
        if (any([BCS_MAX, BCS_BOTH] == ibc_loc)) then
            ! f(nx) contains the boundary condition
            f(:, nx - 3) = u(:, nx - 5)*r1_t(0) + u(:, nx - 4)*r2_t(0) + u(:, nx - 3)*r3_t(0) + u(:, nx - 2)*r4_t(0) + u(:, nx - 1)*r5_t(0) + f(:, nx)*r6_t(0)
            f(:, nx - 2) = u(:, nx - 4)*r1_t(1) + u(:, nx - 3)*r2_t(1) + u(:, nx - 2)*r3_t(1) + u(:, nx - 1)*r4_t(1) + f(:, nx)*r5_t(1)
            f(:, nx - 1) = u(:, nx - 3)*r1_t(2) + u(:, nx - 2)*r2_t(2) + u(:, nx - 1)*r3_t(2) + f(:, nx)*r4_t(2)
            if (present(bcs_t)) bcs_t(:) = u(:, nx - 3)*r5_t(3) + u(:, nx - 2)*r1_t(3) + u(:, nx - 1)*r2_t(3) + f(:, nx)*r3_t(3) ! r5(nx) contains extended stencil
        else
        f(:, nx - 3) = u(:, nx - 5)*r1_i(nx - 3) + u(:, nx - 4)*r2_i(nx - 3) + u(:, nx - 3)*r3_i(nx - 3) + u(:, nx - 2)*r4_i(nx - 3) + u(:, nx - 1)*r5_i(nx - 3)
            f(:, nx - 2) = u(:, nx - 4)*r1_i(nx - 2) + u(:, nx - 3)*r2_i(nx - 2) + u(:, nx - 2)*r3_i(nx - 2) + u(:, nx - 1)*r4_i(nx - 2) + u(:, nx)*r5_i(nx - 2)
            f(:, nx - 1) = u(:, nx - 3)*r1_i(nx - 1) + u(:, nx - 2)*r2_i(nx - 1) + u(:, nx - 1)*r3_i(nx - 1) + u(:, nx)*r4_i(nx - 1)
            f(:, nx) = u(:, nx - 3)*r5_i(nx) + u(:, nx - 2)*r1_i(nx) + u(:, nx - 1)*r2_i(nx) + u(:, nx)*r3_i(nx) ! r5(nx) contains extended stencil
        end if

        return
    end subroutine MatMul_5d

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = f + B u, assuming B is pentadiagonal
    subroutine MatMul_5d_add(nx, nlines, r1, r2, r3, r4, r5, u, f)
        integer(wi), intent(in) :: nx, nlines       ! nlines linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx)
        real(wp), intent(in) :: u(nlines, nx)       ! function u
        real(wp), intent(inout) :: f(nlines, nx)    ! RHS, f = B u

        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, 1) = f(:, 1) + u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1) + u(:, 4)*r1(1)   ! r1(1) contains extended stencil
        f(:, 2) = f(:, 2) + u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 3, nx - 2
            f(:, n) = f(:, n) + u(:, n - 2)*r1(n) + u(:, n - 1)*r2(n) + u(:, n)*r3(n) + u(:, n + 1)*r4(n) + u(:, n + 2)*r5(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, nx - 1) = f(:, nx - 1) + u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1)*r3(nx - 1) + u(:, nx)*r4(nx - 1)
        f(:, nx) = f(:, nx) + u(:, nx - 3)*r5(nx) + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)*r3(nx) ! r5(nx) contains extended stencil

        return
    end subroutine MatMul_5d_add

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric penta-diagonal with 1. upper-diagonal equal to 1
    ! It also assumes equal coefficients in the 2. upper-diagonal for the interior points
    subroutine MatMul_5d_antisym(nx, nlines, r1, r2, r3, r4, r5, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        use TLab_Time, only: mat5dantisym_time

        integer(wi), intent(in) :: nx, nlines          ! nlines linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(nlines, nx)          ! function u
        real(wp), intent(inout) :: f(nlines, nx)       ! RHS, f = B u; f_1 and f_n can contain neumann bcs
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc
        real(wp), intent(in), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(out), optional :: bcs_b(nlines), bcs_t(nlines)

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r5_loc     ! 2. upper-diagonal
        integer ibc_loc

! APU offloading 
#ifdef USE_APU
        integer(wi) l
#endif
        integer clock_0, clock_1, clock_cycle
        ! -------------------------------------------------------------------
        r5_loc = r5(4)      ! The first 3 equations, last 3 equations, can be normalized differently

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! -----------------------------------------------------------------------
        ! Profiling
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_0,clock_cycle) 

#ifndef USE_APU

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx) + r5_loc*(u(:, 3) - u(:, nx - 1))
            f(:, 2) = u(:, 3) - u(:, 1) + r5_loc*(u(:, 4) - u(:, nx))
            f(:, 3) = u(:, 4) - u(:, 2) + r5_loc*(u(:, 5) - u(:, 1))

        else
            if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc_loc)) then
                ! f(1) contains boundary condition
                if (present(bcs_b)) bcs_b(:) = f(:, 1)*r3_b(1) + u(:, 2)*r4_b(1) + u(:, 3)*r5_b(1) + u(:, 4)*r1_b(1) ! r1(1) with extended stencil
                f(:, 2) = f(:, 1)*r2_b(2) + u(:, 2)*r3_b(2) + u(:, 3)*r4_b(2) + u(:, 4)*r5_b(2)
                f(:, 3) = f(:, 1)*r1_b(3) + u(:, 2)*r2_b(3) + u(:, 3)*r3_b(3) + u(:, 4)*r4_b(3) + u(:, 5)*r5_b(3)

            else
                f(:, 1) = u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1) + u(:, 4)*r1(1)   ! r1(1) with extended stencil
                f(:, 2) = u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)
                f(:, 3) = u(:, 1)*r1(3) + u(:, 2)*r2(3) + u(:, 3)*r3(3) + u(:, 4)*r4(3) + u(:, 5)*r5(3)

            end if

        end if

        ! Interior points
        do n = 4, nx - 3
            f(:, n) = u(:, n + 1) - u(:, n - 1) + r5_loc*(u(:, n + 2) - u(:, n - 2))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 2) = u(:, nx - 1) - u(:, nx - 3) + r5_loc*(u(:, nx) - u(:, nx - 4))
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2) + r5_loc*(u(:, 1) - u(:, nx - 3))
            f(:, nx) = u(:, 1) - u(:, nx - 1) + r5_loc*(u(:, 2) - u(:, nx - 2))

        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(n) contains boundary condition
                f(:, nx - 2) = u(:, nx - 4)*r1_t(1) + u(:, nx - 3)*r2_t(1) + u(:, nx - 2)*r3_t(1) + u(:, nx - 1)*r4_t(1) + f(:, nx)*r5_t(1)
                f(:, nx - 1) = u(:, nx - 3)*r1_t(2) + u(:, nx - 2)*r2_t(2) + u(:, nx - 1)*r3_t(2) + f(:, nx)*r4_t(2)
                if (present(bcs_b)) bcs_t(:) = u(:, nx - 3)*r5_t(3) + u(:, nx - 2)*r1_t(3) + u(:, nx - 1)*r2_t(3) + f(:, nx)*r3_t(3) ! r5(nx) with extended stencil

            else
                f(:, nx - 2) = u(:, nx - 4)*r1(nx - 2) + u(:, nx - 3)*r2(nx - 2) + u(:, nx - 2)*r3(nx - 2) + u(:, nx - 1)*r4(nx - 2) + u(:, nx)*r5(nx - 2)
                f(:, nx - 1) = u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1)*r3(nx - 1) + u(:, nx)*r4(nx - 1)
                f(:, nx) = u(:, nx - 3)*r5(nx) + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)*r3(nx)! r5(nx) with extended stencil
            end if

        end if
        ! -------------------------------------------------------------------
        ! -----------------------------------------------------------------------
        ! With APU ACCELERATION 
        ! not possible to use rx_b preprocessors here!
        ! -----------------------------------------------------------------------
#else
        ! Boundary
        if (periodic) then
            !$omp target teams distribute parallel do default(none) &
            !$omp private(l) &
            !$omp shared(nlines,f,u,nx,r5_loc)
            do l = 1, nlines
                f(l, 1) = u(l, 2) - u(l, nx) + r5_loc*(u(l, 3) - u(l, nx - 1))
                f(l, 2) = u(l, 3) - u(l, 1 ) + r5_loc*(u(l, 4) - u(l, nx))
                f(l, 3) = u(l, 4) - u(l, 2 ) + r5_loc*(u(l, 5) - u(l, 1))
            end do
            !$omp end target teams distribute parallel do

        else
            if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc_loc)) then
                ! f(1) contains boundary condition
                if (present(bcs_b)) then
                    !$omp target teams distribute parallel do default(none) &
                    !$omp private(l) &
                    !$omp shared(nlines,bcs_b,f,u,rhs_b)
                    do l = 1, nlines
                        bcs_b(l) = f(l, 1)*rhs_b(1,3) + u(l, 2)*rhs_b(1,4) + u(l, 3)*rhs_b(1,5) + u(l, 4)*rhs_b(1,1) ! r1(1) with extended stencil
                        f(l, 2)  = f(l, 1)*rhs_b(2,2) + u(l, 2)*rhs_b(2,3) + u(l, 3)*rhs_b(2,4) + u(l, 4)*rhs_b(2,5)
                        f(l, 3)  = f(l, 1)*rhs_b(3,1) + u(l, 2)*rhs_b(3,2) + u(l, 3)*rhs_b(3,3) + u(l, 4)*rhs_b(3,4) + u(l, 5)*rhs_b(3,5)
                    end do
                    !$omp end target teams distribute parallel do
                end if

            else
                !$omp target teams distribute parallel do default(none) &
                !$omp private(l) &
                !$omp shared(nlines,f,u,r1,r2,r3,r4,r5)
                do l = 1, nlines
                    f(l, 1) = u(l, 1)*r3(1) + u(l, 2)*r4(1) + u(l, 3)*r5(1) + u(l, 4)*r1(1)   ! r1(1) with extended stencil
                    f(l, 2) = u(l, 1)*r2(2) + u(l, 2)*r3(2) + u(l, 3)*r4(2) + u(l, 4)*r5(2)
                    f(l, 3) = u(l, 1)*r1(3) + u(l, 2)*r2(3) + u(l, 3)*r3(3) + u(l, 4)*r4(3) + u(l, 5)*r5(3)
                end do
                !$omp end target teams distribute parallel do

            end if

        end if

        ! Interior points
        !$omp target teams distribute parallel do collapse(2) default(none) &
        !$omp private(n,l) &
        !$omp shared(nlines,f,u,nx,r5_loc)
        do l = 1, nlines
            do n = 4, nx - 3
                f(l, n) = u(l, n + 1) - u(l, n - 1) + r5_loc*(u(l, n + 2) - u(l, n - 2))
            end do
        end do
        !$omp end target teams distribute parallel do

        ! Boundary
        if (periodic) then
            !$omp target teams distribute parallel do default(none) &
            !$omp private(l) &
            !$omp shared(nlines,f,nx,u,r5_loc)
            do l = 1, nlines
                f(l, nx - 2) = u(l, nx - 1) - u(l, nx - 3) + r5_loc*(u(l, nx) - u(l, nx - 4))
                f(l, nx - 1) = u(l, nx    ) - u(l, nx - 2) + r5_loc*(u(l, 1 ) - u(l, nx - 3))
                f(l, nx    ) = u(l, 1     ) - u(l, nx - 1) + r5_loc*(u(l, 2 ) - u(l, nx - 2))
            end do
            !$omp end target teams distribute parallel do
        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(n) contains boundary condition
                !$omp target teams distribute parallel do default(none) &
                !$omp private(l) &
                !$omp shared(nlines,f,nx,u,rhs_t)
                do l = 1, nlines
                    f(l, nx - 2) = u(l, nx - 4)*rhs_t(1,1) + u(l, nx - 3)*rhs_t(1,2) + u(l, nx - 2)*rhs_t(1,3) + u(l, nx - 1)*rhs_t(1,4) + f(l, nx)*rhs_t(1,5)
                    f(l, nx - 1) = u(l, nx - 3)*rhs_t(2,1) + u(l, nx - 2)*rhs_t(2,2) + u(l, nx - 1)*rhs_t(2,3) + f(l, nx    )*rhs_t(2,4)
                end do
                !$omp end target teams distribute parallel do
                
                if (present(bcs_b)) then
                    !$omp target teams distribute parallel do default(none) &
                    !$omp private(l) &
                    !$omp shared(nlines,bcs_t,f,u,nx,rhs_t)
                    do l = 1, nlines
                        bcs_t(l) = u(l, nx - 3)*rhs_t(3,5) + u(l, nx - 2)*rhs_t(3,1) + u(l, nx - 1)*rhs_t(3,2) + f(l, nx)*rhs_t(3,3) ! r5(nx) with extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                end if

            else
                !$omp target teams distribute parallel do default(none) &
                !$omp private(l) &
                !$omp shared(nlines,f,u,nx,r1,r2,r3,r4,r5)
                do l = 1, nlines
                    f(l, nx - 2) = u(l, nx - 4)*r1(nx - 2) + u(l, nx - 3)*r2(nx - 2) + u(l, nx - 2)*r3(nx - 2) + u(l, nx - 1)*r4(nx - 2) + u(l, nx)*r5(nx - 2)
                    f(l, nx - 1) = u(l, nx - 3)*r1(nx - 1) + u(l, nx - 2)*r2(nx - 1) + u(l, nx - 1)*r3(nx - 1) + u(l, nx    )*r4(nx - 1)
                    f(l, nx    ) = u(l, nx - 3)*r5(nx    ) + u(l, nx - 2)*r1(nx    ) + u(l, nx - 1)*r2(nx    ) + u(l, nx    )*r3(nx    )! r5(nx) with extended stencil
                end do
                !$omp end target teams distribute parallel do
            end if

        end if

#endif

        ! -----------------------------------------------------------------------
        ! Profiling
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_1)
        mat5dantisym_time = mat5dantisym_time + real(clock_1 - clock_0)/ clock_cycle 

        return
    end subroutine MatMul_5d_antisym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_5d_sym(nx, nlines, r1, r2, r3, r4, r5, u, f, periodic, ibc)
        integer(wi), intent(in) :: nx, nlines       ! nlines linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(nlines, nx)       ! function u
        real(wp), intent(out) :: f(nlines, nx)      ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r3_loc     ! center diagonal
        real(wp) r5_loc     ! 2. upper-diagonal
        integer ibc_loc

        ! APU offloading 
#ifdef USE_APU
        integer(wi) l
#endif
        integer clock_0, clock_1, clock_cycle

        ! -------------------------------------------------------------------
        r5_loc = r5(3)      ! The first 2 equations, last 2 equations, are normalized differently
        r3_loc = r3(3)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! Boundary
        if (periodic) then
            f(:, 1) = r3_loc*u(:, 1) + u(:, 2) + u(:, nx) &
                      + r5_loc*(u(:, 3) + u(:, nx - 1))

            f(:, 2) = r3_loc*u(:, 2) + u(:, 3) + u(:, 1) &
                      + r5_loc*(u(:, 4) + u(:, nx))

        else
            f(:, 1) = u(:, 1)*r3(1) + u(:, 2)*r4(1) + u(:, 3)*r5(1) &
                      + u(:, 4)*r1(1)   ! r1(1) contains 3. upper-diagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r2(2) + u(:, 2)*r3(2) + u(:, 3)*r4(2) + u(:, 4)*r5(2)

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! Interior points
        do n = 3, nx - 2
            f(:, n) = r3_loc*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) + u(:, n - 2))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 1) = r3_loc*u(:, nx - 1) + u(:, nx) + u(:, nx - 2) &
                           + r5_loc*(u(:, 1) + u(:, nx - 3))

            f(:, nx) = r3_loc*u(:, nx) + u(:, 1) + u(:, nx - 1) &
                       + r5_loc*(u(:, 2) + u(:, nx - 2))

        else
            f(:, nx - 1) = u(:, nx - 3)*r1(nx - 1) + u(:, nx - 2)*r2(nx - 1) + u(:, nx - 1)*r3(nx - 1) &
                           + u(:, nx)*r4(nx - 1)

            f(:, nx) = u(:, nx - 3)*r5(nx) & ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 2)*r1(nx) + u(:, nx - 1)*r2(nx) + u(:, nx)*r3(nx)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_5d_sym

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric hepta-diagonal with 1. upper-diagonal equal to 1
    ! It also assumes equal coefficients in the 2. and 3. upper-diagonals for the interior points
    subroutine MatMul_7d_antisym(nx, nlines, r1, r2, r3, r4, r5, r6, r7, u, f, periodic, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        integer(wi), intent(in) :: nx, nlines          ! nlines linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx), r6(nx), r7(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(nlines, nx)          ! function u
        real(wp), intent(inout) :: f(nlines, nx)       ! RHS, f = B u; f_1 and f_n can contain neumann bcs
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc
        real(wp), intent(in), optional :: rhs_b(:, :), rhs_t(:, :)
        real(wp), intent(out), optional :: bcs_b(nlines), bcs_t(nlines)

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r6_loc, r7_loc     ! 2. and 3. upper-diagonal
        integer ibc_loc

        ! -------------------------------------------------------------------
        r6_loc = r6(5)      ! The first 4 equations, last 4 equations, can be normalized differently
        r7_loc = r7(5)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! -------------------------------------------------------------------
        ! Boundary
        if (periodic) then
            f(:, 1) = u(:, 2) - u(:, nx) + r6_loc*(u(:, 3) - u(:, nx - 1)) + r7_loc*(u(:, 4) - u(:, nx - 2))
            f(:, 2) = u(:, 3) - u(:, 1) + r6_loc*(u(:, 4) - u(:, nx)) + r7_loc*(u(:, 5) - u(:, nx - 1))
            f(:, 3) = u(:, 4) - u(:, 2) + r6_loc*(u(:, 5) - u(:, 1)) + r7_loc*(u(:, 6) - u(:, nx))
            f(:, 4) = u(:, 5) - u(:, 3) + r6_loc*(u(:, 6) - u(:, 2)) + r7_loc*(u(:, 7) - u(:, 1))

        else
            if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc_loc)) then
                ! f(1) contains boundary condition
                if (present(bcs_b)) bcs_b(:) = f(:, 1)*r4_b(1) + u(:, 2)*r5_b(1) + u(:, 3)*r6_b(1) + u(:, 4)*r7_b(1) + u(:, 5)*r1_b(1)   ! r1(1) with extended stencil
                f(:, 2) = f(:, 1)*r3_b(2) + u(:, 2)*r4_b(2) + u(:, 3)*r5_b(2) + u(:, 4)*r6_b(2) + u(:, 5)*r7_b(2)
                f(:, 3) = f(:, 1)*r2_b(3) + u(:, 2)*r3_b(3) + u(:, 3)*r4_b(3) + u(:, 4)*r5_b(3) + u(:, 5)*r6_b(3) + u(:, 6)*r7_b(3)
                f(:, 4) = f(:, 1)*r1_b(4) + u(:, 2)*r2_b(4) + u(:, 3)*r3_b(4) + u(:, 4)*r4_b(4) + u(:, 5)*r5_b(4) + u(:, 6)*r6_b(4) + u(:, 7)*r7_b(4)

            else
                f(:, 1) = u(:, 1)*r4(1) + u(:, 2)*r5(1) + u(:, 3)*r6(1) + u(:, 4)*r7(1) + u(:, 5)*r1(1)   ! r1(1) with extended stencil
                f(:, 2) = u(:, 1)*r3(2) + u(:, 2)*r4(2) + u(:, 3)*r5(2) + u(:, 4)*r6(2) + u(:, 5)*r7(2)
                f(:, 3) = u(:, 1)*r2(3) + u(:, 2)*r3(3) + u(:, 3)*r4(3) + u(:, 4)*r5(3) + u(:, 5)*r6(3) + u(:, 6)*r7(3)
                f(:, 4) = u(:, 1)*r1(4) + u(:, 2)*r2(4) + u(:, 3)*r3(4) + u(:, 4)*r4(4) + u(:, 5)*r5(4) + u(:, 6)*r6(4) + u(:, 7)*r7(4)

            end if

        end if

        ! Interior points
        do n = 5, nx - 4
            f(:, n) = u(:, n + 1) - u(:, n - 1) + r6_loc*(u(:, n + 2) - u(:, n - 2)) + r7_loc*(u(:, n + 3) - u(:, n - 3))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 3) = u(:, nx - 2) - u(:, nx - 4) + r6_loc*(u(:, nx - 1) - u(:, nx - 5)) + r7_loc*(u(:, nx) - u(:, nx - 6))
            f(:, nx - 2) = u(:, nx - 1) - u(:, nx - 3) + r6_loc*(u(:, nx) - u(:, nx - 4)) + r7_loc*(u(:, 1) - u(:, nx - 5))
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2) + r6_loc*(u(:, 1) - u(:, nx - 3)) + r7_loc*(u(:, 2) - u(:, nx - 4))
            f(:, nx) = u(:, 1) - u(:, nx - 1) + r6_loc*(u(:, 2) - u(:, nx - 2)) + r7_loc*(u(:, 3) - u(:, nx - 3))

        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc_loc)) then
                ! f(nx) contains boundary condition
                f(:, nx - 3) = u(:, nx - 6)*r1_t(1) + u(:, nx - 5)*r2_t(1) + u(:, nx - 4)*r3_t(1) + u(:, nx - 3)*r4_t(1) + u(:, nx - 2)*r5_t(1) + u(:, nx - 1)*r6_t(1)+ f(:, nx)*r7_t(1)
              f(:, nx - 2) = u(:, nx - 5)*r1_t(2) + u(:, nx - 4)*r2_t(2) + u(:, nx - 3)*r3_t(2) + u(:, nx - 2)*r4_t(2) + u(:, nx - 1)*r5_t(2) + f(:, nx)*r6_t(2)
                f(:, nx - 1) = u(:, nx - 4)*r1_t(3) + u(:, nx - 3)*r2_t(3) + u(:, nx - 2)*r3_t(3) + u(:, nx - 1)*r4_t(3) + f(:, nx)*r5_t(3)
                if (present(bcs_b)) bcs_t(:) = u(:, nx - 4)*r7_t(4) + u(:, nx - 3)*r1_t(4) + u(:, nx - 2)*r2_t(4) + u(:, nx - 1)*r3_t(4) + f(:, nx)*r4_t(4) ! r7(nx) with extended stencil

            else
                f(:, nx - 3) = u(:, nx - 6)*r1(nx - 3) + u(:, nx - 5)*r2(nx - 3) + u(:, nx - 4)*r3(nx - 3) + u(:, nx - 3)*r4(nx - 3) + u(:, nx - 2)*r5(nx - 3) + u(:, nx - 1)*r6(nx - 3)+ u(:, nx)*r7(nx - 3)
                f(:, nx - 2) = u(:, nx - 5)*r1(nx - 2) + u(:, nx - 4)*r2(nx - 2) + u(:, nx - 3)*r3(nx - 2) + u(:, nx - 2)*r4(nx - 2) + u(:, nx - 1)*r5(nx - 2) + u(:, nx)*r6(nx - 2)
                f(:, nx - 1) = u(:, nx - 4)*r1(nx - 1) + u(:, nx - 3)*r2(nx - 1) + u(:, nx - 2)*r3(nx - 1) + u(:, nx - 1)*r4(nx - 1) + u(:, nx)*r5(nx - 1)
                f(:, nx) = u(:, nx - 4)*r7(nx) + u(:, nx - 3)*r1(nx) + u(:, nx - 2)*r2(nx) + u(:, nx - 1)*r3(nx) + u(:, nx)*r4(nx) ! r7(nx) with extended stencil
            end if

        end if

        return
    end subroutine MatMul_7d_antisym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_7d_sym(nx, nlines, r1, r2, r3, r4, r5, r6, r7, u, f, periodic, ibc)
        integer(wi), intent(in) :: nx, nlines       ! nlines linear systems or size n
        real(wp), intent(in) :: r1(nx), r2(nx), r3(nx), r4(nx), r5(nx), r6(nx), r7(nx)  ! RHS diagonals
        real(wp), intent(in) :: u(nlines, nx)       ! function u
        real(wp), intent(out) :: f(nlines, nx)      ! RHS, f = B u
        logical, intent(in) :: periodic
        integer, intent(in), optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n
        real(wp) r4_loc     ! center diagonal
        real(wp) r6_loc     ! 2. upper-diagonal
        real(wp) r7_loc     ! 3. upper-diagonal
        integer ibc_loc

        ! -------------------------------------------------------------------
        r7_loc = r7(4)
        r6_loc = r6(4)      ! The first 3 equations, last 3 equations, are normalized differently
        r4_loc = r4(4)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_DD
        end if

        ! Boundary
        if (periodic) then
            f(:, 1) = r4_loc*u(:, 1) + u(:, 2) + u(:, nx) &
                      + r6_loc*(u(:, 3) + u(:, nx - 1)) &
                      + r7_loc*(u(:, 4) + u(:, nx - 2))

            f(:, 2) = r4_loc*u(:, 2) + u(:, 3) + u(:, 1) &
                      + r6_loc*(u(:, 4) + u(:, nx)) &
                      + r7_loc*(u(:, 5) + u(:, nx - 1))

            f(:, 3) = r4_loc*u(:, 3) + u(:, 4) + u(:, 2) &
                      + r6_loc*(u(:, 5) + u(:, 1)) &
                      + r7_loc*(u(:, 6) + u(:, nx))
        else
            f(:, 1) = u(:, 1)*r4(1) + u(:, 2)*r5(1) + u(:, 3)*r6(1) + u(:, 4)*r7(1) &
                      + u(:, 5)*r1(1)   ! r1(1) contains 4. upper-diagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r3(2) + u(:, 2)*r4(2) + u(:, 3)*r5(2) + u(:, 4)*r6(2) + u(:, 5)*r7(2)

            f(:, 3) = u(:, 1)*r2(3) + u(:, 2)*r3(3) + u(:, 3)*r4(3) + u(:, 4)*r5(3) + u(:, 5)*r6(3) + u(:, 6)*r7(3)

            if (any([BCS_ND, BCS_NN] == ibc_loc)) f(:, 1) = 0.0_wp

        end if

        ! Interior points
        do n = 4, nx - 3
            f(:, n) = r4_loc*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                      + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                      + r7_loc*(u(:, n + 3) + u(:, n - 3))
        end do

        ! Boundary
        if (periodic) then
            f(:, nx - 2) = r4_loc*u(:, nx - 2) + u(:, nx - 1) + u(:, nx - 3) &
                           + r6_loc*(u(:, nx) + u(:, nx - 4)) &
                           + r7_loc*(u(:, 1) + u(:, nx - 5))

            f(:, nx - 1) = r4_loc*u(:, nx - 1) + u(:, nx) + u(:, nx - 2) &
                           + r6_loc*(u(:, 1) + u(:, nx - 3)) &
                           + r7_loc*(u(:, 2) + u(:, nx - 4))

            f(:, nx) = r4_loc*u(:, nx) + u(:, 1) + u(:, nx - 1) &
                       + r6_loc*(u(:, 2) + u(:, nx - 2)) &
                       + r7_loc*(u(:, 3) + u(:, nx - 3))
        else
            f(:, nx - 2) = u(:, nx - 5)*r1(nx - 2) + u(:, nx - 4)*r2(nx - 2) + u(:, nx - 3)*r3(nx - 2) + u(:, nx - 2)*r4(nx - 2) + u(:, nx - 1)*r5(nx - 2) &
                           + u(:, nx)*r6(nx - 2)

            f(:, nx - 1) = u(:, nx - 4)*r1(nx - 1) + u(:, nx - 3)*r2(nx - 1) + u(:, nx - 2)*r3(nx - 1) + u(:, nx - 1)*r4(nx - 1) &
                           + u(:, nx)*r5(nx - 1)

            f(:, nx) = u(:, nx - 4)*r7(nx) & ! r7(nx) contains 4. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 3)*r1(nx) + u(:, nx - 2)*r2(nx) + u(:, nx - 1)*r3(nx) + u(:, nx)*r4(nx)

            if (any([BCS_DN, BCS_NN] == ibc_loc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_7d_sym

end module FDM_MatMul
