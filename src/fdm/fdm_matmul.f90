!########################################################################
! Matrix multiplication of n-diagonal matrix with a vector
! The boundary conditions can extend over n/2 + 2 points
! The 1. upper-diagonal in interior points is equal to 1
! This allows use to handle systems A y = B x in which A amd B differ by up to 2 diagonals (see notes)
!########################################################################
module FDM_MatMul
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN
    use TLab_Constants, only: BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Constants, only: BCS_PERIODIC
    use Tlab_Type, only: fdm_integral_dt, fdm_integral_dt2
    implicit none
    private

    ! generic cases
    public MatMul_3d            ! Calculate f = B u, assuming B is tridiagonal
    public MatMul_3d_APU        ! Calculate f = B u, assuming B is tridiagonal, APU offloading
    public MatMul_3d_add        ! Calculate f = f + B u, assuming B is tridiagonal
    public MatMul_3d_add_APU
    ! special cases where coefficients are constant in the interior points
    public MatMul_3d_antisym    ! Calculate f = B u, assuming B is tridiagonal, antisymmetric
    public MatMul_3d_antisym_APU
    public MatMul_3d_sym        ! Calculate f = B u, assuming B is tridiagonal, symmetric
    public MatMul_3d_sym_APU

    ! generic cases
    public MatMul_5d            ! Calculate f = B u, assuming B is pentadiagonal
    public MatMul_5d_APU        ! Calculate f = B u, assuming B is pentadiagonal, APU offloading
    public MatMul_5d_add        ! Calculate f = f + B u, assuming B is pentadiagonal
    public MatMul_5d_add_APU
    ! special cases where coefficients are constant in the interior points
    public MatMul_5d_antisym    ! Calculate f = B u, assuming B is pentadiagonal, antisymmetric
    public MatMul_5d_antisym_APU
    public MatMul_5d_sym        ! Calculate f = B u, assuming B is pentadiagonal, symmetric
    public MatMul_5d_sym_APU

    ! generic cases
    ! tbd when needed
    ! special cases where coefficients are constant in the interior points
    public MatMul_7d_antisym    ! Calculate f = B u, assuming B is heptadiagonal, antisymmetric
    public MatMul_7d_antisym_APU
    public MatMul_7d_sym        ! Calculate f = B u, assuming B is heptadiagonal, symmetric
    public MatMul_7d_sym_APU

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

#define r1_i(j) rhs(j,1)
#define r2_i(j) rhs(j,2)
#define r3_i(j) rhs(j,3)
#define r4_i(j) rhs(j,4)
#define r5_i(j) rhs(j,5)
#define r6_i(j) rhs(j,6)
#define r7_i(j) rhs(j,7)

#define r0b(k,i,j) fdmi%rhs_b(j,k,i,0)
#define r1b(k,i,j) fdmi%rhs_b(j,k,i,1)
#define r2b(k,i,j) fdmi%rhs_b(j,k,i,2)
#define r3b(k,i,j) fdmi%rhs_b(j,k,i,3)
#define r4b(k,i,j) fdmi%rhs_b(j,k,i,4)
#define r5b(k,i,j) fdmi%rhs_b(j,k,i,5)
#define r6b(k,i,j) fdmi%rhs_b(j,k,i,6)
#define r7b(k,i,j) fdmi%rhs_b(j,k,i,7)

#define r1t(k,i,j) fdmi%rhs_t(j,k,i,1)
#define r2t(k,i,j) fdmi%rhs_t(j,k,i,2)
#define r3t(k,i,j) fdmi%rhs_t(j,k,i,3)
#define r4t(k,i,j) fdmi%rhs_t(j,k,i,4)
#define r5t(k,i,j) fdmi%rhs_t(j,k,i,5)
#define r6t(k,i,j) fdmi%rhs_t(j,k,i,6)
#define r7t(k,i,j) fdmi%rhs_t(j,k,i,7)

contains
    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is tridiagonal and 1. upper-diagonal in interior points is equal to 1
    ! Special boundary conditions restricted to 3 points:
    ! r_11 r_12 r_13
    !      r_21 r_22 r_23
    !      r_30 r_31 r_32 r_33
    !                r_41  1.  r_43         <- interior points start here
    !                     ...  ...  ...
    subroutine MatMul_3d(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        ! real(wp), intent(in), optional :: rhs_b(1:3, 0:3), rhs_t(0:2, 1:4)
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx

        ! #######################################################################
        nx = size(rhs, 1)
        ! print *, nd = size(rhs, 2)    ! # diagonals, should be 3
        ! size(u,2) and size(f,2) should be nx
        ! size(u,1) and size(f,1) and size(bcs_b) and size(bcs_t) should be the same (number of equations to solve)

        ! -------------------------------------------------------------------
        ! Boundary; the first 3/2+1+1=3 rows might be different
        if (any([BCS_MIN, BCS_BOTH] == ibc)) then
            if (present(bcs_b)) bcs_b(:) = f(:, 1)*r2_b(1) + u(:, 2)*r3_b(1) + u(:, 3)*r1_b(1) ! r1(1) contains extended stencil
            ! f(1) contains the boundary condition
            f(:, 2) = f(:, 1)*r1_b(2) + u(:, 2)*r2_b(2) + u(:, 3)*r3_b(2)
            f(:, 3) = f(:, 1)*r0_b(3) + u(:, 2)*r1_b(3) + u(:, 3)*r2_b(3) + u(:, 4)*r3_b(3)
        else
            f(:, 1) = u(:, 1)*r2_i(1) + u(:, 2)*r3_i(1) + u(:, 3)*r1_i(1)   ! r1(1) contains extended stencil
            f(:, 2) = u(:, 1)*r1_i(2) + u(:, 2)*r2_i(2) + u(:, 3)*r3_i(2)
            f(:, 3) = u(:, 2)*r1_i(3) + u(:, 3)*r2_i(3) + u(:, 4)*r3_i(3)
        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 4, nx - 3
            f(:, n) = u(:, n - 1)*r1_i(n) + u(:, n)*r2_i(n) + u(:, n + 1)
        end do

        ! -------------------------------------------------------------------
        ! Boundary; the last 3/2+1+1=3 rows might be different
        if (any([BCS_MAX, BCS_BOTH] == ibc)) then
            ! f(nx) contains the boundary condition
            f(:, nx - 2) = u(:, nx - 3)*r1_t(0) + u(:, nx - 2)*r2_t(0) + u(:, nx - 1)*r3_t(0) + f(:, nx)*r4_t(0)
            f(:, nx - 1) = u(:, nx - 2)*r1_t(1) + u(:, nx - 1)*r2_t(1) + f(:, nx)*r3_t(1)
            if (present(bcs_t)) bcs_t(:) = u(:, nx - 2)*r3_t(2) + u(:, nx - 1)*r1_t(2) + f(:, nx)*r2_t(2) ! r3(nx) contains extended stencil
        else
            f(:, nx - 2) = u(:, nx - 3)*r1_i(nx - 2) + u(:, nx - 2)*r2_i(nx - 2) + u(:, nx - 1)*r3_i(nx - 2)
            f(:, nx - 1) = u(:, nx - 2)*r1_i(nx - 1) + u(:, nx - 1)*r2_i(nx - 1) + u(:, nx)*r3_i(nx - 1)
            f(:, nx) = u(:, nx - 2)*r3_i(nx) + u(:, nx - 1)*r1_i(nx) + u(:, nx)*r2_i(nx) ! r3(nx) contains extended stencil
        end if

        return
    end subroutine MatMul_3d

    subroutine MatMul_3d_APU(nlines, klines, ilines, nx, fdmi, rhs, u, f, ibc, bcs_b, bcs_t)
        use TLab_Time, only: mat3d_time
        integer(wi) nlines, ilines, klines, nx
        type(fdm_integral_dt2), intent(in) :: fdmi                          ! rhs_b(1:3, 0:3), rhs_t(0:2, 1:4)  ! Special bcs at bottom and top
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(1:2*nx, 1:klines, 1:ilines)                                     ! vector u
        real(wp), intent(out) :: f(1:2, 1:nx, 1:klines, 1:ilines)                                 ! vector f = B u
        integer, intent(in), optional :: ibc
        real(wp), intent(out), optional :: bcs_b(:, :, :), bcs_t(:, :, :)

        ! -------------------------------------------------------------------
        integer(wi) n, len, i, k, l
        integer(wi) pa, pb, pc, pd, pe, pf
        integer(wi) lp0, lp1, lp2, lp3, lp4, lp5, lp6, lp7
        ! APU offloading 

        integer clock_0, clock_1, clock_2, clock_cycle, clock_3
        ! -----------------------------------------------------------------------
        ! Profiling
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_0,clock_cycle)  ! delete after testing

        ! #######################################################################
        len = size(f,1)

        lp0 = 2*nx; lp1 = 2*nx - 1; lp2 = 2*nx - 2; lp3 = 2*nx - 3 
        lp5 = 2*nx - 5; lp4 = 2*nx - 4; lp7 = 2*nx - 7; lp6 = 2*nx - 6 

            ! -------------------------------------------------------------------
            ! Boundary; the first 3/2+1+1=3 rows might be different
            if (any([BCS_MIN, BCS_BOTH] == ibc)) then
                if (any([BCS_MAX, BCS_BOTH] == ibc)) then
                    if (present(bcs_b) .and. present(bcs_t)) then
                        !$omp target teams distribute parallel do collapse(2) defaultmap(present) &
                        !$omp private(i,k,n,pa,pb,pc,pd,pe,pf) &
                        !$omp shared(ilines,klines,nx,bcs_b,bcs_t,f,u,lp0,lp1,lp2,lp3,lp4,lp5,lp6,lp7,fdmi)
                        do i = 1, ilines
                            do k = 1, klines
                                bcs_b(:, k, i) = f(:, 1, k, i)*r2b(k, i, 1) + u(3:4, k, i)*r3b(k, i, 1) + u(5:6, k, i)*r1b(k, i, 1) ! r1(1) contains extended stencil
                                ! f(1) contains the boundary condition
                                f(:, 2, k, i) = f(:, 1, k, i)*r1b(k, i, 2) + u(3:4, k, i)*r2b(k, i, 2) + u(5:6, k, i)*r3b(k, i, 2)
                                f(:, 3, k, i) = f(:, 1, k, i)*r0b(k, i, 3) + u(3:4, k, i)*r1b(k, i, 3) + u(5:6, k, i)*r2b(k, i, 3) + u(7:8, k, i)*r3b(k, i, 3)
                                do n = 4, nx - 3
                                    pa = 2*n - 3; pb = 2*n - 2; pc = 2*n - 1; pd = 2*n; pe = 2*n + 1; pf = 2*n + 2
                                    f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)
                                end do
                                f(:, nx - 2, k, i) = u(lp7:lp6, k, i)*r1t(k, i, 0) + u(lp5:lp4, k, i)*r2t(k, i, 0) + u(lp3:lp2, k, i)*r3t(k, i, 0) + f(:, nx, k, i)*r4t(k, i, 0)
                                f(:, nx - 1, k, i) = u(lp5:lp4, k, i)*r1t(k, i, 1) + u(lp3:lp2, k, i)*r2t(k, i, 1) + f(:, nx, k, i)*r3t(k, i, 1)
                                bcs_t(:, k, i) = u(lp5:lp4, k, i)*r3t(k, i, 2) + u(lp3:lp2, k, i)*r1t(k, i, 2) + f(:, nx, k, i)*r2t(k, i, 2) ! r3(nx) contains extended stencil
                            end do
                        end do
                        !$omp end target teams distribute parallel do
                    else if (present(bcs_b)) then
                        !$omp target teams distribute parallel do collapse(2) &
                        !$omp private(i,k,n,pa,pb,pc,pd,pe,pf) &
                        !$omp shared(ilines,klines,nx,bcs_b,f,u,lp0,lp1,lp2,lp3,lp4,lp5,lp6,lp7,fdmi)
                        do i = 1, ilines
                            do k = 1, klines
                                bcs_b(:, k, i) = f(:, 1, k, i)*r2b(k, i, 1) + u(3:4, k, i)*r3b(k, i, 1) + u(5:6, k, i)*r1b(k, i, 1) ! r1(1) contains extended stencil
                                ! f(1) contains the boundary condition
                                f(:, 2, k, i) = f(:, 1, k, i)*r1b(k, i, 2) + u(3:4, k, i)*r2b(k, i, 2) + u(5:6, k, i)*r3b(k, i, 2)
                                f(:, 3, k, i) = f(:, 1, k, i)*r0b(k, i, 3) + u(3:4, k, i)*r1b(k, i, 3) + u(5:6, k, i)*r2b(k, i, 3) + u(7:8, k, i)*r3b(k, i, 3)
                                do n = 4, nx - 3
                                    pa = 2*n - 3; pb = 2*n - 2; pc = 2*n - 1; pd = 2*n; pe = 2*n + 1; pf = 2*n + 2
                                    f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)
                                end do
                                f(:, nx - 2, k, i) = u(lp7:lp6, k, i)*r1t(k, i, 0) + u(lp5:lp4, k, i)*r2t(k, i, 0) + u(lp3:lp2, k, i)*r3t(k, i, 0) + f(:, nx, k, i)*r4t(k, i, 0)
                                f(:, nx - 1, k, i) = u(lp5:lp4, k, i)*r1t(k, i, 1) + u(lp3:lp2, k, i)*r2t(k, i, 1) + f(:, nx, k, i)*r3t(k, i, 1)
                            end do
                        end do
                        !$omp end target teams distribute parallel do
                    else if (present(bcs_t)) then
                        ! f(1) contains the boundary condition
                        !$omp target teams distribute parallel do collapse(2) &
                        !$omp private(i,k,n,pa,pb,pc,pd,pe,pf) &
                        !$omp shared(ilines,klines,nx,bcs_t,f,u,lp0,lp1,lp2,lp3,lp4,lp5,lp6,lp7,fdmi)
                        do i = 1, ilines
                            do k = 1, klines
                                f(:, 2, k, i) = f(:, 1, k, i)*r1b(k, i, 2) + u(3:4, k, i)*r2b(k, i, 2) + u(5:6, k, i)*r3b(k, i, 2)
                                f(:, 3, k, i) = f(:, 1, k, i)*r0b(k, i, 3) + u(3:4, k, i)*r1b(k, i, 3) + u(5:6, k, i)*r2b(k, i, 3) + u(7:8, k, i)*r3b(k, i, 3)
                                do n = 4, nx - 3
                                    pa = 2*n - 3; pb = 2*n - 2; pc = 2*n - 1; pd = 2*n; pe = 2*n + 1; pf = 2*n + 2
                                    f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)
                                end do
                                f(:, nx - 2, k, i) = u(lp7:lp6, k, i)*r1t(k, i, 0) + u(lp5:lp4, k, i)*r2t(k, i, 0) + u(lp3:lp2, k, i)*r3t(k, i, 0) + f(:, nx, k, i)*r4t(k, i, 0)
                                f(:, nx - 1, k, i) = u(lp5:lp4, k, i)*r1t(k, i, 1) + u(lp3:lp2, k, i)*r2t(k, i, 1) + f(:, nx, k, i)*r3t(k, i, 1)
                                bcs_t(:, k, i) = u(lp5:lp4, k, i)*r3t(k, i, 2) + u(lp3:lp2, k, i)*r1t(k, i, 2) + f(:, nx, k, i)*r2t(k, i, 2) ! r3(nx) contains extended stencil
                            end do
                        end do
                        !$omp end target teams distribute parallel do
                    else
                        !$omp target teams distribute parallel do collapse(2) &
                        !$omp private(i,k,n,pa,pb,pc,pd,pe,pf) &
                        !$omp shared(ilines,klines,nx,f,u,lp0,lp1,lp2,lp3,lp4,lp5,lp6,lp7,fdmi)
                        do i = 1, ilines
                            do k = 1, klines
                                ! f(1) contains the boundary condition
                                f(:, 2, k, i) = f(:, 1, k, i)*r1b(k, i, 2) + u(3:4, k, i)*r2b(k, i, 2) + u(5:6, k, i)*r3b(k, i, 2)
                                f(:, 3, k, i) = f(:, 1, k, i)*r0b(k, i, 3) + u(3:4, k, i)*r1b(k, i, 3) + u(5:6, k, i)*r2b(k, i, 3) + u(7:8, k, i)*r3b(k, i, 3)
                                do n = 4, nx - 3
                                    pa = 2*n - 3; pb = 2*n - 2; pc = 2*n - 1; pd = 2*n; pe = 2*n + 1; pf = 2*n + 2
                                    f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)
                                end do
                                f(:, nx - 2, k, i) = u(lp7:lp6, k, i)*r1t(k, i, 0) + u(lp5:lp4, k, i)*r2t(k, i, 0) + u(lp3:lp2, k, i)*r3t(k, i, 0) + f(:, nx, k, i)*r4t(k, i, 0)
                                f(:, nx - 1, k, i) = u(lp5:lp4, k, i)*r1t(k, i, 1) + u(lp3:lp2, k, i)*r2t(k, i, 1) + f(:, nx, k, i)*r3t(k, i, 1)
                            end do
                        end do
                        !$omp end target teams distribute parallel do
                    end if
                else
                    if (present(bcs_b)) then
                        !$omp target teams distribute parallel do collapse(2) &
                        !$omp private(i,k,n,pa,pb,pc,pd,pe,pf) &
                        !$omp shared(ilines,klines,nx,bcs_b,f,u,lp0,lp1,lp2,lp3,lp4,lp5,lp6,lp7,fdmi)
                        do i = 1, ilines
                            do k = 1, klines
                                bcs_b(:, k, i) = f(:, 1, k, i)*r2b(k, i, 1) + u(3:4, k, i)*r3b(k, i, 1) + u(5:6, k, i)*r1b(k, i, 1) ! r1(1) contains extended stencil
                                ! f(1) contains the boundary condition
                                f(:, 2, k, i) = f(:, 1, k, i)*r1b(k, i, 2) + u(3:4, k, i)*r2b(k, i, 2) + u(5:6, k, i)*r3b(k, i, 2)
                                f(:, 3, k, i) = f(:, 1, k, i)*r0b(k, i, 3) + u(3:4, k, i)*r1b(k, i, 3) + u(5:6, k, i)*r2b(k, i, 3) + u(7:8, k, i)*r3b(k, i, 3)
                                do n = 4, nx - 3
                                    pa = 2*n - 3; pb = 2*n - 2; pc = 2*n - 1; pd = 2*n; pe = 2*n + 1; pf = 2*n + 2
                                    f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)
                                end do
                                f(:, nx - 2, k, i) = u(lp7:lp6, k, i)*r1_i(nx - 2) + u(lp5:lp4, k, i)*r2_i(nx - 2) + u(lp3:lp2, k, i)*r3_i(nx - 2)
                                f(:, nx - 1, k, i) = u(lp5:lp4, k, i)*r1_i(nx - 1) + u(lp3:lp2, k, i)*r2_i(nx - 1) + u(lp1:lp0, k, i)*r3_i(nx - 1)
                                f(:, nx, k, i) = u(lp5:lp4, k, i)*r3_i(nx) + u(lp3:lp2, k, i)*r1_i(nx) + u(lp1:lp0, k, i)*r2_i(nx) ! r3(nx) contains extended stencil
                            end do
                        end do
                        !$omp end target teams distribute parallel do
                    else
                        !$omp target teams distribute parallel do collapse(2) &
                        !$omp private(i,k,n,pa,pb,pc,pd,pe,pf) &
                        !$omp shared(ilines,klines,nx,f,u,lp0,lp1,lp2,lp3,lp4,lp5,lp6,lp7,rhs,fdmi)
                        do i = 1, ilines
                            do k = 1, klines
                            ! f(1) contains the boundary condition
                                f(:, 2, k, i) = f(:, 1, k, i)*r1b(k, i, 2) + u(3:4, k, i)*r2b(k, i, 2) + u(5:6, k, i)*r3b(k, i, 2)
                                f(:, 3, k, i) = f(:, 1, k, i)*r0b(k, i, 3) + u(3:4, k, i)*r1b(k, i, 3) + u(5:6, k, i)*r2b(k, i, 3) + u(7:8, k, i)*r3b(k, i, 3)
                                do n = 4, nx - 3
                                    pa = 2*n - 3; pb = 2*n - 2; pc = 2*n - 1; pd = 2*n; pe = 2*n + 1; pf = 2*n + 2
                                    f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)
                                end do
                                f(:, nx - 2, k, i) = u(lp7:lp6, k, i)*r1_i(nx - 2) + u(lp5:lp4, k, i)*r2_i(nx - 2) + u(lp3:lp2, k, i)*r3_i(nx - 2)
                                f(:, nx - 1, k, i) = u(lp5:lp4, k, i)*r1_i(nx - 1) + u(lp3:lp2, k, i)*r2_i(nx - 1) + u(lp1:lp0, k, i)*r3_i(nx - 1)
                                f(:, nx, k, i) = u(lp5:lp4, k, i)*r3_i(nx) + u(lp3:lp2, k, i)*r1_i(nx) + u(lp1:lp0, k, i)*r2_i(nx) ! r3(nx) contains extended stencil
                            end do
                        end do
                        !$omp end target teams distribute parallel do
                    end if
                end if
            else
                if (any([BCS_MAX, BCS_BOTH] == ibc)) then
                    if (present(bcs_t)) then
                        !$omp target teams distribute parallel do collapse(2) &
                        !$omp private(i,k,n,pa,pb,pc,pd,pe,pf) &
                        !$omp shared(ilines,klines,nx,bcs_t,f,u,lp0,lp1,lp2,lp3,lp4,lp5,lp6,lp7,rhs,fdmi)
                        do i = 1, ilines
                            do k = 1, klines
                                f(:, 1, k, i) = u(1:2, k, i)*r2_i(1) + u(3:4, k, i)*r3_i(1) + u(5:6, k, i)*r1_i(1)   ! r1(1) contains extended stencil
                                f(:, 2, k, i) = u(1:2, k, i)*r1_i(2) + u(3:4, k, i)*r2_i(2) + u(5:6, k, i)*r3_i(2)
                                f(:, 3, k, i) = u(3:4, k, i)*r1_i(3) + u(5:6, k, i)*r2_i(3) + u(7:8, k, i)*r3_i(3)
                                do n = 4, nx - 3
                                    pa = 2*n - 3; pb = 2*n - 2; pc = 2*n - 1; pd = 2*n; pe = 2*n + 1; pf = 2*n + 2
                                    f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)
                                end do
                                f(:, nx - 2, k, i) = u(lp7:lp6, k, i)*r1t(k, i, 0) + u(lp5:lp4, k, i)*r2t(k, i, 0) + u(lp3:lp2, k, i)*r3t(k, i, 0) + f(:, nx, k, i)*r4t(k, i, 0)
                                f(:, nx - 1, k, i) = u(lp5:lp4, k, i)*r1t(k, i, 1) + u(lp3:lp2, k, i)*r2t(k, i, 1) + f(:, nx, k, i)*r3t(k, i, 1)
                                bcs_t(:, k, i) = u(lp5:lp4, k, i)*r3t(k, i, 2) + u(lp3:lp2, k, i)*r1t(k, i, 2) + f(:, nx, k, i)*r2t(k, i, 2) ! r3(nx) contains extended stencil
                            end do
                        end do
                        !$omp end target teams distribute parallel do
                    else
                        !$omp target teams distribute parallel do collapse(2) &
                        !$omp private(i,k,n,pa,pb,pc,pd,pe,pf) &
                        !$omp shared(ilines,klines,rhs,nx,f,u,lp0,lp1,lp2,lp3,lp4,lp5,lp6,lp7,fdmi)
                        do i = 1, ilines
                            do k = 1, klines
                                f(:, 1, k, i) = u(1:2, k, i)*r2_i(1) + u(3:4, k, i)*r3_i(1) + u(5:6, k, i)*r1_i(1)   ! r1(1) contains extended stencil
                                f(:, 2, k, i) = u(1:2, k, i)*r1_i(2) + u(3:4, k, i)*r2_i(2) + u(5:6, k, i)*r3_i(2)
                                f(:, 3, k, i) = u(3:4, k, i)*r1_i(3) + u(5:6, k, i)*r2_i(3) + u(7:8, k, i)*r3_i(3)
                                do n = 4, nx - 3
                                    pa = 2*n - 3; pb = 2*n - 2; pc = 2*n - 1; pd = 2*n; pe = 2*n + 1; pf = 2*n + 2
                                    f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)
                                end do
                                f(:, nx - 2, k, i) = u(lp7:lp6, k, i)*r1t(k, i, 0) + u(lp5:lp4, k, i)*r2t(k, i, 0) + u(lp3:lp2, k, i)*r3t(k, i, 0) + f(:, nx, k, i)*r4t(k, i, 0)
                                f(:, nx - 1, k, i) = u(lp5:lp4, k, i)*r1t(k, i, 1) + u(lp3:lp2, k, i)*r2t(k, i, 1) + f(:, nx, k, i)*r3t(k, i, 1)
                            end do
                        end do
                        !$omp end target teams distribute parallel do
                    end if
                else
                    !$omp target teams distribute parallel do collapse(2) &
                    !$omp private(i,k,n,pa,pb,pc,pd,pe,pf) &
                    !$omp shared(ilines,klines,rhs,nx,f,u,lp0,lp1,lp2,lp3,lp4,lp5,lp6)
                    do i = 1, ilines
                        do k = 1, klines
                            f(:, 1, k, i) = u(1:2, k, i)*r2_i(1) + u(3:4, k, i)*r3_i(1) + u(5:6, k, i)*r1_i(1)   ! r1(1) contains extended stencil
                            f(:, 2, k, i) = u(1:2, k, i)*r1_i(2) + u(3:4, k, i)*r2_i(2) + u(5:6, k, i)*r3_i(2)
                            f(:, 3, k, i) = u(3:4, k, i)*r1_i(3) + u(5:6, k, i)*r2_i(3) + u(7:8, k, i)*r3_i(3)
                            do n = 4, nx - 3
                                pa = 2*n - 3; pb = 2*n - 2; pc = 2*n - 1; pd = 2*n; pe = 2*n + 1; pf = 2*n + 2
                                f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)
                            end do
                            f(:, nx - 2, k, i) = u(lp7:lp6, k, i)*r1_i(nx - 2) + u(lp5:lp4, k, i)*r2_i(nx - 2) + u(lp3:lp2, k, i)*r3_i(nx - 2)
                            f(:, nx - 1, k, i) = u(lp5:lp4, k, i)*r1_i(nx - 1) + u(lp3:lp2, k, i)*r2_i(nx - 1) + u(lp1:lp0, k, i)*r3_i(nx - 1)
                            f(:, nx, k, i) = u(lp5:lp4, k, i)*r3_i(nx) + u(lp3:lp2, k, i)*r1_i(nx) + u(lp1:lp0, k, i)*r2_i(nx) ! r3(nx) contains extended stencil
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end if
        ! -----------------------------------------------------------------------
        ! Profiling
        ! -----------------------------------------------------------------------
        call SYSTEM_CLOCK(clock_1, clock_cycle)
        mat3d_time = mat3d_time + real(clock_1 - clock_0)/ real(clock_cycle) 
        return
    end subroutine MatMul_3d_APU

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = f + B u, assuming B is tridiagonal
    subroutine MatMul_3d_add(rhs, u, f)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u

        ! -------------------------------------------------------------------
        integer(wi) n, nx

        ! -------------------------------------------------------------------
        nx = size(rhs, 1)

        ! Boundary
        n = 1
        f(:, n) = f(:, n) + u(:, n)*r2_i(n) + u(:, n + 1)*r3_i(n) + u(:, n + 2)*r1_i(n)   ! r1(1) contains extended stencil

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = f(:, n) + u(:, n - 1)*r1_i(n) + u(:, n)*r2_i(n) + u(:, n + 1)*r3_i(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        n = nx
        f(:, n) = f(:, n) + u(:, n - 2)*r3_i(n) + u(:, n - 1)*r1_i(n) + u(:, n)*r2_i(n)   ! r3(nx) contains extended stencil

        return
    end subroutine MatMul_3d_add

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = f + B u, assuming B is tridiagonal
    subroutine MatMul_3d_add_APU(rhs, u, f)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u

        ! -------------------------------------------------------------------
        integer(wi) n, nx, my, m
        real(wp) sum_f, sum_u
        sum_u = 0.0_wp
        sum_f = 0.0_wp
        ! -------------------------------------------------------------------
        nx = size(rhs, 1)
        my = size(f, 1)
        !$omp target teams distribute parallel do &
        !$omp private(m,n) &
        !$omp shared(my,nx,f,u,rhs)
        do m = 1, my
            ! Boundary
            n = 1
            f(m, n) = f(m, n) + u(m, n)*r2_i(n) + u(m, n + 1)*r3_i(n) + u(m, n + 2)*r1_i(n)   ! r1(1) contains extended stencil
            ! -------------------------------------------------------------------
            ! Interior points; accelerate
            do n = 2, nx - 1
                f(m, n) = f(m, n) + u(m, n - 1)*r1_i(n) + u(m, n)*r2_i(n) + u(m, n + 1)*r3_i(n)
            end do
            ! -------------------------------------------------------------------
            ! Boundary
            n = nx
            f(m, n) = f(m, n) + u(m, n - 2)*r3_i(n) + u(m, n - 1)*r1_i(n) + u(m, n)*r2_i(n)   ! r3(nx) contains extended stencil
        end do
        !$omp end target teams distribute parallel do
        return
    end subroutine MatMul_3d_add_APU

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_3d_antisym_APU(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx, m, my
        real(wp) sum_f, sum_u
        ! -------------------------------------------------------------------
        nx = size(rhs, 1)
        my = size(f, 1)
        sum_f = 0.0_wp
        sum_u = 0.0_wp
        ! -------------------------------------------------------------------
        ! Boundary
        if (ibc == BCS_PERIODIC) then
            !$omp target teams distribute parallel do &
            !$omp private(m,n) &
            !$omp shared(my,nx,u,f)
            do m = 1, my
                f(m, 1) = u(m, 2) - u(m, nx)
                f(m, 2) = u(m, 3) - u(m, 1)
                do n = 3, nx - 2
                    f(m, n) = u(m, n + 1) - u(m, n - 1)
                end do
                f(m, nx - 1) = u(m, nx) - u(m, nx - 2)
                f(m, nx) = u(m, 1) - u(m, nx - 1)
            end do
            !$omp end target teams distribute parallel do
        else if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc)) then
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc)) then
                ! f(1) contains the boundary condition
                if (present(bcs_b)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs_b,rhs_t,bcs_b,bcs_t)
                    do m = 1, my
                        bcs_b(m) = f(m, 1)*r2_b(1) + u(m, 2)*r3_b(1) + u(m, 3)*r1_b(1) ! r1(1) contains extended stencil
                        f(m, 2) = f(m, 1)*r1_b(2) + u(m, 2)*r2_b(2) + u(m, 3)*r3_b(2)
                        do n = 3, nx - 2
                            f(m, n) = u(m, n + 1) - u(m, n - 1)
                        end do
                        f(m, nx - 1) = u(m, nx - 2)*r1_t(1) + u(m, nx - 1)*r2_t(1) + f(m, nx)*r3_t(1)
                        bcs_t(m) = u(m, nx - 2)*r3_t(2) + u(m, nx - 1)*r1_t(2) + f(m, nx)*r2_t(2) ! r3(nx) contains extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs_b,rhs_t)
                    do m = 1, my
                        f(m, 2) = f(m, 1)*r1_b(2) + u(m, 2)*r2_b(2) + u(m, 3)*r3_b(2)
                        do n = 3, nx - 2
                            f(m, n) = u(m, n + 1) - u(m, n - 1)
                        end do
                        f(m, nx - 1) = u(m, nx - 2)*r1_t(1) + u(m, nx - 1)*r2_t(1) + f(m, nx)*r3_t(1)
                    end do
                    !$omp end target teams distribute parallel do
                end if
            else
                if (present(bcs_b)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs_b,rhs,bcs_b)
                    do m = 1, my
                        bcs_b(m) = f(m, 1)*r2_b(1) + u(m, 2)*r3_b(1) + u(m, 3)*r1_b(1) ! r1(1) contains extended stencil
                        f(m, 2) = f(m, 1)*r1_b(2) + u(m, 2)*r2_b(2) + u(m, 3)*r3_b(2)
                        do n = 3, nx - 2
                            f(m, n) = u(m, n + 1) - u(m, n - 1)
                        end do
                        f(m, nx - 1) = u(m, nx - 2)*r1_i(nx - 1) + u(m, nx - 1)*r2_i(nx - 1) + u(m, nx)*r3_i(nx - 1)
                        f(m, nx) = u(m, nx - 2)*r3_i(nx) + u(m, nx - 1)*r1_i(nx) + u(m, nx)*r2_i(nx)  ! r3(nx) contains extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs_b,rhs)
                    do m = 1, my
                        f(m, 2) = f(m, 1)*r1_b(2) + u(m, 2)*r2_b(2) + u(m, 3)*r3_b(2)
                        do n = 3, nx - 2
                            f(m, n) = u(m, n + 1) - u(m, n - 1)
                        end do
                        f(m, nx - 1) = u(m, nx - 2)*r1_i(nx - 1) + u(m, nx - 1)*r2_i(nx - 1) + u(m, nx)*r3_i(nx - 1)
                        f(m, nx) = u(m, nx - 2)*r3_i(nx) + u(m, nx - 1)*r1_i(nx) + u(m, nx)*r2_i(nx)  ! r3(nx) contains extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end if
        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc)) then
                if (present(bcs_t)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs_t,rhs,bcs_t)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r2_i(1) + u(m, 2)*r3_i(1) + u(m, 3)*r1_i(1) ! r1(1) contains extended stencil
                        f(m, 2) = u(m, 1)*r1_i(2) + u(m, 2)*r2_i(2) + u(m, 3)*r3_i(2)
                        do n = 3, nx - 2
                            f(m, n) = u(m, n + 1) - u(m, n - 1)
                        end do
                        f(m, nx - 1) = u(m, nx - 2)*r1_t(1) + u(m, nx - 1)*r2_t(1) + f(m, nx)*r3_t(1)
                        bcs_t(m) = u(m, nx - 2)*r3_t(2) + u(m, nx - 1)*r1_t(2) + f(m, nx)*r2_t(2) ! r3(nx) contains extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs_t,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r2_i(1) + u(m, 2)*r3_i(1) + u(m, 3)*r1_i(1) ! r1(1) contains extended stencil
                        f(m, 2) = u(m, 1)*r1_i(2) + u(m, 2)*r2_i(2) + u(m, 3)*r3_i(2)
                        do n = 3, nx - 2
                            f(m, n) = u(m, n + 1) - u(m, n - 1)
                        end do
                        f(m, nx - 1) = u(m, nx - 2)*r1_t(1) + u(m, nx - 1)*r2_t(1) + f(m, nx)*r3_t(1)
                    end do
                    !$omp end target teams distribute parallel do
                end if
            else
                !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs)
                do m = 1, my
                    f(m, 1) = u(m, 1)*r2_i(1) + u(m, 2)*r3_i(1) + u(m, 3)*r1_i(1) ! r1(1) contains extended stencil
                    f(m, 2) = u(m, 1)*r1_i(2) + u(m, 2)*r2_i(2) + u(m, 3)*r3_i(2)
                    do n = 3, nx - 2
                        f(m, n) = u(m, n + 1) - u(m, n - 1)
                    end do
                    f(m, nx - 1) = u(m, nx - 2)*r1_i(nx - 1) + u(m, nx - 1)*r2_i(nx - 1) + u(m, nx)*r3_i(nx - 1)
                    f(m, nx) = u(m, nx - 2)*r3_i(nx) + u(m, nx - 1)*r1_i(nx) + u(m, nx)*r2_i(nx)  ! r3(nx) contains extended stencil
                end do
                !$omp end target teams distribute parallel do
            end if
        end if
        return
    end subroutine MatMul_3d_antisym_APU

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_3d_antisym(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx

        ! -------------------------------------------------------------------
        nx = size(rhs, 1)

        ! -------------------------------------------------------------------
        ! Boundary
        if (ibc == BCS_PERIODIC) then
            f(:, 1) = u(:, 2) - u(:, nx)
            f(:, 2) = u(:, 3) - u(:, 1)

        else if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc)) then
            ! f(1) contains the boundary condition
            if (present(bcs_b)) bcs_b(:) = f(:, 1)*r2_b(1) + u(:, 2)*r3_b(1) + u(:, 3)*r1_b(1) ! r1(1) contains extended stencil
            f(:, 2) = f(:, 1)*r1_b(2) + u(:, 2)*r2_b(2) + u(:, 3)*r3_b(2)

        else
            f(:, 1) = u(:, 1)*r2_i(1) + u(:, 2)*r3_i(1) + u(:, 3)*r1_i(1) ! r1(1) contains extended stencil
            f(:, 2) = u(:, 1)*r1_i(2) + u(:, 2)*r2_i(2) + u(:, 3)*r3_i(2)

        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 3, nx - 2
            f(:, n) = u(:, n + 1) - u(:, n - 1)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        if (ibc == BCS_PERIODIC) then
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2)
            f(:, nx) = u(:, 1) - u(:, nx - 1)

        else if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc)) then
            ! f(nx) contains the boundary condition
            f(:, nx - 1) = u(:, nx - 2)*r1_t(1) + u(:, nx - 1)*r2_t(1) + f(:, nx)*r3_t(1)
            if (present(bcs_t)) bcs_t(:) = u(:, nx - 2)*r3_t(2) + u(:, nx - 1)*r1_t(2) + f(:, nx)*r2_t(2) ! r3(nx) contains extended stencil

        else
            f(:, nx - 1) = u(:, nx - 2)*r1_i(nx - 1) + u(:, nx - 1)*r2_i(nx - 1) + u(:, nx)*r3_i(nx - 1)
            f(:, nx) = u(:, nx - 2)*r3_i(nx) + u(:, nx - 1)*r1_i(nx) + u(:, nx)*r2_i(nx)  ! r3(nx) contains extended stencil

        end if

        return
    end subroutine MatMul_3d_antisym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_3d_sym(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx

        ! -------------------------------------------------------------------
        nx = size(rhs, 1)

        ! -------------------------------------------------------------------
        ! Boundary
        if (ibc == BCS_PERIODIC) then
            f(:, 1) = u(:, 2) + u(:, nx) + u(:, 1)*r2_i(1)

        else
            f(:, 1) = u(:, 1)*r2_i(1) + u(:, 2)*r3_i(1) + u(:, 3)*r1_i(1) ! r1(1) contains extended stencil

            if (any([BCS_ND, BCS_NN] == ibc)) f(:, 1) = 0.0_wp

        end if

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 2, nx - 1
            f(:, n) = u(:, n + 1) + u(:, n - 1) + u(:, n)*r2_i(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        if (ibc == BCS_PERIODIC) then
            f(:, nx) = u(:, 1) + u(:, nx - 1) + u(:, nx)*r2_i(nx)

        else
            f(:, nx) = u(:, nx - 2)*r3_i(nx) + u(:, nx - 1)*r1_i(nx) + u(:, nx)*r2_i(nx) ! r3(nx) contains extended stencil

            if (any([BCS_DN, BCS_NN] == ibc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_3d_sym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_3d_sym_APU(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx, m, my
        ! -------------------------------------------------------------------
        nx = size(rhs, 1)
        my = size(f, 1)
   
        ! -------------------------------------------------------------------
        ! Boundary
        if (ibc == BCS_PERIODIC) then
            !$omp target teams distribute parallel do &
            !$omp private(m,n) &
            !$omp shared(my,nx,u,f,rhs)
            do m = 1, my
                f(m, 1) = u(m, 2) + u(m, nx) + u(m, 1)*r2_i(1)
                do n = 2, nx - 1
                    f(m, n) = u(m, n + 1) + u(m, n - 1) + u(m, n)*r2_i(n)
                end do
                f(m, nx) = u(m, 1) + u(m, nx - 1) + u(m, nx)*r2_i(nx)
            end do
            !$omp end target teams distribute parallel do
        else
            if (any([BCS_ND, BCS_NN] == ibc)) then
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r2_i(1) + u(m, 2)*r3_i(1) + u(m, 3)*r1_i(1) ! r1(1) contains extended stencil
                        f(m, 1) = 0.0_wp
                        do n = 2, nx - 1
                            f(m, n) = u(m, n + 1) + u(m, n - 1) + u(m, n)*r2_i(n)
                        end do
                        f(m, nx) = u(m, nx - 2)*r3_i(nx) + u(m, nx - 1)*r1_i(nx) + u(m, nx)*r2_i(nx) ! r3(nx) contains extended stencil
                        f(m, nx) = 0.0_wp
                    end do
                    !$omp end target teams distribute parallel do
                else 
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r2_i(1) + u(m, 2)*r3_i(1) + u(m, 3)*r1_i(1) ! r1(1) contains extended stencil
                        f(m, 1) = 0.0_wp
                        do n = 2, nx - 1
                            f(m, n) = u(m, n + 1) + u(m, n - 1) + u(m, n)*r2_i(n)
                        end do
                        f(m, nx) = u(m, nx - 2)*r3_i(nx) + u(m, nx - 1)*r1_i(nx) + u(m, nx)*r2_i(nx) ! r3(nx) contains extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                end if
            else
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r2_i(1) + u(m, 2)*r3_i(1) + u(m, 3)*r1_i(1) ! r1(1) contains extended stencil
                        do n = 2, nx - 1
                            f(m, n) = u(m, n + 1) + u(m, n - 1) + u(m, n)*r2_i(n)
                        end do
                        f(m, nx) = u(m, nx - 2)*r3_i(nx) + u(m, nx - 1)*r1_i(nx) + u(m, nx)*r2_i(nx) ! r3(nx) contains extended stencil
                        f(m, nx) = 0.0_wp
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,u,f,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r2_i(1) + u(m, 2)*r3_i(1) + u(m, 3)*r1_i(1) ! r1(1) contains extended stencil
                        do n = 2, nx - 1
                            f(m, n) = u(m, n + 1) + u(m, n - 1) + u(m, n)*r2_i(n)
                        end do
                        f(m, nx) = u(m, nx - 2)*r3_i(nx) + u(m, nx - 1)*r1_i(nx) + u(m, nx)*r2_i(nx) ! r3(nx) contains extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end if
        end if
        return
    end subroutine MatMul_3d_sym_APU

    ! #######################################################################
    ! Calculate f = B u, assuming B is pentadiagonal and 1. upper-diagonal in interior points is equal to 1
    subroutine MatMul_5d_APU(nlines, klines, ilines, nx, fdmi, rhs, u, f, ibc, bcs_b, bcs_t)
        integer(wi) nlines, ilines, klines, nx
        type(fdm_integral_dt2), intent(in) :: fdmi   !rhs_b(1:3, 0:3), rhs_t(0:2, 1:4)  ! Special bcs at bottom and top
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(in) :: u(1:2*nx, 1:klines, 1:ilines)           ! vector u
        real(wp), intent(out) :: f(1:2, 1:nx, 1:klines, 1:ilines)       ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(out), optional :: bcs_b(:,:,:), bcs_t(:,:,:)

        ! -------------------------------------------------------------------
        integer(wi) n, len, i, k
        integer(wi) pa, pb, pc, pd, pe, pf, pg, ph, pi, pj
        integer(wi) lp0, lp1, lp2, lp3, lp4, lp5, lp6, lp7, lp8, lp9, lp10, lp11, lp12
        integer ibc_loc
        ! #######################################################################
        ! print *, nd = size(rhs, 2)    ! # diagonals, should be 5
        ! size(u,2) and size(f,2) should be nx
        ! size(u,1) and size(f,1) and size(bcs_b) and size(bcs_t) should be the same (number of equations to solve)
        lp0 = 2*nx; lp1 = 2*nx - 1; lp2 = 2*nx - 2; lp3 = 2*nx - 3; lp4 = 2*nx - 4; lp5 = 2*nx - 5; lp6 = 2*nx - 6; 
        lp7 = 2*nx - 7; lp8 = 2*nx - 8; lp9 = 2*nx - 9; lp10 = 2*nx - 10; lp11 = 2*nx - 11; lp12 = 2*nx - 12
        ! -------------------------------------------------------------------
        ! Boundary; the first 5/2+1+1=4 rows might be different
        if (any([BCS_MIN, BCS_BOTH] == ibc)) then
            ! f(1) contains the boundary condition
            if (present(bcs_b)) then
                !$omp target teams distribute parallel do collapse(2) &
                !$omp private(i,k,n,pa,pb,pc,pd,pe,pf,pg,ph,pi,pj) &
                !$omp shared(ilines,klines,nx,bcs_b,f,u,fdmi,rhs)
                do i = 1, ilines
                    do k = 1, klines
                        bcs_b(:, k, i) = f(:, 1, k, i)*r3b(k, i, 1) + u(3:4, k, i)*r4b(k, i, 1) + u(5:6, k, i)*r5b(k, i, 1) + u(7:8, k, i)*r1b(k, i, 1) ! r1(1) contains extended stencil
                        f(:, 2, k, i) = f(:, 1, k, i)*r2b(k, i, 2) + u(3:4, k, i)*r3b(k, i, 2) + u(5:6, k, i)*r4b(k, i, 2) + u(7:8, k, i)*r5b(k, i, 2)
                        f(:, 3, k, i) = f(:, 1, k, i)*r1b(k, i, 3) + u(3:4, k, i)*r2b(k, i, 3) + u(5:6, k, i)*r3b(k, i, 3) + u(7:8, k, i)*r4b(k, i, 3) + u(9:10, k, i)*r5b(k, i, 3)
                        f(:, 4, k, i) = f(:, 1, k, i)*r0b(k, i, 4)  + u(3:4, k, i)*r1b(k, i, 4) + u(5:6, k, i)*r2b(k, i, 4) + u(7:8, k, i)*r3b(k, i, 4) + u(9:10, k, i)*r4b(k, i, 4) + u(11:12, k, i)*r5b(k, i, 4)
                        ! -------------------------------------------------------------------
                        ! Interior points; accelerate
                        do n = 5, nx - 4
                            pa = 2*n - 5; pb = 2*n - 4; pc = 2*n - 3; pd = 2*n - 2; pe = 2*n - 1; pf = 2*n; pg = 2*n + 1; ph = 2*n + 2; pi = 2*n + 3; pj = 2*n + 4
                            f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)*r3_i(n) + u(pg:ph, k, i) + u(pi:pj, k, i)*r5_i(n)
                        end do
                    end do
                end do  
                !$omp end target teams distribute parallel do  
                ! -------------------------------------------------------------------
            else
                !$omp target teams distribute parallel do collapse(2) &
                !$omp private(i,k,n,pa,pb,pc,pd,pe,pf,pg,ph,pi,pj) &
                !$omp shared(ilines,klines,nx,f,u,fdmi,rhs)
                do i = 1, ilines
                    do k = 1, klines
                        f(:, 2, k, i) = f(:, 1, k, i)*r2b(k, i, 2) + u(3:4, k, i)*r3b(k, i, 2) + u(5:6, k, i)*r4b(k, i, 2) + u(7:8, k, i)*r5b(k, i, 2)
                        f(:, 3, k, i) = f(:, 1, k, i)*r1b(k, i, 3) + u(3:4, k, i)*r2b(k, i, 3) + u(5:6, k, i)*r3b(k, i, 3) + u(7:8, k, i)*r4b(k, i, 3) + u(9:10, k, i)*r5b(k, i, 3)
                        f(:, 4, k, i) = f(:, 1, k, i)*r0b(k, i, 4)  + u(3:4, k, i)*r1b(k, i, 4) + u(5:6, k, i)*r2b(k, i, 4) + u(7:8, k, i)*r3b(k, i, 4) + u(9:10, k, i)*r4b(k, i, 4) + u(11:12, k, i)*r5b(k, i, 4)
                        ! -------------------------------------------------------------------
                        ! Interior points; accelerate
                        do n = 5, nx - 4
                            f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)*r3_i(n) + u(pg:ph, k, i) + u(pi:pj, k, i)*r5_i(n)
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                    ! -------------------------------------------------------------------
            end if
        else
            !$omp target teams distribute parallel do collapse(2) &
            !$omp private(i,k,n,pa,pb,pc,pd,pe,pf,pg,ph,pi,pj) &
            !$omp shared(ilines,klines,nx,f,u,rhs)
            do i = 1, ilines
                do k = 1, klines
                    f(:, 1, k, i) = u(1:2, k, i)*r3_i(1) + u(3:4, k, i)*r4_i(1) + u(5:6, k, i)*r5_i(1) + u(7:8, k, i)*r1_i(1)   ! r1(1) contains extended stencil
                    f(:, 2, k, i) = u(1:2, k, i)*r2_i(2) + u(3:4, k, i)*r3_i(2) + u(5:6, k, i)*r4_i(2) + u(7:8, k, i)*r5_i(2)
                    f(:, 3, k, i) = u(1:2, k, i)*r1_i(3) + u(3:4, k, i)*r2_i(3) + u(5:6, k, i)*r3_i(3) + u(7:8, k, i)*r4_i(3) + u(9:10, k, i)*r5_i(3)
                    f(:, 4, k, i) = u(3:4, k, i)*r1_i(4) + u(5:6, k, i)*r2_i(4) + u(7:8, k, i)*r3_i(4) + u(9:10, k, i)*r4_i(4) + u(11:12, k, i)*r5_i(4)
                    ! -------------------------------------------------------------------
                    ! Interior points; accelerate
                    do n = 5, nx - 4
                        f(:, n, k, i) = u(pa:pb, k, i)*r1_i(n) + u(pc:pd, k, i)*r2_i(n) + u(pe:pf, k, i)*r3_i(n) + u(pg:ph, k, i) + u(pi:pj, k, i)*r5_i(n)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            ! -------------------------------------------------------------------
        end if       
        ! Boundary; the last 5/2+1+1=4 rows might be different
        if (any([BCS_MAX, BCS_BOTH] == ibc)) then
            ! f(nx) contains the boundary condition
            if (present(bcs_t)) then
                !$omp target teams distribute parallel do collapse(2) &
                !$omp private(i,k) &
                !$omp shared(ilines,klines,nx,lp1,lp2,lp3,lp4,lp5,lp6,lp7,lp8,lp9,lp10,lp11,bcs_t,f,u,fdmi)
                do i = 1, ilines
                    do k = 1, klines
                        f(:, nx-3, k, i) = u(lp11:lp10, k, i)*r1t(i, k, 0) + u(lp9:lp8, k, i)*r2t(i, k, 0) + u(lp7:lp6, k, i)*r3t(i, k, 0) + u(lp5:lp4, k, i)*r4t(i, k, 0) + u(lp3:lp2, k, i)*r5t(i, k, 0) + f(:, nx, k, i)*r6t(i, k, 0)
                        f(:, nx-2, k, i) = u(lp9:lp8, k, i)*r1t(i, k, 1) + u(lp7:lp6, k, i)*r2t(i, k, 1) + u(lp5:lp4, k, i)*r3t(i, k, 1) + u(lp3:lp2, k, i)*r4t(i, k, 1) + f(:, nx, k, i)*r5t(i, k, 1)
                        f(:, nx-1, k, i) = u(lp7:lp6, k, i)*r1t(i, k, 2) + u(lp5:lp4, k, i)*r2t(i, k, 2) + u(lp3:lp2, k, i)*r3t(i, k, 2) + f(:, nx, k, i)*r4t(i, k, 2)
                        bcs_t(:, k, i) = u(lp7:lp6, k, i)*r5t(i, k, 3) + u(lp5:lp4, k, i)*r1t(i, k, 3) + u(lp3:lp2, k, i)*r2t(i, k, 3) + f(:, nx, k, i)*r3t(i, k, 3) ! r5(nx) contains extended stencil
                    end do
                end do
                !$omp end target teams distribute parallel do
            else
                !$omp target teams distribute parallel do collapse(2) &
                !$omp private(i,k) &
                !$omp shared(ilines,klines,nx,lp1,lp2,lp3,lp4,lp5,lp6,lp7,lp8,lp9,lp10,lp11,f,u,fdmi)
                do i = 1, ilines
                    do k = 1, klines
                        f(:, nx-3, k, i) = u(lp11:lp10, k, i)*r1t(i, k, 0) + u(lp9:lp8, k, i)*r2t(i, k, 0) + u(lp7:lp6, k, i)*r3t(i, k, 0) + u(lp5:lp4, k, i)*r4t(i, k, 0) + u(lp3:lp2, k, i)*r5t(i, k, 0) + f(:, nx, k, i)*r6t(i, k, 0)
                        f(:, nx-2, k, i) = u(lp9:lp8, k, i)*r1t(i, k, 1) + u(lp7:lp6, k, i)*r2t(i, k, 1) + u(lp5:lp4, k, i)*r3t(i, k, 1) + u(lp3:lp2, k, i)*r4t(i, k, 1) + f(:, nx, k, i)*r5t(i, k, 1)
                        f(:, nx, k, i) = u(lp7:lp6, k, i)*r1t(i, k, 2) + u(lp5:lp4, k, i)*r2t(i, k, 2) + u(lp3:lp2, k, i)*r3t(i, k, 2) + f(:, nx, k, i)*r4t(i, k, 2)
                    end do
                end do
                !$omp end target teams distribute parallel do
            end if
        else
            !$omp target teams distribute parallel do collapse(2) &
            !$omp private(i,k) &
            !$omp shared(ilines,klines,nx,lp1,lp2,lp3,lp4,lp5,lp6,lp7,lp8,lp9,lp10,lp11,f,u,rhs)
            do i = 1, ilines
                do k = 1, klines
                    f(:, nx-3, k, i) = u(lp11:lp10, k, i)*r1_i(nx - 3) + u(lp9:lp8, k, i)*r2_i(nx - 3) + u(lp7:lp6, k, i)*r3_i(nx - 3) + u(lp5:lp4, k, i)*r4_i(nx - 3) + u(lp3:lp2, k, i)*r5_i(nx - 3)
                    f(:, nx-2, k, i) = u(lp9:lp8, k, i)*r1_i(nx - 2) + u(lp7:lp6, k, i)*r2_i(nx - 2) + u(lp5:lp4, k, i)*r3_i(nx - 2) + u(lp3:lp2, k, i)*r4_i(nx - 2) + u(lp1:lp0, k, nx)*r5_i(nx - 2)
                    f(:, nx-1, k, i) = u(lp7:lp6, k, i)*r1_i(nx - 1) + u(lp5:lp4, k, i)*r2_i(nx - 1) + u(lp3:lp2, k, i)*r3_i(nx - 1) + u(lp1:lp0, k, nx)*r4_i(nx - 1)
                    f(:, nx, k, i) = u(lp7:lp6, k, i)*r5_i(nx) + u(lp5:lp4, k, i)*r1_i(nx) + u(lp3:lp2, k, i)*r2_i(nx) + u(lp1:lp0, k, nx)*r3_i(nx) ! r5(nx) contains extended stencil
                end do
            end do
            !$omp end target teams distribute parallel do
        end if

        return
    end subroutine MatMul_5d_APU
    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is pentadiagonal and 1. upper-diagonal in interior points is equal to 1
    subroutine MatMul_5d(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) sum_f, sum_u

        ! #######################################################################
        nx = size(rhs, 1)
        ! print *, nd = size(rhs, 2)    ! # diagonals, should be 5
        ! size(u,2) and size(f,2) should be nx
        ! size(u,1) and size(f,1) and size(bcs_b) and size(bcs_t) should be the same (number of equations to solve)

        ! -------------------------------------------------------------------
        ! Boundary; the first 5/2+1+1=4 rows might be different
        if (any([BCS_MIN, BCS_BOTH] == ibc)) then
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
        if (any([BCS_MAX, BCS_BOTH] == ibc)) then
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
    subroutine MatMul_5d_add(rhs, u, f)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u

        ! -------------------------------------------------------------------
        integer(wi) n, nx

        ! #######################################################################
        nx = size(rhs, 1)

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, 1) = f(:, 1) + u(:, 1)*r3_i(1) + u(:, 2)*r4_i(1) + u(:, 3)*r5_i(1) + u(:, 4)*r1_i(1)   ! r1(1) contains extended stencil
        f(:, 2) = f(:, 2) + u(:, 1)*r2_i(2) + u(:, 2)*r3_i(2) + u(:, 3)*r4_i(2) + u(:, 4)*r5_i(2)

        ! -------------------------------------------------------------------
        ! Interior points; accelerate
        do n = 3, nx - 2
            f(:, n) = f(:, n) + u(:, n - 2)*r1_i(n) + u(:, n - 1)*r2_i(n) + u(:, n)*r3_i(n) + u(:, n + 1)*r4_i(n) + u(:, n + 2)*r5_i(n)
        end do

        ! -------------------------------------------------------------------
        ! Boundary
        f(:, nx - 1) = f(:, nx - 1) + u(:, nx - 3)*r1_i(nx - 1) + u(:, nx - 2)*r2_i(nx - 1) + u(:, nx - 1)*r3_i(nx - 1) + u(:, nx)*r4_i(nx - 1)
        f(:, nx) = f(:, nx) + u(:, nx - 3)*r5_i(nx) + u(:, nx - 2)*r1_i(nx) + u(:, nx - 1)*r2_i(nx) + u(:, nx)*r3_i(nx) ! r5(nx) contains extended stencil

        return
    end subroutine MatMul_5d_add

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = f + B u, assuming B is pentadiagonal
    subroutine MatMul_5d_add_APU(rhs, u, f)
        real(wp), intent(in) :: rhs(:, :)                                   ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                     ! vector u
        real(wp), intent(out) :: f(:, :)                                    ! vector f = B u

        ! -------------------------------------------------------------------
        integer(wi) n, nx, m, my
        ! #######################################################################
        nx = size(rhs, 1)
        my = size(f, 1)
        !$omp target teams distribute parallel do &
        !$omp private(m,n) &
        !$omp shared(my,nx,f,u,rhs)
        do m = 1, my
            ! -------------------------------------------------------------------
            ! Boundary
            f(m, 1) = f(m, 1) + u(m, 1)*r3_i(1) + u(m, 2)*r4_i(1) + u(m, 3)*r5_i(1) + u(m, 4)*r1_i(1)   ! r1(1) contains extended stencil
            f(m, 2) = f(m, 2) + u(m, 1)*r2_i(2) + u(m, 2)*r3_i(2) + u(m, 3)*r4_i(2) + u(m, 4)*r5_i(2)

            ! -------------------------------------------------------------------
            ! Interior points; accelerate
            do n = 3, nx - 2
                f(m, n) = f(m, n) + u(m, n - 2)*r1_i(n) + u(m, n - 1)*r2_i(n) + u(m, n)*r3_i(n) + u(m, n + 1)*r4_i(n) + u(m, n + 2)*r5_i(n)
            end do

            ! -------------------------------------------------------------------
            ! Boundary
            f(m, nx - 1) = f(m, nx - 1) + u(m, nx - 3)*r1_i(nx - 1) + u(m, nx - 2)*r2_i(nx - 1) + u(m, nx - 1)*r3_i(nx - 1) + u(m, nx)*r4_i(nx - 1)
            f(m, nx) = f(m, nx) + u(m, nx - 3)*r5_i(nx) + u(m, nx - 2)*r1_i(nx) + u(m, nx - 1)*r2_i(nx) + u(m, nx)*r3_i(nx) ! r5(nx) contains extended stencil
        end do
        !$omp end target teams distribute parallel do
        return
    end subroutine MatMul_5d_add_APU

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric pentadiagonal with 1. upper-diagonal equal to 1
    ! It also assumes equal coefficients in the 2. upper-diagonal for the interior points
    subroutine MatMul_5d_antisym(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r5_loc     ! 2. upper-diagonal

        ! #######################################################################
        nx = size(rhs, 1)
        r5_loc = r5_i(4)      ! The first 3 equations, last 3 equations, can be normalized differently

        ! -------------------------------------------------------------------
        ! Boundary
        if (ibc == BCS_PERIODIC) then
            f(:, 1) = u(:, 2) - u(:, nx) + r5_loc*(u(:, 3) - u(:, nx - 1))
            f(:, 2) = u(:, 3) - u(:, 1) + r5_loc*(u(:, 4) - u(:, nx))
            f(:, 3) = u(:, 4) - u(:, 2) + r5_loc*(u(:, 5) - u(:, 1))

        else if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc)) then
            ! f(1) contains boundary condition
            if (present(bcs_b)) bcs_b(:) = f(:, 1)*r3_b(1) + u(:, 2)*r4_b(1) + u(:, 3)*r5_b(1) + u(:, 4)*r1_b(1) ! r1(1) with extended stencil
            f(:, 2) = f(:, 1)*r2_b(2) + u(:, 2)*r3_b(2) + u(:, 3)*r4_b(2) + u(:, 4)*r5_b(2)
            f(:, 3) = f(:, 1)*r1_b(3) + u(:, 2)*r2_b(3) + u(:, 3)*r3_b(3) + u(:, 4)*r4_b(3) + u(:, 5)*r5_b(3)

        else
            f(:, 1) = u(:, 1)*r3_i(1) + u(:, 2)*r4_i(1) + u(:, 3)*r5_i(1) + u(:, 4)*r1_i(1)   ! r1(1) with extended stencil
            f(:, 2) = u(:, 1)*r2_i(2) + u(:, 2)*r3_i(2) + u(:, 3)*r4_i(2) + u(:, 4)*r5_i(2)
            f(:, 3) = u(:, 1)*r1_i(3) + u(:, 2)*r2_i(3) + u(:, 3)*r3_i(3) + u(:, 4)*r4_i(3) + u(:, 5)*r5_i(3)

        end if

        ! Interior points
        do n = 4, nx - 3
            f(:, n) = u(:, n + 1) - u(:, n - 1) + r5_loc*(u(:, n + 2) - u(:, n - 2))
        end do

        ! Boundary
        if (ibc == BCS_PERIODIC) then
            f(:, nx - 2) = u(:, nx - 1) - u(:, nx - 3) + r5_loc*(u(:, nx) - u(:, nx - 4))
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2) + r5_loc*(u(:, 1) - u(:, nx - 3))
            f(:, nx) = u(:, 1) - u(:, nx - 1) + r5_loc*(u(:, 2) - u(:, nx - 2))

        else if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc)) then
            ! f(n) contains boundary condition
            f(:, nx - 2) = u(:, nx - 4)*r1_t(1) + u(:, nx - 3)*r2_t(1) + u(:, nx - 2)*r3_t(1) + u(:, nx - 1)*r4_t(1) + f(:, nx)*r5_t(1)
            f(:, nx - 1) = u(:, nx - 3)*r1_t(2) + u(:, nx - 2)*r2_t(2) + u(:, nx - 1)*r3_t(2) + f(:, nx)*r4_t(2)
            if (present(bcs_b)) bcs_t(:) = u(:, nx - 3)*r5_t(3) + u(:, nx - 2)*r1_t(3) + u(:, nx - 1)*r2_t(3) + f(:, nx)*r3_t(3) ! r5(nx) with extended stencil

        else
            f(:, nx - 2) = u(:, nx - 4)*r1_i(nx - 2) + u(:, nx - 3)*r2_i(nx - 2) + u(:, nx - 2)*r3_i(nx - 2) + u(:, nx - 1)*r4_i(nx - 2) + u(:, nx)*r5_i(nx - 2)
            f(:, nx - 1) = u(:, nx - 3)*r1_i(nx - 1) + u(:, nx - 2)*r2_i(nx - 1) + u(:, nx - 1)*r3_i(nx - 1) + u(:, nx)*r4_i(nx - 1)
            f(:, nx) = u(:, nx - 3)*r5_i(nx) + u(:, nx - 2)*r1_i(nx) + u(:, nx - 1)*r2_i(nx) + u(:, nx)*r3_i(nx)! r5(nx) with extended stencil
        end if

        return
    end subroutine MatMul_5d_antisym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_5d_antisym_APU(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx, m, my
        real(wp) r5_loc     ! 2. upper-diagonal

        ! #######################################################################
        nx = size(rhs, 1)
        my = size(f, 1)
        r5_loc = r5_i(4)      ! The first 3 equations, last 3 equations, can be normalized differently

        ! -------------------------------------------------------------------
        ! Boundary
        if (ibc == BCS_PERIODIC) then
            !$omp target teams distribute parallel do &
            !$omp private(m,n) &
            !$omp shared(my,nx,r5_loc,f,u)
            do m = 1, my
                f(m, 1) = u(m, 2) - u(m, nx) + r5_loc*(u(m, 3) - u(m, nx - 1))
                f(m, 2) = u(m, 3) - u(m, 1) + r5_loc*(u(m, 4) - u(m, nx))
                f(m, 3) = u(m, 4) - u(m, 2) + r5_loc*(u(m, 5) - u(m, 1))
                ! Interior points
                do n = 4, nx - 3
                    f(m, n) = u(m, n + 1) - u(m, n - 1) + r5_loc*(u(m, n + 2) - u(m, n - 2))
                end do
                f(m, nx - 2) = u(m, nx - 1) - u(m, nx - 3) + r5_loc*(u(m, nx) - u(m, nx - 4))
                f(m, nx - 1) = u(m, nx) - u(m, nx - 2) + r5_loc*(u(m, 1) - u(m, nx - 3))
                f(m, nx) = u(m, 1) - u(m, nx - 1) + r5_loc*(u(m, 2) - u(m, nx - 2))
            end do
            !$omp end target teams distribute parallel do
        else if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc)) then
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc)) then
                ! f(1) contains boundary condition
                if (present(bcs_b) .and. present(bcs_t)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,rhs_b,rhs_t,f,u,bcs_b,bcs_t,r5_loc)
                    do m = 1, my
                        bcs_b(m) = f(m, 1)*r3_b(1) + u(m, 2)*r4_b(1) + u(m, 3)*r5_b(1) + u(m, 4)*r1_b(1) ! r1(1) with extended stencil
                        f(m, 2) = f(m, 1)*r2_b(2) + u(m, 2)*r3_b(2) + u(m, 3)*r4_b(2) + u(m, 4)*r5_b(2)
                        f(m, 3) = f(m, 1)*r1_b(3) + u(m, 2)*r2_b(3) + u(m, 3)*r3_b(3) + u(m, 4)*r4_b(3) + u(m, 5)*r5_b(3)
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r5_loc*(u(m, n + 2) - u(m, n - 2))
                        end do
                        f(m, nx - 2) = u(m, nx - 4)*r1_t(1) + u(m, nx - 3)*r2_t(1) + u(m, nx - 2)*r3_t(1) + u(m, nx - 1)*r4_t(1) + f(m, nx)*r5_t(1)
                        f(m, nx - 1) = u(m, nx - 3)*r1_t(2) + u(m, nx - 2)*r2_t(2) + u(m, nx - 1)*r3_t(2) + f(m, nx)*r4_t(2)
                        bcs_t(m) = u(m, nx - 3)*r5_t(3) + u(m, nx - 2)*r1_t(3) + u(m, nx - 1)*r2_t(3) + f(m, nx)*r3_t(3) ! r5(nx) with extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                else if (present(bcs_b)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,rhs_b,rhs_t,f,u,bcs_b,r5_loc)
                    do m = 1, my
                        bcs_b(m) = f(m, 1)*r3_b(1) + u(m, 2)*r4_b(1) + u(m, 3)*r5_b(1) + u(m, 4)*r1_b(1) ! r1(1) with extended stencil
                        f(m, 2) = f(m, 1)*r2_b(2) + u(m, 2)*r3_b(2) + u(m, 3)*r4_b(2) + u(m, 4)*r5_b(2)
                        f(m, 3) = f(m, 1)*r1_b(3) + u(m, 2)*r2_b(3) + u(m, 3)*r3_b(3) + u(m, 4)*r4_b(3) + u(m, 5)*r5_b(3)
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r5_loc*(u(m, n + 2) - u(m, n - 2))
                        end do
                        f(m, nx - 2) = u(m, nx - 4)*r1_t(1) + u(m, nx - 3)*r2_t(1) + u(m, nx - 2)*r3_t(1) + u(m, nx - 1)*r4_t(1) + f(m, nx)*r5_t(1)
                        f(m, nx - 1) = u(m, nx - 3)*r1_t(2) + u(m, nx - 2)*r2_t(2) + u(m, nx - 1)*r3_t(2) + f(m, nx)*r4_t(2)
                    end do
                    !$omp end target teams distribute parallel do
                else if (present(bcs_t)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,rhs_b,rhs_t,f,u,r5_loc)
                    do m = 1, my
                        f(m, 2) = f(m, 1)*r2_b(2) + u(m, 2)*r3_b(2) + u(m, 3)*r4_b(2) + u(m, 4)*r5_b(2)
                        f(m, 3) = f(m, 1)*r1_b(3) + u(m, 2)*r2_b(3) + u(m, 3)*r3_b(3) + u(m, 4)*r4_b(3) + u(m, 5)*r5_b(3)
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r5_loc*(u(m, n + 2) - u(m, n - 2))
                        end do
                        f(m, nx - 2) = u(m, nx - 4)*r1_t(1) + u(m, nx - 3)*r2_t(1) + u(m, nx - 2)*r3_t(1) + u(m, nx - 1)*r4_t(1) + f(m, nx)*r5_t(1)
                        f(m, nx - 1) = u(m, nx - 3)*r1_t(2) + u(m, nx - 2)*r2_t(2) + u(m, nx - 1)*r3_t(2) + f(m, nx)*r4_t(2)
                        bcs_t(m) = u(m, nx - 3)*r5_t(3) + u(m, nx - 2)*r1_t(3) + u(m, nx - 1)*r2_t(3) + f(m, nx)*r3_t(3) ! r5(nx) with extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,rhs_b,rhs_t,f,u,r5_loc)
                    do m = 1, my
                        f(m, 2) = f(m, 1)*r2_b(2) + u(m, 2)*r3_b(2) + u(m, 3)*r4_b(2) + u(m, 4)*r5_b(2)
                        f(m, 3) = f(m, 1)*r1_b(3) + u(m, 2)*r2_b(3) + u(m, 3)*r3_b(3) + u(m, 4)*r4_b(3) + u(m, 5)*r5_b(3)
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r5_loc*(u(m, n + 2) - u(m, n - 2))
                        end do
                        f(m, nx - 2) = u(m, nx - 4)*r1_t(1) + u(m, nx - 3)*r2_t(1) + u(m, nx - 2)*r3_t(1) + u(m, nx - 1)*r4_t(1) + f(m, nx)*r5_t(1)
                        f(m, nx - 1) = u(m, nx - 3)*r1_t(2) + u(m, nx - 2)*r2_t(2) + u(m, nx - 1)*r3_t(2) + f(m, nx)*r4_t(2)
                    end do
                    !$omp end target teams distribute parallel do
                end if 
            else
                if (present(bcs_b)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,rhs_b,rhs_t,f,u,bcs_b,r5_loc)
                    do m = 1, my
                        bcs_b(m) = f(m, 1)*r3_b(1) + u(m, 2)*r4_b(1) + u(m, 3)*r5_b(1) + u(m, 4)*r1_b(1) ! r1(1) with extended stencil
                        f(m, 2) = f(m, 1)*r2_b(2) + u(m, 2)*r3_b(2) + u(m, 3)*r4_b(2) + u(m, 4)*r5_b(2)
                        f(m, 3) = f(m, 1)*r1_b(3) + u(m, 2)*r2_b(3) + u(m, 3)*r3_b(3) + u(m, 4)*r4_b(3) + u(m, 5)*r5_b(3)
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r5_loc*(u(m, n + 2) - u(m, n - 2))
                        end do
                        ! top Boundary
                        f(m, nx - 2) = u(m, nx - 4)*r1_i(nx - 2) + u(m, nx - 3)*r2_i(nx - 2) + u(m, nx - 2)*r3_i(nx - 2) + u(m, nx - 1)*r4_i(nx - 2) + u(m, nx)*r5_i(nx - 2)
                        f(m, nx - 1) = u(m, nx - 3)*r1_i(nx - 1) + u(m, nx - 2)*r2_i(nx - 1) + u(m, nx - 1)*r3_i(nx - 1) + u(m, nx)*r4_i(nx - 1)
                        f(m, nx) = u(m, nx - 3)*r5_i(nx) + u(m, nx - 2)*r1_i(nx) + u(m, nx - 1)*r2_i(nx) + u(m, nx)*r3_i(nx)! r5(nx) with extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,rhs_b,rhs,f,u,r5_loc)
                    do m = 1, my
                        f(m, 2) = f(m, 1)*r2_b(2) + u(m, 2)*r3_b(2) + u(m, 3)*r4_b(2) + u(m, 4)*r5_b(2)
                        f(m, 3) = f(m, 1)*r1_b(3) + u(m, 2)*r2_b(3) + u(m, 3)*r3_b(3) + u(m, 4)*r4_b(3) + u(m, 5)*r5_b(3)
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r5_loc*(u(m, n + 2) - u(m, n - 2))
                        end do
                        ! top Boundary
                        f(m, nx - 2) = u(m, nx - 4)*r1_i(nx - 2) + u(m, nx - 3)*r2_i(nx - 2) + u(m, nx - 2)*r3_i(nx - 2) + u(m, nx - 1)*r4_i(nx - 2) + u(m, nx)*r5_i(nx - 2)
                        f(m, nx - 1) = u(m, nx - 3)*r1_i(nx - 1) + u(m, nx - 2)*r2_i(nx - 1) + u(m, nx - 1)*r3_i(nx - 1) + u(m, nx)*r4_i(nx - 1)
                        f(m, nx) = u(m, nx - 3)*r5_i(nx) + u(m, nx - 2)*r1_i(nx) + u(m, nx - 1)*r2_i(nx) + u(m, nx)*r3_i(nx)! r5(nx) with extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end if
        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc)) then
                if (present(bcs_t)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,rhs,rhs_t,f,u,bcs_t,r5_loc)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r3_i(1) + u(m, 2)*r4_i(1) + u(m, 3)*r5_i(1) + u(m, 4)*r1_i(1)   ! r1(1) with extended stencil
                        f(m, 2) = u(m, 1)*r2_i(2) + u(m, 2)*r3_i(2) + u(m, 3)*r4_i(2) + u(m, 4)*r5_i(2)
                        f(m, 3) = u(m, 1)*r1_i(3) + u(m, 2)*r2_i(3) + u(m, 3)*r3_i(3) + u(m, 4)*r4_i(3) + u(m, 5)*r5_i(3)
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r5_loc*(u(m, n + 2) - u(m, n - 2))
                        end do
                        f(m, nx - 2) = u(m, nx - 4)*r1_t(1) + u(m, nx - 3)*r2_t(1) + u(m, nx - 2)*r3_t(1) + u(m, nx - 1)*r4_t(1) + f(m, nx)*r5_t(1)
                        f(m, nx - 1) = u(m, nx - 3)*r1_t(2) + u(m, nx - 2)*r2_t(2) + u(m, nx - 1)*r3_t(2) + f(m, nx)*r4_t(2)
                        bcs_t(m) = u(m, nx - 3)*r5_t(3) + u(m, nx - 2)*r1_t(3) + u(m, nx - 1)*r2_t(3) + f(m, nx)*r3_t(3) ! r5(nx) with extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,rhs,rhs_t,f,u,r5_loc)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r3_i(1) + u(m, 2)*r4_i(1) + u(m, 3)*r5_i(1) + u(m, 4)*r1_i(1)   ! r1(1) with extended stencil
                        f(m, 2) = u(m, 1)*r2_i(2) + u(m, 2)*r3_i(2) + u(m, 3)*r4_i(2) + u(m, 4)*r5_i(2)
                        f(m, 3) = u(m, 1)*r1_i(3) + u(m, 2)*r2_i(3) + u(m, 3)*r3_i(3) + u(m, 4)*r4_i(3) + u(m, 5)*r5_i(3)
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r5_loc*(u(m, n + 2) - u(m, n - 2))
                        end do
                        f(m, nx - 2) = u(m, nx - 4)*r1_t(1) + u(m, nx - 3)*r2_t(1) + u(m, nx - 2)*r3_t(1) + u(m, nx - 1)*r4_t(1) + f(m, nx)*r5_t(1)
                        f(m, nx - 1) = u(m, nx - 3)*r1_t(2) + u(m, nx - 2)*r2_t(2) + u(m, nx - 1)*r3_t(2) + f(m, nx)*r4_t(2)
                    end do
                    !$omp end target teams distribute parallel do
                end if
            else
                !$omp target teams distribute parallel do &
                !$omp private(m,n) &
                !$omp shared(my,nx,rhs,f,u,bcs_b,r5_loc)
                do m = 1, my
                    f(m, 1) = u(m, 1)*r3_i(1) + u(m, 2)*r4_i(1) + u(m, 3)*r5_i(1) + u(m, 4)*r1_i(1)   ! r1(1) with extended stencil
                    f(m, 2) = u(m, 1)*r2_i(2) + u(m, 2)*r3_i(2) + u(m, 3)*r4_i(2) + u(m, 4)*r5_i(2)
                    f(m, 3) = u(m, 1)*r1_i(3) + u(m, 2)*r2_i(3) + u(m, 3)*r3_i(3) + u(m, 4)*r4_i(3) + u(m, 5)*r5_i(3)
                    ! Interior points
                    do n = 4, nx - 3
                        f(m, n) = u(m, n + 1) - u(m, n - 1) + r5_loc*(u(m, n + 2) - u(m, n - 2))
                    end do
                    ! top Boundary
                    f(m, nx - 2) = u(m, nx - 4)*r1_i(nx - 2) + u(m, nx - 3)*r2_i(nx - 2) + u(m, nx - 2)*r3_i(nx - 2) + u(m, nx - 1)*r4_i(nx - 2) + u(m, nx)*r5_i(nx - 2)
                    f(m, nx - 1) = u(m, nx - 3)*r1_i(nx - 1) + u(m, nx - 2)*r2_i(nx - 1) + u(m, nx - 1)*r3_i(nx - 1) + u(m, nx)*r4_i(nx - 1)
                    f(m, nx) = u(m, nx - 3)*r5_i(nx) + u(m, nx - 2)*r1_i(nx) + u(m, nx - 1)*r2_i(nx) + u(m, nx)*r3_i(nx)! r5(nx) with extended stencil
                end do
                !$omp end target teams distribute parallel do
            end if
        end if

        return
    end subroutine MatMul_5d_antisym_APU

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_5d_sym(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r3_loc     ! center diagonal
        real(wp) r5_loc     ! 2. upper-diagonal

        ! #######################################################################
        nx = size(rhs, 1)
        r5_loc = r5_i(3)      ! The first 2 equations, last 2 equations, are normalized differently
        r3_loc = r3_i(3)

        ! Boundary
        if (ibc == BCS_PERIODIC) then
            f(:, 1) = r3_loc*u(:, 1) + u(:, 2) + u(:, nx) &
                      + r5_loc*(u(:, 3) + u(:, nx - 1))

            f(:, 2) = r3_loc*u(:, 2) + u(:, 3) + u(:, 1) &
                      + r5_loc*(u(:, 4) + u(:, nx))

        else
            f(:, 1) = u(:, 1)*r3_i(1) + u(:, 2)*r4_i(1) + u(:, 3)*r5_i(1) &
                      + u(:, 4)*r1_i(1)   ! r1(1) contains 3. upper-diagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r2_i(2) + u(:, 2)*r3_i(2) + u(:, 3)*r4_i(2) + u(:, 4)*r5_i(2)

            if (any([BCS_ND, BCS_NN] == ibc)) f(:, 1) = 0.0_wp

        end if

        ! Interior points
        do n = 3, nx - 2
            f(:, n) = r3_loc*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                      + r5_loc*(u(:, n + 2) + u(:, n - 2))
        end do

        ! Boundary
        if (ibc == BCS_PERIODIC) then
            f(:, nx - 1) = r3_loc*u(:, nx - 1) + u(:, nx) + u(:, nx - 2) &
                           + r5_loc*(u(:, 1) + u(:, nx - 3))

            f(:, nx) = r3_loc*u(:, nx) + u(:, 1) + u(:, nx - 1) &
                       + r5_loc*(u(:, 2) + u(:, nx - 2))

        else
            f(:, nx - 1) = u(:, nx - 3)*r1_i(nx - 1) + u(:, nx - 2)*r2_i(nx - 1) + u(:, nx - 1)*r3_i(nx - 1) &
                           + u(:, nx)*r4_i(nx - 1)

            f(:, nx) = u(:, nx - 3)*r5_i(nx) & ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 2)*r1_i(nx) + u(:, nx - 1)*r2_i(nx) + u(:, nx)*r3_i(nx)

            if (any([BCS_DN, BCS_NN] == ibc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_5d_sym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_5d_sym_APU(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx, m, my
        real(wp) r3_loc     ! center diagonal
        real(wp) r5_loc     ! 2. upper-diagonal
        real(wp) sum_f, sum_u
        ! #######################################################################
        nx = size(rhs, 1)
        my = size(f, 1)
        r5_loc = r5_i(3)      ! The first 2 equations, last 2 equations, are normalized differently
        r3_loc = r3_i(3)
        sum_f = 0.0_wp
        sum_u = 0.0_wp
        ! Boundary
        !$omp target exit data map(delete: bcs_b, bcs_t)
        if (ibc == BCS_PERIODIC) then
            !$omp target teams distribute parallel do &
            !$omp private(m,n) &
            !$omp shared(my,nx,f,u,r5_loc,r3_loc)
            do m = 1, my
                f(m, 1) = r3_loc*u(m, 1) + u(m, 2) + u(m, nx) + r5_loc*(u(m, 3) + u(m, nx - 1))
                f(m, 2) = r3_loc*u(m, 2) + u(m, 3) + u(m, 1) + r5_loc*(u(m, 4) + u(m, nx))
                do n = 3, nx - 2
                    f(m, n) = r3_loc*u(m, n) + u(m, n + 1) + u(m, n - 1) + r5_loc*(u(m, n + 2) + u(m, n - 2))
                end do
                f(m, nx - 1) = r3_loc*u(m, nx - 1) + u(m, nx) + u(m, nx - 2) + r5_loc*(u(m, 1) + u(m, nx - 3))
                f(m, nx) = r3_loc*u(m, nx) + u(m, 1) + u(m, nx - 1) + r5_loc*(u(m, 2) + u(m, nx - 2))
            end do
            !$omp end target teams distribute parallel do
        else
            if (any([BCS_ND, BCS_NN] == ibc)) then
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r5_loc,r3_loc,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r3_i(1) + u(m, 2)*r4_i(1) + u(m, 3)*r5_i(1) + u(m, 4)*r1_i(1)   
                        ! r1(1) contains 3. upper-diagonal to allow for longer stencil at boundary
                        f(m, 2) = u(m, 1)*r2_i(2) + u(m, 2)*r3_i(2) + u(m, 3)*r4_i(2) + u(m, 4)*r5_i(2)
                        f(m, 1) = 0.0_wp
                        do n = 3, nx - 2
                            f(m, n) = r3_loc*u(m, n) + u(m, n + 1) + u(m, n - 1) + r5_loc*(u(m, n + 2) + u(m, n - 2))
                        end do
                        f(m, nx - 1) = u(m, nx - 3)*r1_i(nx - 1) + u(m, nx - 2)*r2_i(nx - 1) + u(m, nx - 1)*r3_i(nx - 1) + u(m, nx)*r4_i(nx - 1)
                        f(m, nx) = u(m, nx - 3)*r5_i(nx) + u(m, nx - 2)*r1_i(nx) + u(m, nx - 1)*r2_i(nx) + u(m, nx)*r3_i(nx)
                        ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
                        f(m, nx) = 0.0_wp
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r5_loc,r3_loc,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r3_i(1) + u(m, 2)*r4_i(1) + u(m, 3)*r5_i(1) + u(m, 4)*r1_i(1)   
                        ! r1(1) contains 3. upper-diagonal to allow for longer stencil at boundary
                        f(m, 2) = u(m, 1)*r2_i(2) + u(m, 2)*r3_i(2) + u(m, 3)*r4_i(2) + u(m, 4)*r5_i(2)
                        f(m, 1) = 0.0_wp
                        do n = 3, nx - 2
                            f(m, n) = r3_loc*u(m, n) + u(m, n + 1) + u(m, n - 1) + r5_loc*(u(m, n + 2) + u(m, n - 2))
                        end do
                        f(m, nx - 1) = u(m, nx - 3)*r1_i(nx - 1) + u(m, nx - 2)*r2_i(nx - 1) + u(m, nx - 1)*r3_i(nx - 1) + u(m, nx)*r4_i(nx - 1)
                        f(m, nx) = u(m, nx - 3)*r5_i(nx) + u(m, nx - 2)*r1_i(nx) + u(m, nx - 1)*r2_i(nx) + u(m, nx)*r3_i(nx)
                        ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
                    end do
                    !$omp end target teams distribute parallel do
                end if
            else
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r5_loc,r3_loc,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r3_i(1) + u(m, 2)*r4_i(1) + u(m, 3)*r5_i(1) + u(m, 4)*r1_i(1)   
                        ! r1(1) contains 3. upper-diagonal to allow for longer stencil at boundary
                        f(m, 2) = u(m, 1)*r2_i(2) + u(m, 2)*r3_i(2) + u(m, 3)*r4_i(2) + u(m, 4)*r5_i(2)
                        do n = 3, nx - 2
                            f(m, n) = r3_loc*u(m, n) + u(m, n + 1) + u(m, n - 1) + r5_loc*(u(m, n + 2) + u(m, n - 2))
                        end do
                        f(m, nx - 1) = u(m, nx - 3)*r1_i(nx - 1) + u(m, nx - 2)*r2_i(nx - 1) + u(m, nx - 1)*r3_i(nx - 1) + u(m, nx)*r4_i(nx - 1)
                        f(m, nx) = u(m, nx - 3)*r5_i(nx) + u(m, nx - 2)*r1_i(nx) + u(m, nx - 1)*r2_i(nx) + u(m, nx)*r3_i(nx)
                        ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
                        f(m, nx) = 0.0_wp
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r5_loc,r3_loc,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r3_i(1) + u(m, 2)*r4_i(1) + u(m, 3)*r5_i(1) + u(m, 4)*r1_i(1)   
                        ! r1(1) contains 3. upper-diagonal to allow for longer stencil at boundary
                        f(m, 2) = u(m, 1)*r2_i(2) + u(m, 2)*r3_i(2) + u(m, 3)*r4_i(2) + u(m, 4)*r5_i(2)
                        do n = 3, nx - 2
                            f(m, n) = r3_loc*u(m, n) + u(m, n + 1) + u(m, n - 1) + r5_loc*(u(m, n + 2) + u(m, n - 2))
                        end do
                        f(m, nx - 1) = u(m, nx - 3)*r1_i(nx - 1) + u(m, nx - 2)*r2_i(nx - 1) + u(m, nx - 1)*r3_i(nx - 1) + u(m, nx)*r4_i(nx - 1)
                        f(m, nx) = u(m, nx - 3)*r5_i(nx) + u(m, nx - 2)*r1_i(nx) + u(m, nx - 1)*r2_i(nx) + u(m, nx)*r3_i(nx)
                        ! r5(nx) contains 3. subdiagonal to allow for longer stencil at boundary
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end if
        end if

        return
    end subroutine MatMul_5d_sym_APU

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric heptadiagonal with 1. upper-diagonal equal to 1
    ! It also assumes equal coefficients in the 2. and 3. upper-diagonals for the interior points
    subroutine MatMul_7d_antisym(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r6_loc, r7_loc     ! 2. and 3. upper-diagonal

        ! #######################################################################
        nx = size(rhs, 1)
        r6_loc = r6_i(5)      ! The first 4 equations, last 4 equations, can be normalized differently
        r7_loc = r7_i(5)

        ! -------------------------------------------------------------------
        ! Boundary
        if (ibc == BCS_PERIODIC) then
            f(:, 1) = u(:, 2) - u(:, nx) + r6_loc*(u(:, 3) - u(:, nx - 1)) + r7_loc*(u(:, 4) - u(:, nx - 2))
            f(:, 2) = u(:, 3) - u(:, 1) + r6_loc*(u(:, 4) - u(:, nx)) + r7_loc*(u(:, 5) - u(:, nx - 1))
            f(:, 3) = u(:, 4) - u(:, 2) + r6_loc*(u(:, 5) - u(:, 1)) + r7_loc*(u(:, 6) - u(:, nx))
            f(:, 4) = u(:, 5) - u(:, 3) + r6_loc*(u(:, 6) - u(:, 2)) + r7_loc*(u(:, 7) - u(:, 1))

        else if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc)) then
            ! f(1) contains boundary condition
            if (present(bcs_b)) bcs_b(:) = f(:, 1)*r4_b(1) + u(:, 2)*r5_b(1) + u(:, 3)*r6_b(1) + u(:, 4)*r7_b(1) + u(:, 5)*r1_b(1)   ! r1(1) with extended stencil
            f(:, 2) = f(:, 1)*r3_b(2) + u(:, 2)*r4_b(2) + u(:, 3)*r5_b(2) + u(:, 4)*r6_b(2) + u(:, 5)*r7_b(2)
            f(:, 3) = f(:, 1)*r2_b(3) + u(:, 2)*r3_b(3) + u(:, 3)*r4_b(3) + u(:, 4)*r5_b(3) + u(:, 5)*r6_b(3) + u(:, 6)*r7_b(3)
            f(:, 4) = f(:, 1)*r1_b(4) + u(:, 2)*r2_b(4) + u(:, 3)*r3_b(4) + u(:, 4)*r4_b(4) + u(:, 5)*r5_b(4) + u(:, 6)*r6_b(4) + u(:, 7)*r7_b(4)

        else
            f(:, 1) = u(:, 1)*r4_i(1) + u(:, 2)*r5_i(1) + u(:, 3)*r6_i(1) + u(:, 4)*r7_i(1) + u(:, 5)*r1_i(1)   ! r1(1) with extended stencil
            f(:, 2) = u(:, 1)*r3_i(2) + u(:, 2)*r4_i(2) + u(:, 3)*r5_i(2) + u(:, 4)*r6_i(2) + u(:, 5)*r7_i(2)
            f(:, 3) = u(:, 1)*r2_i(3) + u(:, 2)*r3_i(3) + u(:, 3)*r4_i(3) + u(:, 4)*r5_i(3) + u(:, 5)*r6_i(3) + u(:, 6)*r7_i(3)
            f(:, 4) = u(:, 1)*r1_i(4) + u(:, 2)*r2_i(4) + u(:, 3)*r3_i(4) + u(:, 4)*r4_i(4) + u(:, 5)*r5_i(4) + u(:, 6)*r6_i(4) + u(:, 7)*r7_i(4)

        end if

        ! Interior points
        do n = 5, nx - 4
            f(:, n) = u(:, n + 1) - u(:, n - 1) + r6_loc*(u(:, n + 2) - u(:, n - 2)) + r7_loc*(u(:, n + 3) - u(:, n - 3))
        end do

        ! Boundary
        if (ibc == BCS_PERIODIC) then
            f(:, nx - 3) = u(:, nx - 2) - u(:, nx - 4) + r6_loc*(u(:, nx - 1) - u(:, nx - 5)) + r7_loc*(u(:, nx) - u(:, nx - 6))
            f(:, nx - 2) = u(:, nx - 1) - u(:, nx - 3) + r6_loc*(u(:, nx) - u(:, nx - 4)) + r7_loc*(u(:, 1) - u(:, nx - 5))
            f(:, nx - 1) = u(:, nx) - u(:, nx - 2) + r6_loc*(u(:, 1) - u(:, nx - 3)) + r7_loc*(u(:, 2) - u(:, nx - 4))
            f(:, nx) = u(:, 1) - u(:, nx - 1) + r6_loc*(u(:, 2) - u(:, nx - 2)) + r7_loc*(u(:, 3) - u(:, nx - 3))

        else if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc)) then
            ! f(nx) contains boundary condition
            f(:, nx - 3) = u(:, nx - 6)*r1_t(1) + u(:, nx - 5)*r2_t(1) + u(:, nx - 4)*r3_t(1) + u(:, nx - 3)*r4_t(1) + u(:, nx - 2)*r5_t(1) + u(:, nx - 1)*r6_t(1)+ f(:, nx)*r7_t(1)
            f(:, nx - 2) = u(:, nx - 5)*r1_t(2) + u(:, nx - 4)*r2_t(2) + u(:, nx - 3)*r3_t(2) + u(:, nx - 2)*r4_t(2) + u(:, nx - 1)*r5_t(2) + f(:, nx)*r6_t(2)
            f(:, nx - 1) = u(:, nx - 4)*r1_t(3) + u(:, nx - 3)*r2_t(3) + u(:, nx - 2)*r3_t(3) + u(:, nx - 1)*r4_t(3) + f(:, nx)*r5_t(3)
            if (present(bcs_b)) bcs_t(:) = u(:, nx - 4)*r7_t(4) + u(:, nx - 3)*r1_t(4) + u(:, nx - 2)*r2_t(4) + u(:, nx - 1)*r3_t(4) + f(:, nx)*r4_t(4) ! r7(nx) with extended stencil

        else
            f(:, nx - 3) = u(:, nx - 6)*r1_i(nx - 3) + u(:, nx - 5)*r2_i(nx - 3) + u(:, nx - 4)*r3_i(nx - 3) + u(:, nx - 3)*r4_i(nx - 3) + u(:, nx - 2)*r5_i(nx - 3) + u(:, nx - 1)*r6_i(nx - 3)+ u(:, nx)*r7_i(nx - 3)
            f(:, nx - 2) = u(:, nx - 5)*r1_i(nx - 2) + u(:, nx - 4)*r2_i(nx - 2) + u(:, nx - 3)*r3_i(nx - 2) + u(:, nx - 2)*r4_i(nx - 2) + u(:, nx - 1)*r5_i(nx - 2) + u(:, nx)*r6_i(nx - 2)
            f(:, nx - 1) = u(:, nx - 4)*r1_i(nx - 1) + u(:, nx - 3)*r2_i(nx - 1) + u(:, nx - 2)*r3_i(nx - 1) + u(:, nx - 1)*r4_i(nx - 1) + u(:, nx)*r5_i(nx - 1)
            f(:, nx) = u(:, nx - 4)*r7_i(nx) + u(:, nx - 3)*r1_i(nx) + u(:, nx - 2)*r2_i(nx) + u(:, nx - 1)*r3_i(nx) + u(:, nx)*r4_i(nx) ! r7(nx) with extended stencil
        end if

        return
    end subroutine MatMul_7d_antisym

    ! #######################################################################
    ! #######################################################################
    ! Calculate f = B u, assuming B is antisymmetric heptadiagonal with 1. upper-diagonal equal to 1
    ! It also assumes equal coefficients in the 2. and 3. upper-diagonals for the interior points
    subroutine MatMul_7d_antisym_APU(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx, m, my
        real(wp) r6_loc, r7_loc     ! 2. and 3. upper-diagonal
        real(wp) sum_f, sum_u
        ! #######################################################################
        nx = size(rhs, 1)
        my = size(f, 1)
        r6_loc = r6_i(5)      ! The first 4 equations, last 4 equations, can be normalized differently
        r7_loc = r7_i(5)
        ! -------------------------------------------------------------------
        ! Boundary
        !$omp target exit data map(delete: bcs_b, bcs_t)
        if (ibc == BCS_PERIODIC) then
            !$omp target teams distribute parallel do & 
            !$omp private(m,n) &
            !$omp shared(my,nx,f,u,r6_loc,r7_loc)
            do m = 1, my
                f(m, 1) = u(m, 2) - u(m, nx) + r6_loc*(u(m, 3) - u(m, nx - 1)) + r7_loc*(u(m, 4) - u(m, nx - 2))
                f(m, 2) = u(m, 3) - u(m, 1) + r6_loc*(u(m, 4) - u(m, nx)) + r7_loc*(u(m, 5) - u(m, nx - 1))
                f(m, 3) = u(m, 4) - u(m, 2) + r6_loc*(u(m, 5) - u(m, 1)) + r7_loc*(u(m, 6) - u(m, nx))
                f(m, 4) = u(m, 5) - u(m, 3) + r6_loc*(u(m, 6) - u(m, 2)) + r7_loc*(u(m, 7) - u(m, 1))
                ! Interior points
                do n = 5, nx - 4
                    f(m, n) = u(m, n + 1) - u(m, n - 1) + r6_loc*(u(m, n + 2) - u(m, n - 2)) + r7_loc*(u(m, n + 3) - u(m, n - 3))
                end do
                ! Boundary
                f(m, nx - 3) = u(m, nx - 2) - u(m, nx - 4) + r6_loc*(u(m, nx - 1) - u(m, nx - 5)) + r7_loc*(u(m, nx) - u(m, nx - 6))
                f(m, nx - 2) = u(m, nx - 1) - u(m, nx - 3) + r6_loc*(u(m, nx) - u(m, nx - 4)) + r7_loc*(u(m, 1) - u(m, nx - 5))
                f(m, nx - 1) = u(m, nx) - u(m, nx - 2) + r6_loc*(u(m, 1) - u(m, nx - 3)) + r7_loc*(u(m, 2) - u(m, nx - 4))
                f(m, nx) = u(m, 1) - u(m, nx - 1) + r6_loc*(u(m, 2) - u(m, nx - 2)) + r7_loc*(u(m, 3) - u(m, nx - 3))
            end do
            !$omp end target teams distribute parallel do
        else if (any([BCS_ND, BCS_NN, BCS_MIN, BCS_BOTH] == ibc)) then
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc)) then
                ! f(1) contains boundary condition
                if (present(bcs_b)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r6_loc,r7_loc,rhs_b,rhs_t,bcs_b,bcs_t)
                    do m = 1, my
                        bcs_b(m) = f(m, 1)*r4_b(1) + u(m, 2)*r5_b(1) + u(m, 3)*r6_b(1) + u(m, 4)*r7_b(1) + u(m, 5)*r1_b(1)   ! r1(1) with extended stencil
                        f(m, 2) = f(m, 1)*r3_b(2) + u(m, 2)*r4_b(2) + u(m, 3)*r5_b(2) + u(m, 4)*r6_b(2) + u(m, 5)*r7_b(2)
                        f(m, 3) = f(m, 1)*r2_b(3) + u(m, 2)*r3_b(3) + u(m, 3)*r4_b(3) + u(m, 4)*r5_b(3) + u(m, 5)*r6_b(3) + u(m, 6)*r7_b(3)
                        f(m, 4) = f(m, 1)*r1_b(4) + u(m, 2)*r2_b(4) + u(m, 3)*r3_b(4) + u(m, 4)*r4_b(4) + u(m, 5)*r5_b(4) + u(m, 6)*r6_b(4) + u(m, 7)*r7_b(4)
                        ! Interior points
                        do n = 5, nx - 4
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r6_loc*(u(m, n + 2) - u(m, n - 2)) + r7_loc*(u(m, n + 3) - u(m, n - 3))
                        end do
                        ! f(nx) contains boundary condition
                        f(m, nx - 3) = u(m, nx - 6)*r1_t(1) + u(m, nx - 5)*r2_t(1) + u(m, nx - 4)*r3_t(1) + u(m, nx - 3)*r4_t(1) + u(m, nx - 2)*r5_t(1) + u(m, nx - 1)*r6_t(1)+ f(m, nx)*r7_t(1)
                        f(m, nx - 2) = u(m, nx - 5)*r1_t(2) + u(m, nx - 4)*r2_t(2) + u(m, nx - 3)*r3_t(2) + u(m, nx - 2)*r4_t(2) + u(m, nx - 1)*r5_t(2) + f(m, nx)*r6_t(2)
                        f(m, nx - 1) = u(m, nx - 4)*r1_t(3) + u(m, nx - 3)*r2_t(3) + u(m, nx - 2)*r3_t(3) + u(m, nx - 1)*r4_t(3) + f(m, nx)*r5_t(3)
                        bcs_t(m) = u(m, nx - 4)*r7_t(4) + u(m, nx - 3)*r1_t(4) + u(m, nx - 2)*r2_t(4) + u(m, nx - 1)*r3_t(4) + f(m, nx)*r4_t(4) ! r7(nx) with extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r6_loc,r7_loc,rhs_b,rhs_t)
                    do m = 1, my
                        f(m, 2) = f(m, 1)*r3_b(2) + u(m, 2)*r4_b(2) + u(m, 3)*r5_b(2) + u(m, 4)*r6_b(2) + u(m, 5)*r7_b(2)
                        f(m, 3) = f(m, 1)*r2_b(3) + u(m, 2)*r3_b(3) + u(m, 3)*r4_b(3) + u(m, 4)*r5_b(3) + u(m, 5)*r6_b(3) + u(m, 6)*r7_b(3)
                        f(m, 4) = f(m, 1)*r1_b(4) + u(m, 2)*r2_b(4) + u(m, 3)*r3_b(4) + u(m, 4)*r4_b(4) + u(m, 5)*r5_b(4) + u(m, 6)*r6_b(4) + u(m, 7)*r7_b(4)
                        ! Interior points
                        do n = 5, nx - 4
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r6_loc*(u(m, n + 2) - u(m, n - 2)) + r7_loc*(u(m, n + 3) - u(m, n - 3))
                        end do
                        ! f(nx) contains boundary condition
                        f(m, nx - 3) = u(m, nx - 6)*r1_t(1) + u(m, nx - 5)*r2_t(1) + u(m, nx - 4)*r3_t(1) + u(m, nx - 3)*r4_t(1) + u(m, nx - 2)*r5_t(1) + u(m, nx - 1)*r6_t(1)+ f(m, nx)*r7_t(1)
                        f(m, nx - 2) = u(m, nx - 5)*r1_t(2) + u(m, nx - 4)*r2_t(2) + u(m, nx - 3)*r3_t(2) + u(m, nx - 2)*r4_t(2) + u(m, nx - 1)*r5_t(2) + f(m, nx)*r6_t(2)
                        f(m, nx - 1) = u(m, nx - 4)*r1_t(3) + u(m, nx - 3)*r2_t(3) + u(m, nx - 2)*r3_t(3) + u(m, nx - 1)*r4_t(3) + f(m, nx)*r5_t(3)
                    end do
                    !$omp end target teams distribute parallel do
                end if
            else
                if (present(bcs_b)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r6_loc,r7_loc,rhs_b,rhs,bcs_b)
                    do m = 1, my
                        bcs_b(m) = f(m, 1)*r4_b(1) + u(m, 2)*r5_b(1) + u(m, 3)*r6_b(1) + u(m, 4)*r7_b(1) + u(m, 5)*r1_b(1)   ! r1(1) with extended stencil
                        f(m, 2) = f(m, 1)*r3_b(2) + u(m, 2)*r4_b(2) + u(m, 3)*r5_b(2) + u(m, 4)*r6_b(2) + u(m, 5)*r7_b(2)
                        f(m, 3) = f(m, 1)*r2_b(3) + u(m, 2)*r3_b(3) + u(m, 3)*r4_b(3) + u(m, 4)*r5_b(3) + u(m, 5)*r6_b(3) + u(m, 6)*r7_b(3)
                        f(m, 4) = f(m, 1)*r1_b(4) + u(m, 2)*r2_b(4) + u(m, 3)*r3_b(4) + u(m, 4)*r4_b(4) + u(m, 5)*r5_b(4) + u(m, 6)*r6_b(4) + u(m, 7)*r7_b(4)
                        ! Interior points
                        do n = 5, nx - 4
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r6_loc*(u(m, n + 2) - u(m, n - 2)) + r7_loc*(u(m, n + 3) - u(m, n - 3))
                        end do
                        f(m, nx - 3) = u(m, nx - 6)*r1_i(nx - 3) + u(m, nx - 5)*r2_i(nx - 3) + u(m, nx - 4)*r3_i(nx - 3) + u(m, nx - 3)*r4_i(nx - 3) + u(m, nx - 2)*r5_i(nx - 3) + u(m, nx - 1)*r6_i(nx - 3)+ u(m, nx)*r7_i(nx - 3)
                        f(m, nx - 2) = u(m, nx - 5)*r1_i(nx - 2) + u(m, nx - 4)*r2_i(nx - 2) + u(m, nx - 3)*r3_i(nx - 2) + u(m, nx - 2)*r4_i(nx - 2) + u(m, nx - 1)*r5_i(nx - 2) + u(m, nx)*r6_i(nx - 2)
                        f(m, nx - 1) = u(m, nx - 4)*r1_i(nx - 1) + u(m, nx - 3)*r2_i(nx - 1) + u(m, nx - 2)*r3_i(nx - 1) + u(m, nx - 1)*r4_i(nx - 1) + u(m, nx)*r5_i(nx - 1)
                        f(m, nx) = u(m, nx - 4)*r7_i(nx) + u(m, nx - 3)*r1_i(nx) + u(m, nx - 2)*r2_i(nx) + u(m, nx - 1)*r3_i(nx) + u(m, nx)*r4_i(nx) ! r7(nx) with extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r6_loc,r7_loc,rhs_b,rhs)
                    do m = 1, my
                        f(m, 2) = f(m, 1)*r3_b(2) + u(m, 2)*r4_b(2) + u(m, 3)*r5_b(2) + u(m, 4)*r6_b(2) + u(m, 5)*r7_b(2)
                        f(m, 3) = f(m, 1)*r2_b(3) + u(m, 2)*r3_b(3) + u(m, 3)*r4_b(3) + u(m, 4)*r5_b(3) + u(m, 5)*r6_b(3) + u(m, 6)*r7_b(3)
                        f(m, 4) = f(m, 1)*r1_b(4) + u(m, 2)*r2_b(4) + u(m, 3)*r3_b(4) + u(m, 4)*r4_b(4) + u(m, 5)*r5_b(4) + u(m, 6)*r6_b(4) + u(m, 7)*r7_b(4)
                        ! Interior points
                        do n = 5, nx - 4
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r6_loc*(u(m, n + 2) - u(m, n - 2)) + r7_loc*(u(m, n + 3) - u(m, n - 3))
                        end do
                        f(m, nx - 3) = u(m, nx - 6)*r1_i(nx - 3) + u(m, nx - 5)*r2_i(nx - 3) + u(m, nx - 4)*r3_i(nx - 3) + u(m, nx - 3)*r4_i(nx - 3) + u(m, nx - 2)*r5_i(nx - 3) + u(m, nx - 1)*r6_i(nx - 3)+ u(m, nx)*r7_i(nx - 3)
                        f(m, nx - 2) = u(m, nx - 5)*r1_i(nx - 2) + u(m, nx - 4)*r2_i(nx - 2) + u(m, nx - 3)*r3_i(nx - 2) + u(m, nx - 2)*r4_i(nx - 2) + u(m, nx - 1)*r5_i(nx - 2) + u(m, nx)*r6_i(nx - 2)
                        f(m, nx - 1) = u(m, nx - 4)*r1_i(nx - 1) + u(m, nx - 3)*r2_i(nx - 1) + u(m, nx - 2)*r3_i(nx - 1) + u(m, nx - 1)*r4_i(nx - 1) + u(m, nx)*r5_i(nx - 1)
                        f(m, nx) = u(m, nx - 4)*r7_i(nx) + u(m, nx - 3)*r1_i(nx) + u(m, nx - 2)*r2_i(nx) + u(m, nx - 1)*r3_i(nx) + u(m, nx)*r4_i(nx) ! r7(nx) with extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end if
        else
            if (any([BCS_DN, BCS_NN, BCS_MAX, BCS_BOTH] == ibc)) then
                if (present(bcs_t)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r6_loc,r7_loc,rhs,rhs_t,bcs_t)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r4_i(1) + u(m, 2)*r5_i(1) + u(m, 3)*r6_i(1) + u(m, 4)*r7_i(1) + u(m, 5)*r1_i(1)   ! r1(1) with extended stencil
                        f(m, 2) = u(m, 1)*r3_i(2) + u(m, 2)*r4_i(2) + u(m, 3)*r5_i(2) + u(m, 4)*r6_i(2) + u(m, 5)*r7_i(2)
                        f(m, 3) = u(m, 1)*r2_i(3) + u(m, 2)*r3_i(3) + u(m, 3)*r4_i(3) + u(m, 4)*r5_i(3) + u(m, 5)*r6_i(3) + u(m, 6)*r7_i(3)
                        f(m, 4) = u(m, 1)*r1_i(4) + u(m, 2)*r2_i(4) + u(m, 3)*r3_i(4) + u(m, 4)*r4_i(4) + u(m, 5)*r5_i(4) + u(m, 6)*r6_i(4) + u(m, 7)*r7_i(4)
                        ! Interior points
                        do n = 5, nx - 4
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r6_loc*(u(m, n + 2) - u(m, n - 2)) + r7_loc*(u(m, n + 3) - u(m, n - 3))
                        end do
                        ! f(nx) contains boundary condition
                        f(m, nx - 3) = u(m, nx - 6)*r1_t(1) + u(m, nx - 5)*r2_t(1) + u(m, nx - 4)*r3_t(1) + u(m, nx - 3)*r4_t(1) + u(m, nx - 2)*r5_t(1) + u(m, nx - 1)*r6_t(1)+ f(m, nx)*r7_t(1)
                        f(m, nx - 2) = u(m, nx - 5)*r1_t(2) + u(m, nx - 4)*r2_t(2) + u(m, nx - 3)*r3_t(2) + u(m, nx - 2)*r4_t(2) + u(m, nx - 1)*r5_t(2) + f(m, nx)*r6_t(2)
                        f(m, nx - 1) = u(m, nx - 4)*r1_t(3) + u(m, nx - 3)*r2_t(3) + u(m, nx - 2)*r3_t(3) + u(m, nx - 1)*r4_t(3) + f(m, nx)*r5_t(3)
                        bcs_t(m) = u(m, nx - 4)*r7_t(4) + u(m, nx - 3)*r1_t(4) + u(m, nx - 2)*r2_t(4) + u(m, nx - 1)*r3_t(4) + f(m, nx)*r4_t(4) ! r7(nx) with extended stencil
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r6_loc,r7_loc,rhs,rhs_t)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r4_i(1) + u(m, 2)*r5_i(1) + u(m, 3)*r6_i(1) + u(m, 4)*r7_i(1) + u(m, 5)*r1_i(1)   ! r1(1) with extended stencil
                        f(m, 2) = u(m, 1)*r3_i(2) + u(m, 2)*r4_i(2) + u(m, 3)*r5_i(2) + u(m, 4)*r6_i(2) + u(m, 5)*r7_i(2)
                        f(m, 3) = u(m, 1)*r2_i(3) + u(m, 2)*r3_i(3) + u(m, 3)*r4_i(3) + u(m, 4)*r5_i(3) + u(m, 5)*r6_i(3) + u(m, 6)*r7_i(3)
                        f(m, 4) = u(m, 1)*r1_i(4) + u(m, 2)*r2_i(4) + u(m, 3)*r3_i(4) + u(m, 4)*r4_i(4) + u(m, 5)*r5_i(4) + u(m, 6)*r6_i(4) + u(m, 7)*r7_i(4)
                        ! Interior points
                        do n = 5, nx - 4
                            f(m, n) = u(m, n + 1) - u(m, n - 1) + r6_loc*(u(m, n + 2) - u(m, n - 2)) + r7_loc*(u(m, n + 3) - u(m, n - 3))
                        end do
                        ! f(nx) contains boundary condition
                        f(m, nx - 3) = u(m, nx - 6)*r1_t(1) + u(m, nx - 5)*r2_t(1) + u(m, nx - 4)*r3_t(1) + u(m, nx - 3)*r4_t(1) + u(m, nx - 2)*r5_t(1) + u(m, nx - 1)*r6_t(1)+ f(m, nx)*r7_t(1)
                        f(m, nx - 2) = u(m, nx - 5)*r1_t(2) + u(m, nx - 4)*r2_t(2) + u(m, nx - 3)*r3_t(2) + u(m, nx - 2)*r4_t(2) + u(m, nx - 1)*r5_t(2) + f(m, nx)*r6_t(2)
                        f(m, nx - 1) = u(m, nx - 4)*r1_t(3) + u(m, nx - 3)*r2_t(3) + u(m, nx - 2)*r3_t(3) + u(m, nx - 1)*r4_t(3) + f(m, nx)*r5_t(3)
                    end do
                    !$omp end target teams distribute parallel do
                end if
            else
                !$omp target teams distribute parallel do &
                !$omp private(m,n) &
                !$omp shared(my,nx,f,u,r6_loc,r7_loc,rhs)
                do m = 1, my
                    f(m, 1) = u(m, 1)*r4_i(1) + u(m, 2)*r5_i(1) + u(m, 3)*r6_i(1) + u(m, 4)*r7_i(1) + u(m, 5)*r1_i(1)   ! r1(1) with extended stencil
                    f(m, 2) = u(m, 1)*r3_i(2) + u(m, 2)*r4_i(2) + u(m, 3)*r5_i(2) + u(m, 4)*r6_i(2) + u(m, 5)*r7_i(2)
                    f(m, 3) = u(m, 1)*r2_i(3) + u(m, 2)*r3_i(3) + u(m, 3)*r4_i(3) + u(m, 4)*r5_i(3) + u(m, 5)*r6_i(3) + u(m, 6)*r7_i(3)
                    f(m, 4) = u(m, 1)*r1_i(4) + u(m, 2)*r2_i(4) + u(m, 3)*r3_i(4) + u(m, 4)*r4_i(4) + u(m, 5)*r5_i(4) + u(m, 6)*r6_i(4) + u(m, 7)*r7_i(4)
                    ! Interior points
                    do n = 5, nx - 4
                        f(m, n) = u(m, n + 1) - u(m, n - 1) + r6_loc*(u(m, n + 2) - u(m, n - 2)) + r7_loc*(u(m, n + 3) - u(m, n - 3))
                    end do
                    f(m, nx - 3) = u(m, nx - 6)*r1_i(nx - 3) + u(m, nx - 5)*r2_i(nx - 3) + u(m, nx - 4)*r3_i(nx - 3) + u(m, nx - 3)*r4_i(nx - 3) + u(m, nx - 2)*r5_i(nx - 3) + u(m, nx - 1)*r6_i(nx - 3)+ u(m, nx)*r7_i(nx - 3)
                    f(m, nx - 2) = u(m, nx - 5)*r1_i(nx - 2) + u(m, nx - 4)*r2_i(nx - 2) + u(m, nx - 3)*r3_i(nx - 2) + u(m, nx - 2)*r4_i(nx - 2) + u(m, nx - 1)*r5_i(nx - 2) + u(m, nx)*r6_i(nx - 2)
                    f(m, nx - 1) = u(m, nx - 4)*r1_i(nx - 1) + u(m, nx - 3)*r2_i(nx - 1) + u(m, nx - 2)*r3_i(nx - 1) + u(m, nx - 1)*r4_i(nx - 1) + u(m, nx)*r5_i(nx - 1)
                    f(m, nx) = u(m, nx - 4)*r7_i(nx) + u(m, nx - 3)*r1_i(nx) + u(m, nx - 2)*r2_i(nx) + u(m, nx - 1)*r3_i(nx) + u(m, nx)*r4_i(nx) ! r7(nx) with extended stencil
                end do
                !$omp end target teams distribute parallel do
            end if
        end if
        return
    end subroutine MatMul_7d_antisym_APU

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_7d_sym(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx
        real(wp) r4_loc     ! center diagonal
        real(wp) r6_loc     ! 2. upper-diagonal
        real(wp) r7_loc     ! 3. upper-diagonal

        ! #######################################################################
        nx = size(rhs, 1)
        r7_loc = r7_i(4)
        r6_loc = r6_i(4)      ! The first 3 equations, last 3 equations, are normalized differently
        r4_loc = r4_i(4)

        ! Boundary
        if (ibc == BCS_PERIODIC) then
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
            f(:, 1) = u(:, 1)*r4_i(1) + u(:, 2)*r5_i(1) + u(:, 3)*r6_i(1) + u(:, 4)*r7_i(1) &
                      + u(:, 5)*r1_i(1)   ! r1(1) contains 4. upper-diagonal to allow for longer stencil at boundary

            f(:, 2) = u(:, 1)*r3_i(2) + u(:, 2)*r4_i(2) + u(:, 3)*r5_i(2) + u(:, 4)*r6_i(2) + u(:, 5)*r7_i(2)

            f(:, 3) = u(:, 1)*r2_i(3) + u(:, 2)*r3_i(3) + u(:, 3)*r4_i(3) + u(:, 4)*r5_i(3) + u(:, 5)*r6_i(3) + u(:, 6)*r7_i(3)

            if (any([BCS_ND, BCS_NN] == ibc)) f(:, 1) = 0.0_wp

        end if

        ! Interior points
        do n = 4, nx - 3
            f(:, n) = r4_loc*u(:, n) + u(:, n + 1) + u(:, n - 1) &
                      + r6_loc*(u(:, n + 2) + u(:, n - 2)) &
                      + r7_loc*(u(:, n + 3) + u(:, n - 3))
        end do

        ! Boundary
        if (ibc == BCS_PERIODIC) then
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
      f(:, nx - 2) = u(:, nx - 5)*r1_i(nx - 2) + u(:, nx - 4)*r2_i(nx - 2) + u(:, nx - 3)*r3_i(nx - 2) + u(:, nx - 2)*r4_i(nx - 2) + u(:, nx - 1)*r5_i(nx - 2) &
                           + u(:, nx)*r6_i(nx - 2)

            f(:, nx - 1) = u(:, nx - 4)*r1_i(nx - 1) + u(:, nx - 3)*r2_i(nx - 1) + u(:, nx - 2)*r3_i(nx - 1) + u(:, nx - 1)*r4_i(nx - 1) &
                           + u(:, nx)*r5_i(nx - 1)

            f(:, nx) = u(:, nx - 4)*r7_i(nx) & ! r7(nx) contains 4. subdiagonal to allow for longer stencil at boundary
                       + u(:, nx - 3)*r1_i(nx) + u(:, nx - 2)*r2_i(nx) + u(:, nx - 1)*r3_i(nx) + u(:, nx)*r4_i(nx)

            if (any([BCS_DN, BCS_NN] == ibc)) f(:, nx) = 0.0_wp

        end if

        return
    end subroutine MatMul_7d_sym

    ! #######################################################################
    ! #######################################################################
    subroutine MatMul_7d_sym_APU(rhs, u, f, ibc, rhs_b, rhs_t, bcs_b, bcs_t)
        real(wp), intent(in) :: rhs(:, :)                               ! diagonals of B
        real(wp), intent(in) :: u(:, :)                                 ! vector u
        real(wp), intent(out) :: f(:, :)                                ! vector f = B u
        integer, intent(in) :: ibc
        real(wp), intent(in), optional :: rhs_b(:, 0:), rhs_t(0:, :)    ! Special bcs at bottom, top
        real(wp), intent(out), optional :: bcs_b(:), bcs_t(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nx, m, my
        real(wp) r4_loc     ! center diagonal
        real(wp) r6_loc     ! 2. upper-diagonal
        real(wp) r7_loc     ! 3. upper-diagonal
        ! #######################################################################
        nx = size(rhs, 1)
        my = size(f, 1)
        r7_loc = r7_i(4)
        r6_loc = r6_i(4)      ! The first 3 equations, last 3 equations, are normalized differently
        r4_loc = r4_i(4)

        ! Boundary
        if (ibc == BCS_PERIODIC) then
            !$omp target teams distribute parallel do &
            !$omp private(m,n) &
            !$omp shared(my,nx,f,u,r4_loc,r6_loc,r7_loc)
            do m = 1, my
                f(m, 1) = r4_loc*u(m, 1) + u(m, 2) + u(m, nx) + r6_loc*(u(m, 3) + u(m, nx - 1)) + r7_loc*(u(m, 4) + u(m, nx - 2))
                f(m, 2) = r4_loc*u(m, 2) + u(m, 3) + u(m, 1) + r6_loc*(u(m, 4) + u(m, nx)) + r7_loc*(u(m, 5) + u(m, nx - 1))
                f(m, 3) = r4_loc*u(m, 3) + u(m, 4) + u(m, 2) + r6_loc*(u(m, 5) + u(m, 1)) + r7_loc*(u(m, 6) + u(m, nx))
                do n = 4, nx - 3
                    f(m, n) = r4_loc*u(m, n) + u(m, n + 1) + u(m, n - 1) + r6_loc*(u(m, n + 2) + u(m, n - 2)) + r7_loc*(u(m, n + 3) + u(m, n - 3))
                end do
                f(m, nx - 2) = r4_loc*u(m, nx - 2) + u(m, nx - 1) + u(m, nx - 3) + r6_loc*(u(m, nx) + u(m, nx - 4)) + r7_loc*(u(m, 1) + u(m, nx - 5))
                f(m, nx - 1) = r4_loc*u(m, nx - 1) + u(m, nx) + u(m, nx - 2) + r6_loc*(u(m, 1) + u(m, nx - 3)) + r7_loc*(u(m, 2) + u(m, nx - 4))
                f(m, nx) = r4_loc*u(m, nx) + u(m, 1) + u(m, nx - 1) + r6_loc*(u(m, 2) + u(m, nx - 2)) + r7_loc*(u(m, 3) + u(m, nx - 3))
            end do
            !$omp end target teams distribute parallel do
        else
            if (any([BCS_ND, BCS_NN] == ibc)) then
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r4_loc,r6_loc,r7_loc,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r4_i(1) + u(m, 2)*r5_i(1) + u(m, 3)*r6_i(1) + u(m, 4)*r7_i(1) + u(m, 5)*r1_i(1)   
                        ! r1(1) contains 4. upper-diagonal to allow for longer stencil at boundary
                        f(m, 2) = u(m, 1)*r3_i(2) + u(m, 2)*r4_i(2) + u(m, 3)*r5_i(2) + u(m, 4)*r6_i(2) + u(m, 5)*r7_i(2)
                        f(m, 3) = u(m, 1)*r2_i(3) + u(m, 2)*r3_i(3) + u(m, 3)*r4_i(3) + u(m, 4)*r5_i(3) + u(m, 5)*r6_i(3) + u(m, 6)*r7_i(3)
                        f(m, 1) = 0.0_wp
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = r4_loc*u(m, n) + u(m, n + 1) + u(m, n - 1) + r6_loc*(u(m, n + 2) + u(m, n - 2)) + r7_loc*(u(m, n + 3) + u(m, n - 3))
                        end do
                        f(m, nx - 2) = u(m, nx - 5)*r1_i(nx - 2) + u(m, nx - 4)*r2_i(nx - 2) + u(m, nx - 3)*r3_i(nx - 2) + u(m, nx - 2)*r4_i(nx - 2) + u(m, nx - 1)*r5_i(nx - 2) + u(m, nx)*r6_i(nx - 2)
                        f(m, nx - 1) = u(m, nx - 4)*r1_i(nx - 1) + u(m, nx - 3)*r2_i(nx - 1) + u(m, nx - 2)*r3_i(nx - 1) + u(m, nx - 1)*r4_i(nx - 1) + u(m, nx)*r5_i(nx - 1)
                        f(m, nx) = u(m, nx - 4)*r7_i(nx) + u(m, nx - 3)*r1_i(nx) + u(m, nx - 2)*r2_i(nx) + u(m, nx - 1)*r3_i(nx) + u(m, nx)*r4_i(nx)
                        ! r7(nx) contains 4. subdiagonal to allow for longer stencil at boundary
                        f(m, nx) = 0.0_wp   
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r4_loc,r6_loc,r7_loc,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r4_i(1) + u(m, 2)*r5_i(1) + u(m, 3)*r6_i(1) + u(m, 4)*r7_i(1) + u(m, 5)*r1_i(1)   
                        ! r1(1) contains 4. upper-diagonal to allow for longer stencil at boundary
                        f(m, 2) = u(m, 1)*r3_i(2) + u(m, 2)*r4_i(2) + u(m, 3)*r5_i(2) + u(m, 4)*r6_i(2) + u(m, 5)*r7_i(2)
                        f(m, 3) = u(m, 1)*r2_i(3) + u(m, 2)*r3_i(3) + u(m, 3)*r4_i(3) + u(m, 4)*r5_i(3) + u(m, 5)*r6_i(3) + u(m, 6)*r7_i(3)
                        if (any([BCS_ND, BCS_NN] == ibc)) f(m, 1) = 0.0_wp
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = r4_loc*u(m, n) + u(m, n + 1) + u(m, n - 1) + r6_loc*(u(m, n + 2) + u(m, n - 2)) + r7_loc*(u(m, n + 3) + u(m, n - 3))
                        end do
                        f(m, nx - 2) = u(m, nx - 5)*r1_i(nx - 2) + u(m, nx - 4)*r2_i(nx - 2) + u(m, nx - 3)*r3_i(nx - 2) + u(m, nx - 2)*r4_i(nx - 2) + u(m, nx - 1)*r5_i(nx - 2) + u(m, nx)*r6_i(nx - 2)
                        f(m, nx - 1) = u(m, nx - 4)*r1_i(nx - 1) + u(m, nx - 3)*r2_i(nx - 1) + u(m, nx - 2)*r3_i(nx - 1) + u(m, nx - 1)*r4_i(nx - 1) + u(m, nx)*r5_i(nx - 1)
                        f(m, nx) = u(m, nx - 4)*r7_i(nx) + u(m, nx - 3)*r1_i(nx) + u(m, nx - 2)*r2_i(nx) + u(m, nx - 1)*r3_i(nx) + u(m, nx)*r4_i(nx)
                        ! r7(nx) contains 4. subdiagonal to allow for longer stencil at boundary
                    end do
                    !$omp end target teams distribute parallel do
                end if
            else
                if (any([BCS_DN, BCS_NN] == ibc)) then
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r4_loc,r6_loc,r7_loc,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r4_i(1) + u(m, 2)*r5_i(1) + u(m, 3)*r6_i(1) + u(m, 4)*r7_i(1) + u(m, 5)*r1_i(1)   
                        ! r1(1) contains 4. upper-diagonal to allow for longer stencil at boundary
                        f(m, 2) = u(m, 1)*r3_i(2) + u(m, 2)*r4_i(2) + u(m, 3)*r5_i(2) + u(m, 4)*r6_i(2) + u(m, 5)*r7_i(2)
                        f(m, 3) = u(m, 1)*r2_i(3) + u(m, 2)*r3_i(3) + u(m, 3)*r4_i(3) + u(m, 4)*r5_i(3) + u(m, 5)*r6_i(3) + u(m, 6)*r7_i(3)
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = r4_loc*u(m, n) + u(m, n + 1) + u(m, n - 1) + r6_loc*(u(m, n + 2) + u(m, n - 2)) + r7_loc*(u(m, n + 3) + u(m, n - 3))
                        end do
                        f(m, nx - 2) = u(m, nx - 5)*r1_i(nx - 2) + u(m, nx - 4)*r2_i(nx - 2) + u(m, nx - 3)*r3_i(nx - 2) + u(m, nx - 2)*r4_i(nx - 2) + u(m, nx - 1)*r5_i(nx - 2) + u(m, nx)*r6_i(nx - 2)
                        f(m, nx - 1) = u(m, nx - 4)*r1_i(nx - 1) + u(m, nx - 3)*r2_i(nx - 1) + u(m, nx - 2)*r3_i(nx - 1) + u(m, nx - 1)*r4_i(nx - 1) + u(m, nx)*r5_i(nx - 1)
                        f(m, nx) = u(m, nx - 4)*r7_i(nx) + u(m, nx - 3)*r1_i(nx) + u(m, nx - 2)*r2_i(nx) + u(m, nx - 1)*r3_i(nx) + u(m, nx)*r4_i(nx)
                        ! r7(nx) contains 4. subdiagonal to allow for longer stencil at boundary
                        f(m, nx) = 0.0_wp
                    end do
                    !$omp end target teams distribute parallel do
                else
                    !$omp target teams distribute parallel do &
                    !$omp private(m,n) &
                    !$omp shared(my,nx,f,u,r4_loc,r6_loc,r7_loc,rhs)
                    do m = 1, my
                        f(m, 1) = u(m, 1)*r4_i(1) + u(m, 2)*r5_i(1) + u(m, 3)*r6_i(1) + u(m, 4)*r7_i(1) + u(m, 5)*r1_i(1)   
                        ! r1(1) contains 4. upper-diagonal to allow for longer stencil at boundary
                        f(m, 2) = u(m, 1)*r3_i(2) + u(m, 2)*r4_i(2) + u(m, 3)*r5_i(2) + u(m, 4)*r6_i(2) + u(m, 5)*r7_i(2)
                        f(m, 3) = u(m, 1)*r2_i(3) + u(m, 2)*r3_i(3) + u(m, 3)*r4_i(3) + u(m, 4)*r5_i(3) + u(m, 5)*r6_i(3) + u(m, 6)*r7_i(3)
                        ! Interior points
                        do n = 4, nx - 3
                            f(m, n) = r4_loc*u(m, n) + u(m, n + 1) + u(m, n - 1) + r6_loc*(u(m, n + 2) + u(m, n - 2)) + r7_loc*(u(m, n + 3) + u(m, n - 3))
                        end do
                        f(m, nx - 2) = u(m, nx - 5)*r1_i(nx - 2) + u(m, nx - 4)*r2_i(nx - 2) + u(m, nx - 3)*r3_i(nx - 2) + u(m, nx - 2)*r4_i(nx - 2) + u(m, nx - 1)*r5_i(nx - 2) + u(m, nx)*r6_i(nx - 2)
                        f(m, nx - 1) = u(m, nx - 4)*r1_i(nx - 1) + u(m, nx - 3)*r2_i(nx - 1) + u(m, nx - 2)*r3_i(nx - 1) + u(m, nx - 1)*r4_i(nx - 1) + u(m, nx)*r5_i(nx - 1)
                        f(m, nx) = u(m, nx - 4)*r7_i(nx) + u(m, nx - 3)*r1_i(nx) + u(m, nx - 2)*r2_i(nx) + u(m, nx - 1)*r3_i(nx) + u(m, nx)*r4_i(nx)
                        ! r7(nx) contains 4. subdiagonal to allow for longer stencil at boundary
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end if

        end if

        return
    end subroutine MatMul_7d_sym_APU

end module FDM_MatMul
