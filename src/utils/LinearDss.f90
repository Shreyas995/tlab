!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/05/21 - J.P. Mellado
!#              Created
!# 2011/11/01 - C. Ansorge
!#              OpenMP Optimization
!#
!########################################################################
!# DESCRIPTION
!#
!# Solve LU of pentadiagonal system
!#
!########################################################################
!# ARGUMENTS
!#
!# nmax    In     Size of the pentadiagonal system
!# len     In     Number of simultaneous systems to be solved
!# f       In     Array with forcing term
!#         Out    Array with solution
!#
!########################################################################

! #######################################################################
! LU factorization stage
! #######################################################################
#include "dns_error.h"

module LinearDss
    use TLab_Constants, only: wp, wi
    use Tlab_Type, only: fdm_integral_dt

    PUBLIC :: TRIDSS_APU
    PUBLIC :: PENTADSS_APU
    PUBLIC :: HEPTADSS_APU
contains

!########################################################################

    subroutine TRIDSS_APU(len, nmax, klen, ilen, fdmi, f)
        implicit none
        integer(wi), intent(IN) :: nmax  ! dimension of tridiagonal systems
        integer(wi), intent(IN) :: len   ! number of systems to be solved
        integer(wi), intent(IN) :: klen  ! number of equations in each system
        integer(wi), intent(IN) :: ilen  ! number of equations in each system
        type(fdm_integral_dt), intent(in) :: fdmi(:,:)
        real(wp), dimension(1:len, 1:nmax, 1:klen, 1:ilen), intent(INOUT) :: f    ! RHS and solution

    ! -------------------------------------------------------------------
        integer(wi) ::  i, k, n, l
        integer(wi) :: srt, end, siz
        real(wp) :: dummy1, dummy2

    ! ###################################################################
    ! -----------------------------------------------------------------------
    ! Forward sweep
    ! -----------------------------------------------------------------------

        if (len <= 0) then
            goto 999
        end if

        do i = 1, ilen
            do k = 1, klen
                do n = 3, nmax
                    do l = 1, len
                        f(l, n, k, i) = f(l, n, k, i) + fdmi(k,i)%lhs(n, 1)*f(l, n-1, k, i)
                    end do
                end do
            end do
        end do

    ! -----------------------------------------------------------------------
    ! Backward sweep
    ! -----------------------------------------------------------------------
        n = nmax - 1
        do i = 1, ilen
            do k = 1, klen
                do l = 1, len
                    f(l, nmax, k, i) = f(l, nmax, k, i)*fdmi(k,i)%lhs(nmax, 2)
                end do
            end do
        end do

        do i = 1, ilen
            do k = 1, klen
                do n = nmax - 2, 1, -1
                    do l = 1, len
                        f(l, n, k, i) = f(l, n, k, i) + fdmi(k,i)%lhs(n, 3)*f(l, n+1, k, i)*fdmi(k,i)%lhs(n, 2)
                    end do
                end do
            end do
        end do
    999 continue

        return
    end subroutine TRIDSS_APU
!########################################################################


    subroutine PENTADSS_APU(len, nmax, klen, ilen, fdmi, f) !a = fdmi%lhs(2:, 1), b =fdmi%lhs(2:, 2), c =fdmi%lhs(2:, 3), d =fdmi%lhs(2:, 4), e =fdmi%lhs(2:, 5), result(:, 2:))        
        implicit none

        integer(wi) nmax, len, klen, ilen
        type(fdm_integral_dt), intent(in) :: fdmi(:,:) !real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e
        real(wp), dimension(1:len, 1:nmax, 1:klen, 1:ilen), intent(INOUT) :: f
        ! -----------------------------------------------------------------------
        integer(wi) n, i, k, l

        ! #######################################################################
        ! -----------------------------------------------------------------------
        ! Solve Ly=f, forward
        ! -----------------------------------------------------------------------
        do i = 1, ilen
            do k = 1, klen
                n = 3
                do l = 1, len
                    f(l, n, k, i) = f(l, n, k, i) + f(l, n - 1, k, i)*fdmi(k,i)%lhs(3, 2)
                end do

                do n = 4, nmax - 1
                    do l = 1, len
                        f(l, n, k, i) = f(l, n, k, i) + f(l, n - 1, k, i)*fdmi(k,i)%lhs(n, 2) + f(l, n - 2, k, i)*fdmi(k,i)%lhs(n, 1)
                    end do
                end do

                ! -----------------------------------------------------------------------
                ! Solve Ux=y, backward
                ! -----------------------------------------------------------------------
                n = nmax - 1
                do l = 1, len
                    f(l, n, k, i) = f(l, n, k, i)*fdmi(k,i)%lhs(n, 3)
                end do

                n = nmax - 1 - 1
                do l = 1, len
                    f(l, n, k, i) = (f(l, n, k, i) + f(l, n + 1, k, i)*fdmi(k,i)%lhs(n, 4))*fdmi(k,i)%lhs(n, 3)
                end do

                do n = nmax - 3, 2, -1
                    do l = 1, len
                        f(l, n, k,i) = (f(l, n, k, i) + f(l, n + 1, k, i)*fdmi(k,i)%lhs(n, 4) + f(l, n + 2, k, i)*fdmi(k,i)%lhs(n, 5))*fdmi(k,i)%lhs(n, 3)
                    end do
                end do
            end do
        end do
        !! !$omp end parallel
        return
    end subroutine PENTADSS_APU
!########################################################################

    subroutine HEPTADSS_APU(len, nmax, klen, ilen, fdmi, frc)
        implicit none
        
        integer(wi) len, nmax, klen, ilen
        type(fdm_integral_dt), intent(in) :: fdmi(:,:) !real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e
        real(wp), dimension(1:len, 1:nmax, 1:klen, 1:ilen), intent(INOUT) :: frc

    ! -----------------------------------------------------------------------
        integer(wi) n, ij, i, k

    ! #######################################################################
    ! -----------------------------------------------------------------------
    ! Solve Ly=frc, forward
    ! -----------------------------------------------------------------------
        do i = 1, ilen
            do k = 1, klen
                do ij = 1, len
                    frc(ij, 2, k, i) = frc(ij, 2, k, i)*fdmi(k,i)%lhs(1, 3) ! Normalize first eqn. See HEPTADFS
                    frc(ij, 3, k, i) = frc(ij, 3, k, i) - frc(ij, 2, k, i)*fdmi(k,i)%lhs(2, 3)
                    frc(ij, 4, k, i) = frc(ij, 4, k, i) - frc(ij, 3, k, i)*fdmi(k,i)%lhs(3, 3) - frc(ij, 2, k, i)*fdmi(k,i)%lhs(3, 2)
                end do
            end do
        end do

        do i = 1, ilen
            do k = 1, klen
                do n = 5, nmax-1
                    do ij = 1, len
                        frc(ij, n, k, i) = frc(ij, n, k, i) - frc(ij, n-1, k, i)*fdmi(k,i)%lhs(n, 3) - frc(ij, n-2, k, i)*fdmi(k,i)%lhs(n, 2) - frc(ij, n - 3, k, i)*fdmi(k,i)%lhs(n, 1)
                    end do
                end do
            end do
        end do

    ! -----------------------------------------------------------------------
    ! Solve Ux=y, backward
    ! -----------------------------------------------------------------------
        do i = 1, ilen
            do k = 1, klen
                n = nmax - 1
                do ij = 1, len
                    frc(ij, n, k, i) = frc(ij, n, k, i)/fdmi(k,i)%lhs(n, 4)
                    frc(ij, n-1, k, i) = (frc(ij, n-1, k, i) - frc(ij, n, k, i)*fdmi(k,i)%lhs(n-1, 5))/fdmi(k,i)%lhs(n-1, 4)
                    frc(ij, n-2, k, i) = (frc(ij, n-2, k, i) - frc(ij, n-1, k, i)*fdmi(k,i)%lhs(n-2, 5) - frc(ij, n, k, i)*fdmi(k,i)%lhs(n-2, 6))/fdmi(k,i)%lhs(n-2, 4)
                end do
            end do
        end do

        do i = 1, ilen
            do k = 1, klen
                do n = nmax - 4, 1, -1
                    do ij = 1, len
                        frc(ij, n, k, i) = (frc(ij, n, k, i) - frc(ij, n + 1, k, i)*fdmi(k,i)%lhs(n, 5) - frc(ij, n + 2, k ,i)*fdmi(k,i)%lhs(n, 6) - frc(ij, n + 3, k, i)*fdmi(k,i)%lhs(n, 7))/fdmi(k,i)%lhs(n, 4)
                    end do
                end do
            end do
        end do

        return
    end subroutine HEPTADSS_APU
end module LinearDss