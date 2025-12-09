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
!# 2025/12/08 - S. Deshpande
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
    use TLab_Constants, only: wp, wi, mas
    use Tlab_Type, only: fdm_integral_dt, fdm_integral_dt2
    use TLab_Time, only: pentadss_time, heptadss_time, tridss_time
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
        type(fdm_integral_dt2), intent(in) :: fdmi
        real(wp), dimension(1:len, 1:nmax, 1:klen, 1:ilen), intent(INOUT) :: f    ! RHS and solution

    ! -------------------------------------------------------------------
        integer(wi) ::  i, k, n, l    
    ! ###################################################################
    ! -----------------------------------------------------------------------
    ! Forward sweep
    ! -----------------------------------------------------------------------
         if (len <= 0) then
            goto 999
        end if
#ifdef USE_APU
        !$omp target teams distribute parallel do collapse(2) &
        !$omp private(i, k, n) &
        !$omp shared(ilen,klen,len,nmax,f,fdmi) &
        !$omp if (ilen*klen*len > mas)
#endif
        do i = 1, ilen
            do k = 1, klen
                do l = 1, len
                    do n = 3, nmax
                        f(l, n, k, i) = f(l, n, k, i) + fdmi%lhs(k, i ,n ,1)*f(l, n-1, k, i)
                    end do
                end do
        ! -----------------------------------------------------------------------
        ! Backward sweep
        ! -----------------------------------------------------------------------
                n = nmax - 1
                do l = 1, len
                    f(l, nmax, k, i) = f(l, nmax, k, i)*fdmi%lhs(k, i, nmax, 2)
                end do

                do l = 1, len
                    do n = nmax - 2, 1, -1
                        f(l, n, k, i) = f(l, n, k, i) + fdmi%lhs(k, i, n, 3)*f(l, n+1, k, i)*fdmi%lhs(k, i, n, 2)
                    end do
                end do
            end do
        end do
#ifdef USE_APU
        !$omp end target teams distribute parallel do
#endif
        999 continue
        return
    end subroutine TRIDSS_APU
!########################################################################


    subroutine PENTADSS_APU(len, nmax, klen, ilen, fdmi, f) !a = fdmi%lhs(2:, 1), b =fdmi%lhs(2:, 2), c =fdmi%lhs(2:, 3), d =fdmi%lhs(2:, 4), e =fdmi%lhs(2:, 5), result(:, 2:))        
        implicit none

        integer(wi) nmax, len, klen, ilen
        type(fdm_integral_dt2), intent(in) :: fdmi !real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e
        real(wp), dimension(1:len, 1:nmax, 1:klen, 1:ilen), intent(INOUT) :: f
        ! -----------------------------------------------------------------------
        integer(wi) n, i, k, l
        ! -----------------------------------------------------------------------
        ! Solve Ly=f, forward
        ! -----------------------------------------------------------------------
#ifdef USE_APU
        !$omp target teams distribute parallel do collapse(2) &
        !$omp private(i, k, n) &
        !$omp shared(ilen,klen,len,nmax,f,fdmi) &
        !$omp if (ilen*klen*len > mas)
#endif
        do i = 1, ilen
            do k = 1, klen
                do l = 1, len
                    f(l, 3, k, i) = f(l, 3, k, i) + f(l, 3 - 1, k, i)*fdmi%lhs(k, i, 3, 2)
                    do n = 4, nmax - 1
                        f(l, n, k, i) = f(l, n, k, i) + f(l, n - 1, k, i)*fdmi%lhs(k, i, n, 2) + f(l, n - 2, k, i)*fdmi%lhs(k, i, n, 1)
                    end do
                ! end do

                ! -----------------------------------------------------------------------
                ! Solve Ux=y, backward
                ! -----------------------------------------------------------------------
                ! do l = 1, len
                    f(l, nmax - 1, k, i) = f(l, nmax - 1, k, i)*fdmi%lhs(k, i, nmax - 1, 3)
                    f(l, nmax - 2, k, i) = (f(l, nmax - 2, k, i) + f(l, nmax - 2 + 1, k, i)*fdmi%lhs(k, i, nmax - 2, 4))*fdmi%lhs(k, i, nmax - 2, 3)
                    do n = nmax - 3, 2, -1
                        f(l, n, k, i) = (f(l, n, k, i) + f(l, n + 1, k, i)*fdmi%lhs(k, i, n, 4) + f(l, n + 2, k, i)*fdmi%lhs(k, i, n, 5))*fdmi%lhs(k, i, n, 3)
                    end do
                end do
            end do
        end do
#ifdef USE_APU
        !$omp end target teams distribute parallel do
#endif
        return
    end subroutine PENTADSS_APU
!########################################################################

    subroutine HEPTADSS_APU(len, nmax, klen, ilen, fdmi, frc)
        implicit none
        
        integer(wi) len, nmax, klen, ilen
        type(fdm_integral_dt2), intent(in) :: fdmi   !real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e
        real(wp), dimension(1:len, 1:nmax, 1:klen, 1:ilen), intent(INOUT) :: frc

    ! -----------------------------------------------------------------------
        integer(wi) n, ij, i, k
    ! #######################################################################
    ! -----------------------------------------------------------------------
    ! Solve Ly=frc, forward
    ! -----------------------------------------------------------------------
#ifdef USE_APU
        !$omp target teams distribute parallel do collapse(2) private(i, k, n) &
        !$omp if (ilen*klen*len > mas)
#endif
        do i = 1, ilen
            do k = 1, klen
                do ij = 1, len
                    frc(ij, 2, k, i) = frc(ij, 2, k, i)*fdmi%lhs(k, i, 1, 3) ! Normalize first eqn. See HEPTADFS
                    frc(ij, 3, k, i) = frc(ij, 3, k, i) - frc(ij, 2, k, i)*fdmi%lhs(k, i, 2, 3)
                    frc(ij, 4, k, i) = frc(ij, 4, k, i) - frc(ij, 3, k, i)*fdmi%lhs(k, i, 3, 3) - frc(ij, 2, k, i)*fdmi%lhs(k, i, 3, 2)
                end do

                do n = 5, nmax-1
                    do ij = 1, len
                        frc(ij, n, k, i) = frc(ij, n, k, i) - frc(ij, n-1, k, i)*fdmi%lhs(k, i, n, 3) - frc(ij, n-2, k, i)*fdmi%lhs(k, i, n, 2) - frc(ij, n - 3, k, i)*fdmi%lhs(k, i, n, 1)
                    end do
                end do
                ! -----------------------------------------------------------------------
                ! Solve Ux=y, backward
                ! -----------------------------------------------------------------------
                n = nmax - 1
                do ij = 1, len
                    frc(ij, n, k, i) = frc(ij, n, k, i)/fdmi%lhs(k, i, n, 4)
                    frc(ij, n-1, k, i) = (frc(ij, n-1, k, i) - frc(ij, n, k, i)*fdmi%lhs(k, i, n-1, 5))/fdmi%lhs(k, i, n-1, 4)
                    frc(ij, n-2, k, i) = (frc(ij, n-2, k, i) - frc(ij, n-1, k, i)*fdmi%lhs(k, i, n-2, 5) - frc(ij, n, k, i)*fdmi%lhs(k, i, n-2, 6))/fdmi%lhs(k, i, n-2, 4)
                end do

                do n = nmax - 4, 1, -1
#ifdef USE_APU
                    !$omp simd
#endif
                    do ij = 1, len
                        frc(ij, n, k, i) = (frc(ij, n, k, i) - frc(ij, n + 1, k, i)*fdmi%lhs(k, i, n, 5) - frc(ij, n + 2, k ,i)*fdmi%lhs(k, i, n, 6) - frc(ij, n + 3, k, i)*fdmi%lhs(k, i, n, 7))/fdmi%lhs(k, i, n, 4)
                    end do
                end do
            end do
        end do
#ifdef USE_APU
        !$omp end target teams distribute parallel do
#endif
        return
    end subroutine HEPTADSS_APU
end module LinearDss