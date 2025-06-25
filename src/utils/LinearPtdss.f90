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

module LinearPtdss
    PUBLIC :: PENTADSS_test
contains

!########################################################################

    subroutine PENTADSS_test(nx, len, ilen, fdmi, f) !a = fdmi%lhs(2:, 1), b =fdmi%lhs(2:, 2), c =fdmi%lhs(2:, 3), d =fdmi%lhs(2:, 4), e =fdmi%lhs(2:, 5), result(:, 2:))
        use TLab_Constants, only: wp, wi
        use Tlab_Type, only: fdm_integral_dt
        
        implicit none

        integer(wi) nx, len, ilen
        type(fdm_integral_dt), intent(in) :: fdmi(:) !real(wp), dimension(nmax), intent(IN) :: a, b, c, d, e
        real(wp), dimension(len, nx, ilen), intent(INOUT) :: f
        ! -----------------------------------------------------------------------
        integer(wi) n, i, l, omp_srt, omp_end, omp_siz

        !! !$omp parallel default(none) &
        !! !$omp shared(f,a,b,c,d,e,nmax,len) &
        !! !$omp private(l,n,omp_srt,omp_end,omp_siz)

        !  CALL TLab_OMP_PARTITION(len,omp_srt,omp_end,omp_siz)
        omp_srt = 1
        omp_end = len
        omp_siz = len

        ! #######################################################################
        ! -----------------------------------------------------------------------
        ! Solve Ly=f, forward
        ! -----------------------------------------------------------------------
        do i = 1, ilen
            n = 3
            do l = omp_srt, omp_end
                f(l, n, i) = f(l, n, i) + f(l, n - 1, i)*fdmi(i)%lhs(3, 2)
            end do

            do n = 4, nx - 1
                do l = omp_srt, omp_end
                    f(l, n, i) = f(l, n, i) + f(l, n - 1, i)*fdmi(i)%lhs(n, 2) + f(l, n - 2, i)*fdmi(i)%lhs(n, 1)
                end do
            end do

            ! -----------------------------------------------------------------------
            ! Solve Ux=y, backward
            ! -----------------------------------------------------------------------
            n = nx - 1
            do l = omp_srt, omp_end
                f(l, n, i) = f(l, n, i)*fdmi(i)%lhs(n, 3)
            end do

            n = nx - 1 - 1
            do l = omp_srt, omp_end
                f(l, n, i) = (f(l, n, i) + f(l, n + 1, i)*fdmi(i)%lhs(n, 4))*fdmi(i)%lhs(n, 3)
            end do

            do n = nx - 3, 2, -1
                do l = omp_srt, omp_end
                    f(l, n, i) = (f(l, n, i) + f(l, n + 1, i)*fdmi(i)%lhs(n, 4) + f(l, n + 2, i)*fdmi(i)%lhs(n, 5))*fdmi(i)%lhs(n, 3)
                end do
            end do
        end do
        !! !$omp end parallel
        return
    end subroutine PENTADSS_test
!########################################################################
end module LinearPtdss