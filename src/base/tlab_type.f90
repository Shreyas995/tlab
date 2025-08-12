!########################################################################
!# HISTORY / AUTHORS
!#
!# 2025/05/06 - S. Deshpande
!#              Created
!#
!########################################################################
!# DESCRIPTION OF MODLE
!#   Add FDM realted data types here
!#                    
!#
!########################################################################

module Tlab_Type

    use  TLab_Constants, only : wp

    implicit none 

    save 
    
    type, public :: fdm_integral_dt
        sequence
        integer mode_fdm                            ! original finite-difference method; only informative
        real(wp) :: lambda                          ! constant of the equation
        integer :: bc                               ! type of boundary condition, [ BCS_MIN, BCS_MAX ]
        real(wp) :: rhs_b(1:5, 0:7), rhs_t(0:4, 8)  ! # of diagonals is 7, # rows is 7/2+1
        real(wp), allocatable :: lhs(:, :)          ! Often overwritten to LU decomposition.
        real(wp), allocatable :: rhs(:, :)
    end type fdm_integral_dt

end module Tlab_Type
