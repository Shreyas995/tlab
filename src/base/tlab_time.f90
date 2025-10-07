module TLab_Time
    use TLab_Constants, only: wp, wi, sp
    implicit none

    integer(wi) :: itime                ! iteration number
    real(wp) :: rtime                   ! physical time

! ###################################################################
! APU profiling timings
! ###################################################################
    real(wp) :: add_time, pps_time
    real(wp) :: trans_time, tridss_time, tridssadd_time, tridpss_time, tridpssadd_time
    real(wp) :: mat5dantisym_time, mat5dsym_time, mat3dadd_time, mat3d_time
    real(wp) :: elliptic_time, elliptic_setup_time
    real(wp) :: fdm_solve2_time, pentadss_time, heptadss_time, burger1D_time, source_time, rhs_incomp_time, courant_time
end module TLab_Time

!########################################################################