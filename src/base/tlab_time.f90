module TLab_Time
    use TLab_Constants, only: wp, wi, sp
    implicit none

    integer(wi) :: itime                ! iteration number
    real(wp) :: rtime                   ! physical time

! ###################################################################
! APU profiling timings
! ###################################################################    
    real(wp) :: trans_time, tridss_time, tridpss_time
    real(wp) :: mat5dantisym_time,mat5dsym_time,mat3dadd_time,mat3d_time, t_map_in, t_compute, t_map_out
    real(wp) :: elliptic_time, elliptic_setup_time
    real(wp) :: fdm_solve2_time
end module TLab_Time
