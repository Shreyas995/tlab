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
    real(wp), DIMENSION(:), ALLOCATABLE ::  runtime
    
    ! -------------------------------------------------------------------
    
    PUBLIC :: TLab_Time_Data
    PUBLIC :: TLab_Time_Initialize
    ! -------------------------------------------------------------------
    contains

    ! ###################################################################
    subroutine TLab_Time_Initialize()
        ALLOCATE(runtime(10))
        add_time = 0.0_wp
        pps_time = 0.0_wp
        trans_time = 0.0_wp
        tridss_time = 0.0_wp
        tridssadd_time = 0.0_wp
        tridpss_time = 0.0_wp
        tridpssadd_time = 0.0_wp
        mat5dantisym_time = 0.0_wp
        mat5dsym_time = 0.0_wp
        mat3dadd_time = 0.0_wp
        mat3d_time = 0.0_wp
        fdm_solve2_time = 0.0_wp
        elliptic_time = 0.0_wp
        elliptic_setup_time = 0.0_wp
        pentadss_time = 0.0_wp
        heptadss_time = 0.0_wp
        burger1D_time = 0.0_wp
        source_time = 0.0_wp
        rhs_incomp_time = 0.0_wp
        courant_time = 0.0_wp
    end subroutine TLab_Time_Initialize

    ! ###################################################################
    subroutine TLab_Time_Data()
        PRINT 100,SUM(runtime)/10, MINVAL(runtime),MAXVAL(runtime)
        PRINT 101, 'Transpose      ',trans_time/10, 100*trans_time/SUM(runtime)
        PRINT 101, 'Addition      ',add_time/10, 100*add_time/SUM(runtime) 
        PRINT 101, 'PressurePoiss ',pps_time/10, 100*pps_time/SUM(runtime) 
        PRINT 101, 'TRIDSS        ',tridss_time/10, 100*tridss_time/SUM(runtime)
        PRINT 101, 'TRIDSSADD     ',tridssadd_time/10, 100*tridssadd_time/SUM(runtime)
        PRINT 101, 'TRIDPSS       ',tridpss_time/10, 100*tridpss_time/SUM(runtime)
        PRINT 101, 'TRIDPSSADD    ',tridpssadd_time/10, 100*tridpssadd_time/SUM(runtime)
        PRINT 101, 'MATMUL5D_ANTI ',mat5dantisym_time/10, 100*mat5dantisym_time/SUM(runtime)
        PRINT 101, 'MATMUL5D_SYM  ',mat5dsym_time/10, 100*mat5dsym_time/SUM(runtime)
        PRINT 101, 'MATMUL3D_ADD  ',mat3dadd_time/10, 100*mat3dadd_time/SUM(runtime)
        PRINT 101, 'MATMUL3D      ',mat3d_time/10, 100*mat3d_time/SUM(runtime)
        PRINT 101, 'FDM_SOLVE2    ',fdm_solve2_time/10, 100*fdm_solve2_time/SUM(runtime)
        PRINT 101, 'ELLIPTIC      ',elliptic_time/10, 100*elliptic_time/SUM(runtime)
        PRINT 101, 'ELLIPTIC_SETUP',elliptic_setup_time/10, 100*elliptic_setup_time/SUM(runtime)
        PRINT 101, 'PENTADSS      ',pentadss_time/10, 100*pentadss_time/SUM(runtime)
        PRINT 101, 'HEPTADSS      ',heptadss_time/10, 100*heptadss_time/SUM(runtime)
        PRINT 101, 'BURGER1D      ',burger1D_time/10, 100*burger1D_time/SUM(runtime)
        PRINT 101, 'SOURCE        ',source_time/10, 100*source_time/SUM(runtime)
        PRINT 101, 'RHS_INCOMP    ',rhs_incomp_time/10, 100*rhs_incomp_time/SUM(runtime)
        PRINT 101, 'COURANT       ',courant_time/10, 100*courant_time/SUM(runtime)
        100 FORMAT('T MEAN|MIN|MAX                 [s]:', F7.3, 1x, F7.3, 1x , F7.3)
        101 FORMAT('Time per run in ',A15,'[s]:', F9.5,'s (', F7.4,'%)') 
    end subroutine TLab_Time_Data


end module TLab_Time

!########################################################################