#include "dns_const.h"
#define C_FILE_LOC "VBURGERS_MEASURE"

   program VBURGERS

   use TLab_Constants, only: wp, wi, big_wp, gfile, ifile
   use TLab_Constants, only: BCS_NN, wfile
   use TLab_Time, only: itime
   use TLab_Time, only: mat5dantisym_time,mat5dsym_time,mat3dadd_time,mat3d_time
   use TLab_Time, only: trans_time, tridss_time, tridpss_time
   use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start, stagger_on
   use TLab_Memory, only: TLab_Initialize_Memory
   use TLab_Memory, only: imax, jmax, kmax, inb_txc
   use TLab_Arrays
   use TLab_Pointers_3D, only: tmp1
   use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
   use MPI
   use TLabMPI_PROCS, only: TLabMPI_Initialize
   use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
   use TLabMPI_VARS
#endif
   use FDM, only: g, FDM_Initialize
    use FDM_Integral, only: fdm_Int0
    use NavierStokes, only: NavierStokes_Initialize_Parameters, visc
    use TLab_Grid
    use IO_Fields
    use OPR_PARTIAL
    use OPR_Burgers
    use OPR_FILTERS
    use OPR_ELLIPTIC
    use OPR_FOURIER
    use TLab_Background, only: TLab_Initialize_Background

   implicit none
  
   !$omp requires unified_shared_memory
  
#ifndef USE_MPI
   integer(wi), parameter :: ims_pro = 0
#endif

   real(wp), dimension(:, :, :), pointer :: a, b, c, d, e
   real(wp), dimension(:, :), allocatable :: bcs_hb, bcs_ht

   integer(wi) i, j, k, ig, bcs(2, 2)

   integer ibc

   integer irun,nrun,stat
   integer clock_0, clock_1,clock_cycle
   integer clock_add0, clock_add1
   CHARACTER(len=64) :: nrun_string 
   real(wp), DIMENSION(:), ALLOCATABLE ::  runtime 
   real(wp) add_time, pps_time
   real(wp) dummy, error, params(0)

   trans_time = 0.0_wp
   tridss_time = 0.0_wp
   mat5dantisym_time = 0.0_wp
   mat5dsym_time = 0.0_wp
   mat3dadd_time = 0.0_wp
   mat3d_time = 0.0_wp
   tridpss_time = 0.0_wp
   add_time = 0.0_wp
   pps_time = 0.0_wp

   call SYSTEM_CLOCK(clock_0,clock_cycle)
   IF ( COMMAND_ARGUMENT_COUNT() .GE. 1 ) THEN 
      call GetArg(1,nrun_String)
      read(nrun_string,*) nrun
   ELSE
      nrun = 1 
   ENDIF

   PRINT *,'EXECUTING ',nrun, ' RUNS for Performance Measurement'

   ALLOCATE(runtime(nrun))


   ! ###################################################################
   call TLab_Start()

   call TLab_Initialize_Parameters(ifile) 
#ifdef USE_MPI
   call TLabMPI_Initialize(ifile)
   call TLabMPI_Transpose_Initialize(ifile)
#endif
   call NavierStokes_Initialize_Parameters(ifile) 

   inb_txc = 6

   call TLab_Initialize_Memory(__FILE__)

   allocate (bcs_ht(imax, kmax), bcs_hb(imax, kmax))
   a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 2)
   b(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
   c(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
   d(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 5)
   e(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 6)

   visc = 1.0_wp/big_wp    ! inviscid

   call TLab_Grid_Read(gfile, x, y, z, [g(1)%size, g(2)%size, g(3)%size])

   call FDM_Initialize(x, g(1))
   call FDM_Initialize(y, g(2), fdm_Int0)
   call FDM_Initialize(z, g(3))

   call TLab_Initialize_Background(ifile) 
   call OPR_Burgers_Initialize(ifile)
   call OPR_Elliptic_Initialize(ifile)

   ! Staggering of the pressure grid not implemented here
   if (stagger_on) then
      call TLab_Write_ASCII(wfile, C_FILE_LOC//'. Staggering of the pressure grid not implemented here.')
      stagger_on = .false. ! turn staggering off for OPR_Poisson_FourierXZ_Factorize(...)
   end if

   if (any(PressureFilter%type /= DNS_FILTER_NONE)) then
      call TLab_Write_ASCII(wfile, C_FILE_LOC//'. Pressure and dpdy Filter not implemented here.')
   end if

   call OPR_FOURIER_INITIALIZE()
   call OPR_CHECK()
   ! call BOUNDARY_BCS_INITIALIZE()

   bcs = 0

   ! ###################################################################
   ! Define forcing term
   ! ###################################################################
   ! call IO_READ_FIELDS('field.inp', IO_SCAL, imax, jmax, kmax, 1, 0, a)
   call IO_Read_Fields('field.inp', imax, jmax, kmax, itime, 1, 0, a, params)

   visc = 1.0_wp/big_wp

   ! For the Pressure Poisson solver
   ibc = BCS_NN
   bcs_hb(:, :) = a(:, 1, :)
   bcs_ht(:, :) = a(:, jmax, :)
   ! bcs_hb = 0.0_wp; bcs_ht = 0.0_wp

   if (ims_pro == 0) then
      write (*, *) '-------------------------------------------------------- '
   end if

   ! Call everything once to get the memory initialized on the APUs
   ! (avoids measuring the first-touch penalty) 
   call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), a, b, c)
   call output_sum(b, 'OPR_PAR_X')
   call output_sum(c, 'OPR_PAR_X')
   if (ims_pro == 0) write (*, *) '------------------- '

   call OPR_BURGERS_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, a, a, c, tmp1)
   call output_sum(c, 'OPR_BUR_X')
   call output_sum(tmp1, 'OPR_BUR_X')
   if (ims_pro == 0) write (*, *) '------------------- '
   
   call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), a, b, c)
   call output_sum(b, 'OPR_PAR_Y')
   call output_sum(c, 'OPR_PAR_Y')
   if (ims_pro == 0) write (*, *) '------------------- '
   
   call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, a, a, c, tmp1)
   call output_sum(c, 'OPR_BUR_Y')
   call output_sum(tmp1, 'OPR_BUR_Y')
   if (ims_pro == 0) write (*, *) '------------------- '
   
   call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), a, b, c)
   call output_sum(b, 'OPR_PAR_Z')
   call output_sum(c, 'OPR_PAR_Z')
   if (ims_pro == 0) write (*, *) '------------------- '
   
   call OPR_BURGERS_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, a, a, c, tmp1)
   call output_sum(c, 'OPR_BUR_Z')
   call output_sum(tmp1, 'OPR_BUR_Z')
   if (ims_pro == 0) write (*, *) '------------------- '

   e = a
   call OPR_Poisson(imax, jmax, kmax, g, ibc, e, b, d, bcs_hb, bcs_ht, c)
   call output_sum(c, 'OPR_ELLIP') ! dpdy
   call output_sum(e, 'OPR_ELLIP') ! p

   if (ims_pro == 0) then
      write (*, *) '-------------------------------------------------------- '
   end if

   ! reset time for transpositions to avoid measuring first-touch penalty
   ! (else nrun needs to be incremented by one in normalization) 
   
   trans_time = 0.0_wp
   tridss_time = 0.0_wp
   tridpss_time = 0.0_wp
   mat5dantisym_time = 0.0_wp
   mat5dsym_time = 0.0_wp
   mat3dadd_time = 0.0_wp
   mat3d_time = 0.0_wp
   MatMul_3d_time = 0.0_wp

   DO irun=1,nrun


      call SYSTEM_CLOCK(clock_0) 
      ! ###################################################################
      call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), a, b, c)
      do k = 1, kmax
         do j = 1, jmax
            do i = 1, imax
               b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
            end do
         end do
      end do

      call OPR_BURGERS_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, a, a, c, tmp1)
      call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), a, b, c)
      CALL SYSTEM_CLOCK(clock_add0) 
      !$omp target teams distribute parallel do collapse(3) default (none) &
      !$omp private(i,j,k) &
      !$omp shared(a,b,c,imax,jmax,kmax,visc) 
      do k = 1, kmax
         do j = 1, jmax
            do i = 1, imax
               b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
            end do
         end do
      end do
      !$omp end target teams distribute parallel do
      CALL SYSTEM_CLOCK(clock_add1)
      add_time = add_time + real(clock_add1 - clock_add0) / clock_cycle 
      call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, a, a, c, tmp1)

      ! ###################################################################
      if (g(3)%size > 1) then
         CALL SYSTEM_CLOCK(clock_add0)
         !$omp target teams distribute parallel do collapse(3) default(none) &
         !$omp private(i,j,k) &
         !$omp shared(a,b,c,imax,jmax,kmax,visc) 
         do k = 1, kmax
            do j = 1, jmax
               do i = 1, imax
                  b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
               end do
            end do
         end do
         !$omp end target teams distribute parallel do
         CALL SYSTEM_CLOCK(clock_add1)
         add_time = add_time + real(clock_add1 - clock_add0) / clock_cycle 
         call OPR_BURGERS_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, a, a, c, tmp1)

         ! ------------------------------------------
         
         ! call Pressure Poisson solver
         CALL SYSTEM_CLOCK(clock_add0)
         call OPR_Poisson(imax, jmax, kmax, g, ibc, e, b, d, bcs_hb, bcs_ht, c)
         CALL SYSTEM_CLOCK(clock_add1)
         pps_time = pps_time + real(clock_add1 - clock_add0) / clock_cycle 
      
      end if

      call SYSTEM_CLOCK(clock_1)
      runtime(irun) = real(clock_1-clock_0)/clock_cycle
   
   end do
   
   PRINT 100,SUM(runtime)/nrun, MINVAL(runtime),MAXVAL(runtime)
   PRINT 101, 'Transpos      ',trans_time/nrun, 100*trans_time/SUM(runtime)
   PRINT 101, 'Addition      ',add_time/nrun, 100*add_time/SUM(runtime) 
   PRINT 101, 'PressurePoiss ',pps_time/nrun, 100*pps_time/SUM(runtime) 
   PRINT 101, 'TRIDSS        ',tridss_time/nrun, 100*tridss_time/SUM(runtime)
   PRINT 101, 'TRIDPSS       ',tridpss_time/nrun, 100*tridpss_time/SUM(runtime)
   PRINT 101, 'MATMUL5D_ANTI ',mat5dantisym_time/nrun, 100*mat5dantisym_time/SUM(runtime)
   PRINT 101, 'MATMUL5D_SYM  ',mat5dsym_time/nrun, 100*mat5dsym_time/SUM(runtime)
   PRINT 101, 'MATMUL3D_ADD  ',mat3dadd_time/nrun, 100*mat3dadd_time/SUM(runtime)
   PRINT 101, 'MATMUL3D      ',mat3d_time/nrun, 100*mat3d_time/SUM(runtime)
100 FORMAT('T MEAN|MIN|MAX                 [s]:', F7.3, 1x, F7.3, 1x , F7.3)
101 FORMAT('Time per run in ',A15,'[s]:', F9.5,'s (', F7.4,'%)') 

   call TLab_STOP(0)

   ! ----------------------------------------------------

   contains 
      subroutine output_sum(fld, name)

         implicit none
      
         real(wp), intent(in) :: fld(imax,jmax,kmax)
         character(len=*), intent(in) :: name
         
         real(wp) dummy
#ifdef USE_MPI
         real(wp) dummy2
#endif
         dummy = 0.0_wp
         
         dummy = sum(b)
#ifdef USE_MPI
         dummy2 = 0_wp
         dummy2 = dummy
         call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
         if (ims_pro == 0) then
            write (*, '(1X,A,A10,ES22.14)') 'Sum ', name, dummy
         end if

      end subroutine output_sum


   end program VBURGERS


