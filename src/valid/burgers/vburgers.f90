#include "dns_const.h"

program VBURGERS

    use TLab_Constants, only: wp, wi, big_wp, gfile, ifile
    use TLAB_VARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLab_Arrays
    use TLab_Pointers_3D, only: tmp1
#ifdef USE_MPI
    use MPI
    use TLabMPI_PROCS
    use TLabMPI_VARS
#endif
    use FDM, only: g,  FDM_Initialize
    use IO_FIELDS
    use OPR_PARTIAL
    use OPR_BURGERS
    use OPR_FILTERS
    use TLab_Background, only: TLab_Initialize_Background
    implicit none

#ifdef USE_MPI
    real(wp) error2, dummy2
#else
    integer(wi), parameter :: ims_pro = 0
#endif

    real(wp), dimension(:, :, :), pointer :: a, b, c

    integer(wi) i, j, k, ig, bcs(2, 2)
    real(wp) dummy, error

! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
#endif
    call NavierStokes_Initialize_Parameters(ifile)

    inb_txc = 4

    call TLab_Initialize_Memory(__FILE__)

    a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 2)
    b(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    c(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)

    visc = 1.0_wp/big_wp    ! inviscid

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:,1), wrk1d(:,2), wrk1d(:,3))
    call FDM_Initialize(x, g(1), wrk1d(:,1), wrk1d(:,4))
    call FDM_Initialize(y, g(2), wrk1d(:,2), wrk1d(:,4))
    call FDM_Initialize(z, g(3), wrk1d(:,3), wrk1d(:,4))

    call TLab_Initialize_Background(ifile)

    bcs = 0

    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig))
    end do

! ###################################################################
! Define forcing term
! ###################################################################
    call IO_READ_FIELDS('field.inp', IO_SCAL, imax, jmax, kmax, 1, 0, a)

    visc = 1.0_wp/big_wp

! ###################################################################
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), a, b, c)
    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
                b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
                ! b(i, j, k) = b(i, j, k)*visc*ribackground(j) - a(i, j, k)*c(i, j, k)
            end do
        end do
    end do
    call IO_WRITE_FIELDS('fieldXdirect.out', IO_SCAL, imax, jmax, kmax, 1, b)

    call OPR_BURGERS_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(1), a, a, c, tmp1)
    call IO_WRITE_FIELDS('fieldXburgers.out', IO_SCAL, imax, jmax, kmax, 1, c)

    c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
    error2 = error; dummy2 = dummy
    call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
    if (ims_pro == 0) then
        write (*, *) 'Relative error X ...........: ', sqrt(error)/sqrt(dummy)
    end if
    call IO_WRITE_FIELDS('fieldX.dif', IO_SCAL, imax, jmax, kmax, 1, c)

! ###################################################################
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), a, b, c)
    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
                b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
                ! b(i, j, k) = b(i, j, k)*visc*ribackground(j) - a(i, j, k)*c(i, j, k)
            end do
        end do
    end do
    call IO_WRITE_FIELDS('fieldYdirect.out', IO_SCAL, imax, jmax, kmax, 1, b)

    call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(2), a, a, c, tmp1)
    call IO_WRITE_FIELDS('fieldYburgers.out', IO_SCAL, imax, jmax, kmax, 1, c)

    c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
    error2 = error; dummy2 = dummy
    call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
    if (ims_pro == 0) then
        write (*, *) 'Relative error Y ...........: ', sqrt(error)/sqrt(dummy)
    end if
    call IO_WRITE_FIELDS('fieldY.dif', IO_SCAL, imax, jmax, kmax, 1, c)

! ###################################################################
    if (g(3)%size > 1) then

        call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), a, b, c)
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
                    ! b(i, j, k) = b(i, j, k)*visc*ribackground(j) - a(i, j, k)*c(i, j, k)
                end do
            end do
        end do
        call IO_WRITE_FIELDS('fieldZdirect.out', IO_SCAL, imax, jmax, kmax, 1, b)

        call OPR_BURGERS_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(3), a, a, c, tmp1)
        call IO_WRITE_FIELDS('fieldZburgers.out', IO_SCAL, imax, jmax, kmax, 1, c)

        c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
        error2 = error; dummy2 = dummy
        call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
        if (ims_pro == 0) then
            write (*, *) 'Relative error Z ...........: ', sqrt(error)/sqrt(dummy)
        end if
        call IO_WRITE_FIELDS('fieldZ.dif', IO_SCAL, imax, jmax, kmax, 1, c)

    end if

    call TLab_Stop(0)
end program VBURGERS
