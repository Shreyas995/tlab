! ###################################################################
! Jet Statistic
! ###################################################################
module Avg_Spatial
    use TLab_Constants, only: wp
    implicit none
    private

    integer, parameter :: MAX_STATS_SPATIAL = 100 ! Running statistics

    integer, public :: nstatavg, statavg(MAX_STATS_SPATIAL) ! Ox planes at which to accumulate statistics
    integer, public :: nstatavg_points                      ! number of accumulated points
    integer, public :: istattimeorg                         ! time at which accumulation of statistics started
    real(wp), public :: rstattimeorg

    public :: IO_WRITE_AVG_SPATIAL, IO_READ_AVG_SPATIAL

contains

#include "dns_error.h"
#include "avgij_map.h"

#define LOC_UNIT_ID 57
#define LOC_STATUS 'unknown'

    !########################################################################
    !#
    !# Note that the array mean1d is duplicated in each processor, and only
    !# that from PE0 is written
    !#
    !########################################################################
    subroutine IO_WRITE_AVG_SPATIAL(name, mean_flow, mean_scal)
        use TLab_Constants, only: lfile, wp, wi
        use TLab_Time, only: itime, rtime
        use TLab_Memory, only: jmax, inb_scal

#ifdef USE_MPI
        use TLabMPI_VARS, only: ims_pro
#endif

        implicit none

        character*(*) name
        real(wp) mean_flow(nstatavg, jmax, MA_MOMENTUM_SIZE)
        real(wp) mean_scal(nstatavg, jmax, MS_SCALAR_SIZE, inb_scal)

        ! -------------------------------------------------------------------
        character*128 :: line
        integer(wi) nstat
        integer, parameter :: i0 = 0

        ! ###################################################################
        line = 'Writing field '//trim(adjustl(name))//'...'

#ifdef USE_MPI
        if (ims_pro == 0) then
#endif
#include "dns_open_file.h"
            nstat = MA_MOMENTUM_SIZE + MS_SCALAR_SIZE*inb_scal
            call WRT_STHD(LOC_UNIT_ID, i0, &
                          itime, rtime, istattimeorg, rstattimeorg, &
                          nstatavg, jmax, nstat, nstatavg_points, statavg)
            write (LOC_UNIT_ID) mean_flow
            write (LOC_UNIT_ID) mean_scal
            close (LOC_UNIT_ID)

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine IO_WRITE_AVG_SPATIAL

    ! ###################################################################
    ! ###################################################################
    subroutine WRT_STHD(unit, irec, &
                        iter, rtime, iterorg, rtimeorg, &
                        nstatavg, jmax, nstat, nstatavg_points, statavg)

        use TLab_Constants, only: wp, wi, sizeofint, sizeofreal
        implicit none

        integer(wi) unit, irec
        integer(wi) iter, iterorg
        integer(wi) nstatavg, jmax, nstat, nstatavg_points
        real(wp) rtime, rtimeorg
        integer(wi) statavg(nstatavg)

        real(wp) tmp(1)
        integer(wi) reclen

        tmp(1) = rtime
        reclen = SIZEOFINT + SIZEOFREAL
        if (irec == 1) then
            write (unit) reclen
            write (unit) iter
            write (unit) rtime
            write (unit) reclen
        else
            write (unit) iter, tmp
        end if

        tmp(1) = rtimeorg
        reclen = SIZEOFINT + SIZEOFREAL
        if (irec == 1) then
            write (unit) reclen
            write (unit) iterorg
            write (unit) rtimeorg
            write (unit) reclen
        else
            write (unit) iterorg, tmp
        end if

        reclen = 4*SIZEOFINT
        if (irec == 1) then
            write (unit) reclen
            write (unit) nstatavg
            write (unit) jmax
            write (unit) nstat
            write (unit) nstatavg_points
            write (unit) reclen
        else
            write (unit) nstatavg, jmax, nstat, nstatavg_points
        end if

        reclen = nstatavg*SIZEOFINT
        if (irec == 1) then
            write (unit) reclen
            write (unit) statavg
            write (unit) reclen
        else
            write (unit) statavg
        end if

        return
    end subroutine WRT_STHD

#undef LOC_STATUS

    !########################################################################
    !########################################################################
#define LOC_STATUS 'old'

    subroutine IO_READ_AVG_SPATIAL(name, mean_flow, mean_scal)
        use TLab_Constants, only: lfile, wp, wi
        use TLab_Time, only: itime, rtime
        use TLab_Memory, only: jmax, inb_scal
        use TLab_WorkFlow, only: TLab_Write_ASCII

#ifdef USE_MPI
        use TLabMPI_VARS
#endif

        implicit none

        character*(*) name
        real(wp), dimension(nstatavg*jmax*MA_MOMENTUM_SIZE) :: mean_flow
        real(wp), dimension(nstatavg*jmax*MS_SCALAR_SIZE*inb_scal) :: mean_scal

        ! -------------------------------------------------------------------
        character*128 :: line
        logical lfilexist
        integer(wi) nstat
        integer, parameter :: i0 = 0

        ! ###################################################################
#ifdef USE_MPI
        if (ims_pro == 0) then
#endif
            lfilexist = .false.
            inquire (file=name, EXIST=lfilexist)

            ! -------------------------------------------------------------------
            ! Read data
            ! -------------------------------------------------------------------
            if (lfilexist) then
                line = 'Reading field '//trim(adjustl(name))//'...'
                call TLab_Write_ASCII(lfile, line)

#include "dns_open_file.h"
                rewind (LOC_UNIT_ID)
                nstat = MA_MOMENTUM_SIZE + MS_SCALAR_SIZE*inb_scal
                call RD_STHD(LOC_UNIT_ID, i0, &
                             itime, rtime, istattimeorg, rstattimeorg, &
                             nstatavg, jmax, nstat, nstatavg_points, statavg)
                read (LOC_UNIT_ID) mean_flow
                read (LOC_UNIT_ID) mean_scal
                close (LOC_UNIT_ID)

                ! -------------------------------------------------------------------
                ! Initialize data
                ! -------------------------------------------------------------------
            else
                nstatavg_points = 0
                istattimeorg = itime
                rstattimeorg = rtime
                mean_flow = 0.0_wp
                mean_scal = 0.0_wp
                call TLab_Write_ASCII(lfile, 'Statistics have been initialized.')
            end if

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine IO_READ_AVG_SPATIAL

    ! ###################################################################
    ! ###################################################################
    subroutine RD_STHD(unit, irec, iter, rtime, iterorg, rtimeorg, &
                       nstatavg, jmax, nstat, nstatavg_points, statavg)
        use TLab_Constants, only: efile, wp, wi
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop

        implicit none

        integer(wi) unit, irec
        integer(wi) iter, iterorg
        integer(wi) nstatavg, jmax, nstat, nstatavg_points
        real(wp) rtime, rtimeorg
        integer(wi) statavg(nstatavg)

        integer(wi) iterdum
        integer(wi) nstatavgdum, jmaxdum, nstatdum
        real(wp) tmp(1)
        integer(wi) reclen

        if (irec == 1) then
            read (unit) reclen
            read (unit) iterdum
            read (unit) rtime
            read (unit) reclen
        else
            read (unit) iterdum, tmp
            rtime = tmp(1)
        end if

        if (irec == 1) then
            read (unit) reclen
            read (unit) iterorg
            read (unit) rtimeorg
            read (unit) reclen
        else
            read (unit) iterorg, tmp
            rtimeorg = tmp(1)
        end if

        if (irec == 1) then
            read (unit) reclen
            read (unit) nstatavgdum
            read (unit) jmaxdum
            read (unit) nstatdum
            read (unit) nstatavg_points
            read (unit) reclen
        else
            read (unit) nstatavgdum, jmaxdum, nstatdum, nstatavg_points
        end if

        if (irec == 1) then
            read (unit) reclen
            read (unit) statavg
            read (unit) reclen
        else
            read (unit) statavg
        end if

        ! #####################
        ! # Checking
        ! #####################

        if (iterdum /= iter) then
            call TLab_Write_ASCII(efile, 'Stat file error (iter mismatch).')
            call TLab_Stop(DNS_ERROR_STFILE)
        end if

        if (jmaxdum /= jmax) then
            call TLab_Write_ASCII(efile, 'Stat file error (jmax mismatch).')
            call TLab_Stop(DNS_ERROR_STFILE)
        end if

        if (nstatavgdum /= nstatavg) then
            call TLab_Write_ASCII(efile, 'Stat file error (nstatavg mismatch).')
            call TLab_Stop(DNS_ERROR_STFILE)
        end if

        if (nstatdum /= nstat) then
            call TLab_Write_ASCII(efile, 'Stat file error (nstat mismatch).')
            call TLab_Stop(DNS_ERROR_STFILE)
        else
            nstat = nstatdum
        end if

        return
    end subroutine RD_STHD

end module Avg_Spatial
