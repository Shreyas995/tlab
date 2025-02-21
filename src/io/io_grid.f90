#include "dns_error.h"

module IO_Grid
    use TLab_Constants, only: efile, wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    public :: IO_READ_GRID
    public :: IO_WRITE_GRID

contains
!########################################################################
!########################################################################
    subroutine IO_READ_GRID(name, imax, jmax, kmax, scalex, scaley, scalez, x, y, z)
        character*(*) name
        integer(wi), intent(in) :: imax, jmax, kmax
        real(wp) scalex, scaley, scalez
        real(wp), allocatable, intent(out) :: x(:), y(:), z(:)

        ! -----------------------------------------------------------------------
        integer(wi) imax_loc, jmax_loc, kmax_loc
        real(wp) scale(3)
        character*(32) line

        ! #######################################################################
        open (50, file=name, status='old', form='unformatted')
        rewind (50)

        ! -----------------------------------------------------------------------
        read (50) imax_loc, jmax_loc, kmax_loc
        read (50) scale
        scalex = scale(1)
        scaley = scale(2)
        scalez = scale(3)

        ! -----------------------------------------------------------------------
        if (imax_loc /= imax .or. jmax_loc /= jmax .or. kmax_loc /= kmax) then
            close (50)
            write (line, 100) imax_loc, jmax_loc, kmax_loc
            call TLab_Write_ASCII(efile, 'IO_READ_GRID. Dimensions ('//trim(line)//') unmatched.')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if

        ! -----------------------------------------------------------------------
        allocate (x(imax), y(jmax), z(kmax))
        read (50) x
        read (50) y
        read (50) z
        close (50)

        return

100     format(I5, ',', I5, ',', I5)

    end subroutine IO_READ_GRID

!########################################################################
!########################################################################
    subroutine IO_WRITE_GRID(name, imax, jmax, kmax, scalex, scaley, scalez, x, y, z)
        character*(*) name
        integer(wi) imax, jmax, kmax
        real(wp) scalex, scaley, scalez
        real(wp), intent(in) :: x(imax), y(jmax), z(kmax)

        ! -----------------------------------------------------------------------
        real(wp) scale(3)

        !########################################################################
        open (unit=51, file=name, form='unformatted', status='unknown')

        ! -----------------------------------------------------------------------
        scale(1) = scalex
        scale(2) = scaley
        scale(3) = scalez
        write (51) imax, jmax, kmax
        write (51) scale

        ! -----------------------------------------------------------------------
        write (51) x
        write (51) y
        write (51) z
        close (51)

        return
    end subroutine IO_WRITE_GRID

end module IO_Grid
