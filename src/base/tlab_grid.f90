#include "dns_error.h"

module TLab_Grid
    use TLab_Constants, only: efile, wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    type, public :: grid_dt
        sequence
        character*8 name
        integer(wi) size
        ! logical :: uniform = .true.
        logical :: periodic = .false.
        real(wp) scale
        real(wp), allocatable :: nodes(:)
    end type
    type(grid_dt), public :: x, y, z

    public :: TLab_Grid_Read
    public :: TLab_Grid_Write

contains
!########################################################################
!########################################################################
    subroutine TLab_Grid_Read(name, x, y, z, sizes)
        character*(*) name
        type(grid_dt), intent(inout) :: x, y, z
        integer(wi), intent(in), optional :: sizes(3)

        ! -----------------------------------------------------------------------
        character*(32) line

        ! #######################################################################
        write (0, '(*(G0))') '[GRID-1] opening file ', trim(name), ' present(sizes)=', present(sizes)
        open (50, file=name, status='old', form='unformatted')
        rewind (50)
        write (0, '(A)') '[GRID-2] file opened; reading sizes record'

        ! -----------------------------------------------------------------------
        read (50) x%size, y%size, z%size
        write (0, '(*(G0))') '[GRID-3] sizes read: x%size=', x%size, ' y%size=', y%size, ' z%size=', z%size
        if (present(sizes)) write (0, '(*(G0))') '[GRID-3b] expected sizes=', sizes(1), ' x', sizes(2), ' x', sizes(3)

        if (present(sizes)) then        ! check
            if (any([x%size, y%size, z%size] /= sizes)) then
                close (50)
                write (line, 100) x%size, y%size, z%size
                call TLab_Write_ASCII(efile, __FILE__//'. Dimensions ('//trim(line)//') unmatched.')
                call TLab_Stop(DNS_ERROR_DIMGRID)
            end if
        end if

        write (0, '(A)') '[GRID-4] size check passed; reading scale record'
        read (50) x%scale, y%scale, z%scale
        write (0, '(*(G0))') '[GRID-5] scales read: x%scale=', x%scale, ' y%scale=', y%scale, ' z%scale=', z%scale

        if (allocated(x%nodes)) deallocate (x%nodes)
        if (allocated(y%nodes)) deallocate (y%nodes)
        if (allocated(z%nodes)) deallocate (z%nodes)
        write (0, '(A)') '[GRID-6] allocating nodes arrays'
        allocate (x%nodes(x%size), y%nodes(y%size), z%nodes(z%size))
        write (0, '(A)') '[GRID-7] nodes allocated; reading x%nodes'

        read (50) x%nodes(:)
        write (0, '(A)') '[GRID-8] x%nodes read; reading y%nodes'
        read (50) y%nodes(:)
        write (0, '(A)') '[GRID-9] y%nodes read; reading z%nodes'
        read (50) z%nodes(:)
        write (0, '(A)') '[GRID-10] z%nodes read; closing file'

        ! -----------------------------------------------------------------------
        close (50)
        write (0, '(A)') '[GRID-11] file closed; TLab_Grid_Read returning'

        return

100     format(I5, ',', I5, ',', I5)

    end subroutine TLab_Grid_Read

!########################################################################
!########################################################################
    subroutine TLab_Grid_Write(name, x, y, z)
        character*(*) name
        type(grid_dt), intent(in) :: x, y, z

        !########################################################################
        open (unit=51, file=name, form='unformatted', status='unknown')

        write (51) x%size, y%size, z%size
        write (51) x%scale, y%scale, z%scale

        write (51) x%nodes(1:x%size)
        write (51) y%nodes(1:y%size)
        write (51) z%nodes(1:z%size)

        close (51)

        return
    end subroutine TLab_Grid_Write

end module TLab_Grid
