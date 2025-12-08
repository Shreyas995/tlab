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

module Tlab_Debug
    
    use  TLab_Constants, only : wp, wi
#ifdef USE_MPI
    use  TLabMPI_VARS, only : ims_pro
#endif
    implicit none
    save
    public :: TLab_Debug_Initialize, TLab_Debug_Print_1D, TLab_Debug_Print_2D, TLab_Debug_Print_3D, TLab_Debug_Print_4D
    public :: TLab_Debug_Print_2D_c, TLab_Debug_Print_1D_c, TLab_Debug_Print_3D_c
    public :: TLab_Debug_Print_2D_i
    integer(wi), parameter :: FILE_UNIT_BASE = 500
#ifndef USE_MPI
    integer(wi), parameter :: ims_pro = 0
#endif
contains 
    
    subroutine Tlab_Debug_Initialize()
        implicit none
        integer :: unit_num, iostat
        character(len=256) :: filename
        logical :: opened

        ! Calculate thread-specific unit number
        unit_num = FILE_UNIT_BASE + ims_pro

        write(filename, '(a, i0, a)') 'debug_thread_n', ims_pro, '.log'

        inquire(unit=unit_num, opened=opened)

        if (.not. opened) then
            open(unit=unit_num, file=trim(filename), &
                position='append', status='unknown')

            if (iostat /= 0) then
                write(*,*) 'ERROR: Could not open debug file for PE', ims_pro, &
                           ' (Unit:', unit_num, ', File:', trim(filename), ')'
            end if
        end if

    end subroutine Tlab_Debug_Initialize

    subroutine TLab_Debug_Print_1D(msg, var, msg2)
        implicit none
        character(len=*), intent(in) :: msg
        real(wp), intent(in) :: var(:)
        character(len=*), intent(in), optional :: msg2
        integer(wi) :: i
        integer :: unit_num
        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', sum(var)
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', sum(var)
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_1D

    subroutine TLab_Debug_Print_2D(msg, var, msg2)
        character(len=*), intent(in) :: msg
        real(wp), intent(in) :: var(:,:)
        character(len=*), intent(in), optional :: msg2
        integer(wp) :: i
        integer :: unit_num
        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', sum(var)
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', sum(var)
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_2D

    subroutine TLab_Debug_Print_3D(msg, var, msg2)
        character(len=*), intent(in) :: msg
        real(wp), intent(in) :: var(:,:,:)
        character(len=*), intent(in), optional :: msg2
        integer(wp) :: i
        integer :: unit_num
        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', sum(var)
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', sum(var)
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_3D

    subroutine TLab_Debug_Print_2D_c(msg, var, msg2)
        character(len=*), intent(in) :: msg
        complex(wp), intent(in) :: var(:,:)
        character(len=*), intent(in), optional :: msg2

        integer(wp) :: i
        integer :: unit_num

        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', sum(var)
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', sum(var)
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_2D_c

    subroutine TLab_Debug_Print_1D_c(msg, var, msg2)
        implicit none
        character(len=*), intent(in) :: msg
        complex(wp), intent(in) :: var(:)
        character(len=*), intent(in), optional :: msg2

        integer(wi) :: i
        integer :: unit_num

        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', sum(var)
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', sum(var)
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_1D_c

     subroutine TLab_Debug_Print_3D_c(msg, var, msg2)
        character(len=*), intent(in) :: msg
        complex(wp), intent(in) :: var(:,:,:)
        character(len=*), intent(in), optional :: msg2

        integer(wp) :: i
        integer :: unit_num

        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', sum(var)
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', sum(var)
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_3D_c

    subroutine TLab_Debug_Print_4D(msg, var, msg2)
        character(len=*), intent(in) :: msg
        real(wp), intent(in) :: var(:,:,:,:)
        character(len=*), intent(in), optional :: msg2
        integer(wp) :: i
        integer :: unit_num
        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', sum(var)
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', sum(var)
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_4D

    subroutine TLab_Debug_Print_1D_i(msg, var, msg2)
        character(len=*), intent(in) :: msg
        integer(wi), intent(in) :: var(:)
        character(len=*), intent(in), optional :: msg2
        integer(wp) :: i
        integer :: unit_num
        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', sum(var)
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', sum(var)
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_1D_i

    subroutine TLab_Debug_Print_2D_i(msg, var, msg2)
        character(len=*), intent(in) :: msg
        integer(wi), intent(in) :: var(:,:)
        character(len=*), intent(in), optional :: msg2
        integer(wp) :: i
        integer :: unit_num
        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', sum(var)
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', sum(var)
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_2D_i

    subroutine TLab_Debug_Print_int(msg, var, msg2)
        character(len=*), intent(in) :: msg
        integer(wi), intent(in) :: var
        character(len=*), intent(in), optional :: msg2
        integer(wp) :: i
        integer :: unit_num
        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', var
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', var
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_int

    subroutine TLab_Debug_Print_real(msg, var, msg2)
        character(len=*), intent(in) :: msg
        real(wi), intent(in) :: var
        character(len=*), intent(in), optional :: msg2
        integer(wp) :: i
        integer :: unit_num
        unit_num = FILE_UNIT_BASE + ims_pro
        if (present(msg2)) then
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', trim(msg2), ' ', var
        else
            write(unit_num, *) trim(msg), 'DEBUG (PE', ims_pro, '): ', var
        end if
        flush(unit_num)
    end subroutine TLab_Debug_Print_real

end module Tlab_Debug
!########################################################################