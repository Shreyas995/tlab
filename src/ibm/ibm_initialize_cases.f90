#include "dns_error.h"
!########################################################################
!# HISTORY / AUTHORS
!#
!# 2023/10/26 - Shreyas Deshpande
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   Initiliza IBM case for all solids
!#
!#
!########################################################################
!# ARGUMENTS
!#
!#
!########################################################################
!# REQUIREMENTS
!#
!#
!########################################################################

subroutine IBM_INITIALIZE_CASES(g, isize_nob, isize_nob_be, nob, nob_b, nob_e, ibm_case)
    use TLab_Constants, only: efile, wp, wi
    use IBM_VARS, only: nflu
    use FDM, only: fdm_dt
    use TLab_WorkFlow,     only: TLab_Write_ASCII, TLab_Stop
    implicit none

    type(fdm_dt), intent(in) :: g
    integer(wi), intent(in) :: isize_nob, isize_nob_be
    integer(wi), dimension(isize_nob), intent(in) :: nob
    integer(wi), dimension(isize_nob_be), intent(in) :: nob_b, nob_e
    integer(wi), dimension(isize_nob_be), intent(out) :: ibm_case

    integer(wi) :: ii, ip, iob, nlines

    ! ================================================================== !
    ! initialize ibm_case
    ibm_case(:) = 0

    ! cf. ibm_allocate.f90
    nlines = isize_nob

    ! index convention on contiguous lines
    ! ||...-ip_fl-x-(fluid points)-x-ip_il||---(solid points)---||ip_ir-x-(fluid points)-x-ip_fr-...||

    ! index ii (dummy index; for x,y,z: ii == jk,ik,ij)
    do ii = 1, nlines        ! index of ii-plane, loop over plane and check for objects in each line
        if (nob(ii) /= 0) then ! if line contains immersed object(s) --yes-->  spline interpolation
            ip = 0
            do iob = 1, nob(ii)    ! loop over immersed object(s)
                ! select different cases of immersed objects
                if (nob_b(ip + ii) == 1) then
                    ! ================================================================== !
                    if (nob_e(ip + ii) == g%size) then
                        ! 1. case: object over full extend of line
                        ibm_case((iob - 1)*nlines + ii) = 1
                        ! .............................................................. !

                    else if ((nob_e(ip + ii) <= (g%size - nflu)) .and. (g%periodic .eqv. .true.)) then
                        ! 2. case: object is semi-immersed - periodic case
                        ibm_case((iob - 1)*nlines + ii) = 2
                        ! .........................SANITY CHECK......................... !
                        call GEOMETRY_CHK(g, nob_e, nob_b, isize_nob_be, nlines, iob, nob(ii), ibm_case((iob - 1)*nlines + ii), ii, ip)
                        ! .............................................................. !

                    else if ((nob_e(ip + ii) <= (g%size - nflu)) .and. (g%periodic .eqv. .false.) &
                             .and. ((nob_e(ip + ii) /= 1 + nob_b(ip + ii)) .and. (nob_e(ip + ii) /= nob_b(ip + ii)))) then ! e.g. in vertical direction
                        ! 3. case: object is semi-immersed - non-periodic case - lower boundary
                        ibm_case((iob - 1)*nlines + ii) = 3
                        ! .........................SANITY CHECK......................... !
                        call GEOMETRY_CHK(g, nob_e, nob_b, isize_nob_be, nlines, iob, nob(ii), ibm_case((iob - 1)*nlines + ii), ii, ip)
                        ! .............................................................. !

                    else if ((nob_e(ip + ii) == 1 + nob_b(ip + ii)) .and. (g%periodic .eqv. .false.)) then
                        ! 8. case: Two solid point at lower boundary
                        ibm_case((iob - 1)*nlines + ii) = 8
                        ! .............................................................. !

                    else if ((nob_e(ip + ii) == nob_b(ip + ii)) .and. (g%periodic .eqv. .false.)) then
                        ! 9. case: single solid point at lower boundary
                        ibm_case((iob - 1)*nlines + ii) = 9
                        ! .............................................................. !

                    else
                        call TLab_Write_ASCII(efile, 'IBM_INITIALIZE. Not enough fluid points right of the right interface.')
                        call TLab_Stop(DNS_ERROR_IBM_SPLINE)
                    end if
                    ! ================================================================== !
                else if (nob_b(ip + ii) >= (nflu + 1)) then
                    if ((nob_e(ip + ii) <= (g%size - nflu)) .and. (nob_b(ip + ii) < nob_e(ip + ii))) then
                        ! 4. case: object is fully immersed
                        ibm_case((iob - 1)*nlines + ii) = 4
                        ! .........................SANITY CHECK......................... !
                        call GEOMETRY_CHK(g, nob_e, nob_b, isize_nob_be, nlines, iob, nob(ii), ibm_case((iob - 1)*nlines + ii), ii, ip)
                        ! .............................................................. !

                    else if ((nob_e(ip + ii) == g%size) .and. (g%periodic .eqv. .true.)) then
                        ! 5. case: object is semi-immersed - periodic case
                        ibm_case((iob - 1)*nlines + ii) = 5
                        ! .........................SANITY CHECK......................... !
                        call GEOMETRY_CHK(g, nob_e, nob_b, isize_nob_be, nlines, iob, nob(ii), ibm_case((iob - 1)*nlines + ii), ii, ip)
                        ! .............................................................. !

                    else if ((nob_e(ip + ii) == g%size) .and. (g%periodic .eqv. .false.)) then  ! e.g. in vertical direction
                        ! 6. case: object is semi-immersed - non-periodic case - upper boundary
                        ibm_case((iob - 1)*nlines + ii) = 6
                        ! .........................SANITY CHECK......................... !
                        call GEOMETRY_CHK(g, nob_e, nob_b, isize_nob_be, nlines, iob, nob(ii), ibm_case((iob - 1)*nlines + ii), ii, ip)
                        ! .............................................................. !

                    else if (((nob_e(ip + ii)) < (nob_b(ip + ii))) .and. (g%periodic .eqv. .true.)) then
                        ! 7. case: object is semi-immersed at both the domain boundaries - periodic case
                        ibm_case((iob - 1)*nlines + ii) = 7
                        ! .............................................................. !

                    else
                        call TLab_Write_ASCII(efile, 'IBM_INITIALIZE. Case not implemented.')
                        call TLab_Stop(DNS_ERROR_IBM_INITIALIZE)
                    end if
                    ! .............................................................. !

                else
                    call TLab_Write_ASCII(efile, 'IBM_INITIALIZE. Invalid case number check ibm_case.')
                    call TLab_Stop(DNS_ERROR_IBM_INITIALIZE)
                end if
                if ((ibm_case((iob - 1)*nlines + ii) < 1) .or. (ibm_case((iob - 1)*nlines + ii) > 9)) then
                    call TLab_Write_ASCII(efile, 'IBM_INITIALIZE. Case number not between 2 to 8. Array overwritten')
                    call TLab_Stop(DNS_ERROR_IBM_INITIALIZE)
                end if
                if (ibm_case((iob - 1)*nlines + ii) == 0) then
                    call TLab_Write_ASCII(efile, 'The object should not exist. Incorrect assignment of IBM_Cases')
                    call TLab_Stop(DNS_ERROR_IBM_INITIALIZE)
                end if
                ip = ip + nlines
            end do
        end if
    end do

    return
end subroutine IBM_INITIALIZE_CASES

!########################################################################

subroutine GEOMETRY_CHK(g, nob_e, nob_b, isize_nob_be, nlines, iob, nob, ibm_case, ii, ip)

    use FDM, only: fdm_dt
    use TLab_Constants, only: efile, wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop

    implicit none

    type(fdm_dt), intent(in) :: g
    integer(wi), dimension(isize_nob_be), intent(in) :: nob_b, nob_e
    integer(wi), intent(in) :: ibm_case
    integer(wi), intent(in) :: nlines, nob, iob, ii, ip, isize_nob_be

    character(len=128) :: line

    ! ================================================================== !
    ! Check length of the gap vector before code execution
    if (ibm_case < 8) then
        if (abs(nob_e(ip + ii) - nob_b(ip + ii) + 1) < 3) then
            write (line, *) ibm_case
            line = 'IBM_INITIALIZE_CASES: case = '//trim(adjustl(line))//'. Less than 3 solid points. Check geometry.'
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_IBM_INITIALIZE)
        end if
        if (abs(nob_e(ip + ii) - nob_b(ip + ii)) > (g%size - 3)) then
            write (line, *) ibm_case
            line = 'IBM_INITIALIZE_CASES: case = '//trim(adjustl(line))//'. There must be at least 3 fluid points. Check geometry.'
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_IBM_INITIALIZE)
        end if
        if ((iob + 1) <= nob) then
            if (nob_b(ip + nlines + ii) - (nob_e(ip + ii)) < 3) then
                write (line, *) ibm_case
                line = 'IBM_INITIALIZE_CASES: case = '//trim(adjustl(line))//'. Insufficient gap between 2 solid objects. Check geometry.'
                call TLab_Write_ASCII(efile, line)
                call TLab_Stop(DNS_ERROR_IBM_INITIALIZE)
            end if
        end if
    else if (ibm_case == 8) then
        if ((iob + 1) <= nob) then
            if ((nob_e(ip + ii) - nob_b(ip + nlines + ii)) < 3) then
                write (line, *) ibm_case
                line = 'IBM_INITIALIZE_CASES: case = '//trim(adjustl(line))//'. Insufficient gap bewtween 2 solid objects. Check geometry.'
                call TLab_Write_ASCII(efile, line)
                call TLab_Stop(DNS_ERROR_IBM_INITIALIZE)
            end if
        end if
    end if
end subroutine GEOMETRY_CHK
!########################################################################
