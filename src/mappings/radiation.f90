#include "dns_const.h"
#include "dns_error.h"

module Radiation
    use TLAB_CONSTANTS, only: wp, wi, pi_wp, BCS_MAX, BCS_MIN, efile, MAX_PROF
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_VARS, only: imode_eqns, inb_scal_array
    ! use TLAB_VARS, only: infrared
    use TLAB_ARRAYS, only: wrk3d
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    use THERMO_VARS, only: imixture
    use OPR_PARTIAL, only: OPR_PARTIAL_Y
    use OPR_ODES
    implicit none

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_LW_BULK1D_LIQUID = 1
    integer, parameter :: TYPE_LW_BULK1D = 2

    real(wp) :: sigma = 5.67037442e-8_wp ! W /m^2 /K

    public :: Radiation_Initialize
    public :: Radiation_Infrared

contains
    ! ###################################################################
    ! ###################################################################
    subroutine RADIATION_READBLOCK(bakfile, inifile, block, var)
        character(len=*), intent(in) :: bakfile, inifile, block
        type(term_dt), intent(out) :: var

        character(len=512) sRes
        integer(wi) idummy

        ! -------------------------------------------------------------------
        call TLAB_WRITE_ASCII(bakfile, '#')
        call TLAB_WRITE_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLAB_WRITE_ASCII(bakfile, '#Type=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

        if (var%type == EQNS_NONE) then        ! backwards compatibility, to be removed
            call SCANINICHAR(bakfile, inifile, block, 'Type', 'None', sRes)
            if (trim(adjustl(sRes)) == 'none') then; var%type = TYPE_NONE
            else if (trim(adjustl(sRes)) == 'irbulk1dliquid') then; var%type = TYPE_LW_BULK1D_LIQUID
            else if (trim(adjustl(sRes)) == 'irbulk1d') then; var%type = TYPE_LW_BULK1D
            else if (trim(adjustl(sRes)) == 'bulk1dlocal') then; var%type = EQNS_RAD_BULK1D_LOCAL
            else if (trim(adjustl(sRes)) == 'bulk1dglobal') then; var%type = EQNS_RAD_BULK1D_GLOBAL
            else
                call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Radiation option.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            end if
        end if

        var%active = .false.
        if (var%type /= EQNS_NONE) then
            call SCANINIINT(bakfile, inifile, block, 'Scalar', '1', idummy)
            var%active(idummy) = .true.

            var%parameters(:) = 0.0_wp
            call SCANINICHAR(bakfile, inifile, block, 'Parameters', '1.0', sRes)
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, var%parameters)

        end if

        return
    end subroutine RADIATION_READBLOCK

!########################################################################
!########################################################################
    subroutine Radiation_Initialize()
        use TLAB_VARS, only: infrared

        ! -------------------------------------------------------------------
        ! By default, transport and radiation are caused by last scalar
        infrared%scalar = inb_scal_array

        if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            if (infrared%type /= EQNS_NONE) then
                infrared%active(inb_scal_array) = .true. ! liquid
                infrared%active(inb_scal_array + 1) = .true. ! buoyancy
            end if

        end if

        ! backwards compatibility
        if (any([EQNS_RAD_BULK1D_LOCAL, EQNS_RAD_BULK1D_GLOBAL] == infrared%type)) then
            infrared%parameters(1) = infrared%parameters(1)*infrared%parameters(2)
            infrared%parameters(3) = infrared%parameters(3)*infrared%parameters(2)
            infrared%parameters(2) = 1.0_wp/infrared%parameters(2)
        end if

        ! -------------------------------------------------------------------
        ! in case nondimensional we need to adjust sigma

        ! ! testing
        ! infrared%type = TYPE_LW_BULK1D
        ! infrared%parameters(4) = infrared%parameters(3)
        ! infrared%parameters(3) = 0.0_wp
        ! sigma = 0.0_wp
        !
        return
    end subroutine Radiation_Initialize

!########################################################################
!########################################################################
    subroutine IR_Bulk1D_Liquid(infrared, nx, ny, nz, g, a_source, flux)
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(inout) :: a_source(nx*nz, ny)      ! input as bulk absorption coefficent, output as source
        real(wp), intent(out), optional :: flux(nx*nz, ny)

        target a_source, flux

! -----------------------------------------------------------------------
        integer(wi) j, nxy, nxz
        real(wp) fd, fu
        real(wp), pointer :: p_a(:, :) => null()
        real(wp), pointer :: p_tau(:, :) => null()
        real(wp), pointer :: p_source(:, :) => null()
        real(wp), pointer :: p_flux(:, :) => null()

! #######################################################################
        nxy = nx*ny     ! For transposition to make y direction the last one
        nxz = nx*nz

        if (nz == 1) then
            p_a => a_source
            p_source => a_source
            p_tau(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            if (present(flux)) p_flux => flux
        else
            p_a(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            p_source(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            if (present(flux)) then
                p_tau => flux
                p_flux(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
            else
                p_tau => a_source
            end if

#ifdef USE_ESSL
            call DGETMO(a_source, nxy, nxy, nz, p_a, nz)
#else
            call DNS_TRANSPOSE(a_source, nxy, nz, nxy, p_a, nz)
#endif
        end if

! ###################################################################
        ! Calculate (negative) optical thickness; integrating from the top, to be checked
        p_tau(:, ny) = 0.0_wp   ! boundary condition
        call OPR_Integral1(nxz, g, p_a, p_tau, BCS_MAX)
        ! p_tau(:, 1) = 0.0_wp     ! boundary condition
        ! call OPR_Integral1(nxz, g, p_org, p_tau, BCS_MIN)

! Calculate exp(-tau(z, zmax)/\mu)
        do j = ny, 1, -1
            p_tau(:, j) = exp(p_tau(:, j))
        end do
        !  p_tau = dexp(p_tau)         seg-fault; need ulimit -u unlimited

! ###################################################################
! Calculate heating rate
        fd = infrared%parameters(1)
        fu = infrared%parameters(3)
        if (abs(infrared%parameters(3)) > 0.0_wp) then
            do j = ny, 1, -1
                p_source(:, j) = p_a(:, j)*(p_tau(:, j)*fd &                       ! downward flux
                                            + p_tau(:, 1)/p_tau(:, j)*fu)       ! upward flux
            end do
        else
            p_source = p_a*p_tau*fd
        end if

        if (nz > 1) then ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
            call DGETMO(p_source, nz, nz, nxy, a_source, nxy)
#else
            call DNS_TRANSPOSE(p_source, nz, nxy, nz, a_source, nxy)
#endif
        end if

! ###################################################################
! Calculate radiative flux, if necessary
        if (present(flux)) then
            fd = -infrared%parameters(1)
            fu = infrared%parameters(3)
            if (abs(infrared%parameters(3)) > 0.0_wp) then
                do j = ny, 1, -1
                    p_flux(:, j) = fd*p_tau(:, j) &                       ! downward flux
                                   + fu*p_tau(:, 1)/p_tau(:, j)       ! upward flux
                end do
            else
                p_flux = p_tau*fd
            end if

            if (nz > 1) then ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
                call DGETMO(p_flux, nz, nz, nxy, flux, nxy)
#else
                call DNS_TRANSPOSE(p_flux, nz, nxy, nz, flux, nxy)
#endif
            end if

        end if

! ###################################################################
        nullify (p_a, p_tau, p_source, p_flux)

        return
    end subroutine IR_Bulk1D_Liquid

!########################################################################
!########################################################################
    subroutine IR_Bulk1D(infrared, nx, ny, nz, g, a, b, source, tmp1, flux)
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(inout) :: a(nx*nz, ny)                ! bulk absorption coefficent
        real(wp), intent(inout) :: b(nx*nz, ny)                ! emission function
        real(wp), intent(out) :: source(nx*nz, ny)
        real(wp), intent(inout) :: tmp1(nx*nz, ny)                ! emission function
        real(wp), intent(out), optional :: flux(nx*nz, ny)

        target a, b, source, flux, tmp1

! -----------------------------------------------------------------------
        integer(wi) j, nxy, nxz
        real(wp) fd, fu, mu, dummy
        real(wp), pointer :: p_a(:, :) => null()
        real(wp), pointer :: p_tau(:, :) => null()
        real(wp), pointer :: p_ab(:, :) => null()
        real(wp), pointer :: p_source(:, :) => null()
        real(wp), pointer :: p_flux(:, :) => null()
        real(wp), pointer :: p_wrk(:, :) => null()

! #######################################################################
        nxy = nx*ny     ! For transposition to make y direction the last one
        nxz = nx*nz

        p_a => source
        p_ab => a
        p_tau(1:nx*nz, 1:ny) => b
        p_wrk(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
        p_source(1:nx*nz, 1:ny) => wrk3d(1:nx*ny*nz)
        p_flux(1:nx*nz, 1:ny) => tmp1

#ifdef USE_ESSL
        call DGETMO(a, nxy, nxy, nz, p_a, nz)
        call DGETMO(b, nxy, nxy, nz, p_ab, nz)
#else
        call DNS_TRANSPOSE(a, nxy, nz, nxy, p_a, nz)
        call DNS_TRANSPOSE(b, nxy, nz, nxy, p_ab, nz)
#endif

! ###################################################################
        mu = 1.0_wp     ! spherical angle correction

        dummy = 1.0_wp/mu
        p_a = p_a*dummy
        ! Calculate (negative) optical thickness; integrating from the top, to be checked
        p_tau(:, ny) = 0.0_wp   ! boundary condition
        call OPR_Integral1(nxz, g, p_a, p_tau, BCS_MAX)
        ! p_tau(:, 1) = 0.0_wp     ! boundary condition
        ! call OPR_Integral1(nxz, g, p_org, p_tau, BCS_MIN)

        ! Calculate exp(-tau(z, zmax)/\mu)
        do j = ny, 1, -1
            p_tau(:, j) = exp(p_tau(:, j))
        end do
        !  p_tau = dexp(p_tau)         seg-fault; need ulimit -u unlimited

        ! product of absorption coefficient times emission function
        dummy = 2.0_wp*pi_wp*mu
        p_ab = p_a*p_ab*dummy

! ###################################################################
        ! downward flux
        fd = -infrared%parameters(1)
        do j = ny, 1, -1
            p_wrk(:, j) = p_ab(:, j)*p_tau(:, 1)/p_tau(:, j)
        end do
        p_flux(:, ny) = 0.0_wp                                  ! boundary condition
        call OPR_Integral1(nxz, g, p_wrk, p_flux, BCS_MAX)     ! recall this gives the negative of the integral
        p_flux = fd*p_tau - p_flux

        ! calculate upward flux at the surface
        fu = infrared%parameters(4)                             ! to be updated

        ! upward flux
        if (present(flux)) then
            p_wrk = p_ab*p_tau
            flux(:, ny) = 0.0_wp                                  ! boundary condition
            call OPR_Integral1(nxz, g, p_wrk, flux, BCS_MAX)     ! recall this gives the negative of the integral
            do j = ny, 1, -1
                p_flux(:, j) = fu*p_tau(:, 1)/p_tau(:, j) + flux(:, j) - flux(:, 1) &
                               + p_flux(:, j)
            end do

#ifdef USE_ESSL
            call DGETMO(p_flux, nz, nz, nxy, flux, nxy)
#else
            call DNS_TRANSPOSE(p_flux, nz, nxy, nz, flux, nxy)
#endif
        end if

! ###################################################################
! Calculate heating rate
        fd = infrared%parameters(1)
        fu = infrared%parameters(4)
        do j = ny, 1, -1
            p_source(:, j) = (fd*p_a(:, j) - p_ab(:, j))*p_tau(:, j) &
                             + (fu*p_a(:, j) - p_ab(:, j))*p_tau(:, 1)/p_tau(:, j)
        end do

#ifdef USE_ESSL
        call DGETMO(p_source, nz, nz, nxy, source, nxy)
#else
        call DNS_TRANSPOSE(p_source, nz, nxy, nz, source, nxy)
#endif

! ###################################################################
        nullify (p_a, p_tau, p_source, p_flux)

        return
    end subroutine IR_Bulk1D

!########################################################################
!########################################################################
    subroutine Radiation_Infrared(infrared, nx, ny, nz, g, s, source, a, b, tmp1, flux)
        use THERMO_ANELASTIC
        type(term_dt), intent(in) :: infrared
        integer(wi), intent(in) :: nx, ny, nz
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: s(nx*ny*nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny*nz)
        real(wp), intent(inout) :: a(nx*ny*nz)      ! space for bulk absorption coefficient
        real(wp), intent(inout) :: b(nx*ny*nz)      ! space for emission function
        real(wp), intent(inout) :: tmp1(nx*ny*nz)
        real(wp), intent(out), optional :: flux(nx*ny*nz)

        ! -----------------------------------------------------------------------
        real(wp) kappa, kappal, kappav

        !########################################################################
        select case (infrared%type)
        case (TYPE_LW_BULK1D_LIQUID, EQNS_RAD_BULK1D_LOCAL)
            ! bulk absorption coefficient; we save it in array source to save memory
            kappa = infrared%parameters(2)
            source = kappa*s(:, infrared%scalar(1))
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, rbackground, source)
            end if

            if (present(flux)) then
                call IR_Bulk1D_Liquid(infrared, nx, ny, nz, g, source, flux)
            else
                call IR_Bulk1D_Liquid(infrared, nx, ny, nz, g, source)
            end if

        case (TYPE_LW_BULK1D)
            ! bulk absorption coefficients and emission function
            kappal = infrared%parameters(2)
            kappav = infrared%parameters(3)
            a = kappal*s(:, infrared%scalar(1)) + kappav*(s(:, 2) - s(:, infrared%scalar(1)))
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_INPLACE(nx, ny, nz, rbackground, a)

                call THERMO_ANELASTIC_TEMPERATURE(nx, ny, nz, s, wrk3d)
            else
                ! calculate temperature
            end if
            b = sigma*wrk3d**4./pi_wp

            if (present(flux)) then
                call IR_Bulk1D(infrared, nx, ny, nz, g, a, b, source, tmp1, flux)
            else
                call IR_Bulk1D(infrared, nx, ny, nz, g, a, b, source, tmp1)
            end if

        end select

    end subroutine Radiation_Infrared

end module
