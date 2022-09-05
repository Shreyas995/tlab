#include "types.h"

subroutine RAND_PSD(nx, ny, nz, u)
    use TLAB_VARS, only: isize_txc_dimz
    use TLAB_VARS, only: g
    use RAND_LOCAL
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif

    implicit none

    integer(ci) nx, ny, nz
    real(cp), dimension(isize_txc_dimz, nz), intent(INOUT) :: u

    ! -----------------------------------------------------------------------
    integer(ci) j, k, iglobal, jglobal, kglobal, ip
    real(cp) pow_dst, pow_org, phase
    real(cp) f, f0, f1, fi, fj, fk
    real(cp) RAN0

    ! #######################################################################
    f0 = spc_param(1); f1 = spc_param(4)

    do k = 1, nz
#ifdef USE_MPI
        kglobal = k + ims_offset_k
#else
        kglobal = k
#endif
        if (kglobal <= g(3)%size/2 + 1) then; fk = M_REAL(kglobal - 1)/g(3)%scale
        else; fk = -M_REAL(g(3)%size + 1 - kglobal)/g(3)%scale; end if

        do j = 1, ny
            jglobal = j ! No MPI decomposition along Oy
            if (jglobal <= g(2)%size/2 + 1) then; fj = M_REAL(jglobal - 1)/g(2)%scale
            else; fj = -M_REAL(g(2)%size + 1 - jglobal)/g(2)%scale; end if

            do i = 1, nx/2 + 1
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif
                if (iglobal <= g(1)%size/2 + 1) then; fi = M_REAL(iglobal - 1)/g(1)%scale
                else; fi = -M_REAL(g(1)%size + 1 - iglobal)/g(1)%scale; end if

                f = SQRT(fi**2 + fj**2 + fk**2)

                ! -----------------------------------------------------------------------
                ! target psd
                ! -----------------------------------------------------------------------
                if (ispectrum == 1) then; pow_dst = C_1_R
                elseif (ispectrum == 2) then; pow_dst = (f/f0)**4/(1.+12./5.*(f/f0)**2)**(17./6.)
                elseif (ispectrum == 3) then; pow_dst = f**4*EXP(-C_2_R*(f/f0)**2)
                elseif (ispectrum == 4) then; pow_dst = f**2*EXP(-C_2_R*f/f0)
                elseif (ispectrum == 5) then; pow_dst = f**4/16.*EXP(-C_2_R*(f/f0)**2)
                elseif (ispectrum == 6) then; pow_dst = EXP(-C_05_R*((f - f0)/f1)**2)/(f1*SQRT(C_2_R*C_PI_R))
                end if

                if ((f - spc_param(2))*(spc_param(3) - f) > 0.0_cp) pow_dst = 0.0_cp ! Clip

                if (f == 0.0_cp) then
                    pow_dst = 0.0_cp
                else
                    if (g(2)%size == 1 .or. g(3)%size == 1) then ! 2D spectrum
                        pow_dst = pow_dst/(C_PI_R*f)
                    else
                        pow_dst = pow_dst/(2*C_PI_R*f**2)
                    end if
                end if

                ! -----------------------------------------------------------------------
                ! phase and scaling of complex data
                ! -----------------------------------------------------------------------
                ip = (nx + 2)*(j - 1) + 2*i

                if (ipdf == 0) then
                    if (iglobal == 1 .or. iglobal == g(1)%size/2 + 1) then; phase = 0.0_cp
                    else; phase = (RAN0(seed) - C_05_R)*C_2_R*C_PI_R; end if
                    u(ip - 1, k) = SQRT(pow_dst)*COS(phase)
                    u(ip, k) = SQRT(pow_dst)*SIN(phase)

                else
                    pow_org = u(ip - 1, k)**2 + u(ip, k)**2

                    if (pow_org > 0.0_cp) pow_dst = SQRT(pow_dst/pow_org)

                    u(ip - 1, k) = u(ip - 1, k)*pow_dst
                    u(ip, k) = u(ip, k)*pow_dst

                end if

            end do
        end do
    end do

    return
end subroutine RAND_PSD
