#include "dns_const.h"

!########################################################################
!#
!# Implementation of skewsymmetric formulation and viscous explicit
!# in one routine, to avoid duplication of derivatives. Internal energy
!# formulation only.
!# 30 first derivative operations and 9 second derivative operations.
!# Viscosity needs to be homogeneous, because of the explicit treatment
!# of diffusion terms.
!#
!########################################################################
subroutine RHS_FLOW_GLOBAL_2()
    use TLab_Constants, only: efile, wp, wi
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
    use TLab_WorkFlow, only: TLab_Write_ASCII
#endif
    use TLab_Memory, only: imax, jmax, kmax, inb_scal
    use FDM, only: g
    use NavierStokes, only: nse_diffusion, EQNS_NONE
    use NavierStokes, only: visc, prandtl
    use Gravity, only: buoyancy
    use TLab_Arrays, only: s
    use TLab_Pointers
    use Thermodynamics, only: CRATIO_INV
    use THERMO_CALORIC
    use DNS_ARRAYS
    use BOUNDARY_BCS
    use OPR_Partial

#ifdef USE_OPENMP
    use OMP_LIB
#endif

    implicit none

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1), i, is
    real(wp) g1, g2, g3, cond, dummy, c13, dum1, dum2, dum3

! ###################################################################
#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'ENTERING RHS_FLOW_GLOBAL_2')
#endif

    bcs = 0

    g1 = buoyancy%vector(1)
    g2 = buoyancy%vector(2)
    g3 = buoyancy%vector(3)
    c13 = 1.0_wp/3.0_wp

! ###################################################################
! Terms \rho u in mass, u-momentum equations
! ###################################################################
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared(imax,jmax,kmax,tmp4,tmp3,tmp2,tmp1,u,v,w,p,rho)
#endif
    do i = 1, imax*jmax*kmax
        tmp4(i) = 0.5_wp*rho(i)*u(i)
        tmp3(i) = tmp4(i)*w(i)
        tmp2(i) = tmp4(i)*v(i)
        tmp1(i) = tmp4(i)*u(i) + p(i)
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif
! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp6)

#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,hq,tmp6 )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 5) = hq(i, 5) - 2.0_wp*tmp6(i)
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif
! -------------------------------------------------------------------
! u-momentum
! -------------------------------------------------------------------
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,hq,tmp2,tmp3,tmp4,rho,g1 )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 1) = hq(i, 1) - (tmp2(i) + tmp3(i) + tmp4(i)) + g1*rho(i)
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

! ###################################################################
! Terms \rho v in mass, v-momentum equations
! ###################################################################
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp1,tmp2,tmp3,tmp4,rho,u,v,w,p )
#endif
    do i = 1, imax*jmax*kmax
        tmp4(i) = 0.5_wp*rho(i)*v(i)
        tmp3(i) = tmp4(i)*w(i)
        tmp2(i) = tmp4(i)*v(i) + p(i)
        tmp1(i) = tmp4(i)*u(i)
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp5)
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp5,tmp6,hq )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 5) = hq(i, 5) - 2.0_wp*tmp5(i)
        tmp6(i) = tmp6(i) + tmp5(i)
    end do
#ifdef USE_APU
    !$omp end parallel do
#endif

! -------------------------------------------------------------------
! v-momentum
! -------------------------------------------------------------------
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp2,tmp3,tmp4,g2,rho,hq )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 2) = hq(i, 2) - (tmp2(i) + tmp3(i) + tmp4(i)) + g2*rho(i)
    end do
#ifdef USE_APU
    !$omp end parallel do
#endif
! ###################################################################
! Terms \rho w in mass, w-momentum equations
! ###################################################################
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp1,tmp2,tmp3,tmp4,u,v,w,p,rho )
#endif
    do i = 1, imax*jmax*kmax
        tmp4(i) = 0.5_wp*rho(i)*w(i)
        tmp3(i) = tmp4(i)*w(i) + p(i)
        tmp2(i) = tmp4(i)*v(i)
        tmp1(i) = tmp4(i)*u(i)
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif
! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp5)
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp5,tmp6,hq )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 5) = hq(i, 5) - 2.0_wp*tmp5(i)
        tmp6(i) = tmp6(i) + tmp5(i)
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

! -------------------------------------------------------------------
! w-momentum
! -------------------------------------------------------------------
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp2,tmp3,tmp4,g3,rho,hq )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 3) = hq(i, 3) - (tmp2(i) + tmp3(i) + tmp4(i)) + g3*rho(i)
    end do
#ifdef USE_APU
    !$omp end parallel do
#endif
! ###################################################################
! Term \rho e in energy equation
! ###################################################################
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i,dummy ) &
    !$omp shared( imax,jmax,kmax,tmp1,tmp2,tmp3,u,v,w,e,rho )
#endif

    do i = 1, imax*jmax*kmax
        dummy = 0.5_wp*rho(i)*e(i)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)

#ifdef USE_APU
    !$omp target teams distribute parallel do &
    !$omp private( i ) &
    !$omp shared( imax,jmax,kmax,hq,tmp2,tmp3,2mp4)
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 4) = hq(i, 4) - (tmp2(i) + tmp3(i) + tmp4(i))
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! Terms u_i d(\rho u_k)/dx_k
! ###################################################################
#ifdef USE_APU
    !$omp target teams distribute parallel do &
    !$omp private( i ) &
    !$omp shared( imax,jmax,kmax,hq,u,v,w,e,tmp6)
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 1) = hq(i, 1) - u(i)*tmp6(i)
        hq(i, 2) = hq(i, 2) - v(i)*tmp6(i)
        hq(i, 3) = hq(i, 3) - w(i)*tmp6(i)
        hq(i, 4) = hq(i, 4) - e(i)*tmp6(i)
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

#ifdef USE_APU
    !$omp target teams distribute parallel do collaspe (2) &
    !$omp private( i,is ) &
    !$omp shared( imax,jmax,kmax,inb_scal,hs,s,tmp6 )
#endif
    do is = 1, inb_scal
        do i = 1, imax*jmax*kmax
            hs(i, is) = hs(i, is) - s(i, is)*tmp6(i)
        end do
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! Terms \rho u_k du_i/dx_k
! First derivatives in internal energy equation from dissipation
! Second derivative terms in the momentum equation
! Note that second derivative routines give also first derivatives !
! ###################################################################
! -------------------------------------------------------------------
! Cross derivatives
! -------------------------------------------------------------------
! energy equation
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), e, tmp4)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), e, tmp3)
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), e, tmp2)
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp2,tmp3,tmp4,rho,hq,u,v,w )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 4) = hq(i, 4) - 0.5_wp*rho(i)*(u(i)*tmp2(i) + v(i)*tmp3(i) + w(i)*tmp4(i))
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

! momentum equations
    dummy = CRATIO_INV*visc

    call OPR_Partial_Y(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), u, tmp5, tmp1)
    call OPR_Partial_Z(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), u, tmp6, tmp2)
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp1,tmp2,tmp5,tmp6,hq,rho,visc,v,w )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 1) = hq(i, 1) + visc*(tmp5(i) + tmp6(i)) - 0.5_wp*rho(i)*(v(i)*tmp1(i) + w(i)*tmp2(i))
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

    call OPR_Partial_X(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), v, tmp5, tmp3)
    call OPR_Partial_Z(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), v, tmp6, tmp4)

#ifdef USE_APU
    !$omp target teams distribute parallel do private( i, dum1 ) &
    !$omp shared( imax,jmax,kmax,tmp1,tmp3,tmp4,tmp5,tmp6,hq,rho,visc,u,w,dummy )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 2) = hq(i, 2) + visc*(tmp5(i) + tmp6(i)) - 0.5_wp*rho(i)*(u(i)*tmp3(i) + w(i)*tmp4(i))
        dum1 = tmp1(i) + tmp3(i)
        hq(i, 4) = hq(i, 4) + dummy*dum1*dum1
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

! arrays tmp1 and tmp3 can be reused
    call OPR_Partial_X(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), w, tmp5, tmp1)
    call OPR_Partial_Y(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), w, tmp6, tmp3)

#ifdef USE_APU
    !$omp target teams distribute parallel do private( i,dum2,dum3 ) &
    !$omp shared( imax,jmax,kmax,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,hq,rho,visc,u,v,dummy )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 3) = hq(i, 3) + visc*(tmp5(i) + tmp6(i)) - 0.5_wp*rho(i)*(u(i)*tmp1(i) + v(i)*tmp3(i))
        dum2 = tmp4(i) + tmp3(i)
        dum3 = tmp1(i) + tmp2(i)
        hq(i, 4) = hq(i, 4) + dummy*(dum2*dum2 + dum3*dum3)
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

! -------------------------------------------------------------------
! Dilatation part
! -------------------------------------------------------------------
    call OPR_Partial_Z(OPR_P2_P1, imax, jmax, kmax, bcs_inf(:, :, 3), g(3), w, tmp6, tmp3)
    call OPR_Partial_Y(OPR_P2_P1, imax, jmax, kmax, bcs_inf(:, :, 2), g(2), v, tmp5, tmp2)
    call OPR_Partial_X(OPR_P2_P1, imax, jmax, kmax, bcs_inf(:, :, 1), g(1), u, tmp4, tmp1)
    dummy = 2.0_wp*visc

#ifdef USE_APU
    !$omp target teams distribute parallel do private( i,dum1 ) &
    !$omp shared( imax,jmax,kmax,tmp1,tmp2,tmp3,hq,rho,u,v,w,p,CRATIO_INV,dummy,c13 )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 1) = hq(i, 1) - 0.5_wp*rho(i)*u(i)*tmp1(i)
        hq(i, 2) = hq(i, 2) - 0.5_wp*rho(i)*v(i)*tmp2(i)
        hq(i, 3) = hq(i, 3) - 0.5_wp*rho(i)*w(i)*tmp3(i)

        dum1 = tmp1(i) + tmp2(i) + tmp3(i)
        hq(i, 4) = hq(i, 4) + CRATIO_INV*( &
                   dummy*(tmp3(i)*tmp3(i) + tmp2(i)*tmp2(i) + tmp1(i)*tmp1(i) - c13*dum1*dum1) - p(i)*dum1) ! array tmp1 no longer needed
        tmp1(i) = c13*dum1
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

! Second derivative terms in the momentum equation
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 1), g(1), tmp1, tmp2)
#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp2,tmp4,hq,visc )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 1) = hq(i, 1) + visc*(tmp4(i) + tmp2(i))
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 2), g(2), tmp1, tmp2)

#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp2,tmp5,hq,visc )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 2) = hq(i, 2) + visc*(tmp5(i) + tmp2(i))
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 3), g(3), tmp1, tmp2)

#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp2,tmp6,hq,visc )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 3) = hq(i, 3) + visc*(tmp6(i) + tmp2(i))
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

! ###################################################################
! Enthalpy diffusion in energy equation
! ###################################################################
    if (nse_diffusion == EQNS_NONE) then; cond = 0.0_wp
    else; cond = visc/prandtl; end if

! calculate the enthalpy
    call THERMO_CALORIC_ENTHALPY(imax*jmax*kmax, s, T, tmp4)

! total flux
    call OPR_Partial_Z(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tmp4, tmp3, tmp5)
    call OPR_Partial_Y(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tmp4, tmp2, tmp5)
    call OPR_Partial_X(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tmp4, tmp1, tmp5)

#ifdef USE_APU
    !$omp target teams distribute parallel do private( i ) &
    !$omp shared( imax,jmax,kmax,tmp1,tmp2,tmp3,hq,cond )
#endif
    do i = 1, imax*jmax*kmax
        hq(i, 4) = hq(i, 4) + cond*(tmp1(i) + tmp2(i) + tmp3(i))
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'LEAVING RHS_FLOW_GLOBAL_2')
#endif

    return
end subroutine RHS_FLOW_GLOBAL_2
