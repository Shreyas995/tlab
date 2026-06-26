#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!#
!# Evolution equations, nonlinear term in convective form and the
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!#
!# It is written such that u and w transposes are calculated first for the
!# Ox and Oz momentum equations, stored in tmp5 and tmp6 and then used as needed.
!# This saves 2 MPI transpositions.
!# Includes the scalar to benefit from the same reduction
!#
!########################################################################
subroutine RHS_GLOBAL_INCOMPRESSIBLE_1()
#ifdef USE_OPENMP
    use OMP_LIB
#endif
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
#endif
    use TLab_Constants, only: wp, wi, BCS_NN, mas
    use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use FDM, only: g
    use TLab_WorkFlow, only: stagger_on
    use TLab_Time, only: itime
#ifdef USE_APU
    use TLabMPI_Transpose, only: trp_dbg_fft   ! gate node-window I-transpose probes around self-burgX
#endif
    use TLab_Arrays
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
    use TLab_Pointers_3D, only: p_tmp2
    use Thermo_Anelastic
    use TLab_OpenMP
    use DNS_ARRAYS
    use DNS_LOCAL, only: remove_divergence
    use DNS_LOCAL, only: use_tower
    use DNS_LOCAL, only: nitera_first, nitera_save
    use TIME, only: rkm_substep, rkm_endstep, dte
    use DNS_TOWER
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    use IBM_VARS, only: imode_ibm, imode_ibm_scal, ibm_burgers
    use OPR_Partial
    use OPR_Burgers
    use OPR_Elliptic
    use OPR_FILTERS
    use AVG_PHASE
    use TLab_Debug
    use TLab_Arrays, only: wrk3d

    implicit none

#ifdef USE_APU
    interface
        function hipDeviceSynchronize() bind(C, name='hipDeviceSynchronize') result(ierr)
            integer :: ierr
        end function hipDeviceSynchronize
    end interface
    integer :: hip_sync_err
#endif

    ! -----------------------------------------------------------------------
    integer(wi) iq, is, ij
    integer ibc, bcs(2, 2)
    real(wp) dummy

    integer(wi) siz, srt, end    !  Variables for OpenMP Partitioning

    real(wp), dimension(:, :, :), pointer :: p_bcs

#ifdef USE_ESSL
    integer ilen
#endif

#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'ENTERING SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1')
#endif

    ! #######################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

#ifdef USE_ESSL
    ilen = isize_field
#endif

    ! #######################################################################
    ! Preliminaries for Scalar BC
    ! (flow BCs initialized below as they are used for pressure in between)
    ! #######################################################################
    ! Default is zero
    BcsScalJmin%ref(:, :, :) = 0.0_wp
    BcsScalJmax%ref(:, :, :) = 0.0_wp

    ! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
    if (any(BcsScalJmin%SfcType(1:inb_scal) == DNS_SFC_LINEAR) .or. any(BcsScalJmax%SfcType(1:inb_scal) == DNS_SFC_LINEAR)) then
        do is = 1, inb_scal
            p_bcs(1:imax, 1:jmax, 1:kmax) => hs(1:imax*jmax*kmax, is)
            if (BcsScalJmin%SfcType(is) == DNS_SFC_LINEAR) BcsScalJmin%ref(:, :, is) = p_bcs(:, 1, :)
            if (BcsScalJmax%SfcType(is) == DNS_SFC_LINEAR) BcsScalJmax%ref(:, :, is) = p_bcs(:, jmax, :)
        end do
    end if


    ! COARSE crash-localization: per-component velocity ENTERING the RHS (u=q1, v=q2, w=q3). If one of these
    ! is already corrupt here the seed is upstream (RK update / sources / post-step output), not the RHS.
    DNS_PROBE('RHS1:in-u', u(1), isize_field, rkm_substep)
    DNS_PROBE('RHS1:in-v', v(1), isize_field, rkm_substep)
    DNS_PROBE('RHS1:in-w', w(1), isize_field, rkm_substep)

    ! #######################################################################
    ! Diffusion and advection terms
    ! #######################################################################
    ! Preliminaries for IBM use
    ! (OPR_Burgers_X/Y/Z uses modified fields for derivatives)
    if (imode_ibm == 1) ibm_burgers = .true.

    ! Diagonal terms and transposed velocity arrays
#ifdef DNS_DEBUG_PROBES
    trp_dbg_fft = .true.    ! emit the NWFR/NWBR node-window I-transpose probes for THIS call (the 2026-06-24 seed)
#endif
    call OPR_Burgers_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, u, u, tmp1, tmp4) ! store u transposed in tmp4
#ifdef DNS_DEBUG_PROBES
    trp_dbg_fft = .false.
#endif

    call OPR_Burgers_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, v, v, tmp2, tmp5) ! store v transposed in tmp5
    call OPR_Burgers_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, w, w, tmp3, tmp6) ! store w transposed in tmp6

    ! COARSE crash-localization: diagonal self-advection-diffusion terms (Burgers X/Y/Z SELF). tmp3 uses
    ! Burgers_Z = K-transpose, tmp1 uses Burgers_X = I-transpose; a jump here pins the faulting Burgers/transpose.
    DNS_PROBE('RHS1:self-burgX-u', tmp1, isize_field, rkm_substep)
    DNS_PROBE('RHS1:self-burgY-v', tmp2, isize_field, rkm_substep)
    DNS_PROBE('RHS1:self-burgZ-w', tmp3, isize_field, rkm_substep)

    ! Ox momentum equation
    call OPR_Burgers_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, u, v, tmp7, tmp9, tmp5) ! tmp5 contains v transposed
    call OPR_Burgers_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, u, w, tmp8, tmp9, tmp6) ! tmp6 contains w transposed

    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
#ifdef USE_APU
    !$omp target teams distribute parallel do &
    !$omp private( ij ) &
    !$omp shared( srt,end,hq,tmp1,tmp7,tmp8 ) &
    !$omp if(end > mas)
#endif
    do ij = srt, end
        hq(ij, 1) = hq(ij, 1) + tmp1(ij) + tmp7(ij) + tmp8(ij)
    end do

#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

    ! COARSE crash-localization: Ox momentum tendency after its advection/diffusion (Burgers X/Y/Z of u).
    DNS_PROBE('RHS1:hq1-Ox', hq(1, 1), isize_field, rkm_substep)

    ! Oy momentum equation
    call OPR_Burgers_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, v, u, tmp7, tmp9, tmp4) ! tmp4 contains u transposed
    call OPR_Burgers_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, v, w, tmp8, tmp9, tmp6) ! tmp6 contains w transposed

    call TLab_OMP_PARTITION(isize_field, srt, end, siz)

#ifdef USE_APU
    !$omp target teams distribute parallel do &
    !$omp private( ij ) &
    !$omp shared( srt,end,hq,tmp2,tmp7,tmp8 ) &
    !$omp if(end > mas)
#endif

    do ij = srt, end
        hq(ij, 2) = hq(ij, 2) + tmp2(ij) + tmp7(ij) + tmp8(ij)
    end do

#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

    ! COARSE crash-localization: Oy momentum tendency after its advection/diffusion.
    DNS_PROBE('RHS1:hq2-Oy', hq(1, 2), isize_field, rkm_substep)

    ! Oz momentum equation
    call OPR_Burgers_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, w, u, tmp7, tmp9, tmp4) ! tmp4 contains u transposed
    call OPR_Burgers_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, w, v, tmp8, tmp9, tmp5) ! tmp5 contains v transposed
    call TLab_OMP_PARTITION(isize_field, srt, end, siz)

#ifdef USE_APU
    !$omp target teams distribute parallel do &
    !$omp private( ij ) &
    !$omp shared( srt,end,hq,tmp2,tmp7,tmp8 ) &
    !$omp if(end > mas)
#endif

    do ij = srt, end
        hq(ij, 3) = hq(ij, 3) + tmp3(ij) + tmp7(ij) + tmp8(ij)
    end do

#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif

    ! COARSE crash-localization: Oz momentum tendency after its advection/diffusion (Burgers Z = K-transpose).
    DNS_PROBE('RHS1:hq3-Oz', hq(1, 3), isize_field, rkm_substep)

    ! Sentinel: catch pollution from the advection/diffusion (OPR_Burgers) stage
    DNS_PROBE('RHS1:post-adv', hq(1, 1), isize_field*inb_flow, rkm_substep)

    ! IBM
    if (imode_ibm == 1) then
        ibm_burgers = .false. ! until here, IBM is used for flow fields
        if (imode_ibm_scal == 1) then ! IBM usage for scalar field
            ! (requirenments: only possible with objects on bottom boundary
            !  with homogeneous temperature in solid regions)
            ibm_burgers = .true.
        end if
    end if

    ! Scalar equations
    do is = 1, inb_scal
        call OPR_Burgers_X(OPR_B_U_IN, is, imax, jmax, kmax, bcs, s(1, is), u, tmp1, tmp9, tmp4) ! tmp4 contains u transposed
        call OPR_Burgers_Y(OPR_B_U_IN, is, imax, jmax, kmax, bcs, s(1, is), v, tmp2, tmp9, tmp5) ! tmp5 contains v transposed
        call OPR_Burgers_Z(OPR_B_U_IN, is, imax, jmax, kmax, bcs, s(1, is), w, tmp3, tmp9, tmp6) ! tmp6 contains w transposed

        call TLab_OMP_PARTITION(isize_field, srt, end, siz)
#ifdef USE_APU
    !$omp target teams distribute parallel do &
    !$omp private( ij ) &
    !$omp shared( srt,end,hs,tmp1,tmp2,tmp3 ) &
    !$omp if(end > mas)
#endif
        do ij = srt, end ! offload to APU
            hs(ij, is) = hs(ij, is) + tmp1(ij) + tmp2(ij) + tmp3(ij)
        end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif
    end do

    ! IBM usage for scalar field, done
    if (imode_ibm_scal == 1) ibm_burgers = .false.

    ! COARSE crash-localization: scalar tendency after its advection/diffusion.
    if (inb_scal > 0) DNS_PROBE('RHS1:scalar-hs', hs(1, 1), isize_field, rkm_substep)

    ! #######################################################################
    ! Impose buffer zone as relaxation terms
    ! #######################################################################
    if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
        call BOUNDARY_BUFFER_RELAX_FLOW()
    end if

    ! COARSE crash-localization: flow tendency after the buffer-relaxation stage (before the pressure step).
    DNS_PROBE('RHS1:post-buffer', hq(1, 1), isize_field*inb_flow, rkm_substep)


    ! #######################################################################
    ! Pressure term
    ! #######################################################################
    if (remove_divergence) then ! remove residual divergence

        call TLab_OMP_PARTITION(isize_field, srt, end, siz)
        dummy = 1.0_wp/dte

#ifdef USE_APU
        !$omp target teams distribute parallel do &
        !$omp private( ij ) &
        !$omp shared( srt,end,tmp2,tmp3,tmp4,hq,u,v,w,dummy ) &
        !$omp if(end > mas)
        do ij = srt, end ! offload to APU
            tmp2(ij) = hq(ij, 2) + v(ij)*dummy
            tmp3(ij) = hq(ij, 1) + u(ij)*dummy
            tmp4(ij) = hq(ij, 3) + w(ij)*dummy
        end do
        !$omp end target teams distribute parallel do
#elif defined(USE_ESSL)
        ilen = siz
        call DZAXPY(ilen, dummy, v(srt), 1, hq(srt, 2), 1, tmp2(srt), 1)
        call DZAXPY(ilen, dummy, u(srt), 1, hq(srt, 1), 1, tmp3(srt), 1)
        call DZAXPY(ilen, dummy, w(srt), 1, hq(srt, 3), 1, tmp4(srt), 1)

#else
        do ij = srt, end 
            tmp2(ij) = hq(ij, 2) + v(ij)*dummy
            tmp3(ij) = hq(ij, 1) + u(ij)*dummy
            tmp4(ij) = hq(ij, 3) + w(ij)*dummy
        end do
#endif

        if (imode_ibm == 1) then
            call IBM_BCS_FIELD(tmp2)
            call IBM_BCS_FIELD(tmp3)
            call IBM_BCS_FIELD(tmp4)
#ifdef USE_APU
            hip_sync_err = hipDeviceSynchronize()  ! flush GPU L2 before CPU-side MPI in OPR_Partial_X/Z
#endif
        end if
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp2)
            call Thermo_Anelastic_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp3)
            call Thermo_Anelastic_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp4)
        end if

        ! CRASH-LOC: the provisional-velocity divergence (div-forcing) is the first all-NaN field this run.
        ! Split it: provvel = input (should be clean); postX = after OPR_Partial_X (I-transpose = all-intra
        ! node-window for npro_i=6 = the COMMON mechanism, also in apudirect); postZ = after OPR_Partial_Z
        ! (K-transpose = inter-node MPI = fabricdirect-specific). Whichever first goes NaN pins the leg.
        DNS_PROBE('RHS1:dvg-provvel', tmp2, isize_field, rkm_substep)
        if (stagger_on) then ! staggering on horizontal pressure nodes
            !  Oy derivative
            call OPR_Partial_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp2, tmp5)
            DNS_PROBE('RHS1:dvg-postX', tmp5, isize_field, rkm_substep)
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp5, tmp2)
            call OPR_Partial_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp2, tmp1)
            DNS_PROBE('RHS1:dvg-postZ', tmp1, isize_field, rkm_substep)
            !  Ox derivative
            call OPR_Partial_X(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(1), tmp3, tmp5)
            call OPR_Partial_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp2)
            !  Oz derivative
            call OPR_Partial_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp4, tmp5)
            call OPR_Partial_Z(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp3)
        else
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp1)
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp3, tmp2)
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)
        end if

    else
        if (imode_ibm == 1) then
            call IBM_BCS_FIELD(hq(:, 2))
            call IBM_BCS_FIELD(hq(:, 1))
            call IBM_BCS_FIELD(hq(:, 3))
        end if
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call Thermo_Anelastic_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 2), tmp2)
            call Thermo_Anelastic_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 1), tmp3)
            call Thermo_Anelastic_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 3), tmp4)
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp1)
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp3, tmp2)
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)
        else
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), hq(:, 2), tmp1)
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), hq(:, 1), tmp2)
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), hq(:, 3), tmp3)
        end if

    end if
    ! -----------------------------------------------------------------------
    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
#ifdef USE_APU
    !$omp target teams distribute parallel do &
    !$omp private( ij ) &
    !$omp shared( srt,end,tmp1,tmp2,tmp3 ) &
    !$omp if(end > mas)
#endif
    do ij = srt, end
        tmp1(ij) = tmp1(ij) + tmp2(ij) + tmp3(ij) ! forcing term in tmp1
    end do
#ifdef USE_APU
    !$omp end target teams distribute parallel do
#endif
    ! DIAG: provisional-velocity divergence (= Poisson forcing). If THIS stays flat while [DIL] grows, the
    ! projection is failing to clean the mode; if THIS grows, the divergence is injected upstream.
    DNS_PROBE('RHS1:div-forcing', tmp1, isize_field, rkm_substep)



    ! -----------------------------------------------------------------------
    ! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    ! Stagger also Bcs
    if (imode_ibm == 1) then
        call IBM_BCS_FIELD(hq(:, 2))
#ifdef USE_APU
        if (stagger_on) hip_sync_err = hipDeviceSynchronize()  ! flush GPU L2 before CPU-side MPI in OPR_Partial_X/Z
#endif
    end if
    if (stagger_on) then ! todo: only need to stagger upper/lower boundary plane, not full h2-array
        call OPR_Partial_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), hq(:, 2), tmp5)
        call OPR_Partial_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp4)
        if (imode_ibm == 1) call IBM_BCS_FIELD_STAGGER(tmp4)
        p_bcs(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    else
        p_bcs(1:imax, 1:jmax, 1:kmax) => hq(1:imax*jmax*kmax, 2)
    end if

    if (nse_eqns == DNS_EQNS_ANELASTIC) then
        BcsFlowJmin%ref(:, :, 2) = p_bcs(:, 1, :)*rbackground(1)
        BcsFlowJmax%ref(:, :, 2) = p_bcs(:, jmax, :)*rbackground(g(2)%size)
    else
        BcsFlowJmin%ref(:, :, 2) = p_bcs(:, 1, :)
        BcsFlowJmax%ref(:, :, 2) = p_bcs(:, jmax, :)
    end if

    ! pressure in tmp1, Oy derivative in tmp3
    call OPR_Poisson(imax, jmax, kmax, BCS_NN, tmp1, p_tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)
    ! DIAG: pressure out of the Poisson solve (tmp1) on the fabricdirect path.
    DNS_PROBE('RHS1:p-poisson', tmp1, isize_field, rkm_substep)


    ! filter pressure p and its vertical gradient dpdy
    if (any(PressureFilter(:)%type /= DNS_FILTER_NONE)) then
        call OPR_FILTER(imax, jmax, kmax, PressureFilter, tmp1, txc(1:isize_field,4:6))
        call OPR_FILTER(imax, jmax, kmax, PressureFilter, tmp3, txc(1:isize_field,4:6))
    end if

    ! Saving pressure for towers to tmp array
    if (rkm_substep == rkm_endstep) then
        if (stagger_on .and. ( use_tower .or. PhAvg%active )) then ! Stagger pressure field back on velocity grid (only for towers)
            call OPR_Partial_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), tmp1, tmp5)
            call OPR_Partial_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp4)
        endif
        if ( use_tower ) &
            call DNS_TOWER_ACCUMULATE(tmp4, 4, wrk1d)
        if ( PhAvg%active) then   
            if (mod((itime+1),PhAvg%stride) == 0)  then
                call AvgPhaseSpace(wrk2d, 1, (itime+1)/PhAvg%stride, nitera_first, nitera_save/PhAvg%stride, tmp4)
            end if
        end if
    end if

    if (stagger_on) then
        !  vertical pressure derivative   dpdy - back on horizontal velocity nodes
        call OPR_Partial_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), tmp3, tmp5)
        ! print *, 'rhs_global_incompressible1 tmp3 cointains pressure ', sum(tmp3)
        call OPR_Partial_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp3)
        ! print *, 'rhs_global_incompressible1 tmp5 cointains pressure ', sum(tmp5)

        !  horizontal pressure derivative dpdz - back on horizontal velocity nodes
        call OPR_Partial_Z(OPR_P1_INT_PV, imax, jmax, kmax, bcs, g(3), tmp1, tmp5)
        ! print *, 'rhs_global_incompressible1 tmp1 cointains pressure ', sum(tmp1)

        call OPR_Partial_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp4)
        ! print *, 'rhs_global_incompressible1 tmp5 cointains pressure ', sum(tmp5)

        !  horizontal pressure derivative dpdx - back on horizontal velocity nodes
        call OPR_Partial_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), tmp1, tmp5)
        ! print *, 'rhs_global_incompressible1 tmp1 cointains pressure ', sum(tmp1)

        call OPR_Partial_X(OPR_P1_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp2)
        ! print *, 'rhs_global_incompressible1 tmp5 cointains pressure ', sum(tmp5)
    else
        !  horizontal pressure derivatives
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp1, tmp4)
    end if

    ! 'RHS1:post-pgrad' (below) trips, the host !$omp parallel do that subtracts them
    ! into hq is the culprit (the suspected GPU-written -> CPU-read coherency hazard).

    ! -----------------------------------------------------------------------
    ! Add pressure gradient
    ! -----------------------------------------------------------------------
    if (nse_eqns == DNS_EQNS_ANELASTIC) then
        call Thermo_Anelastic_WEIGHT_SUBTRACT(imax, jmax, kmax, ribackground, tmp2, hq(:, 1))
        call Thermo_Anelastic_WEIGHT_SUBTRACT(imax, jmax, kmax, ribackground, tmp3, hq(:, 2))
        call Thermo_Anelastic_WEIGHT_SUBTRACT(imax, jmax, kmax, ribackground, tmp4, hq(:, 3))

    else

    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
#ifdef USE_APU
        ! DIAG: pressure-gradient components (dpdx/dpdz/dpdy) just before they are subtracted from hq.
        DNS_PROBE('RHS1:pgrad-x', tmp2, isize_field, rkm_substep)
        DNS_PROBE('RHS1:pgrad-y', tmp3, isize_field, rkm_substep)
        DNS_PROBE('RHS1:pgrad-z', tmp4, isize_field, rkm_substep)
        !$omp parallel do default( shared ) private ( ij ) &
        !$omp if(end > mas)
        do ij = srt, end
            hq(ij, 1) = hq(ij, 1) - tmp2(ij)
            hq(ij, 2) = hq(ij, 2) - tmp3(ij)
            hq(ij, 3) = hq(ij, 3) - tmp4(ij)
        end do
        !$omp end parallel do
        ! DIAG: corrected velocity tendency after the pressure-gradient subtraction (the projected hq).
        DNS_PROBE('RHS1:hq-corrected', hq(1, 1), isize_field*inb_flow, rkm_substep)
#elif defined(USE_ESSL)
        ilen = siz
        dummy = -1.0_wp
        call DAXPY(ilen, dummy, tmp2(srt), 1, hq(srt, 1), 1)
        call DAXPY(ilen, dummy, tmp3(srt), 1, hq(srt, 2), 1)
        call DAXPY(ilen, dummy, tmp4(srt), 1, hq(srt, 3), 1)
#else
        do ij = srt, end
            hq(ij, 1) = hq(ij, 1) - tmp2(ij)
            hq(ij, 2) = hq(ij, 2) - tmp3(ij)
            hq(ij, 3) = hq(ij, 3) - tmp4(ij)
        end do
#endif
    end if

    ! (in particular the host !$omp parallel do at ~L388 reading GPU-written tmp2/3/4)

    ! #######################################################################
    ! Boundary conditions
    ! #######################################################################
    BcsFlowJmin%ref = 0.0_wp ! default is no-slip (dirichlet)
    BcsFlowJmax%ref = 0.0_wp ! Scalar BCs initialized at start of routine

    do iq = 1, inb_flow
        ibc = 0
        if (BcsFlowJmin%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsFlowJmax%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), hq(1, iq), &
                                        BcsFlowJmin%ref(1, 1, iq), BcsFlowJmax%ref(1, 1, iq), tmp1)
        end if
        if (imode_ibm == 1) call IBM_BCS_FIELD(hq(1, iq)) ! set tendency in solid to zero

        p_bcs(1:imax, 1:jmax, 1:kmax) => hq(1:imax*jmax*kmax, iq)
        p_bcs(:, 1, :) = BcsFlowJmin%ref(:, :, iq)
        p_bcs(:, jmax, :) = BcsFlowJmax%ref(:, :, iq)

    end do
    do is = 1, inb_scal
        ibc = 0
        if (BcsScalJmin%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsScalJmax%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), hs(1, is), &
                                        BcsScalJmin%ref(1, 1, is), BcsScalJmax%ref(1, 1, is), tmp1)
        end if

        if (BcsScalJmin%type(is) /= DNS_SFC_STATIC .or. &
            BcsScalJmax%type(is) /= DNS_SFC_STATIC) then
            call BOUNDARY_BCS_SURFACE_Y(is, bcs, s, hs, tmp1, tmp2)
        end if
        if (imode_ibm == 1) call IBM_BCS_FIELD(hs(1, is)) ! set tendency in solid to zero

        p_bcs(1:imax, 1:jmax, 1:kmax) => hs(1:imax*jmax*kmax, is)
        p_bcs(:, 1, :) = BcsScalJmin%ref(:, :, is)
        p_bcs(:, jmax, :) = BcsScalJmax%ref(:, :, is)

    end do


#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'LEAVING SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1')
#endif

    return
end subroutine RHS_GLOBAL_INCOMPRESSIBLE_1
