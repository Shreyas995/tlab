#include "types.h"

#include "dns_error.h"
#include "dns_const.h"

module AVG_PHASE

    use TLab_WorkFlow
    use TLab_Constants, only: wp, wi, longi, efile, mas
    use TLAB_CONSTANTS, only: sizeofint, sizeofreal
    use FDM, only: g
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Time, only: rtime
    use NavierStokes, only: visc, froude, rossby, prandtl
    use Thermodynamics, only: mach
    use NavierStokes, only: nse_eqns, DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL
    use TLab_Memory, only: inb_flow, inb_scal
    use TLAB_ARRAYS, only: q, s
    use TLab_Arrays, only: wrk2d, wrk3d
    use Thermodynamics, only: gama0
    use TLab_Memory, only: Tlab_Allocate_Real_LONG
    ! NOTE: c_f_pointer/c_loc are deliberately NOT imported. The plane reductions
    ! read the source fields by name; aliasing them would reintroduce the MI300A
    ! "different device address for an aliased pointer" fault (see CLAUDE.md).

    implicit none
#ifdef USE_APU
    ! MI300A: required so this unit's !$omp target regions share host memory
    ! coherently (matches the program-scope declaration in dns_main.f90). Per
    ! OpenMP the directive must appear in EVERY compilation unit containing
    ! device constructs; a mixed USM/non-USM binary is UB on Cray CCE.
    !$omp requires unified_shared_memory

    interface
        ! Device-scope flush. The plane reductions below are written by the GPU and
        ! then read by the CPU for MPI_Reduce; on the APU `target update`/`map(from:)`
        ! are no-ops, so hipDeviceSynchronize is the only reliable handoff.
        function hipDeviceSynchronize() bind(C, name='hipDeviceSynchronize') result(ierr)
            use, intrinsic :: iso_c_binding, only: c_int
            integer(c_int) :: ierr
        end function hipDeviceSynchronize
    end interface
    ! This module has no `private` default, so without this the interface is exported
    ! and collides with the local hipDeviceSynchronize interface in every unit that
    ! does an unrestricted `use AVG_PHASE` (Cray ftn-613: "defined with more than one
    ! explicit interface"). The codebase convention is a local interface per scoping
    ! unit; keep this one internal.
    private :: hipDeviceSynchronize
#endif
    type phaseavg_dt
        sequence
        logical :: active
        integer(wi) :: stride
        character(32) :: type
    end type phaseavg_dt

    interface AvgPhaseSpace
        module procedure AvgPhaseSpaceFieldPtr, AvgPhaseSpaceIndex
    end interface AvgPhaseSpace

    type(phaseavg_dt) :: PhAvg
    real(wp), dimension(:), allocatable, target :: avg_flow, avg_stress, avg_p, avg_scal, avg_flux
    integer(wi) :: nxy, nxz, nyz, nz_total
    integer(wi) :: avg_planes
    character(len=32), parameter :: avgu_name = 'avg_flow'
    character(len=32), parameter :: avgstr_name = 'avg_stress'
    character(len=32), parameter :: avgp_name = 'avg_p'
    character(len=32), parameter :: avgs_name = 'avg_scal'
    character(len=32), parameter :: avgflux_name = 'avg_flux'

    integer, parameter, public :: IO_SCAL = 1       ! Header of scalar field
    integer, parameter, public :: IO_FLOW = 2       ! Header of flow field

    ! Destination selector for the shared AvgPhaseCalcProduct worker. The target
    ! arrays are allocated only on ims_pro_k == 0, so they are resolved to a pointer
    ! inside the worker (as AvgPhaseSpaceExec already does) rather than passed in.
    integer(wi), parameter :: AVGPH_STRESS = 1
    integer(wi), parameter :: AVGPH_FLUX = 2

    public :: AvgPhaseSpace
    public :: avg_flow, avg_p, avg_scal, avg_stress, avg_flux, avg_planes
    public :: PhAvg
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AvgPhaseInitializeMemory(C_FILE_LOC, restart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_pro, ims_npro_i, ims_pro_k
#endif

        implicit none
        character(len=*), intent(in) :: C_FILE_LOC
        integer(wi), intent(in) :: restart
        integer(longi) :: alloc_size
        ! ================================================================== !

        nxy = imax*jmax
        nyz = jmax*kmax
        nxz = imax*kmax

        if (restart == -1) then ! used for calls from outside dns_main / dns.x
            avg_planes = 0
        elseif (mod(restart, PhAvg%stride) == 0) then
            avg_planes = restart/PhAvg%stride
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Number of average planes not an integer. Change stride.')
            call TLAB_STOP(DNS_ERROR_AVG_PHASE)
        end if

#ifdef USE_MPI
        if (ims_pro_k == 0) then
#endif
            alloc_size = imax*jmax*(avg_planes + 1)
            call Tlab_Allocate_Real_LONG(C_FILE_LOC, avg_flow, [alloc_size*inb_flow], 'avgflow.')
            call Tlab_Allocate_Real_LONG(C_FILE_LOC, avg_stress, [alloc_size*6], 'avgstr.') ! allocated not yet coded
            call Tlab_Allocate_Real_LONG(C_FILE_LOC, avg_p, [alloc_size*1], 'avgp.')
            call Tlab_Allocate_Real_LONG(C_FILE_LOC, avg_scal, [alloc_size*inb_scal], 'avgscal.')
            call Tlab_Allocate_Real_LONG(C_FILE_LOC, avg_flux, [alloc_size*3], 'avgflux.') ! velocity-scalar flux u_i*s1 (3 components)

            avg_flow(:) = 0.0_wp
            avg_stress(:) = 0.0_wp
            avg_p(:) = 0.0_wp
            avg_scal(:) = 0.0_wp
            avg_flux(:) = 0.0_wp
#ifdef USE_MPI
        end if
#endif
        return
    end subroutine AvgPhaseInitializeMemory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Plane reductions: average over the LOCAL k-extent into an (imax*jmax) plane.
!
! Both workers are FUSED -- the pointwise term is consumed directly by the
! k-reduction, so the isize_field staging buffer the old code materialized in
! wrk3d (zero it, fill it, read it straight back) is never created at all. That
! removes ~3 of the ~5 full-field memory passes per call.
!
! The loop ORDER differs by target, and the difference matters a lot:
!
!   GPU: ij outer / k inner, one private accumulator per ij. Neighbouring threads
!        (ij, ij+1) touch neighbouring addresses => fully coalesced, and each
!        output element is owned by one thread => NO atomics (unlike AVG_IK_V).
!
!   CPU: k outer / ij inner, accumulating into the plane. The GPU order would walk
!        the inner loop with a stride of nxy*8 bytes (~288 KB at 4 ranks here) --
!        a new cache line and often a new page every step, which measured ~18%
!        SLOWER overall than the staged version it replaced. Streaming ij-inner
!        keeps both source reads sequential and the whole plane hot in L2.
!
! Both orders sum the same terms in the same sequence per output element, so they
! agree bit-for-bit with each other.
!
! The fields (q, s, and the pressure field) are already GPU-resident in unified
! memory; they are passed as explicit-shape dummies and read BY NAME inside the
! target region. Do NOT reintroduce a c_f_pointer/c_loc alias here -- on MI300A an
! aliased pointer and its backing array have different device addresses inside a
! target region (see CLAUDE.md, "same handle" rule).
!
! Division by g(3)%size (the GLOBAL z size) happens once at the end instead of per
! k-term; the MPI_Reduce over ims_comm_z at the call site completes the average.
! Fewer divides and slightly less rounding than the old per-term form.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AvgPhasePlaneAvg(f, plane)
        real(wp), intent(in) :: f(isize_field)
        real(wp), intent(out) :: plane(nxy)

        integer(wi) :: ij, k, base
        real(wp) :: acc, znorm

        znorm = real(g(3)%size, wp)     ! hoisted: no derived-type access on device

#ifdef USE_APU
        !$omp target teams distribute parallel do default(shared) private(ij,k,acc) &
        !$omp if (nxy*kmax > mas)
        do ij = 1, nxy
            acc = 0.0_wp
            do k = 1, kmax
                acc = acc + f(nxy*(k - 1) + ij)
            end do
            plane(ij) = acc/znorm
        end do
        !$omp end target teams distribute parallel do
#else
        do ij = 1, nxy
            plane(ij) = 0.0_wp
        end do
        do k = 1, kmax
            base = nxy*(k - 1)
            do ij = 1, nxy
                plane(ij) = plane(ij) + f(base + ij)
            end do
        end do
        do ij = 1, nxy
            plane(ij) = plane(ij)/znorm
        end do
#endif
    end subroutine AvgPhasePlaneAvg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AvgPhasePlaneProd(f1, f2, plane)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(wp), intent(in) :: f1(isize_field), f2(isize_field)
        real(wp), intent(out) :: plane(nxy)

        integer(wi) :: ij, k, idx, base
        real(wp) :: acc, znorm

        znorm = real(g(3)%size, wp)

#ifdef USE_APU
        !$omp target teams distribute parallel do default(shared) private(ij,k,idx,acc) &
        !$omp if (nxy*kmax > mas)
        do ij = 1, nxy
            acc = 0.0_wp
            do k = 1, kmax
                idx = nxy*(k - 1) + ij
                acc = acc + f1(idx)*f2(idx)
            end do
            plane(ij) = acc/znorm
        end do
        !$omp end target teams distribute parallel do
#else
        do ij = 1, nxy
            plane(ij) = 0.0_wp
        end do
        do k = 1, kmax
            base = nxy*(k - 1)
            do ij = 1, nxy
                plane(ij) = plane(ij) + f1(base + ij)*f2(base + ij)
            end do
        end do
        do ij = 1, nxy
            plane(ij) = plane(ij)/znorm
        end do
#endif
    end subroutine AvgPhasePlaneProd

    subroutine AvgPhaseSpaceFieldPtr(localsum, nfield, itr, it_first, it_save, field)
        implicit none
        real(wp), dimension(imax, jmax), intent(inout) :: localsum
        integer(wi), intent(in) :: nfield
        integer(wi), intent(in) :: itr, it_first, it_save
        real(wp), pointer, intent(in) :: field(:)

        integer :: index_loc = 4 ! Index needs to be set to appropriate value for pressure. Needed later for if statement

        call AvgPhaseSpaceExec(localsum, nfield, itr, it_first, it_save, index_loc, field)
    end subroutine AvgPhaseSpaceFieldPtr

    subroutine AvgPhaseSpaceIndex(localsum, nfield, itr, it_first, it_save, index)
        implicit none
        real(wp), dimension(imax, jmax), intent(inout) :: localsum
        integer(wi), intent(in) :: nfield
        integer(wi), intent(in) :: itr, it_first, it_save, index
        real(wp), pointer :: field_loc(:) => null()

        call AvgPhaseSpaceExec(localsum, nfield, itr, it_first, it_save, index, field_loc)
    end subroutine AvgPhaseSpaceIndex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AvgPhaseSpaceExec(localsum, nfield, itr, it_first, it_save, index, field)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif

        implicit none
        real(wp), dimension(imax*jmax), intent(inout) :: localsum
        integer(wi), intent(in) :: nfield
        integer(wi), intent(in) :: itr, it_first, it_save, index
        ! Assumed-size in the field index so ifld > 1 stays in bounds (the old code
        ! reached the same elements through a c_f_pointer shaped [isize_field*nfield]).
        real(wp), dimension(isize_field, *), target, intent(in) :: field

        integer(wi) :: ifld, plane_id
        real(wp), dimension(:), pointer :: avg_ptr
        integer(wi) :: iavg_srt, iavg_end, lpl_srt, lpl_end
#ifdef USE_APU
        integer :: hip_err
#endif
        ! ================================================================== !
        ! Calculation of the plane id to write the spatial average
        plane_id = 1
        if (it_save /= 0) plane_id = mod((itr - 1) - (it_first), it_save) + 1

        ! Determing the tendency to be written. The source field is NOT aliased via
        ! c_f_pointer any more: each branch below hands the plane reduction the real
        ! array by name (see the "same handle" note on AvgPhasePlaneAvg).
        if (index == 1) then
            avg_ptr => avg_flow
        elseif (index == 2) then
            avg_ptr => avg_scal
        elseif (index == 4) then
            avg_ptr => avg_p
        elseif (index == 5) then
            avg_ptr => avg_stress
            ! Not yet coded
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Unassigned case type check the index of the field in AvgPhaseSpaceExec')
            call TLAB_STOP(DNS_ERROR_AVG_PHASE)
        end if

        if ((index == 1) .or. (index == 2) .or. (index == 4)) then
            do ifld = 1, nfield
                ! Fused plane reduction directly off the source field.
                if (index == 1) then
                    call AvgPhasePlaneAvg(q(1, ifld), localsum)
                elseif (index == 2) then
                    call AvgPhasePlaneAvg(s(1, ifld), localsum)
                else
                    call AvgPhasePlaneAvg(field(1, ifld), localsum)
                end if
#ifdef USE_APU
                hip_err = hipDeviceSynchronize()    ! GPU-written localsum -> CPU MPI
#endif

                ! Computing the local sum from start and end of the field for accumulating the space averages
                iavg_srt = (ifld - 1)*nxy*(avg_planes + 1) + nxy*(plane_id - 1) + 1
                iavg_end = (ifld - 1)*nxy*(avg_planes + 1) + nxy*plane_id
#ifdef USE_MPI
                if (ims_pro_k == 0) then
                    call MPI_Reduce(localsum, avg_ptr(iavg_srt:iavg_end), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err) ! avg_ptr(imax*jmax*restarts*fld)
                else
                    ! Non-root: recvbuf is not significant, but it must NOT be MPI_IN_PLACE
                    ! (that sentinel is only legal in the SEND buffer at the root). OpenMPI
                    ! (Curta) rejects MPI_IN_PLACE here; pass a real, distinct scratch (wrk3d).
                    call MPI_Reduce(localsum, wrk3d, nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
                end if
#else
                avg_ptr(iavg_srt:iavg_end) = localsum
#endif

                lpl_srt = (ifld - 1)*nxy*(avg_planes + 1) + nxy*avg_planes + 1
                lpl_end = ifld*nxy*(avg_planes + 1)

#ifdef USE_MPI
                if (ims_pro_k == 0) then
#endif
                    avg_ptr(lpl_srt:lpl_end) = avg_ptr(lpl_srt:lpl_end) + avg_ptr(iavg_srt:iavg_end)/avg_planes
#ifdef USE_MPI
                end if
#endif
            end do
        end if
        return
    end subroutine AvgPhaseSpaceExec

    subroutine AvgPhaseStress(q, itr, it_first, it_save)
        ! Assumed-size in the component index: q is handed through as a base address,
        ! so the q(1,iq) actual arguments below need no descriptor and no copy. (An
        ! assumed-shape dummy would let the compiler emit a contiguity check and, in
        ! the worst case, a full isize_field temporary per component.)
        real(wp), dimension(isize_field, *), intent(in) :: q
        integer(wi), intent(in) :: itr
        integer(wi), intent(in) :: it_first
        integer(wi), intent(in) :: it_save

        integer(wi) :: plane_id

        plane_id = 1
        if (it_save /= 0) plane_id = mod((itr - 1) - (it_first), it_save) + 1

        ! Component slots in avg_stress: uu 1, uv 2, uw 3, vv 4, vw 5, ww 6
        ! (q components: 1 = u, 2 = v, 3 = w). Order chosen to reuse cache.
        call AvgPhaseCalcProduct(q(1, 1), q(1, 1), AVGPH_STRESS, 1, plane_id)
        call AvgPhaseCalcProduct(q(1, 1), q(1, 2), AVGPH_STRESS, 2, plane_id)
        call AvgPhaseCalcProduct(q(1, 2), q(1, 2), AVGPH_STRESS, 4, plane_id)
        call AvgPhaseCalcProduct(q(1, 2), q(1, 3), AVGPH_STRESS, 5, plane_id)
        call AvgPhaseCalcProduct(q(1, 1), q(1, 3), AVGPH_STRESS, 3, plane_id)
        call AvgPhaseCalcProduct(q(1, 3), q(1, 3), AVGPH_STRESS, 6, plane_id)

    end subroutine AvgPhaseStress

    subroutine AvgPhaseFlux(q, s, itr, it_first, it_save)
        real(wp), dimension(isize_field, *), intent(in) :: q
        real(wp), dimension(isize_field, *), intent(in) :: s
        integer(wi), intent(in) :: itr
        integer(wi), intent(in) :: it_first
        integer(wi), intent(in) :: it_save

        integer(wi) :: plane_id

        ! Velocity-scalar flux uses the FIRST scalar only; needs at least one scalar.
        if (inb_scal < 1) return

        plane_id = 1
        if (it_save /= 0) plane_id = mod((itr - 1) - (it_first), it_save) + 1

        ! Component slots in avg_flux: u*s1 1, v*s1 2, w*s1 3
        call AvgPhaseCalcProduct(q(1, 1), s(1, 1), AVGPH_FLUX, 1, plane_id)
        call AvgPhaseCalcProduct(q(1, 2), s(1, 1), AVGPH_FLUX, 2, plane_id)
        call AvgPhaseCalcProduct(q(1, 3), s(1, 1), AVGPH_FLUX, 3, plane_id)

    end subroutine AvgPhaseFlux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AvgPhaseCalcProduct(field1, field2, dest_id, comp_id, plane_id)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Single worker for both the Reynolds-stress and the velocity-scalar-flux phase
! averages. The two were byte-identical apart from the destination array, so they
! are merged here and selected by dest_id.
!
! Was: zero all of wrk3d (dead -- the next line overwrote every element the
! k-loop ever reads), stage field1*field2 into wrk3d, then read wrk3d straight
! back in the k-loop. That is ~5 full-field memory passes per call, on one CPU
! core, nine times per iteration.
!
! Now: one fused device pass (2 reads, no full-field write) into the plane
! accumulator. wrk3d is no longer touched as a compute buffer anywhere in this
! module -- it survives only as the ignored non-root MPI_Reduce recvbuf.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_comm_z, ims_err, ims_pro, ims_pro_k
#endif
        real(wp), intent(in) :: field1(isize_field), field2(isize_field)
        integer(wi), intent(in) :: dest_id      ! AVGPH_STRESS or AVGPH_FLUX
        integer(wi), intent(in) :: comp_id      ! component slot within that array
        integer(wi), intent(in) :: plane_id

        real(wp), dimension(:), pointer :: avg_ptr
        integer(wi) :: iavg_srt, iavg_end, lpl_srt, lpl_end
#ifdef USE_APU
        integer :: hip_err
#endif

        if (dest_id == AVGPH_STRESS) then
            avg_ptr => avg_stress
        else
            avg_ptr => avg_flux
        end if

        ! Fused product + k-reduction. Writes every element of wrk2d(1:nxy,1), so
        ! the old defensive wrk2d/wrk3d zeroing is not needed.
        call AvgPhasePlaneProd(field1, field2, wrk2d(1, 1))
#ifdef USE_APU
        hip_err = hipDeviceSynchronize()        ! GPU-written wrk2d -> CPU-side MPI
#endif

        iavg_srt = (comp_id - 1)*nxy*(avg_planes + 1) + nxy*(plane_id - 1) + 1
        iavg_end = (comp_id - 1)*nxy*(avg_planes + 1) + nxy*(plane_id)

#ifdef USE_MPI
        if (ims_pro_k == 0) then
            call MPI_Reduce(wrk2d, avg_ptr(iavg_srt:iavg_end), nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
        else
            ! Non-root: recvbuf is not significant, but it must NOT be MPI_IN_PLACE
            ! (that sentinel is only legal in the SEND buffer at the root; OpenMPI
            ! rejects it here). wrk3d is entirely free now, so reuse it as scratch.
            call MPI_Reduce(wrk2d, wrk3d, nxy, MPI_REAL8, MPI_SUM, 0, ims_comm_z, ims_err)
        end if
#else
        avg_ptr(iavg_srt:iavg_end) = wrk2d(1:nxy, 1)
#endif

        lpl_srt = (comp_id - 1)*nxy*(avg_planes + 1) + nxy*avg_planes + 1
        lpl_end = (comp_id)*nxy*(avg_planes + 1)

#ifdef USE_MPI
        if (ims_pro_k == 0) then
#endif
            avg_ptr(lpl_srt:lpl_end) = avg_ptr(lpl_srt:lpl_end) + avg_ptr(iavg_srt:iavg_end)/avg_planes
#ifdef USE_MPI
        end if
#endif

    end subroutine AvgPhaseCalcProduct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine IO_WRITE_HEADER(unit, isize, nx, ny, nz, nt, params)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer, intent(in) :: unit, isize
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(in) :: params(isize)

        ! -------------------------------------------------------------------
        integer(wi) offset

        !########################################################################
        offset = 5*SIZEOFINT + isize*SIZEOFREAL

        write (unit) offset, nx, ny, nz, nt

        if (isize > 0) then   ! do not write params to file if there are none
            write (unit) params(1:isize)
        end if

        return
    end subroutine IO_WRITE_HEADER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine IO_Write_AvgPhase(avg_planes, nfield, iheader, it_save, stride, basename, index, avg_ptr, avg_start)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use TLab_Memory, only: imax, jmax
        use FDM, only: g
        use TLab_Time, only: rtime, itime
        use NavierStokes, only: visc, froude, rossby, prandtl
        use Thermodynamics, only: mach
        use NavierStokes, only: nse_eqns
        use TLAB_CONSTANTS, only: sizeofint, sizeofreal
        use Thermodynamics, only: gama0

#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_comm_x, ims_err, ims_npro_i, ims_pro_i, ims_pro, ims_comm_z, ims_pro_k
#endif
        implicit none
        integer(wi), intent(in) :: avg_planes
        integer(wi), intent(in) :: nfield
        integer(wi), intent(in) :: it_save
        integer(wi), intent(in) :: stride
        character(len=*), intent(in) :: basename
        integer(wi), intent(in) :: index
        real(wp), dimension(:), intent(in) :: avg_ptr
        integer(wi), intent(in), optional :: avg_start

        character(len=128) :: name
        character(len=32) :: varname(1)
        integer(wi), parameter :: isize_max = 20
        real(wp) :: params(isize_max)
        integer(wi) :: isize, iheader, ifld, ifld_srt, ifld_end
        character(len=10) :: start, end, fld_id
        integer(wi) :: arr_planes, header_offset, ioffset_local
        integer(wi) :: nxy

#ifdef USE_MPI
        integer(kind=MPI_OFFSET_KIND) :: f_offset
        type(MPI_File) :: f_handle
        TYPE(MPI_Datatype) :: ftype, mtype
        type(MPI_Status) :: status
#endif
        nxy = imax*jmax

        if (index > 9 .or. index == 3 .or. index == 5 .or. index == 6 .or. index == 7) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Unassigned case type check the index of the field in PhaseAvg_Write')
            call TLAB_STOP(DNS_ERROR_AVG_PHASE)
        end if

        nz_total = it_save/stride + 1

        isize = 0
        isize = isize + 1; params(isize) = rtime
        isize = isize + 1; params(isize) = visc ! inverse of reynolds
        if (iheader == IO_SCAL) then
            isize = isize + 1 + 1                     ! prepare space for schmidt and damkohler

        else if (iheader == IO_FLOW) then
            isize = isize + 1; params(isize) = froude
            isize = isize + 1; params(isize) = rossby
            if (nse_eqns == DNS_EQNS_INTERNAL .or. nse_eqns == DNS_EQNS_TOTAL) then
                isize = isize + 1; params(isize) = gama0
                isize = isize + 1; params(isize) = prandtl
                isize = isize + 1; params(isize) = mach
            end if
        end if
        ! INITIALIZATION OF MPI TYPES SHOULD ONLY BE CARRIED OUT ONCE *AND* DURING INITIALIZATION
        header_offset = 5*SIZEOFINT + isize*SIZEOFREAL

        if (present(avg_start)) then
            write (start, '(I10)') (avg_start)
            write (end, '(I10)') avg_start
        else
            write (start, '(I10)') (itime - it_save + 1)
            write (end, '(I10)') itime
        end if

        do ifld = 1, nfield
            ifld_srt = (ifld - 1)*nxy*(avg_planes + 1) + 1
            ifld_end = ifld*nxy*(avg_planes + 1)
            write (fld_id, '(I10)') ifld
            varname(1) = ''
            if (start == end) then ! write single iteration
                name = trim(adjustl(basename))//trim(adjustl(start))//'.'//trim(adjustl(fld_id))
            else ! write multiple iteration including phase average
                name = trim(adjustl(basename))//trim(adjustl(start)) &
                       //'_'//trim(adjustl(end))//'.'//trim(adjustl(fld_id))
            end if

            arr_planes = (jmax*(avg_planes + 1))
            ! Define the array size for planes and file offset
#ifdef USE_MPI
            f_offset = header_offset + ims_pro_i*imax*8
#endif

#ifdef USE_MPI
            if (ims_pro == 0) then
#endif
#define LOC_STATUS "unknown"
#define LOC_UNIT_ID 75
#include "dns_open_file.h"
                call IO_WRITE_HEADER(LOC_UNIT_ID, isize, g(1)%size, g(2)%size, nz_total, itime, params)
                close (LOC_UNIT_ID)
#ifdef USE_MPI
            end if
#endif

#ifdef USE_MPI
            call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
            if (ims_pro_k == 0) then
                ! Create the MPI derived data types for the file view and contiguous blocks
                call MPI_TYPE_VECTOR(arr_planes, imax, imax*ims_npro_i, MPI_REAL8, ftype, ims_err)
                call MPI_TYPE_COMMIT(ftype, ims_err)
                call MPI_TYPE_CONTIGUOUS(imax, MPI_REAL8, mtype, ims_err)
                call MPI_TYPE_COMMIT(mtype, ims_err)

                ! Open the file for writing
                call MPI_FILE_OPEN(ims_comm_x, name, ior(MPI_MODE_CREATE, MPI_MODE_WRONLY), MPI_INFO_NULL, f_handle, ims_err)

                ! Set the file view
                call MPI_File_set_view(f_handle, f_offset, MPI_REAL8, ftype, 'native', MPI_INFO_NULL, ims_err)

                ! Write the data to the file
                call MPI_FILE_WRITE_ALL(f_handle, avg_ptr(ifld_srt), arr_planes, mtype, status, ims_err)

                ! Close the file
                call MPI_FILE_CLOSE(f_handle, ims_err)

                ! Free the MPI derived data types
                call MPI_TYPE_FREE(ftype, ims_err)
                call MPI_TYPE_FREE(mtype, ims_err)
            end if
#else
#include "dns_open_file.h"
            ioffset_local = header_offset + 1
            write (LOC_UNIT_ID, POS=ioffset_local) avg_ptr(ifld_srt:ifld_end)
            close (LOC_UNIT_ID)
#endif
        end do
        return
    end subroutine IO_Write_AvgPhase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AvgPhaseResetVariable()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_comm_x, ims_err, ims_pro, ims_pro_k
#endif

#ifdef USE_MPI
        if (ims_pro_k == 0) then
#endif
            avg_flow(:) = 0.0_wp
            avg_stress(:) = 0.0_wp
            avg_p(:) = 0.0_wp
            avg_scal(:) = 0.0_wp
            avg_flux(:) = 0.0_wp
#ifdef USE_MPI
        end if
#endif
    end subroutine AvgPhaseResetVariable
end module
