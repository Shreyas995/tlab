! =====================================================================================================
! vhipfft.f90 — hipFFT bring-up harness (Phase 0 of the GPU-FFT roadmap step).
!
! PURPOSE: prove hipFFT reproduces CPU FFTW to machine precision, standalone, BEFORE touching the
! production opr_fourier.f90. Runs BOTH backends on the same input and diffs them — the FFT analogue of
! src/valid/mpi/vmpi_nodewin.f90. No MPI (the production FFTs are per-rank/local; the MPI transpose is a
! separate, already-validated piece), so this is a single-rank GPU program.
!
! It replicates the EXACT FFTW plan layouts used by the Poisson solver (src/operators/opr_fourier.f90,
! single-rank / no-MPI branch):
!   Z (complex-to-complex, batched, INTERLEAVED): n=nz, howmany=nlines_z=(nx/2+1)*ny,
!       istride=ostride=nlines_z, idist=odist=1                         [opr_fourier.f90:113-121]
!   X (real-to-complex fwd / complex-to-real bwd, batched, CONTIGUOUS): n=nx, howmany=nlines_x=ny*nz,
!       istride=ostride=1, idist=nx (real) / nxc (cplx), odist=nxc / nx [opr_fourier.f90:164-171]
! Both FFTW and hipFFT are UNNORMALIZED; the roundtrip divides by n to recover the input.
!
! BUILD: standalone (NOT the main CMake), like the src/valid/mpi harnesses. See Makefile.hipfft.
!   - FFTW-only (e.g. local syntax check):  no -DUSE_HIPFFT  -> hipFFT path compiled out.
!   - Full (Hunter):  -DUSE_HIPFFT  ->  links -lhipfft -lamdhip64, runs both backends and diffs.
!
! HUNTER ITERATION POINTS (resolve on first run — these are the documented unknowns):
!   [B1] hipfftHandle is treated here as an OPAQUE POINTER (type(c_ptr)); hipfftPlanMany takes hipfftHandle*
!        (so plan is passed by reference). If the toolchain's hipfftHandle is `int`, switch plan to
!        integer(c_int).  [B2] data pointers: on the MI300A unified memory (HSA_XNACK=1) the device pointer
!        == host pointer, so we pass c_loc(host_array) straight to hipfftExec*. If hipFFT rejects a host
!        pointer, wrap the exec in `!$omp target data use_device_addr(...)` or allocate via omp_target_alloc.
!   [B3] hipfftType enum values below are the cuFFT-compatible ones; verify against <hipfft.h> if a plan fails.
! =====================================================================================================

#ifdef USE_HIPFFT
module hipfft_m
    use iso_c_binding
    implicit none
    ! hipfftType enum (cuFFT-compatible values; [B3]) and transform directions.
    integer(c_int), parameter :: HIPFFT_C2C = int(z'29', c_int)   ! 41
    integer(c_int), parameter :: HIPFFT_D2Z = int(z'6a', c_int)   ! 106 (double real -> double complex)
    integer(c_int), parameter :: HIPFFT_Z2D = int(z'6c', c_int)   ! 108 (double complex -> double real)
    integer(c_int), parameter :: HIPFFT_Z2Z = int(z'69', c_int)   ! 105 (double complex c2c)
    integer(c_int), parameter :: HIPFFT_FORWARD = -1
    integer(c_int), parameter :: HIPFFT_BACKWARD = 1
    interface
        ! hipfftResult hipfftPlanMany(hipfftHandle* plan, int rank, int* n, int* inembed, int istride,
        !     int idist, int* onembed, int ostride, int odist, hipfftType type, int batch)
        integer(c_int) function hipfftPlanMany(plan, rank, n, inembed, istride, idist, &
                                               onembed, ostride, odist, ttype, batch) bind(C, name="hipfftPlanMany")
            import :: c_int, c_ptr
            type(c_ptr)            :: plan                 ! hipfftHandle* (handle written here); [B1]
            integer(c_int), value :: rank
            integer(c_int)        :: n(*), inembed(*), onembed(*)
            integer(c_int), value :: istride, idist, ostride, odist, ttype, batch
        end function
        integer(c_int) function hipfftExecZ2Z(plan, idata, odata, direction) bind(C, name="hipfftExecZ2Z")
            import :: c_int, c_ptr
            type(c_ptr), value    :: plan, idata, odata
            integer(c_int), value :: direction
        end function
        integer(c_int) function hipfftExecD2Z(plan, idata, odata) bind(C, name="hipfftExecD2Z")
            import :: c_int, c_ptr
            type(c_ptr), value :: plan, idata, odata
        end function
        integer(c_int) function hipfftExecZ2D(plan, idata, odata) bind(C, name="hipfftExecZ2D")
            import :: c_int, c_ptr
            type(c_ptr), value :: plan, idata, odata
        end function
        integer(c_int) function hipfftDestroy(plan) bind(C, name="hipfftDestroy")
            import :: c_int, c_ptr
            type(c_ptr), value :: plan
        end function
        integer(c_int) function hipDeviceSynchronize() bind(C, name="hipDeviceSynchronize")
            import :: c_int
        end function
    end interface
end module hipfft_m
#endif

program VHIPFFT
    use iso_c_binding
    use omp_lib, only: omp_is_initial_device
#ifdef USE_HIPFFT
    use hipfft_m
#endif
    implicit none

    ! FFTW3 constants (define locally so the harness is self-contained — no fftw3.f03 include needed).
    integer, parameter :: FFTW_FORWARD = -1, FFTW_BACKWARD = +1, FFTW_ESTIMATE = 64

    integer :: nx, ny, nz, nxc, nlines_z, nlines_x, i
    integer(8) :: pf, pb                                  ! FFTW (F77) plan handles
    complex(c_double_complex), allocatable, target :: cz0(:), czf_ref(:), czf_hip(:), czb(:)
    real(c_double),            allocatable, target :: rx0(:), rxb_ref(:), rxb_hip(:)
    complex(c_double_complex), allocatable, target :: cxf_ref(:), cxf_hip(:)
    integer :: n1(1), ie1(1), oe1(1)
    logical :: on_dev, all_ok
    real(c_double) :: dz, dx
#ifdef USE_HIPFFT
    type(c_ptr) :: hplan
    integer(c_int) :: ierr
#endif

    ! ----- problem size (optionally from argv: nx ny nz). Defaults small for a quick first run. -----
    call get_dims(nx, ny, nz)
    nxc = nx/2 + 1
    nlines_z = nxc*ny            ! Z-FFT batches (single-rank, matches opr_fourier nlines=(imax/2+1)*jmax)
    nlines_x = ny*nz             ! X-FFT batches (single-rank, matches opr_fourier nlines=jmax*kmax)
    all_ok = .true.

    write (*, '(A)') '================ vhipfft : hipFFT vs FFTW bring-up ================'
    write (*, '(A,I5,A,I5,A,I5,A,I8,A,I8)') ' dims nx=', nx, ' ny=', ny, ' nz=', nz, &
        '   nlines_z=', nlines_z, '  nlines_x=', nlines_x
#ifdef USE_HIPFFT
    write (*, '(A)') ' backends: FFTW (reference) + hipFFT (test)'
#else
    write (*, '(A)') ' backends: FFTW only (build with -DUSE_HIPFFT for the hipFFT comparison)'
#endif

    ! ----- offload probe: confirm a target region actually runs on the GPU. -----
    on_dev = .true.
    !$omp target map(tofrom: on_dev)
    on_dev = omp_is_initial_device()
    !$omp end target
    write (*, '(A,L1,A)') ' [ENV] target_on_host=', on_dev, '   (must be F for GPU offload)'

    ! =================================================================================================
    ! STAGE Z : complex-to-complex batched, INTERLEAVED layout (istride=idist swapped vs X).
    ! =================================================================================================
    allocate (cz0(nlines_z*nz), czf_ref(nlines_z*nz), czf_hip(nlines_z*nz), czb(nlines_z*nz))
    call fill_complex(cz0, nlines_z*nz)

    ! FFTW reference: forward then backward(/nz) roundtrip.
    call dfftw_plan_many_dft(pf, 1, nz, nlines_z, cz0, nz, nlines_z, 1, czf_ref, nz, nlines_z, 1, &
                             FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_many_dft(pb, 1, nz, nlines_z, czf_ref, nz, nlines_z, 1, czb, nz, nlines_z, 1, &
                             FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(pf, cz0, czf_ref)
    call dfftw_execute_dft(pb, czf_ref, czb)
    czb = czb/real(nz, c_double)
    call dfftw_destroy_plan(pf); call dfftw_destroy_plan(pb)
    call report_real('Z FFTW roundtrip (fwd->bwd/nz vs input)', cabsmax(czb - cz0), cabsmax(cz0), all_ok)

#ifdef USE_HIPFFT
    n1(1) = nz; ie1(1) = nz; oe1(1) = nz
    ierr = hipfftPlanMany(hplan, 1, n1, ie1, nlines_z, 1, oe1, nlines_z, 1, HIPFFT_Z2Z, nlines_z)
    call chk(ierr, 'hipfftPlanMany Z2Z')
    ierr = hipfftExecZ2Z(hplan, c_loc(cz0), c_loc(czf_hip), HIPFFT_FORWARD)   ! [B2] host ptr via XNACK
    call chk(ierr, 'hipfftExecZ2Z fwd')
    ierr = hipDeviceSynchronize()
    ierr = hipfftDestroy(hplan)
    call report_real('Z hipFFT-fwd vs FFTW-fwd', cabsmax(czf_hip - czf_ref), cabsmax(czf_ref), all_ok)
#endif

    ! =================================================================================================
    ! STAGE X : real-to-complex (fwd) / complex-to-real (bwd) batched, CONTIGUOUS layout.
    ! =================================================================================================
    allocate (rx0(nx*nlines_x), rxb_ref(nx*nlines_x), rxb_hip(nx*nlines_x))
    allocate (cxf_ref(nxc*nlines_x), cxf_hip(nxc*nlines_x))
    call fill_real(rx0, nx*nlines_x)

    ! FFTW reference: r2c forward then c2r backward(/nx) roundtrip.
    call dfftw_plan_many_dft_r2c(pf, 1, nx, nlines_x, rx0, nx, 1, nx, cxf_ref, nxc, 1, nxc, FFTW_ESTIMATE)
    call dfftw_plan_many_dft_c2r(pb, 1, nx, nlines_x, cxf_ref, nxc, 1, nxc, rxb_ref, nx, 1, nx, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(pf, rx0, cxf_ref)
    call dfftw_execute_dft_c2r(pb, cxf_ref, rxb_ref)
    rxb_ref = rxb_ref/real(nx, c_double)
    call dfftw_destroy_plan(pf); call dfftw_destroy_plan(pb)
    call report_real('X FFTW roundtrip (r2c->c2r/nx vs input)', rabsmax(rxb_ref - rx0), rabsmax(rx0), all_ok)

#ifdef USE_HIPFFT
    n1(1) = nx; ie1(1) = nx; oe1(1) = nxc
    ierr = hipfftPlanMany(hplan, 1, n1, ie1, 1, nx, oe1, 1, nxc, HIPFFT_D2Z, nlines_x)
    call chk(ierr, 'hipfftPlanMany D2Z')
    ierr = hipfftExecD2Z(hplan, c_loc(rx0), c_loc(cxf_hip))
    call chk(ierr, 'hipfftExecD2Z')
    ierr = hipDeviceSynchronize()
    ierr = hipfftDestroy(hplan)
    call report_real('X hipFFT-r2c vs FFTW-r2c', cabsmax(cxf_hip - cxf_ref), cabsmax(cxf_ref), all_ok)

    ! c2r roundtrip of the hipFFT forward result -> should recover rx0 after /nx.
    n1(1) = nx; ie1(1) = nxc; oe1(1) = nx
    ierr = hipfftPlanMany(hplan, 1, n1, ie1, 1, nxc, oe1, 1, nx, HIPFFT_Z2D, nlines_x)
    call chk(ierr, 'hipfftPlanMany Z2D')
    ierr = hipfftExecZ2D(hplan, c_loc(cxf_hip), c_loc(rxb_hip))
    call chk(ierr, 'hipfftExecZ2D')
    ierr = hipDeviceSynchronize()
    ierr = hipfftDestroy(hplan)
    rxb_hip = rxb_hip/real(nx, c_double)
    call report_real('X hipFFT roundtrip (r2c->c2r/nx vs input)', rabsmax(rxb_hip - rx0), rabsmax(rx0), all_ok)
#endif

    write (*, '(A)') '------------------------------------------------------------------'
    if (all_ok) then
        write (*, '(A)') ' OVERALL: PASS (all stages at machine precision)'
    else
        write (*, '(A)') ' OVERALL: FAIL (see stages above)'
    end if

contains

    subroutine get_dims(nx, ny, nz)
        integer, intent(out) :: nx, ny, nz
        character(len=32) :: a
        integer :: na
        nx = 32; ny = 16; nz = 32                 ! defaults (small, even); pass "nx ny nz" to override
        na = command_argument_count()
        if (na >= 3) then
            call get_command_argument(1, a); read (a, *) nx
            call get_command_argument(2, a); read (a, *) ny
            call get_command_argument(3, a); read (a, *) nz
        end if
    end subroutine

    subroutine fill_complex(c, n)                  ! deterministic, non-trivial input
        complex(c_double_complex), intent(out) :: c(*)
        integer, intent(in) :: n
        integer :: k
        do k = 1, n
            c(k) = cmplx(sin(0.1_c_double*k) + 0.3_c_double*cos(0.027_c_double*k), &
                         cos(0.05_c_double*k), c_double_complex)
        end do
    end subroutine

    subroutine fill_real(r, n)
        real(c_double), intent(out) :: r(*)
        integer, intent(in) :: n
        integer :: k
        do k = 1, n
            r(k) = sin(0.1_c_double*k) + 0.3_c_double*cos(0.027_c_double*k)
        end do
    end subroutine

    real(c_double) function cabsmax(c) result(m)
        complex(c_double_complex), intent(in) :: c(:)
        m = maxval(abs(c))
    end function

    real(c_double) function rabsmax(r) result(m)
        real(c_double), intent(in) :: r(:)
        m = maxval(abs(r))
    end function

    subroutine report_real(tag, diff, scale, ok)
        character(len=*), intent(in) :: tag
        real(c_double), intent(in) :: diff, scale
        logical, intent(inout) :: ok
        real(c_double) :: norm
        logical :: pass
        norm = diff/max(scale, tiny(scale))
        pass = (norm <= 1.0e-11_c_double)          ! FFT reorder noise floor for double precision
        if (.not. pass) ok = .false.
        write (*, '(A,A40,A,ES11.4,A,ES11.4,A,L1)') ' ', tag, '  max|diff|=', diff, &
            '  norm=', norm, '  PASS=', pass
    end subroutine

#ifdef USE_HIPFFT
    subroutine chk(ierr, tag)
        integer(c_int), intent(in) :: ierr
        character(len=*), intent(in) :: tag
        if (ierr /= 0) write (*, '(A,A,A,I0)') ' !! ', tag, ' returned hipfftResult=', ierr
    end subroutine
#endif

end program VHIPFFT
