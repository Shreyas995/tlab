#include "dns_error.h"

! Reader-side acquire for the node/apu shared-window I-transposes. On MI300A the cross-rank push needs the FULL
! acquire/release protocol on the COARSE-GRAINED (cacheable) MPI_Win_allocate_shared window: writer = release
! (__threadfence_system, via the fused hip_write_with_fence or the separate hip_system_fence) AND reader =
! acquire (hip_invalidate_recv, an L2 invalidate before the GPU unpack). MPI_Win_fence orders the processes on
! the CPU side but does NOT invalidate the reader rank's GPU L2 -> without the acquire the unpack can read a
! stale cached line. So whenever EITHER writer fix is on, the reader invalidate must also be on.
#if defined(TRP_I_SYSFENCE) || defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
#define TRP_I_READER_ACQUIRE 1
#endif

! Circular transposition within directional communicators
module TLabMPI_Transpose
#ifdef USE_APU
    use omp_lib, only: omp_get_default_device, omp_target_alloc, omp_target_free
#endif
    use TLab_Constants, only: lfile, efile, wp, dp, sp, wi, sizeofreal
    use TLab_Memory, only: imax, jmax, kmax, isize_wrk3d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: TLab_Allocate_Real
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc, c_null_ptr, c_associated, c_size_t, c_sizeof, c_int
    ! c_ptr / c_intptr_t are accessible via mpi_f08 (which re-exports iso_c_binding); re-declaring causes ambiguity.
    use TLabMPI_VARS
    implicit none
    private

    public :: TLabMPI_Trp_Initialize
    public :: TLabMPI_Trp_PlanI, TLabMPI_Trp_PlanK
    public :: TLabMPI_Trp_ExecK_Forward, TLabMPI_Trp_ExecK_Backward
    public :: TLabMPI_Trp_ExecI_Forward, TLabMPI_Trp_ExecI_Backward

    type, public :: tmpi_transpose_dt
        ! sequence
        type(MPI_Datatype) :: type_s, type_r                        ! derived send/recv types (kept for ALLTOALL/SENDRECV)
        type(MPI_Datatype) :: base_type                             ! scalar base type for flat MPI_ISEND/IRECV
        integer(wi) :: nlines                                        ! lines per rank per peer
        integer(wi) :: nmax                                          ! elements per line (imax for I-dir, kmax for K-dir)
        integer(wi) :: size3d
        integer(wi), allocatable :: disp_s(:), disp_r(:)            ! send/recv displacements
    end type tmpi_transpose_dt
    type(tmpi_transpose_dt), public :: tmpi_plan_dx                 ! general plans used in derivatives and other operators
    type(tmpi_transpose_dt), public :: tmpi_plan_dz

    ! -----------------------------------------------------------------------
    integer :: trp_mode_i, trp_mode_k                               ! Mode of transposition
    integer, parameter :: TLAB_MPI_TRP_NONE = 0
    integer, parameter :: TLAB_MPI_TRP_ASYNCHRONOUS = 1
    integer, parameter :: TLAB_MPI_TRP_SENDRECV = 2
    integer, parameter :: TLAB_MPI_TRP_ALLTOALL = 3
    integer, parameter :: TLAB_MPI_TRP_APU_DIRECT = 4    ! APU single-node: fused GPU writes into peers' shared windows, MPI_Win_fence barrier
    integer, parameter :: TLAB_MPI_TRP_FABRIC_DIRECT = 5 ! APU multi-node: two-sided MPI (ISEND/IRECV) per peer over a clean MPI_COMM_WORLD-split comm

#ifdef USE_APU
    ! APU_DIRECT state: one shared-memory MPI window per direction. Each rank allocates its
    ! recv buffer as a shared segment (MPI_Win_allocate_shared on the full directional comm,
    ! which is node-local) so every peer GPU-writes directly into it via MPI_Win_shared_query
    ! pointers. Cross-XCD coherency within a node is provided by the MPI_Win_fence barrier
    ! (apudirect is the verified bit-reference on one node).
    integer(wi) :: apu_size_k = 0_wi, apu_size_i = 0_wi
    type(MPI_Win) :: apu_win_k, apu_win_i
    real(dp), pointer :: apu_recv_fptr_k(:) => null()
    real(dp), pointer, contiguous :: apu_recv_fptr_i(:) => null()   ! contiguous: allows whole-array passing to hip_invalidate_recv (complex apudirect reader fix)
    type(c_ptr), allocatable :: apu_peer_cptr_k(:), apu_peer_cptr_i(:) ! per-rank pointers from Win_shared_query
    complex(dp), pointer :: apu_cx_recv_fptr_k(:) => null(), apu_cx_recv_fptr_i(:) => null()  ! complex aliases of the same windows
    ! Contiguous span over all peers' recv segments for the fused single-kernel write.
    ! apu_all_k/i(m*stride + 1 : (m+1)*stride) = rank m's recv buffer (segments are contiguous).
    integer(wi) :: apu_stride_k = 0_wi, apu_stride_i = 0_wi
    real(dp), pointer :: apu_all_k(:) => null()
    real(dp), pointer, contiguous :: apu_all_i(:) => null()   ! contiguous: allows section-passing to hip_write_with_fence (V1 apudirect I-transpose fix)
    ! FABRIC_DIRECT state (multi-node): every peer goes through plain two-sided MPI. Traffic is
    ! routed on a fresh MPI_Comm_split(MPI_COMM_WORLD, ...) comm — NOT a dup of the Cartesian
    ! ims_comm_z/x, because allocating any shared window in the job retroactively taints
    ! Cartesian-lineage comms for two-sided traffic on Cray MPICH. Pack/unpack and the staging
    ! buffer (wrk_mpi_dp) are CPU-side; GPU-written sources are flushed with hipDeviceSynchronize.
    type(MPI_Comm) :: fabric_mpi_comm_k, fabric_mpi_comm_i   ! FABRIC_DIRECT two-sided comm (MPI_COMM_WORLD split)
    ! OPTION 1 (intra-node shared-window K): serve the intra-node K-peers via node-local shared-window
    ! GPU writes (the apudirect mechanism on a true NODE comm) and keep MPI only for inter-node K-peers.
    ! Validated standalone by src/valid/mpi/vmpi_nodewin.f90 (node window contiguous + cross-XCD coherent).
    type(MPI_Comm) :: node_comm_k                          ! NODE comm (split by hostname hash)
    type(MPI_Win)  :: node_win_k                           ! node-local shared recv window (real K)
    integer        :: node_size_k = 0
#ifdef DNS_DEBUG_PROBES
    logical, public :: trp_dbg_fft = .false.               ! DEBUG: gate internal X-FFT-region sentinels (set by OPR_Fourier_X_Backward)
#else
    logical, parameter, public :: trp_dbg_fft = .false.    ! production: PARAMETER -> the compiler dead-code-eliminates every `if (trp_dbg_fft)` probe branch
#endif
    logical        :: use_node_win_k = .false.             ! true iff the node window came up contiguous
    real(dp), pointer :: node_recv_fptr_k(:) => null()     ! our own node-window segment
    real(dp), pointer :: node_all_k(:) => null()           ! fused span over the whole node window
    complex(dp), pointer :: node_cx_recv_fptr_k(:) => null() ! complex view of our segment (Poisson cx-K)
    complex(dp), pointer :: node_cx_all_k(:) => null()       ! complex view of the fused node-window span
    logical, allocatable :: is_intra_k(:)                  ! K-peer m on our node?
    integer, allocatable :: node_lrank_k(:)                ! K-peer m -> node-comm local rank
    ! I-direction node-window (same mechanism; for npro_i=6 all 6 I-peers are intra-node, no inter leg).
    type(MPI_Comm) :: node_comm_i
    type(MPI_Win)  :: node_win_i
    integer        :: node_size_i = 0
    logical        :: use_node_win_i = .false.
    real(dp), pointer, contiguous :: node_recv_fptr_i(:) => null()   ! contiguous: c_f_pointer'd window segment; allows element-passing to hip_invalidate_recv
    real(dp), pointer, contiguous :: node_all_i(:) => null()   ! contiguous: c_f_pointer'd window; allows section-passing to hip_write_with_fence
    complex(dp), pointer :: node_cx_recv_fptr_i(:) => null() ! complex view of our segment (Poisson cx-I)
    complex(dp), pointer :: node_cx_all_i(:) => null()       ! complex view of the fused node-window span
    logical, allocatable :: is_intra_i(:)
    integer, allocatable :: node_lrank_i(:)
#endif

    integer(wi) :: trp_sizBlock_i, trp_sizBlock_k                   ! explicit send/recv: group sizes of send/recv messages
    integer(wi), allocatable :: maps_send_i(:), maps_recv_i(:)      ! PE maps to use explicit send/recv
    integer(wi), allocatable :: maps_send_k(:), maps_recv_k(:)
    type(MPI_Datatype), allocatable :: types_send(:), types_recv(:) ! alltoallw
    integer, allocatable :: counts(:)

    type(MPI_Datatype) :: trp_datatype_i, trp_datatype_k            ! Transposition in double or single precision

    ! sp work buffer: three sections of size3d sp elements each.
    !   [1..size]          = a_wrk: dp→sp input copy
    !   [size+1..2*size]   = b_wrk: sp recv/output buffer
    !   [2*size+1..3*size] = c_wrk: flat staging buffer for strided send/recv
    ! Allocated via omp_target_alloc so it resides in device-accessible memory.
    ! Only allocated when ASYNCHRONOUS mode uses single precision.
    type(c_ptr) :: wrk_mpi_cptr = c_null_ptr
    real(sp), pointer :: wrk_mpi_fptr(:) => null()
    real(sp), pointer :: a_wrk(:) => null(), b_wrk(:) => null(), c_wrk(:) => null()

    ! dp/complex staging buffer for flat MPI send/recv.
    ! Size = imax*jmax*kmax: covers real(dp) and complex(dp) paths.
    ! Allocated with standard Fortran allocate; on APU unified memory this is device-accessible.
    real(dp), allocatable, target :: wrk_mpi_dp(:)
    real(dp), pointer, contiguous    :: c_wrk_dp(:) => null()   ! contiguous: MPI must see the real buffer, not a temporary
    complex(dp), pointer, contiguous :: c_wrk_cx(:) => null()
    complex(dp), pointer, contiguous :: c_wrk_send_cx(:) => null()   ! SEND staging (2nd half of wrk_mpi_dp): GPU-packed send buffer for the complex K-backward inter leg (must be GPU-written, not the CPU-FFTW b, for GPU-aware MPI)
    type(MPI_Status) status(128)
    type(MPI_Request) request(128)

    interface TLabMPI_Trp_ExecK_Forward
        module procedure TLabMPI_Trp_ExecK_Forward_Real, TLabMPI_Trp_ExecK_Forward_Complex
    end interface TLabMPI_Trp_ExecK_Forward
    interface TLabMPI_Trp_ExecK_Backward
        module procedure TLabMPI_Trp_ExecK_Backward_Real, TLabMPI_Trp_ExecK_Backward_Complex
    end interface TLabMPI_Trp_ExecK_Backward

    interface TLabMPI_Trp_ExecI_Forward
        module procedure TLabMPI_Trp_ExecI_Forward_Real, TLabMPI_Trp_ExecI_Forward_Complex
    end interface TLabMPI_Trp_ExecI_Forward
    interface TLabMPI_Trp_ExecI_Backward
        module procedure TLabMPI_Trp_ExecI_Backward_Real, TLabMPI_Trp_ExecI_Backward_Complex
    end interface TLabMPI_Trp_ExecI_Backward

#ifdef USE_APU
    interface
        subroutine hip_write_with_fence(src, dst, n) bind(C, name='hip_write_with_fence')
            use iso_c_binding
            real(c_double), intent(in)  :: src(*)
            real(c_double), intent(out) :: dst(*)
            integer(c_int), value       :: n
        end subroutine

        ! System-scope L2 write-back (__threadfence_system) — commits the preceding node-window push to MALL
        ! so other ranks/XCDs see the fresh data (device-scope hipDeviceSynchronize does NOT). Used by Fix 2.
        subroutine hip_system_fence() bind(C, name='hip_system_fence')
        end subroutine

        ! Reader-side L2 refresh — volatile-loads buf(1:n) from MALL before the GPU unpack reads the recv
        ! window, so a stale L2 copy on THIS rank is not used. Fix 2 reader-side half (the writer fence alone
        ! left a ~6-order residual = COARSE_GRAINED memory; the reader's L2 isn't invalidated by the writer).
        subroutine hip_invalidate_recv(buf, n) bind(C, name='hip_invalidate_recv')
            use iso_c_binding
            real(c_double), intent(inout) :: buf(*)
            integer(c_int), value         :: n
        end subroutine

        ! V2 (TRP_I_MEMCPY): blocking, system-coherent push of n doubles src -> dst via
        ! hipMemcpy(...hipMemcpyDefault). Synchronous + committed to MALL on return, so no writer-side
        ! __threadfence_system is needed (it replaces the fused-fence kernel for the node-window K push).
        ! The reader still acquires via hip_invalidate_recv. dst/src are first elements of contiguous arrays.
        subroutine hip_memcpy_push(dst, src, n) bind(C, name='hip_memcpy_push')
            use iso_c_binding
            real(c_double), intent(out) :: dst(*)
            real(c_double), intent(in)  :: src(*)
            integer(c_int), value       :: n
        end subroutine

        function hipHostRegister(ptr, sz, flags) bind(C, name='hipHostRegister') result(ierr)
            use iso_c_binding
            integer(c_int) :: ierr
            type(c_ptr), value       :: ptr
            integer(c_size_t), value :: sz
            integer(c_int), value    :: flags
        end function

        function hipDeviceSynchronize() bind(C, name='hipDeviceSynchronize') result(ierr)
            use iso_c_binding
            integer(c_int) :: ierr
        end function hipDeviceSynchronize
    end interface
#endif

contains

    ! ######################################################################
    ! ######################################################################
    subroutine TLabMPI_Trp_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        integer(c_size_t) :: bytes_sp
        integer :: dev_id
        integer(wi) :: total_elements
        integer(MPI_ADDRESS_KIND) :: win_query_size
        integer :: win_disp_unit
        type(c_ptr) :: win_baseptr
#ifdef USE_APU
        type(MPI_Info)  :: win_info_contig   ! forces contiguous shared-window segments (alloc_shared_noncontig=false)
#endif
        ! -----------------------------------------------------------------------
        integer(wi) ip, npage, dummy
        character(len=32) bakfile, block
        character(len=512) sRes, line
        character*64 lstr

        ! #######################################################################
        ! Read data
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Parallel'

        call ScanFile_Char(bakfile, inifile, block, 'TransposeModeI', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') &
            call ScanFile_Char(bakfile, inifile, 'Main', 'ComModeITranspose', 'asynchronous', sRes)
        if (trim(adjustl(sRes)) == 'none') then; trp_mode_i = TLAB_MPI_TRP_NONE
        elseif (trim(adjustl(sRes)) == 'asynchronous') then; trp_mode_i = TLAB_MPI_TRP_ASYNCHRONOUS
        elseif (trim(adjustl(sRes)) == 'sendrecv') then; trp_mode_i = TLAB_MPI_TRP_SENDRECV
        elseif (trim(adjustl(sRes)) == 'alltoall') then; trp_mode_i = TLAB_MPI_TRP_ALLTOALL
        elseif (trim(adjustl(sRes)) == 'apudirect') then
#ifdef USE_APU
            trp_mode_i = TLAB_MPI_TRP_APU_DIRECT
#else
            call TLab_Write_ASCII(efile, __FILE__//'. TransposeModeI=apudirect requires USE_APU.')
            call TLab_Stop(DNS_ERROR_OPTION)
#endif
        elseif (trim(adjustl(sRes)) == 'fabricdirect') then
#ifdef USE_APU
            trp_mode_i = TLAB_MPI_TRP_FABRIC_DIRECT
#else
            call TLab_Write_ASCII(efile, __FILE__//'. TransposeModeI=fabricdirect requires USE_APU.')
            call TLab_Stop(DNS_ERROR_OPTION)
#endif
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeModeI option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeModeK', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') &
            call ScanFile_Char(bakfile, inifile, 'Main', 'ComModeKTranspose', 'asynchronous', sRes)
        if (trim(adjustl(sRes)) == 'none') then; trp_mode_k = TLAB_MPI_TRP_NONE
        elseif (trim(adjustl(sRes)) == 'asynchronous') then; trp_mode_k = TLAB_MPI_TRP_ASYNCHRONOUS
        elseif (trim(adjustl(sRes)) == 'sendrecv') then; trp_mode_k = TLAB_MPI_TRP_SENDRECV
        elseif (trim(adjustl(sRes)) == 'alltoall') then; trp_mode_k = TLAB_MPI_TRP_ALLTOALL
        elseif (trim(adjustl(sRes)) == 'apudirect') then
#ifdef USE_APU
            trp_mode_k = TLAB_MPI_TRP_APU_DIRECT
#else
            call TLab_Write_ASCII(efile, __FILE__//'. TransposeModeK=apudirect requires USE_APU.')
            call TLab_Stop(DNS_ERROR_OPTION)
#endif
        elseif (trim(adjustl(sRes)) == 'fabricdirect') then
#ifdef USE_APU
            trp_mode_k = TLAB_MPI_TRP_FABRIC_DIRECT
#else
            call TLab_Write_ASCII(efile, __FILE__//'. TransposeModeK=fabricdirect requires USE_APU.')
            call TLab_Stop(DNS_ERROR_OPTION)
#endif
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeModeK option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeTypeK', 'Double', sRes)
        if (trim(adjustl(sRes)) == 'double') then; trp_datatype_k = MPI_REAL8
        elseif (trim(adjustl(sRes)) == 'single') then; trp_datatype_k = MPI_REAL4
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeTypeK.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeTypeI', 'Double', sRes)
        if (trim(adjustl(sRes)) == 'double') then; trp_datatype_i = MPI_REAL8
        elseif (trim(adjustl(sRes)) == 'single') then; trp_datatype_i = MPI_REAL4
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeTypeI.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        ! #######################################################################
        ! Initialize

        ! Size of communication in explicit send/recv
#ifdef HLRS_HAWK
        ! On hawk, we tested that 192 yields optimum performance;
        ! Blocking will thus only take effect in very large cases
        trp_sizBlock_k = 192
        trp_sizBlock_i = 384
#else
        ! We assume that this will help to release some of the very heavy
        ! network load in transpositions on most systems
        trp_sizBlock_k = 64
        trp_sizBlock_i = 128
        ! trp_sizBlock_k=1e5   -- would essentially switch off the blocking
#endif

        if (ims_npro_i > trp_sizBlock_i) then
            write (line, *) trp_sizBlock_i
            line = 'Using blocking of '//trim(adjustl(line))//' in TLabMPI_TRP<F,B>_I'
            call TLab_Write_ASCII(lfile, line)
        end if

        if (ims_npro_k > trp_sizBlock_k) then
            write (line, *) trp_sizBlock_k
            line = 'Using blocking of '//trim(adjustl(line))//' in TLabMPI_TRP<F,B>_K'
            call TLab_Write_ASCII(lfile, line)
        end if

        ! -----------------------------------------------------------------------
        ! local PE mappings for explicit send/recv (only needed when multiple ranks per direction)
        if (ims_npro_i > 1) then
            allocate (maps_send_i(ims_npro_i))
            allocate (maps_recv_i(ims_npro_i))
            do ip = 0, ims_npro_i - 1
                maps_send_i(ip + 1) = ip
                maps_recv_i(ip + 1) = mod(ims_npro_i - ip, ims_npro_i)
            end do
            maps_send_i = cshift(maps_send_i, ims_pro_i)
            maps_recv_i = cshift(maps_recv_i, -ims_pro_i)
        end if

        if (ims_npro_k > 1) then
            allocate (maps_send_k(ims_npro_k))
            allocate (maps_recv_k(ims_npro_k))
            do ip = 0, ims_npro_k - 1
                maps_send_k(ip + 1) = ip
                maps_recv_k(ip + 1) = mod(ims_npro_k - ip, ims_npro_k)
            end do
            maps_send_k = cshift(maps_send_k, ims_pro_k)
            maps_recv_k = cshift(maps_recv_k, -ims_pro_k)
        end if

        ! -----------------------------------------------------------------------
        ! ALLTOALLW scratch arrays: only allocated when at least one direction uses ALLTOALL mode.
        if (trp_mode_i == TLAB_MPI_TRP_ALLTOALL .or. trp_mode_k == TLAB_MPI_TRP_ALLTOALL) then
            allocate (counts(max(ims_npro_i, ims_npro_j, ims_npro_k)))
            allocate (types_send(max(ims_npro_i, ims_npro_j, ims_npro_k)))
            allocate (types_recv(max(ims_npro_i, ims_npro_j, ims_npro_k)))
            counts(:) = 1
        end if

        ! -----------------------------------------------------------------------
        ! Flat-MPI staging buffers — allocated only for ASYNCHRONOUS mode.
        !
        ! wrk_mpi_fptr (sp, 3 x size3d): dp<->sp conversion buffers + flat staging.
        !   Allocated via omp_target_alloc to guarantee device-accessible memory on APU.
        !   Only needed when ASYNCHRONOUS transport uses single precision (MPI_REAL4).
#ifdef USE_APU
        if ((trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS .and. trp_datatype_i == MPI_REAL4) .or. &
            (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS .and. trp_datatype_k == MPI_REAL4)) then
            dev_id = omp_get_default_device()
            total_elements = int(imax, wi) * int(jmax, wi) * int(kmax, wi)
            bytes_sp = 3_c_size_t * total_elements * c_sizeof(1.0_sp)
            wrk_mpi_cptr = omp_target_alloc(bytes_sp, dev_id)
            if (.not. c_associated(wrk_mpi_cptr)) then
                call TLab_Write_ASCII(efile, __FILE__//'. omp_target_alloc failed for sp staging buffer.')
                call TLab_Stop(DNS_ERROR_ALLOC)
            end if
            call c_f_pointer(wrk_mpi_cptr, wrk_mpi_fptr, [3 * total_elements])
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated sp flat-MPI staging buffer (3 x size3d).')
        end if
#endif

        ! wrk_mpi_dp (dp, 2 x size3d): staging buffer for dp-real and complex ASYNCHRONOUS paths.
        !   2× size3d is required: the complex Fourier X plan uses nmax = imax/2+1, so its
        !   size3d = jmax*kmax*(imax/2+1) complex elements = jmax*kmax*(imax+2) real elements,
        !   which exceeds 1×size3d = imax*jmax*kmax real elements. 2× covers all cases.
        !   Standard Fortran allocate: on APU unified memory this is device-accessible.
        !   Not needed for SENDRECV or ALLTOALL which use MPI_TYPE_VECTOR via the kernels.
        ! APU_DIRECT also needs wrk_mpi_dp: complex exec functions fall back to ASYNC behaviour
        ! (no mode-4 path for complex — complex transposes are infrequent and dp already).
        if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS .or. trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS .or. &
            trp_mode_i == TLAB_MPI_TRP_APU_DIRECT   .or. trp_mode_k == TLAB_MPI_TRP_APU_DIRECT   .or. &
            trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT .or. trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
            allocate (wrk_mpi_dp(4*imax*jmax*kmax))   ! 4x (was 2x): 1st half = recv c_wrk_cx, 2nd half = GPU-packed SEND staging c_wrk_send_cx (complex K-backward inter leg)
            wrk_mpi_dp = 0.0_dp   ! INIT: not zeroing left UNFILLED staging slots as garbage (6.46e208). A cross-rank
                                  ! write that does not land coherently then reads that garbage -> overflow. Zero = bounded.
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated dp flat-MPI staging buffer (4 x size3d).')
        end if

#ifdef USE_APU
        ! Force CONTIGUOUS shared-window segments. APU_DIRECT treats all peer segments as one
        ! contiguous span (apu_all_X(m*stride + ...)). Cray MPICH defaults to a non-contiguous
        ! layout on multi-node (per-XCD VA blocks / 2x-segsize guard regions), which breaks that
        ! span -> GPU SIGSEGV. alloc_shared_noncontig=false requests one back-to-back block.
        call MPI_Info_create(win_info_contig, ims_err)
        call MPI_Info_set(win_info_contig, 'alloc_shared_noncontig', 'false', ims_err)

        ! APU direct mode: allocate per-rank shared-memory MPI windows so every rank in the
        ! communicator gets a pointer into every other rank's recv buffer segment.
        ! MPI_Win_allocate_shared maps all segments into the calling process's address space
        ! (via sysv/posix shared memory), solving the cross-process virtual-address problem
        ! that broke the previous omp_target_alloc + integer(8) address sharing approach.
        ! Layout in recv buffer: slot r*chunk holds data written by rank r.
        if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT .and. ims_npro_k > 1) then
            ! Buffer sized for complex FFT plans: (imax/2+1)*jmax*kmax cx = (imax+2)*jmax*kmax dp.
            apu_size_k = int(imax + 2, wi)*int(jmax, wi)*int(kmax, wi)
            call MPI_Win_allocate_shared(int(apu_size_k, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
                                         int(c_sizeof(1.0_dp)), win_info_contig, ims_comm_z, win_baseptr, apu_win_k, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU K recv buffer.')
                call TLab_Stop(DNS_ERROR_ALLOC)
            end if
            call c_f_pointer(win_baseptr, apu_recv_fptr_k, [apu_size_k])
            call c_f_pointer(win_baseptr, apu_cx_recv_fptr_k, [apu_size_k/2])
            apu_recv_fptr_k = 0.0_dp   ! INIT: a cross-rank push that misses a slot then reads 0 (bounded), not allocation-garbage
            allocate (apu_peer_cptr_k(0:ims_npro_k - 1))
            do ip = 0, ims_npro_k - 1
                call MPI_Win_shared_query(apu_win_k, ip, win_query_size, win_disp_unit, apu_peer_cptr_k(ip), ims_err)
            end do
            ! apu_all_k spans all peers' contiguous segments (MPI-3 shared windows are always contiguous).
            apu_stride_k = apu_size_k
            call c_f_pointer(apu_peer_cptr_k(0), apu_all_k, [apu_stride_k*ims_npro_k])
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated K APU direct recv buffer.')
        end if
        if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT .and. ims_npro_i > 1) then
            ! Buffer sized for complex FFT plans: (imax/2+1)*jmax*kmax cx = (imax+2)*jmax*kmax dp.
            apu_size_i = int(imax + 2, wi)*int(jmax, wi)*int(kmax, wi)
            call MPI_Win_allocate_shared(int(apu_size_i, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
                                         int(c_sizeof(1.0_dp)), win_info_contig, ims_comm_x, win_baseptr, apu_win_i, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU I recv buffer.')
                call TLab_Stop(DNS_ERROR_ALLOC)
            end if
            call c_f_pointer(win_baseptr, apu_recv_fptr_i, [apu_size_i])
            call c_f_pointer(win_baseptr, apu_cx_recv_fptr_i, [apu_size_i/2])
            apu_recv_fptr_i = 0.0_dp   ! INIT (see apu_recv_fptr_k)
            allocate (apu_peer_cptr_i(0:ims_npro_i - 1))
            do ip = 0, ims_npro_i - 1
                call MPI_Win_shared_query(apu_win_i, ip, win_query_size, win_disp_unit, apu_peer_cptr_i(ip), ims_err)
            end do
            ! apu_all_i spans all peers' segments as one contiguous block. Validity now relies on
            ! alloc_shared_noncontig=false (win_info_contig) forcing Cray MPICH to lay the segments
            ! out back-to-back; otherwise the per-XCD/guard-region layout breaks this span on multi-node.
            apu_stride_i = apu_size_i
            call c_f_pointer(apu_peer_cptr_i(0), apu_all_i, [apu_stride_i*ims_npro_i])
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated I APU direct recv buffer.')
        end if
        ! win_info_contig has been copied into any windows that used it; safe to free now.
        call MPI_Info_free(win_info_contig, ims_err)
        ! FABRIC_DIRECT: route ALL two-sided traffic on fresh MPI_Comm_split(MPI_COMM_WORLD,...)
        ! comms. On Cray MPICH, allocating any MPI_Win_allocate_shared in the job retroactively
        ! taints the Cartesian ims_comm_z/x (and any dup of them) for two-sided traffic — confirmed
        ! by standalone commtest on Hunter (2026-05-29). MPI_Comm_split of MPI_COMM_WORLD has no
        ! Cartesian/shmem lineage, so it stays clean.
        ! color=ims_pro_k groups same-K-row ranks (I-direction peers, local rank = ims_pro_i);
        ! color=ims_pro_i groups same-I-column ranks (K-direction peers, local rank = ims_pro_k).
        if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT .and. ims_npro_k > 1) then
            call MPI_Comm_split(MPI_COMM_WORLD, ims_pro_i, ims_pro_k, fabric_mpi_comm_k, ims_err)
            ! Option 1: NODE-local shared window for the intra-node K-peers (apudirect mechanism on a
            ! true NODE comm built by hostname — MPI_Comm_split_type splits per-XCD on MI300A, not per
            ! node). If the window does not come up contiguous (e.g. alongside apu_win_i), use_node_win_k
            ! stays .false. and the exec falls back to the all-MPI path. Validated: src/valid/mpi/vmpi_nodewin.f90.
            block
                character(len=MPI_MAX_PROCESSOR_NAME) :: hname
                integer :: hlen, hh, ic, mpk, my_lrank, ncontig
                integer :: peer_world(0:ims_npro_k - 1)
                type(MPI_Group) :: world_grp, node_grp
                type(MPI_Info)  :: node_info
                type(c_ptr)     :: seg_cptr
                integer(MPI_ADDRESS_KIND) :: va0, va, segb
                call MPI_Get_processor_name(hname, hlen, ims_err)
                hh = 0
                do ic = 1, hlen
                    hh = mod(hh*31 + ichar(hname(ic:ic)), 1000000007)
                end do
                call MPI_Comm_split(MPI_COMM_WORLD, hh, ims_pro, node_comm_k, ims_err)
                call MPI_Comm_rank(node_comm_k, my_lrank, ims_err)
                call MPI_Comm_size(node_comm_k, node_size_k, ims_err)
                apu_size_k = int(imax + 2, wi)*int(jmax, wi)*int(kmax, wi)        ! real-K recv buffer per rank
                segb = int(apu_size_k, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND)
                call MPI_Info_create(node_info, ims_err)
                call MPI_Info_set(node_info, 'alloc_shared_noncontig', 'false', ims_err)
                call MPI_Win_allocate_shared(segb, int(c_sizeof(1.0_dp)), node_info, node_comm_k, &
                                             win_baseptr, node_win_k, ims_err)
                call MPI_Info_free(node_info, ims_err)
                if (ims_err /= MPI_SUCCESS) then
                    call TLab_Write_ASCII(efile, __FILE__//'. node K MPI_Win_allocate_shared failed.')
                    call TLab_Stop(DNS_ERROR_ALLOC)
                end if
                ! own segment + fused span over the whole node window, then a contiguity check.
                call MPI_Win_shared_query(node_win_k, my_lrank, win_query_size, win_disp_unit, seg_cptr, ims_err)
                call c_f_pointer(seg_cptr, node_recv_fptr_k, [apu_size_k])
                call c_f_pointer(seg_cptr, node_cx_recv_fptr_k, [apu_size_k/2])   ! complex view (Poisson cx-K)
                node_recv_fptr_k = 0.0_dp   ! INIT (see apu_recv_fptr_k): un-landed cross-rank push reads 0, not garbage
                call MPI_Win_shared_query(node_win_k, 0, win_query_size, win_disp_unit, win_baseptr, ims_err)
                call c_f_pointer(win_baseptr, node_all_k, [int(apu_size_k, 8)*int(node_size_k, 8)])
                call c_f_pointer(win_baseptr, node_cx_all_k, [int(apu_size_k, 8)*int(node_size_k, 8)/2])
                va0 = transfer(win_baseptr, va0)
                ncontig = 0
                do ic = 0, node_size_k - 1
                    call MPI_Win_shared_query(node_win_k, ic, win_query_size, win_disp_unit, seg_cptr, ims_err)
                    va = transfer(seg_cptr, va)
                    if (va - va0 /= int(ic, MPI_ADDRESS_KIND)*segb) ncontig = ncontig + 1
                end do
                use_node_win_k = (ncontig == 0)
                ! classify K-peers: world rank = m*npro_i + pro_i; node-comm local rank via translate.
                allocate (is_intra_k(0:ims_npro_k - 1), node_lrank_k(0:ims_npro_k - 1))
                do mpk = 0, ims_npro_k - 1
                    peer_world(mpk) = mpk*ims_npro_i + ims_pro_i
                end do
                call MPI_Comm_group(MPI_COMM_WORLD, world_grp, ims_err)
                call MPI_Comm_group(node_comm_k, node_grp, ims_err)
                call MPI_Group_translate_ranks(world_grp, ims_npro_k, peer_world, node_grp, node_lrank_k, ims_err)
                call MPI_Group_free(world_grp, ims_err)
                call MPI_Group_free(node_grp, ims_err)
                do mpk = 0, ims_npro_k - 1
                    is_intra_k(mpk) = (node_lrank_k(mpk) /= MPI_UNDEFINED)
                end do
            end block
            if (use_node_win_k) then
                call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: K FABRIC_DIRECT + node-local intra-window ready (real + complex).')
            else
                call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: K FABRIC_DIRECT (node window non-contiguous; all-MPI fallback).')
            end if
        end if
        if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT .and. ims_npro_i > 1) then
            call MPI_Comm_split(MPI_COMM_WORLD, ims_pro_k, ims_pro_i, fabric_mpi_comm_i, ims_err)
            ! Option 1: NODE-local shared window for the intra-node I-peers (mirror of the K block).
            ! For npro_i=6 all I-peers are intra-node (one XCD) -> all go via the window, no MPI.
            block
                character(len=MPI_MAX_PROCESSOR_NAME) :: hname
                integer :: hlen, hh, ic, mpi_, my_lrank, ncontig
                integer :: peer_world(0:ims_npro_i - 1)
                type(MPI_Group) :: world_grp, node_grp
                type(MPI_Info)  :: node_info
                type(c_ptr)     :: seg_cptr
                integer(MPI_ADDRESS_KIND) :: va0, va, segb
                call MPI_Get_processor_name(hname, hlen, ims_err)
                hh = 0
                do ic = 1, hlen
                    hh = mod(hh*31 + ichar(hname(ic:ic)), 1000000007)
                end do
                call MPI_Comm_split(MPI_COMM_WORLD, hh, ims_pro, node_comm_i, ims_err)
                call MPI_Comm_rank(node_comm_i, my_lrank, ims_err)
                call MPI_Comm_size(node_comm_i, node_size_i, ims_err)
                apu_size_i = int(imax + 2, wi)*int(jmax, wi)*int(kmax, wi)        ! real-I recv buffer per rank
                segb = int(apu_size_i, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND)
                call MPI_Info_create(node_info, ims_err)
                call MPI_Info_set(node_info, 'alloc_shared_noncontig', 'false', ims_err)
                call MPI_Win_allocate_shared(segb, int(c_sizeof(1.0_dp)), node_info, node_comm_i, &
                                             win_baseptr, node_win_i, ims_err)
                call MPI_Info_free(node_info, ims_err)
                if (ims_err /= MPI_SUCCESS) then
                    call TLab_Write_ASCII(efile, __FILE__//'. node I MPI_Win_allocate_shared failed.')
                    call TLab_Stop(DNS_ERROR_ALLOC)
                end if
                call MPI_Win_shared_query(node_win_i, my_lrank, win_query_size, win_disp_unit, seg_cptr, ims_err)
                call c_f_pointer(seg_cptr, node_recv_fptr_i, [apu_size_i])
                call c_f_pointer(seg_cptr, node_cx_recv_fptr_i, [apu_size_i/2])   ! complex view (Poisson cx-I)
                node_recv_fptr_i = 0.0_dp   ! INIT (see apu_recv_fptr_k)
                call MPI_Win_shared_query(node_win_i, 0, win_query_size, win_disp_unit, win_baseptr, ims_err)
                call c_f_pointer(win_baseptr, node_all_i, [int(apu_size_i, 8)*int(node_size_i, 8)])
                call c_f_pointer(win_baseptr, node_cx_all_i, [int(apu_size_i, 8)*int(node_size_i, 8)/2])
                va0 = transfer(win_baseptr, va0)
                ncontig = 0
                do ic = 0, node_size_i - 1
                    call MPI_Win_shared_query(node_win_i, ic, win_query_size, win_disp_unit, seg_cptr, ims_err)
                    va = transfer(seg_cptr, va)
                    if (va - va0 /= int(ic, MPI_ADDRESS_KIND)*segb) ncontig = ncontig + 1
                end do
                use_node_win_i = (ncontig == 0)
#ifdef TRP_I_FORCE_MPI
                ! TEMP STABILIZER: route all 4 I-transposes onto the robust all-MPI fallback. The node-window
                ! cross-rank GPU push delivers stale data (device-scope flush, not system-scope); the fallback
                ! flushes GPU->HBM and moves data via reliable intra-node MPI. K-direction is untouched.
                use_node_win_i = .false.
#endif
                allocate (is_intra_i(0:ims_npro_i - 1), node_lrank_i(0:ims_npro_i - 1))
                do mpi_ = 0, ims_npro_i - 1
                    peer_world(mpi_) = ims_pro_k*ims_npro_i + mpi_        ! I-peer mpi_: pro_k fixed, pro_i = mpi_
                end do
                call MPI_Comm_group(MPI_COMM_WORLD, world_grp, ims_err)
                call MPI_Comm_group(node_comm_i, node_grp, ims_err)
                call MPI_Group_translate_ranks(world_grp, ims_npro_i, peer_world, node_grp, node_lrank_i, ims_err)
                call MPI_Group_free(world_grp, ims_err)
                call MPI_Group_free(node_grp, ims_err)
                do mpi_ = 0, ims_npro_i - 1
                    is_intra_i(mpi_) = (node_lrank_i(mpi_) /= MPI_UNDEFINED)
                end do
            end block
            if (use_node_win_i) then
                call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: I FABRIC_DIRECT + node-local intra-window ready (real + complex).')
            else
                call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: I FABRIC_DIRECT (node window non-contiguous; all-MPI fallback).')
            end if
        end if
#ifdef TRP_CX_MPI
        ! TRP_CX_MPI: route the apudirect COMPLEX transposes through the all-MPI fallback on a CLEAN comm
        ! (the apudirect REAL transposes stay on apu_win). Build fabric_mpi_comm_i/k (fresh MPI_COMM_WORLD
        ! splits -> no Cartesian-comm taint) and set is_intra_*=.false. so the reused fabricdirect complex
        ! path takes its all-MPI branch (use_node_win_* is .false. in apudirect).
        if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT .and. ims_npro_k > 1) then
            call MPI_Comm_split(MPI_COMM_WORLD, ims_pro_i, ims_pro_k, fabric_mpi_comm_k, ims_err)
            if (.not. allocated(is_intra_k)) allocate (is_intra_k(0:ims_npro_k - 1))
            is_intra_k = .false.
        end if
        if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT .and. ims_npro_i > 1) then
            call MPI_Comm_split(MPI_COMM_WORLD, ims_pro_k, ims_pro_i, fabric_mpi_comm_i, ims_err)
            if (.not. allocated(is_intra_i)) allocate (is_intra_i(0:ims_npro_i - 1))
            is_intra_i = .false.
        end if
#endif
#endif

        ! -----------------------------------------------------------------------
        ! Create basic transposition plans used for partial X and partial Z; could be in another module...
        if (ims_npro_i > 1) then
            npage = kmax*jmax
            tmpi_plan_dx = TLabMPI_Trp_PlanI(imax, npage, message='Ox derivatives.')
        end if

        if (ims_npro_k > 1) then
            npage = imax*jmax
            tmpi_plan_dz = TLabMPI_Trp_PlanK(kmax, npage, message='Oz derivatives.')
        end if

        return
    end subroutine TLabMPI_Trp_Initialize

    ! ######################################################################
    ! ######################################################################
    ! Pointers and types for transposition across processors
    function TLabMPI_Trp_PlanI(nmax, npage, locStride, locType, message) result(trp_plan)
        integer(wi), intent(in) :: npage, nmax
        integer(wi), intent(in), optional :: locStride
        type(MPI_Datatype), intent(in), optional :: locType
        character(len=*), intent(in), optional :: message
        type(tmpi_transpose_dt) :: trp_plan

        ! -----------------------------------------------------------------------
        integer(wi) i
        type(MPI_Datatype) :: datatype
        integer block_count, block_length, stride
        integer ims_ss, ims_rs
        character*64 str, line

        ! #######################################################################
        if (present(message)) &
            call TLab_Write_ASCII(lfile, 'Creating derived MPI types for '//trim(adjustl(message)))

        if (mod(npage, ims_npro_i) == 0) then
            trp_plan%nlines = npage/ims_npro_i
            allocate (trp_plan%disp_s(ims_npro_i), trp_plan%disp_r(ims_npro_i))
            trp_plan%size3d = npage*nmax
        else
            call TLab_Write_ASCII(efile, 'TLabMPI_TypeI_Create. Ratio npage/npro not an integer.')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        block_count = trp_plan%nlines
        block_length = nmax

        ! Calculate array displacements in Forward Send/Receive
        trp_plan%disp_s(1) = 0
        trp_plan%disp_r(1) = 0
        do i = 2, ims_npro_i
            trp_plan%disp_s(i) = trp_plan%disp_s(i - 1) + block_length*block_count
            trp_plan%disp_r(i) = trp_plan%disp_r(i - 1) + block_length
        end do

        ! #######################################################################
        if (present(locType)) then
            datatype = locType
        else
            datatype = trp_datatype_i
        end if

        trp_plan%nmax      = nmax
        trp_plan%base_type = datatype

        stride = block_length       ! stride = block_length because things are together
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, trp_plan%type_s, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_s, ims_err)

        stride = nmax*ims_npro_i    ! stride is a multiple of nmax_total=nmax*ims_npro_i
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, trp_plan%type_r, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_r, ims_err)

        ! -----------------------------------------------------------------------
        call MPI_TYPE_SIZE(trp_plan%type_s, ims_ss, ims_err)
        call MPI_TYPE_SIZE(trp_plan%type_r, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            write (str, *) ims_ss; write (line, *) ims_rs
            line = 'Send size '//trim(adjustl(str))//'differs from recv size '//trim(adjustl(line))
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_MPITYPECHECK)
        end if

        return
    end function TLabMPI_Trp_PlanI

    !########################################################################
    !########################################################################
    function TLabMPI_Trp_PlanK(nmax, npage, locStride, locType, message) result(trp_plan)
        integer(wi), intent(in) :: npage, nmax
        integer(wi), intent(in), optional :: locStride
        type(MPI_Datatype), intent(in), optional :: locType
        character(len=*), intent(in), optional :: message
        type(tmpi_transpose_dt) :: trp_plan

        ! -----------------------------------------------------------------------
        integer(wi) i
        type(MPI_Datatype) :: datatype
        integer block_count, block_length, stride
        integer ims_ss, ims_rs
        character*64 str, line

        ! #######################################################################
        if (present(message)) &
            call TLab_Write_ASCII(lfile, 'Creating derived MPI types for '//trim(adjustl(message)))

        if (mod(npage, ims_npro_k) == 0) then
            trp_plan%nlines = npage/ims_npro_k
            allocate (trp_plan%disp_s(ims_npro_k), trp_plan%disp_r(ims_npro_k))
            trp_plan%size3d = npage*nmax
        else
            call TLab_Write_ASCII(efile, 'TLabMPI_TypeI_Create. Ratio npage/npro not an integer.')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        block_count = nmax
        block_length = trp_plan%nlines

        ! Calculate array displacements in Forward Send/Receive
        trp_plan%disp_s(1) = 0
        trp_plan%disp_r(1) = 0
        do i = 2, ims_npro_k
            trp_plan%disp_s(i) = trp_plan%disp_s(i - 1) + block_length
            trp_plan%disp_r(i) = trp_plan%disp_r(i - 1) + block_length*block_count
        end do

        ! #######################################################################
        if (present(locType)) then
            datatype = locType
        else
            datatype = trp_datatype_k                               ! fixed: was trp_datatype_i
        end if

        trp_plan%nmax      = nmax
        trp_plan%base_type = datatype

        stride = npage
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, trp_plan%type_s, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_s, ims_err)

        stride = block_length       ! stride = block_length to put things together
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, trp_plan%type_r, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_r, ims_err)

        ! -----------------------------------------------------------------------
        call MPI_TYPE_SIZE(trp_plan%type_s, ims_ss, ims_err)
        call MPI_TYPE_SIZE(trp_plan%type_r, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            write (str, *) ims_ss; write (line, *) ims_rs
            line = 'Send size '//trim(adjustl(str))//'differs from recv size '//trim(adjustl(line))
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_MPITYPECHECK)
        end if

        return
    end function TLabMPI_Trp_PlanK

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Forward_Real(a, b, trp_plan)
        ! K-Forward: scatter from strided Z-space (a) to flat K-space (b).
        ! Buffer layout: a(m*nlines_p + i*npage + j) for peer m, element i, line j.
        ! After transposition: b(r*chunk + i*nlines_p + j) holds data from rank r.
        real(wp), intent(in) :: a(:)
        real(wp), intent(out) :: b(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr
        integer(wi) :: nmax_p, nlines_p, npage, flat_off, disp_ns, mas
#ifdef USE_APU
        integer(c_int) :: hip_sync_err
        integer(8) :: off                     ! int64 window offset for the fused write+fence real-K push (2026-06-27)
#endif
#ifdef PROFILE_ON
        real(wp) :: time_loc_1, time_loc_2
#endif

#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif
        nmax_p   = trp_plan%nmax
        nlines_p = trp_plan%nlines
        npage    = nlines_p * ims_npro_k   ! total Z-lines across all K ranks
        mas      = nmax_p * nlines_p       ! elements per peer chunk

        ! ==================================================================== !
        ! APU paths — GPU direct writes between peer recv buffers.            !
        ! APU_DIRECT: all ranks on one node share one big window; one fused   !
        !   kernel writes to all peers; MPI_Win_fence is the barrier.         !
        ! FABRIC_DIRECT: multi-node; every peer goes through plain two-sided  !
        !   MPI on fabric_mpi_comm_k (clean MPI_COMM_WORLD split). Strided    !
        !   pack/unpack is CPU-side via wrk_mpi_dp.                           !
        ! ==================================================================== !
#ifdef USE_APU
        if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
            ! -- Push: one fused GPU kernel writes our chunk to ALL peers simultaneously.
            ! Eliminates per-peer kernel-launch overhead (npro separate launches → 1).
            size = trp_plan%size3d
            call MPI_Win_fence(0, apu_win_k, ims_err)
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
            ! Real-K fused-fence push (2026-06-27, the real-K-transpose coherence fix). The source a is STRIDED,
            ! so gather each peer's chunk into contiguous staging (wrk_mpi_dp, free in apudirect K), then
            ! WRITE+per-workgroup __threadfence_system per peer via hip_write_with_fence (write and system fence
            ! in the same wavefront -> reliable cross-XCD visibility; the bare !$omp push + device-scope
            ! hipDeviceSynchronize is NOT a system-scope L2 write-back to MALL).
            !$omp target teams distribute parallel do collapse(3)
            do m = 0, ims_npro_k - 1
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        wrk_mpi_dp(m*mas + i*nlines_p + j + 1) = a(m*nlines_p + i*npage + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            ! No hipDeviceSynchronize: the blocking !$omp end target above already completed the gather, and
            ! hip_write_with_fence/hip_memcpy_push are themselves synchronous, so the push reads committed staging
            ! and the epoch closes after it. (Re-add a sync ONLY if these target regions become nowait.)
            do m = 0, ims_npro_k - 1
                off = int(m, 8)*int(apu_stride_k, 8) + int(ims_pro_k, 8)*int(mas, 8)
#ifdef TRP_I_MEMCPY
                call hip_memcpy_push(apu_all_k(off + 1:off + mas), wrk_mpi_dp(m*mas + 1:m*mas + mas), int(mas, c_int))
#else
                call hip_write_with_fence(wrk_mpi_dp(m*mas + 1:m*mas + mas), apu_all_k(off + 1:off + mas), int(mas, c_int))
#endif
            end do
#else
            !$omp target teams distribute parallel do collapse(3)
            do m = 0, ims_npro_k - 1         ! peer rank (0-based)
                do i = 0, nmax_p - 1          ! element along K axis (kmax total)
                    do j = 0, nlines_p - 1    ! line within peer's chunk
                        ! slot: peer m receives from our rank at offset own_rank*chunk in its buffer
                        apu_all_k(m*apu_stride_k + ims_pro_k*nmax_p*nlines_p + i*nlines_p + j + 1) = &
                            a(m*nlines_p + i*npage + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#endif
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: our recv buffer is now fully populated
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
            call hip_invalidate_recv(apu_recv_fptr_k, int(size, c_int))   ! reader acquire: invalidate L2 -> fresh MALL
#endif
            ! -- Unpack: recv buffer is already in the flat K-space layout; one-to-one copy to b.
            !$omp target teams distribute parallel do
            do i = 1, size
                b(i) = apu_recv_fptr_k(i)
            end do
            !$omp end target teams distribute parallel do

        else if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
            ! K-Forward FABRIC_DIRECT. Option 1 (use_node_win_k): the 4 intra-node K-peers go via
            ! node-local shared-window GPU writes (apudirect mechanism) + MPI_Win_fence; the 4 inter-node
            ! peers via two-sided MPI. Else: all-MPI fallback (the original path). IRECV target is plain
            ! heap (wrk_mpi_dp second half) — shared-window memory is unreliable as an IRECV target.
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
            if (use_node_win_k) then
                ! 1. IRECV inter-node peers into wrk_mpi_dp second half.
                l = 0
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_IRECV(wrk_mpi_dp(size + m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                ! 2. open node-window epoch.
                call MPI_Win_fence(0, node_win_k, ims_err)
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
                ! 3+4. Real-K fused-fence push (2026-06-27, the real-K-transpose coherence fix). ONE kernel
                !      routes each (m,i,j): intra peers gather their STRIDED source chunk into contiguous
                !      staging (wrk_mpi_dp third quarter; offset 2*size clears c_wrk_dp[0:size] inter-pack and
                !      the inter IRECV target [size:2*size]); inter peers pack into c_wrk_dp as before.
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_k - 1
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            if (is_intra_k(m)) then
                                wrk_mpi_dp(2*size + m*mas + i*nlines_p + j + 1) = a(m*nlines_p + i*npage + j + 1)
                            else
                                c_wrk_dp(m*mas + i*nlines_p + j + 1) = a(m*nlines_p + i*npage + j + 1)
                            end if
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                hip_sync_err = hipDeviceSynchronize()   ! gather + inter-pack complete before the HIP push reads staging and before ISEND
                ! inter peers: ISEND the GPU-packed slot (now committed to HBM).
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                ! intra peers: WRITE+per-workgroup __threadfence_system per peer (system-scope visibility in one
                ! wavefront; the bare !$omp push + device-scope hipDeviceSynchronize is NOT a system-scope L2
                ! write-back to MALL -> this is the fix for the cross-XCD real-K push race seen at it=235100).
                do m = 0, ims_npro_k - 1
                    if (.not. is_intra_k(m)) cycle
                    off = int(node_lrank_k(m), 8)*int(apu_size_k, 8) + int(ims_pro_k, 8)*int(mas, 8)
#ifdef TRP_I_MEMCPY
                    call hip_memcpy_push(node_all_k(off + 1:off + mas), wrk_mpi_dp(2*size + m*mas + 1:2*size + m*mas + mas), int(mas, c_int))
#else
                    call hip_write_with_fence(wrk_mpi_dp(2*size + m*mas + 1:2*size + m*mas + mas), &
                                              node_all_k(off + 1:off + mas), int(mas, c_int))
#endif
                end do
                hip_sync_err = hipDeviceSynchronize()   ! wait for the async push+fence kernels before closing the epoch
#else
                ! 3. intra-node peers: ONE fused cross-XCD GPU write over all K-peers (inter masked out) —
                !    collapse(3) over (m,i,j); one kernel launch instead of one per intra peer.
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_k - 1
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            if (is_intra_k(m)) then
                                node_all_k(int(node_lrank_k(m),8)*int(apu_size_k,8) + ims_pro_k*mas + i*nlines_p + j + 1) = &
                                    a(m*nlines_p + i*npage + j + 1)
                            end if
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                ! 4. inter-node peers (Option 2 = GPU-aware MPI): ONE fused GPU pack of all inter peers
                !    a -> c_wrk_dp (intra masked out), then ISEND each GPU-resident slot (no flush; MPICH
                !    orders the stream). m*mas / m*nlines_p inlined — no per-peer scalars in the kernel.
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_k - 1
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            if (.not. is_intra_k(m)) then
                                c_wrk_dp(m*mas + i*nlines_p + j + 1) = a(m*nlines_p + i*npage + j + 1)
                            end if
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#endif
                ! 5. close epoch (intra writes committed) then wait for inter MPI.
                call MPI_Win_fence(0, node_win_k, ims_err)
                if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
                call hip_invalidate_recv(node_recv_fptr_k, int(size, c_int))   ! reader acquire: invalidate L2 -> fresh MALL before assemble
#endif
                ! 6. assemble b per slot m in ONE fused kernel: intra from our node segment, inter from
                !    wrk_mpi_dp — collapse(2) over (m,i), branch inside (was 2 kernels per peer).
                !$omp target teams distribute parallel do collapse(2)
                do m = 0, ims_npro_k - 1
                    do i = 1, mas
                        if (is_intra_k(m)) then
                            b(m*mas + i) = node_recv_fptr_k(m*mas + i)
                        else
                            b(m*mas + i) = wrk_mpi_dp(size + m*mas + i)
                        end if
                    end do
                end do
                !$omp end target teams distribute parallel do
            else
                ! Fallback: all K-peers via two-sided MPI (original fabricdirect path).
                l = 0
                do m = 0, ims_npro_k - 1
                    l = l + 1
                    call MPI_IRECV(wrk_mpi_dp(size + m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                do m = 0, ims_npro_k - 1
                    flat_off = m * mas
                    disp_ns  = m * nlines_p
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            c_wrk_dp(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                        end do
                    end do
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp(flat_off + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
                do i = 1, size
                    b(i) = wrk_mpi_dp(size + i)
                end do
            end if
            nullify (c_wrk_dp)

        else   ! CPU paths: ASYNCHRONOUS, SENDRECV, ALLTOALL
#endif
        ! ==================================================================== !
        ! CPU paths — MPI ISEND/IRECV (ASYNCHRONOUS), SENDRECV, ALLTOALL.     !
        ! When TransposeType=single (MPI_REAL4 with wp=dp) add dp<->sp casts. !
        ! ==================================================================== !
            if (trp_datatype_k == MPI_REAL4 .and. wp == dp) then
                ! Single-precision path: dp→sp before send, sp→dp after recv.
                ! Uses three sp work slots: a_wrk (send copy), b_wrk (recv), c_wrk (flat pack).
                size = trp_plan%size3d
                a_wrk => wrk_mpi_fptr(1:size)
                b_wrk => wrk_mpi_fptr(size + 1:2*size)
                c_wrk => wrk_mpi_fptr(2*size + 1:3*size)
#ifdef USE_APU
                !$omp target teams distribute parallel do
#endif
                do i = 1, size
                    a_wrk(i) = real(a(i), sp)   ! dp→sp
                end do
#ifdef USE_APU
                !$omp end target teams distribute parallel do
#endif
                if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    ! IRECVs posted before GPU pack so network can prepare while GPU runs.
                    ! Step 1: post all IRECVs
                    l = 0
                    do m = 1, ims_npro_k
                        nr = maps_recv_k(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_IRECV(b_wrk(trp_plan%disp_r(nr) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                       ipr, ims_tag, ims_comm_z, request(l), ims_err)
                    end do
                    ! Step 2: GPU gather: a_wrk (strided) → c_wrk (flat)
                    do m = 1, ims_npro_k
                        ns = maps_send_k(m) + 1
                        flat_off = (ns - 1)*nmax_p*nlines_p
                        disp_ns  = trp_plan%disp_s(ns)
#ifdef USE_APU
                        !$omp target teams distribute parallel do collapse(2)
#endif
                        do i = 0, nmax_p - 1
                            do j = 0, nlines_p - 1
                                c_wrk(flat_off + i*nlines_p + j + 1) = a_wrk(disp_ns + i*npage + j + 1)
                            end do
                        end do
#ifdef USE_APU
                        !$omp end target teams distribute parallel do
#endif
                    end do
                    ! Step 3: post all ISENDs from flat c_wrk
                    do m = 1, ims_npro_k
                        ns = maps_send_k(m) + 1; ips = ns - 1
                        l = l + 1
                        call MPI_ISEND(c_wrk((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, MPI_REAL4, &
                                       ips, ims_tag, ims_comm_z, request(l), ims_err)
                    end do
                    ! Step 4: single WAITALL for all sends and receives
                    call MPI_WAITALL(l, request, status, ims_err)
                else
                    call Transpose_Kernel_Single(a_wrk, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                                 b_wrk, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                                 ims_comm_z, trp_sizBlock_k, trp_mode_k)
                end if
#ifdef USE_APU
                !$omp target teams distribute parallel do
#endif
                do i = 1, size
                    b(i) = real(b_wrk(i), dp)   ! sp→dp
                end do
#ifdef USE_APU
                !$omp end target teams distribute parallel do
#endif
                nullify (a_wrk, b_wrk, c_wrk)

            else
                ! Double-precision path: pack directly from dp array a into flat staging buffer.
                if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    ! IRECVs posted before GPU pack so network and GPU overlap.
                    size = trp_plan%size3d
                    call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                    ! Step 1: post all IRECVs
                    l = 0
                    do m = 1, ims_npro_k
                        nr = maps_recv_k(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_IRECV(b(trp_plan%disp_r(nr) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, ims_comm_z, request(l), ims_err)
                    end do
                    ! Step 2: CPU gather: a (strided) → c_wrk_dp (flat)
                    ! NOTE: must be CPU-only. GPU !$omp target write to c_wrk_dp followed by
                    ! CPU MPI_ISEND causes cache-coherency failures on cross-node sends on MI300A
                    ! (GPU L2 not flushed to HBM before MPI reads it over the network fabric).
                    do m = 1, ims_npro_k
                        ns = maps_send_k(m) + 1
                        flat_off = (ns - 1)*nmax_p*nlines_p
                        disp_ns  = trp_plan%disp_s(ns)
                        do i = 0, nmax_p - 1
                            do j = 0, nlines_p - 1
                                c_wrk_dp(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                            end do
                        end do
                    end do
                    ! Step 3: post all ISENDs from flat c_wrk_dp
                    do m = 1, ims_npro_k
                        ns = maps_send_k(m) + 1; ips = ns - 1
                        l = l + 1
                        call MPI_ISEND(c_wrk_dp((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, ims_comm_z, request(l), ims_err)
                    end do
                    ! Step 4: single WAITALL for all sends and receives
                    call MPI_WAITALL(l, request, status, ims_err)
                    nullify (c_wrk_dp)
                else
                    call Transpose_Kernel_Double(a(1:trp_plan%size3d), maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                                 b(1:trp_plan%size3d), maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                                 ims_comm_z, trp_sizBlock_k, trp_mode_k)
                end if
            end if
#ifdef USE_APU
        end if   ! end APU/CPU dispatch
#endif

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine TLabMPI_Trp_ExecK_Forward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Forward_Complex(a, b, trp_plan)
        ! K-Forward (complex): same strided-to-flat scatter as the real version but for complex(dp).
        ! apu_cx_all spans all peers' windows reinterpreted as complex; stride in complex units = apu_stride_k/2.
        ! Complex transposes have no APU_DIRECT path; FABRIC_DIRECT falls back to the MPI path below.
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr
        integer(wi) :: nmax_p, nlines_p, npage, flat_off, disp_ns, mas
        integer :: send_to, recv_from, fbd_tag
#ifdef USE_APU
        complex(dp), pointer :: apu_cx_all(:) => null()   ! complex view of apu_all_k across all peers
        integer(c_int) :: hip_sync_err                    ! flush GPU-assembled b before the CPU FFTW reads it
        integer(8) :: off                                 ! int64 real-view window offset (complex-K fused-fence localization)
        real(dp) :: dbg1(1)                               ! 1-elt CPU scratch for DNS_PROBE_CPU localization markers
#endif
        type(MPI_Comm) :: trp_comm_k   ! fabric_mpi_comm_k (MPI_COMM_WORLD split) or ims_comm_z
#ifdef USE_APU
        if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
            trp_comm_k = fabric_mpi_comm_k   ! clean MPI_COMM_WORLD split; local K-rank = dir rank
            fbd_tag = ims_tag
        else
            trp_comm_k = ims_comm_z
            fbd_tag = ims_tag
        end if
#else
        trp_comm_k = ims_comm_z
        fbd_tag = ims_tag
#endif

        nmax_p   = trp_plan%nmax
        nlines_p = trp_plan%nlines
        npage    = nlines_p * ims_npro_k
        mas      = nmax_p * nlines_p

        ! ==================================================================== !
        ! APU_DIRECT path — fused GPU kernel; same logic as real version but   !
        ! using complex-typed window pointers (apu_stride_k/2 complex elements  !
        ! per segment instead of apu_stride_k real elements).                  !
        ! ==================================================================== !
#ifdef USE_APU
#ifndef TRP_CX_MPI
        if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
            size = trp_plan%size3d
            ! Reinterpret the contiguous real window as complex; size in complex units = apu_stride_k/2
            call c_f_pointer(apu_peer_cptr_k(0), apu_cx_all, [apu_stride_k*ims_npro_k/2])
            call MPI_Win_fence(0, apu_win_k, ims_err)
            ! Push: strided gather from a into all peers' recv buffers in one fused kernel
            !$omp target teams distribute parallel do collapse(3)
            do m = 0, ims_npro_k - 1
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        ! apu_stride_k/2 = segment size in complex units; own rank writes at slot own_rank*chunk
                        apu_cx_all(m*(apu_stride_k/2) + ims_pro_k*nmax_p*nlines_p + i*nlines_p + j + 1) = &
                            a(m*nlines_p + i*npage + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            hip_sync_err = hipDeviceSynchronize()   ! ROOT FIX: flush GPU push to HBM before the fence (cross-rank visibility)
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: recv buffer now fully populated
            ! Unpack: flat copy from complex-typed recv alias to b
            !$omp target teams distribute parallel do
            do i = 1, size
                b(i) = apu_cx_recv_fptr_k(i)
            end do
            !$omp end target teams distribute parallel do
            ! b is GPU-written here but its consumer is the CPU FFTW in OPR_Fourier_Z_Forward/Backward.
            ! MPI_Win_fence orders the RMA epoch but does NOT flush the GPU write to be CPU-coherent on
            ! MI300A (unlike the real K-transpose, whose consumer is the GPU FDM solve and is stream-ordered).
            ! Without this flush the CPU FFTW intermittently reads stale b -> the random one-step blow-up.
            hip_sync_err = hipDeviceSynchronize()
            nullify (apu_cx_all)

        else if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
#else
        if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT .or. trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
#endif
            ! K-Forward complex FABRIC_DIRECT. Option 1 (use_node_win_k): the intra-node K-peers go via
            ! node-local complex shared-window GPU writes + MPI_Win_fence (the apudirect-complex mechanism on
            ! the node comm); the inter-node K-peers via two-sided MPI, IRECV straight into b (b is consumed
            ! by the LOCAL CPU FFTW, so its inter slots are NIC-written = CPU-coherent and its intra slots are
            ! GPU-assembled = CPU-coherent after the synchronous target region, per the apudirect-complex path
            ! above). Only the inter NIC send of GPU-written a needs the flush. Complex view of the window:
            ! apu_size_k/2 complex units per segment. Else: all-MPI complex fallback (original path).
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
            ! Option 2: the inter-node leg is GPU-aware (GPU pack + ISEND of the GPU buffer), so no
            ! hipDeviceSynchronize is needed (MPICH_GPU_SUPPORT orders the GPU stream — vmpi_gpuaware M2,
            ! vmpi_nodewin window+GPU-aware PASS). Frees the CPU pack the old flush+CPU-pack path required.
            if (use_node_win_k) then
                ! 1. IRECV inter-node peers directly into b's slots (disp_r(m+1) = m*mas on fabric comm).
                l = 0
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_IRECV(b(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                ! 2. open node-window epoch.
                call MPI_Win_fence(0, node_win_k, ims_err)
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
                ! ===== INSTRUMENTED complex-K forward fused-fence (LOCALIZATION; last CXKF:* in the dead rank's
                !       fort.5xx = the faulting op; 'max=' = index/bound, 'sub=' = peer m). =====
                dbg1(1) = real(apu_size_k,dp)*real(node_size_k,dp); DNS_PROBE_CPU('CXKF:bnd-nodeall', dbg1, 1, -1)
                dbg1(1) = real(4,dp)*real(imax,dp)*real(jmax,dp)*real(kmax,dp); DNS_PROBE_CPU('CXKF:bnd-wrk', dbg1, 1, -1)
                dbg1(1) = real(2*size + 2*ims_npro_k*mas,dp);       DNS_PROBE_CPU('CXKF:stage-end', dbg1, 1, -1)
                dbg1(1) = real(mas,dp);                             DNS_PROBE_CPU('CXKF:0-enter', dbg1, 1, ims_pro_k)
                ! gather strided a: intra -> wrk_mpi_dp BY NAME (real/imag, one device handle); inter -> c_wrk_cx
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_k - 1
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            if (is_intra_k(m)) then
                                wrk_mpi_dp(2*size + 2*(m*mas + i*nlines_p + j + 1) - 1) = real(a(m*nlines_p + i*npage + j + 1), dp)
                                wrk_mpi_dp(2*size + 2*(m*mas + i*nlines_p + j + 1))     = aimag(a(m*nlines_p + i*npage + j + 1))
                            else
                                c_wrk_cx(m*mas + i*nlines_p + j + 1) = a(m*nlines_p + i*npage + j + 1)
                            end if
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                hip_sync_err = hipDeviceSynchronize()
                DNS_PROBE_CPU('CXKF:1-post-gather', dbg1, 1, ims_pro_k)
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_ISEND(c_wrk_cx(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                DNS_PROBE_CPU('CXKF:2-post-isend', dbg1, 1, l)
                do m = 0, ims_npro_k - 1
                    if (.not. is_intra_k(m)) cycle
                    off = int(node_lrank_k(m),8)*int(apu_size_k,8) + int(ims_pro_k,8)*int(2*mas,8)
                    dbg1(1) = real(off + 2*mas, dp); DNS_PROBE_CPU('CXKF:3-pre-fence-end', dbg1, 1, m)
#ifdef TRP_I_MEMCPY
                    call hip_memcpy_push(node_all_k(off + 1:off + 2*mas), wrk_mpi_dp(2*size + m*2*mas + 1:2*size + m*2*mas + 2*mas), int(2*mas, c_int))
#else
                    call hip_write_with_fence(wrk_mpi_dp(2*size + m*2*mas + 1:2*size + m*2*mas + 2*mas), &
                                              node_all_k(off + 1:off + 2*mas), int(2*mas, c_int))
#endif
                    DNS_PROBE_CPU('CXKF:3-post-fence-peer', dbg1, 1, m)
                end do
                hip_sync_err = hipDeviceSynchronize()
                DNS_PROBE_CPU('CXKF:4-post-fence-all', dbg1, 1, ims_pro_k)
#else
                ! 3. intra-node peers: ONE fused cross-XCD GPU write over all K-peers (inter masked out) —
                !    collapse(3) over (m,i,j); one kernel instead of one per intra peer.
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_k - 1
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            if (is_intra_k(m)) then
                                node_cx_all_k(int(node_lrank_k(m),8)*int(apu_size_k/2,8) + ims_pro_k*mas + i*nlines_p + j + 1) = &
                                    a(m*nlines_p + i*npage + j + 1)
                            end if
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                ! 4. inter-node peers (Option 2 = GPU-aware MPI): ONE fused GPU pack of all inter peers
                !    a -> c_wrk_cx (intra masked out), then ISEND each GPU-resident slot (no flush; MPICH
                !    orders the stream). m*mas / m*nlines_p inlined — no per-peer scalars in the kernel.
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_k - 1
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            if (.not. is_intra_k(m)) then
                                c_wrk_cx(m*mas + i*nlines_p + j + 1) = a(m*nlines_p + i*npage + j + 1)
                            end if
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_ISEND(c_wrk_cx(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#endif
                ! 5. close epoch (intra writes committed) then wait for inter MPI.
                call MPI_Win_fence(0, node_win_k, ims_err)
                if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
                DNS_PROBE_CPU('CXKF:5-pre-inval', dbg1, 1, ims_pro_k)
                call hip_invalidate_recv(node_recv_fptr_k, int(apu_size_k, c_int))   ! reader acquire (real view)
                DNS_PROBE_CPU('CXKF:6-post-inval', dbg1, 1, ims_pro_k)
#endif
                ! 6. assemble intra slots of b in ONE fused kernel (inter slots already in b via IRECV) —
                !    collapse(2) over (m,i), masked is_intra; was a kernel per intra peer.
                !$omp target teams distribute parallel do collapse(2)
                do m = 0, ims_npro_k - 1
                    do i = 1, mas
                        if (is_intra_k(m)) then
                            b(m*mas + i) = node_cx_recv_fptr_k(m*mas + i)
                        end if
                    end do
                end do
                !$omp end target teams distribute parallel do
                ! b's intra slots are GPU-assembled; its consumer is the CPU FFTW. Flush so the GPU writes
                ! are CPU-coherent before the FFTW reads them (the inter slots came via IRECV, already
                ! CPU-coherent). Same hazard the APU_DIRECT path above fixes.
                hip_sync_err = hipDeviceSynchronize()
            else
                ! Fallback: all K-peers via two-sided MPI (original complex fabricdirect path).
                do m = 1, ims_npro_k
                    ns = maps_send_k(m) + 1
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_s(ns)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            c_wrk_cx(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                        end do
                    end do
                end do
                do j = 1, ims_npro_k, trp_sizBlock_k
                    l = 0
                    do m = j, min(j + trp_sizBlock_k - 1, ims_npro_k)
                        ns = maps_send_k(m) + 1; ips = ns - 1
                        nr = maps_recv_k(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_ISEND(c_wrk_cx((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(b(trp_plan%disp_r(nr) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                    end do
                    call MPI_WAITALL(l, request, status, ims_err)
                end do
            end if
            nullify (c_wrk_cx)

        else   ! ASYNCHRONOUS, SENDRECV, ALLTOALL
#endif
        ! ==================================================================== !
        ! MPI path — pack into flat staging buffer c_wrk_cx, then ISEND/IRECV. !
        ! ==================================================================== !
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
                ! Pack: strided a → flat c_wrk_cx
                do m = 1, ims_npro_k
                    ns = maps_send_k(m) + 1
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_s(ns)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            c_wrk_cx(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                        end do
                    end do
                end do
                ! ISEND/IRECV in batches of trp_sizBlock_k peers per WAITALL.
                ! trp_comm_k is ims_comm_z (ASYNCHRONOUS); ips/ipr are local K-comm ranks (0..npro_k-1).
                do j = 1, ims_npro_k, trp_sizBlock_k
                    l = 0
                    do m = j, min(j + trp_sizBlock_k - 1, ims_npro_k)
                        ns = maps_send_k(m) + 1; ips = ns - 1
                        nr = maps_recv_k(m) + 1; ipr = nr - 1
                        send_to   = ips
                        recv_from = ipr
                        l = l + 1
                        call MPI_ISEND(c_wrk_cx((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, send_to, fbd_tag, trp_comm_k, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(b(trp_plan%disp_r(nr) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, recv_from, fbd_tag, trp_comm_k, request(l), ims_err)
                    end do
                    call MPI_WAITALL(l, request, status, ims_err)
                end do
                nullify (c_wrk_cx)
            else
                call Transpose_Kernel_Complex(a, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                              b, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                              trp_comm_k, trp_sizBlock_k, trp_mode_k)
            end if
#ifdef USE_APU
        end if   ! end APU/CPU dispatch
#endif
        return
    end subroutine TLabMPI_Trp_ExecK_Forward_Complex

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Backward_Real(b, a, trp_plan)
        ! K-Backward: reassemble from flat K-space (b) to strided Z-space (a). Reverse of K-Forward.
        ! b layout: b(m*chunk + i) for peer m, element i (flat per peer).
        ! a layout: a(m*nlines_p + i*npage + j) for peer m, element i, line j (strided).
        real(wp), intent(in) :: b(:)
        real(wp), intent(out) :: a(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr
        integer(wi) :: nmax_p, nlines_p, npage, flat_off, disp_nr, mas
#ifdef USE_APU
        integer(c_int) :: hip_sync_err
        integer(8) :: off                     ! int64 window offset for the fused write+fence real-K push (2026-06-27)
#endif
#ifdef PROFILE_ON
        real(wp) :: time_loc_1, time_loc_2
#endif

#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif
        nmax_p   = trp_plan%nmax
        nlines_p = trp_plan%nlines
        npage    = nlines_p * ims_npro_k
        mas      = nmax_p * nlines_p
        ! ==================================================================== !
        ! APU paths — GPU direct writes; inverse of K-Forward.                  !
        ! b is flat (r*chunk layout); each rank pushes its chunk to all peers.  !
        ! After fence, apu_recv_fptr_k is in strided Z-space layout (= a).     !
        ! ==================================================================== !
#ifdef USE_APU
        if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
            ! -- Push: b is flat K-space; push our chunk to all peers' recv buffers in one fused kernel.
            size = trp_plan%size3d
            call MPI_Win_fence(0, apu_win_k, ims_err)
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
            ! Real-K fused-fence push (2026-06-27): b is flat (contiguous mas per peer), so WRITE+per-workgroup
            ! __threadfence_system per peer directly via hip_write_with_fence (no gather needed). System-scope
            ! L2 write-back in one wavefront -> reliable cross-XCD visibility (vs the device-scope bare push).
            ! No hipDeviceSynchronize: b is produced by the caller's blocking !$omp target (the codebase uses no
            ! nowait), and hip_write_with_fence/hip_memcpy_push are themselves synchronous, so the push reads
            ! committed b and the epoch closes after it. (Re-add a sync ONLY if the producer/these targets go nowait.)
            do m = 0, ims_npro_k - 1
                off = int(m, 8)*int(apu_stride_k, 8) + int(ims_pro_k, 8)*int(mas, 8)
#ifdef TRP_I_MEMCPY
                call hip_memcpy_push(apu_all_k(off + 1:off + mas), b(m*mas + 1:m*mas + mas), int(mas, c_int))
#else
                call hip_write_with_fence(b(m*mas + 1:m*mas + mas), apu_all_k(off + 1:off + mas), int(mas, c_int))
#endif
            end do
#else
            !$omp target teams distribute parallel do collapse(2)
            do m = 0, ims_npro_k - 1
                do i = 1, nmax_p * nlines_p
                    ! Write our chunk (b[m*chunk]) into peer m's slot (own_rank*chunk) in their buffer
                    apu_all_k(m*apu_stride_k + ims_pro_k*nmax_p*nlines_p + i) = b(m * nmax_p * nlines_p + i)
                end do
            end do
            !$omp end target teams distribute parallel do
            hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#endif
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: all peers have written to our buffer
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
            call hip_invalidate_recv(apu_recv_fptr_k, int(size, c_int))   ! reader acquire: invalidate L2 -> fresh MALL
#endif
            ! -- Unpack: recv buffer holds sorted chunks; scatter to strided a in one fused kernel.
            !$omp target teams distribute parallel do collapse(3)
            do m = 0, ims_npro_k - 1
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        a(m*nlines_p + i*npage + j + 1) = &
                            apu_recv_fptr_k(m*nmax_p*nlines_p + i*nlines_p + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do

        else if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
            ! K-Backward FABRIC_DIRECT. Option 1: intra-node peers via node-window GPU push + fence;
            ! inter-node peers via two-sided MPI. b is flat per-peer; a is strided. IRECV target = heap.
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
            if (use_node_win_k) then
                ! 1. IRECV inter peers into c_wrk_dp(m*mas).
                l = 0
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_IRECV(c_wrk_dp(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                ! 2. open node-window epoch.
                call MPI_Win_fence(0, node_win_k, ims_err)
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
                ! 3+4. Real-K fused-fence push (2026-06-27, the real-K-transpose coherence fix). b is flat
                !      (contiguous mas per peer). Sync first (b is produced on the OMP stream and is also the
                !      inter ISEND source), then WRITE+per-workgroup __threadfence_system per intra peer via
                !      hip_write_with_fence (system-scope L2 write-back in one wavefront -> fixes the cross-XCD
                !      real-K push race). Inter peers ISEND b directly as before.
                hip_sync_err = hipDeviceSynchronize()
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_ISEND(b(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                do m = 0, ims_npro_k - 1
                    if (.not. is_intra_k(m)) cycle
                    off = int(node_lrank_k(m), 8)*int(apu_size_k, 8) + int(ims_pro_k, 8)*int(mas, 8)
#ifdef TRP_I_MEMCPY
                    call hip_memcpy_push(node_all_k(off + 1:off + mas), b(m*mas + 1:m*mas + mas), int(mas, c_int))
#else
                    call hip_write_with_fence(b(m*mas + 1:m*mas + mas), node_all_k(off + 1:off + mas), int(mas, c_int))
#endif
                end do
                hip_sync_err = hipDeviceSynchronize()   ! wait for the async push+fence kernels before closing the epoch
#else
                ! 3. intra peers: ONE fused GPU push of b into each peer's segment (inter masked out) —
                !    collapse(2) over (m,i); one kernel instead of one per intra peer.
                !$omp target teams distribute parallel do collapse(2)
                do m = 0, ims_npro_k - 1
                    do i = 1, mas
                        if (is_intra_k(m)) then
                            node_all_k(int(node_lrank_k(m),8)*int(apu_size_k,8) + ims_pro_k*mas + i) = b(m*mas + i)
                        end if
                    end do
                end do
                !$omp end target teams distribute parallel do
                ! 4. inter peers: ISEND b[m*chunk] directly (b is flat).
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_ISEND(b(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#endif
                ! 5. close epoch then wait.
                call MPI_Win_fence(0, node_win_k, ims_err)
                if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
                call hip_invalidate_recv(node_recv_fptr_k, int(size, c_int))   ! reader acquire: invalidate L2 -> fresh MALL before scatter
#endif
                ! 6. scatter to strided a in ONE fused kernel: intra from our node segment, inter from
                !    c_wrk_dp — collapse(3) over (m,i,j), branch inside (was a kernel per peer).
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_k - 1
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            if (is_intra_k(m)) then
                                a(m*nlines_p + i*npage + j + 1) = node_recv_fptr_k(m*mas + i*nlines_p + j + 1)
                            else
                                a(m*nlines_p + i*npage + j + 1) = c_wrk_dp(m*mas + i*nlines_p + j + 1)
                            end if
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
            else
                ! Fallback: all K-peers via two-sided MPI (original fabricdirect path).
                l = 0
                do m = 0, ims_npro_k - 1
                    l = l + 1
                    call MPI_IRECV(c_wrk_dp(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                do m = 0, ims_npro_k - 1
                    l = l + 1
                    call MPI_ISEND(b(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
                do m = 0, ims_npro_k - 1
                    flat_off = m * mas
                    disp_nr  = m * nlines_p
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            a(disp_nr + i*npage + j + 1) = c_wrk_dp(flat_off + i*nlines_p + j + 1)
                        end do
                    end do
                end do
            end if
            nullify (c_wrk_dp)

        else   ! CPU paths: ASYNCHRONOUS, SENDRECV, ALLTOALL
#endif
        ! ==================================================================== !
        ! CPU paths — MPI ISEND/IRECV (ASYNCHRONOUS), SENDRECV, ALLTOALL.     !
        ! ==================================================================== !
            if (trp_datatype_k == MPI_REAL4 .and. wp == dp) then
                ! Single-precision path: dp→sp before send, sp→dp merged into unpack.
                ! b_wrk = dp→sp copy of b (send buffer); c_wrk = flat recv staging.
                size = trp_plan%size3d
                if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    b_wrk => wrk_mpi_fptr(1:size)           ! send: dp→sp copy of b
                    c_wrk => wrk_mpi_fptr(size + 1:2*size)  ! recv: flat staging
#ifdef USE_APU
                    !$omp target teams distribute parallel do
#endif
                    do i = 1, size
                        b_wrk(i) = real(b(i), sp)   ! dp→sp
                    end do
#ifdef USE_APU
                    !$omp end target teams distribute parallel do
#endif
                    ! Step 1: IRECVs (note backward map swap: mrecv=maps_send_k, msend=maps_recv_k)
                    l = 0
                    do m = 1, ims_npro_k
                        nr = maps_send_k(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_IRECV(c_wrk((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, MPI_REAL4, &
                                       ipr, ims_tag, ims_comm_z, request(l), ims_err)
                    end do
                    ! Step 2: ISENDs
                    do m = 1, ims_npro_k
                        ns = maps_recv_k(m) + 1; ips = ns - 1
                        l = l + 1
                        call MPI_ISEND(b_wrk(trp_plan%disp_r(ns) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                       ips, ims_tag, ims_comm_z, request(l), ims_err)
                    end do
                    ! Step 3: WAITALL
                    call MPI_WAITALL(l, request, status, ims_err)
                    ! Step 4: GPU unpack: flat c_wrk → strided a (sp→dp merged)
                    do m = 1, ims_npro_k
                        nr = maps_send_k(m) + 1
                        flat_off = (nr - 1)*nmax_p*nlines_p
                        disp_nr  = trp_plan%disp_s(nr)
#ifdef USE_APU
                        !$omp target teams distribute parallel do collapse(2)
#endif
                        do i = 0, nmax_p - 1
                            do j = 0, nlines_p - 1
                                a(disp_nr + i*npage + j + 1) = real(c_wrk(flat_off + i*nlines_p + j + 1), dp)
                            end do
                        end do
#ifdef USE_APU
                        !$omp end target teams distribute parallel do
#endif
                    end do
                    nullify (b_wrk, c_wrk)
                else
                    b_wrk => wrk_mpi_fptr(1:size)
                    a_wrk => wrk_mpi_fptr(size + 1:2*size)
#ifdef USE_APU
                    !$omp target teams distribute parallel do
#endif
                    do i = 1, size
                        b_wrk(i) = real(b(i), sp)   ! dp→sp
                    end do
#ifdef USE_APU
                    !$omp end target teams distribute parallel do
#endif
                    call Transpose_Kernel_Single(b_wrk, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                                 a_wrk, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                                 ims_comm_z, trp_sizBlock_k, trp_mode_k)
#ifdef USE_APU
                    !$omp target teams distribute parallel do
#endif
                    do i = 1, size
                        a(i) = real(a_wrk(i), dp)   ! sp→dp
                    end do
#ifdef USE_APU
                    !$omp end target teams distribute parallel do
#endif
                    nullify (a_wrk, b_wrk)
                end if

            else
                ! Double-precision path
                if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    ! b is already flat K-space; ISEND directly, recv into c_wrk_dp, scatter to a.
                    size = trp_plan%size3d
                    call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                    ! Step 1: IRECVs
                    l = 0
                    do m = 1, ims_npro_k
                        nr = maps_send_k(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_IRECV(c_wrk_dp((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, ims_comm_z, request(l), ims_err)
                    end do
                    ! Step 2: ISENDs from flat b
                    do m = 1, ims_npro_k
                        ns = maps_recv_k(m) + 1; ips = ns - 1
                        l = l + 1
                        call MPI_ISEND(b(trp_plan%disp_r(ns) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, ims_comm_z, request(l), ims_err)
                    end do
                    ! Step 3: WAITALL
                    call MPI_WAITALL(l, request, status, ims_err)
                    ! Step 4: GPU scatter: flat c_wrk_dp → strided a
                    do m = 1, ims_npro_k
                        nr = maps_send_k(m) + 1
                        flat_off = (nr - 1)*nmax_p*nlines_p
                        disp_nr  = trp_plan%disp_s(nr)
#ifdef USE_APU
                        !$omp target teams distribute parallel do collapse(2)
#endif
                        do i = 0, nmax_p - 1
                            do j = 0, nlines_p - 1
                                a(disp_nr + i*npage + j + 1) = c_wrk_dp(flat_off + i*nlines_p + j + 1)
                            end do
                        end do
#ifdef USE_APU
                        !$omp end target teams distribute parallel do
#endif
                    end do
                    nullify (c_wrk_dp)
                else
                    call Transpose_Kernel_Double(b(1:trp_plan%size3d), maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                                 a(1:trp_plan%size3d), maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                                 ims_comm_z, trp_sizBlock_k, trp_mode_k)
                end if
            end if
#ifdef USE_APU
        end if   ! end APU/CPU dispatch
#endif

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine TLabMPI_Trp_ExecK_Backward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Backward_Complex(b, a, trp_plan)
        ! K-Backward (complex): inverse of K-Forward_Complex. b is flat; a is strided.
        complex(wp), intent(in) :: b(:)
        complex(wp), intent(out) :: a(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr
        integer(wi) :: nmax_p, nlines_p, npage, flat_off, disp_nr, mas
        integer :: send_to, recv_from, fbd_tag
#ifdef USE_APU
        complex(dp), pointer :: apu_cx_all(:) => null()
        integer(c_int) :: hip_sync_err                    ! flush GPU-assembled a before the CPU FFTW reads it
        integer(8) :: off                                 ! int64 real-view window offset (complex-K fused-fence localization)
        real(dp) :: dbg1(1)                               ! 1-elt CPU scratch for DNS_PROBE_CPU localization markers
#endif
        type(MPI_Comm) :: trp_comm_k   ! fabric_mpi_comm_k (MPI_COMM_WORLD split) or ims_comm_z
#ifdef USE_APU
        if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
            trp_comm_k = fabric_mpi_comm_k   ! clean MPI_COMM_WORLD split; local K-rank = dir rank m
            fbd_tag = ims_tag
        else
            trp_comm_k = ims_comm_z
            fbd_tag = ims_tag
        end if
#else
        trp_comm_k = ims_comm_z
        fbd_tag = ims_tag
#endif

        nmax_p   = trp_plan%nmax
        nlines_p = trp_plan%nlines
        npage    = nlines_p * ims_npro_k
        mas      = nmax_p * nlines_p

        ! ==================================================================== !
        ! APU_DIRECT path — fused GPU kernels; inverse of K-Forward_Complex.   !
        ! ==================================================================== !
#ifdef USE_APU
#ifndef TRP_CX_MPI
        if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
            size = trp_plan%size3d
            call c_f_pointer(apu_peer_cptr_k(0), apu_cx_all, [apu_stride_k*ims_npro_k/2])
            call MPI_Win_fence(0, apu_win_k, ims_err)
            ! Push: b is flat K-space; write our chunk (b[m*chunk]) to every peer
            !$omp target teams distribute parallel do collapse(2)
            do m = 0, ims_npro_k - 1
                do i = 1, nmax_p * nlines_p
                    apu_cx_all(m*(apu_stride_k/2) + ims_pro_k*nmax_p*nlines_p + i) = &
                        b(m * nmax_p * nlines_p + i)
                end do
            end do
            !$omp end target teams distribute parallel do
            hip_sync_err = hipDeviceSynchronize()   ! ROOT FIX: flush GPU push to HBM before the fence (cross-rank visibility)
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: recv buffer fully populated
            ! CRASH-LOC (apudirect K-bwd-cplx): our recv window AFTER cross-rank push+fence, BEFORE the unpack.
            ! SYMMETRIC to the fabricdirect ZKBC probes — covers the OTHER communication branch (single-node
            ! apudirect window) so a crash on either branch is pinned internally. Input b (=ZFWD:post-fft)
            ! healthy; if this window is blown the apudirect cross-rank push/fence is the seed.
            DNS_PROBE('ZKBC:apuwin', apu_recv_fptr_k(1), 2*size, -1)
            ! Unpack: complex recv buffer → strided Z-space a
            !$omp target teams distribute parallel do collapse(3)
            do m = 0, ims_npro_k - 1
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        a(m*nlines_p + i*npage + j + 1) = &
                            apu_cx_recv_fptr_k(m*nmax_p*nlines_p + i*nlines_p + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            ! a is GPU-written but its consumer is the CPU FFTW (OPR_Fourier_X_Backward). Flush so the GPU
            ! write is CPU-coherent before the FFTW reads it (the missing flush = the random one-step blow-up).
            hip_sync_err = hipDeviceSynchronize()
            nullify (apu_cx_all)

        else if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
#else
        if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT .or. trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
#endif
            ! K-Backward complex FABRIC_DIRECT. Option 1 (use_node_win_k): intra-node K-peers via node-local
            ! complex shared-window GPU push + MPI_Win_fence, then GPU-unpack our recv segment → strided a
            ! (apudirect-complex pattern); inter-node K-peers via two-sided MPI into flat c_wrk_cx, then a
            ! CPU scatter → strided a (the existing fabricdirect-complex pattern — a is strided so we cannot
            ! IRECV into it directly). a is consumed by the GPU elliptic Y-solve: intra slots GPU-written and
            ! inter slots CPU-written are both coherent for it on unified memory (apudirect / fabricdirect
            ! precedents). Else: all-MPI complex fallback (original path).
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi_dp(2*size + 1)), c_wrk_send_cx, shape=[size])   ! GPU-packed SEND staging (2nd half)
            ! BUGFIX [2026-06-27]: the inter-node leg must ISEND a GPU-WRITTEN buffer for GPU-aware MPI to be
            ! coherent. The FORWARD complex K GPU-packs its input into c_wrk_cx then ISENDs that (correct); this
            ! BACKWARD used to ISEND b DIRECTLY, but b is the CPU-FFTW z2z output, NOT GPU-written -> GPU-aware
            ! MPI shipped uninitialized garbage (the 235056 NaN seed: ZKBC:cwrk-inter 6.46e208 from a healthy b).
            ! Fix: GPU-pack b's inter chunks into c_wrk_send_cx (GPU-written) in the SAME kernel as the intra
            ! push, then ISEND c_wrk_send_cx (stream-ordered, exactly like the forward).
            if (use_node_win_k) then
                ! 1. IRECV inter-node peers into flat c_wrk_cx slots (m*mas).
                l = 0
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_IRECV(c_wrk_cx(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                ! 2. open node-window epoch.
                call MPI_Win_fence(0, node_win_k, ims_err)
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
                ! ===== INSTRUMENTED complex-K backward fused-fence (LOCALIZATION: build -DTRP_I_FUSEDFENCE
                !       -DDNS_DEBUG_PROBES; the LAST CXKB:* line in the dead rank's fort.5xx names the faulting op,
                !       and 'max=' carries the index/bound, 'sub=' the peer m). =====
                dbg1(1) = real(apu_size_k,dp)*real(node_size_k,dp); DNS_PROBE_CPU('CXKB:bnd-nodeall', dbg1, 1, -1)
                dbg1(1) = real(4,dp)*real(imax,dp)*real(jmax,dp)*real(kmax,dp); DNS_PROBE_CPU('CXKB:bnd-wrk', dbg1, 1, -1)
                dbg1(1) = real(2*size + 2*ims_npro_k*mas,dp);       DNS_PROBE_CPU('CXKB:stage-end', dbg1, 1, -1)
                dbg1(1) = real(mas,dp);                             DNS_PROBE_CPU('CXKB:0-enter', dbg1, 1, ims_pro_k)
                ! gather: intra -> wrk_mpi_dp BY NAME (real/imag, one device handle); inter -> c_wrk_send_cx
                !$omp target teams distribute parallel do collapse(2)
                do m = 0, ims_npro_k - 1
                    do i = 1, mas
                        if (is_intra_k(m)) then
                            wrk_mpi_dp(2*size + 2*(m*mas + i) - 1) = real(b(m*mas + i), dp)
                            wrk_mpi_dp(2*size + 2*(m*mas + i))     = aimag(b(m*mas + i))
                        else
                            c_wrk_send_cx(m*mas + i) = b(m*mas + i)
                        end if
                    end do
                end do
                !$omp end target teams distribute parallel do
                hip_sync_err = hipDeviceSynchronize()
                DNS_PROBE_CPU('CXKB:1-post-gather', dbg1, 1, ims_pro_k)
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_ISEND(c_wrk_send_cx(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                DNS_PROBE_CPU('CXKB:2-post-isend', dbg1, 1, l)
                do m = 0, ims_npro_k - 1
                    if (.not. is_intra_k(m)) cycle
                    off = int(node_lrank_k(m),8)*int(apu_size_k,8) + int(ims_pro_k,8)*int(2*mas,8)
                    dbg1(1) = real(off + 2*mas, dp); DNS_PROBE_CPU('CXKB:3-pre-fence-end', dbg1, 1, m)
#ifdef TRP_I_MEMCPY
                    call hip_memcpy_push(node_all_k(off + 1:off + 2*mas), wrk_mpi_dp(2*size + m*2*mas + 1:2*size + m*2*mas + 2*mas), int(2*mas, c_int))
#else
                    call hip_write_with_fence(wrk_mpi_dp(2*size + m*2*mas + 1:2*size + m*2*mas + 2*mas), &
                                              node_all_k(off + 1:off + 2*mas), int(2*mas, c_int))
#endif
                    DNS_PROBE_CPU('CXKB:3-post-fence-peer', dbg1, 1, m)
                end do
                hip_sync_err = hipDeviceSynchronize()
                DNS_PROBE_CPU('CXKB:4-post-fence-all', dbg1, 1, ims_pro_k)
#else
                ! 3. ONE fused GPU kernel: intra peers -> node-window segment; inter peers -> c_wrk_send_cx
                !    (GPU-written send staging). Both read b on the GPU (coherent), so the ISEND'd inter buffer
                !    is GPU-written and GPU-aware MPI is valid (the whole point of the bugfix above).
                !$omp target teams distribute parallel do collapse(2)
                do m = 0, ims_npro_k - 1
                    do i = 1, mas
                        if (is_intra_k(m)) then
                            node_cx_all_k(int(node_lrank_k(m),8)*int(apu_size_k/2,8) + ims_pro_k*mas + i) = &
                                b(m*mas + i)
                        else
                            c_wrk_send_cx(m*mas + i) = b(m*mas + i)
                        end if
                    end do
                end do
                !$omp end target teams distribute parallel do
                ! 4. inter-node peers: ISEND our GPU-packed c_wrk_send_cx[m*mas] chunk (overlaps the window writes).
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    l = l + 1
                    call MPI_ISEND(c_wrk_send_cx(m*mas + 1), mas, &
                                   trp_plan%base_type, m, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                end do
                hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#endif
                ! 5. close epoch (intra writes committed) then wait for inter MPI.
                call MPI_Win_fence(0, node_win_k, ims_err)
                if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
#if defined(TRP_I_FUSEDFENCE) || defined(TRP_I_MEMCPY)
                DNS_PROBE_CPU('CXKB:5-pre-inval', dbg1, 1, ims_pro_k)
                call hip_invalidate_recv(node_recv_fptr_k, int(apu_size_k, c_int))   ! reader acquire (real view)
                DNS_PROBE_CPU('CXKB:6-post-inval', dbg1, 1, ims_pro_k)
#endif
                ! CRASH-LOC (fabricdirect K-bwd-cplx, 2-node): split the two recv legs to pin the faulting one.
                ! cwrk-inter = inter-node MPI recv (4 inter peers, c_wrk_cx); nodewin = our node-window segment
                ! (4 intra peers pushed in). Input b (=ZFWD:post-fft) is healthy; whichever recv is blown =
                ! the faulting leg (inter-node GPU-aware MPI vs intra node-window). Apudirect lacks the inter leg.
                DNS_PROBE('ZKBC:cwrk-inter', wrk_mpi_dp(1), 2*size, -1)
                DNS_PROBE('ZKBC:nodewin', node_recv_fptr_k(1), apu_size_k, -1)
                ! 6a. intra peers: ONE fused GPU unpack of our recv segment → strided a (inter masked out) —
                !     collapse(3) over (m,i,j); one kernel instead of one per intra peer. (6b inter stays CPU.)
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_k - 1
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            if (is_intra_k(m)) then
                                a(m*nlines_p + i*npage + j + 1) = node_cx_recv_fptr_k(m*mas + i*nlines_p + j + 1)
                            end if
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                ! a's intra slots are GPU-assembled; its consumer is the CPU FFTW. Flush so the GPU writes
                ! are CPU-coherent before the FFTW reads them. Same hazard the APU_DIRECT path above fixes.
                hip_sync_err = hipDeviceSynchronize()
                ! 6b. inter peers: CPU scatter flat c_wrk_cx → strided a.
                do m = 0, ims_npro_k - 1
                    if (is_intra_k(m)) cycle
                    flat_off = m * mas
                    disp_nr  = m * nlines_p
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            a(disp_nr + i*npage + j + 1) = c_wrk_cx(flat_off + i*nlines_p + j + 1)
                        end do
                    end do
                end do
            else
                ! Fallback: all K-peers via two-sided MPI (original complex fabricdirect path).
                do j = 1, ims_npro_k, trp_sizBlock_k
                    l = 0
                    do m = j, min(j + trp_sizBlock_k - 1, ims_npro_k)
                        ns = maps_recv_k(m) + 1; ips = ns - 1
                        nr = maps_send_k(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_ISEND(b(trp_plan%disp_r(ns) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(c_wrk_cx((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, fabric_mpi_comm_k, request(l), ims_err)
                    end do
                    call MPI_WAITALL(l, request, status, ims_err)
                end do
                do m = 1, ims_npro_k
                    nr = maps_send_k(m) + 1
                    flat_off = (nr - 1)*nmax_p*nlines_p
                    disp_nr  = trp_plan%disp_s(nr)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            a(disp_nr + i*npage + j + 1) = c_wrk_cx(flat_off + i*nlines_p + j + 1)
                        end do
                    end do
                end do
            end if
            nullify (c_wrk_cx)

        else   ! ASYNCHRONOUS, SENDRECV, ALLTOALL
#endif
        ! ==================================================================== !
        ! MPI path — ISEND/IRECV into flat c_wrk_cx, then scatter to a.       !
        ! ==================================================================== !
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
                ! ISEND/IRECV in batches; recv into flat c_wrk_cx.
                ! For FABRIC_DIRECT: ips/ipr are K-ranks (= peer pro_k); global = pro_k*npro_i + ims_pro_i.
                do j = 1, ims_npro_k, trp_sizBlock_k
                    l = 0
                    do m = j, min(j + trp_sizBlock_k - 1, ims_npro_k)
                        ns = maps_recv_k(m) + 1; ips = ns - 1   ! backward: send/recv maps swapped
                        nr = maps_send_k(m) + 1; ipr = nr - 1
                        send_to   = ips
                        recv_from = ipr
                        l = l + 1
                        call MPI_ISEND(b(trp_plan%disp_r(ns) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, send_to, fbd_tag, trp_comm_k, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(c_wrk_cx((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, recv_from, fbd_tag, trp_comm_k, request(l), ims_err)
                    end do
                    call MPI_WAITALL(l, request, status, ims_err)
                end do
                ! Scatter: flat c_wrk_cx → strided a
                do m = 1, ims_npro_k
                    nr = maps_send_k(m) + 1
                    flat_off = (nr - 1)*nmax_p*nlines_p
                    disp_nr  = trp_plan%disp_s(nr)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            a(disp_nr + i*npage + j + 1) = c_wrk_cx(flat_off + i*nlines_p + j + 1)
                        end do
                    end do
                end do
                nullify (c_wrk_cx)
            else
                call Transpose_Kernel_Complex(b, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                              a, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                              trp_comm_k, trp_sizBlock_k, trp_mode_k)
            end if
#ifdef USE_APU
        end if   ! end APU/CPU dispatch
#endif
        return
    end subroutine TLabMPI_Trp_ExecK_Backward_Complex

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Forward_Real(a, b, trp_plan)
        ! I-Forward: scatter from flat X-space (a) to strided I-space (b). Analogous to K-Forward.
        ! a layout: a(m*chunk + i) for peer m (flat, contiguous per peer).
        ! b layout: b(m*nmax_p + i*nmax_full + j) for peer m, element i, line j (strided).
        real(wp), dimension(:), intent(in) :: a
        real(wp), dimension(:), intent(out) :: b
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr
        integer(wi) :: nmax_p, nlines_p, nmax_full, flat_off, disp_nr, mas
#ifdef USE_APU
        integer(c_int) :: hip_sync_err
        integer(8) :: off                     ! int64 window offset for the V1 fused write+fence push
#endif

        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i   ! total X-elements per line (stride in b)
        mas       = nmax_p * nlines_p

        ! ==================================================================== !
        ! APU paths — GPU direct writes between shared-memory windows.         !
        ! a is flat (chunk per peer); each rank pushes a[m*chunk] to peer m's  !
        ! recv buffer at slot own_rank*chunk. After fence, apu_recv_fptr_i has  !
        ! all data; unpack to strided b.                                        !
        ! ==================================================================== !
#ifdef USE_APU
        if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT) then
            ! -- Push: a is flat; one fused kernel writes our chunk to ALL peers simultaneously.
            size = trp_plan%size3d
            call MPI_Win_fence(0, apu_win_i, ims_err)
#ifdef TRP_I_FUSEDFENCE
            ! V1 (apudirect): cross-rank push via the fused hip_write_with_fence (write + __threadfence_system
            ! in ONE wavefront = system-scope L2 write-back), mirroring the fabricdirect node-window fix. The
            ! source a(m*mas+1:..) is contiguous per peer; apu_all_i is contiguous. No hipDeviceSynchronize: a is
            ! produced by the caller's blocking !$omp target (no nowait in the codebase) and hip_write_with_fence
            ! is itself synchronous, so the push reads committed a and the epoch closes after it. (Re-add a sync
            ! ONLY if the producer/these targets go nowait.)
            do m = 0, ims_npro_i - 1
                off = int(m, 8)*int(apu_stride_i, 8) + int(ims_pro_i, 8)*int(mas, 8)
                call hip_write_with_fence(a(m*mas + 1:m*mas + mas), apu_all_i(off + 1:off + mas), int(mas, c_int))
            end do
#else
            !$omp target teams distribute parallel do collapse(2)
            do m = 0, ims_npro_i - 1
                do i = 1, nmax_p * nlines_p
                    ! a[m*chunk] is the flat chunk destined for peer m; write to their recv slot
                    apu_all_i(m*apu_stride_i + ims_pro_i*nmax_p*nlines_p + i) = a(m * nmax_p * nlines_p + i)
                end do
            end do
            !$omp end target teams distribute parallel do
            hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#endif
            call MPI_Win_fence(0, apu_win_i, ims_err)   ! barrier: recv buffer fully populated
#ifdef TRP_I_FUSEDFENCE
            call hip_invalidate_recv(apu_recv_fptr_i, int(size, c_int))   ! reader acquire: invalidate L2 -> fresh MALL (pairs with the fused writer release)
#endif
            ! -- Unpack: recv buffer holds sorted flat chunks; scatter to strided b in one fused kernel.
            !$omp target teams distribute parallel do collapse(3)
            do m = 0, ims_npro_i - 1
                do i = 0, nlines_p - 1
                    do j = 0, nmax_p - 1
                        ! recv slot m*chunk + i*nmax_p + j maps to strided b position
                        b(m*nmax_p + i*nmax_full + j + 1) = &
                            apu_recv_fptr_i(m*nmax_p*nlines_p + i*nmax_p + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do

        else if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
            size = trp_plan%size3d
            if (use_node_win_i .and. all(is_intra_i(0:ims_npro_i - 1))) then
                ! Option 1: all I-peers intra-node -> node-window GPU writes (apudirect-I pattern), no MPI.
                ! ----- FINE LOCALIZATION (forward node-window I-transpose; mirrors NWBR). gated by trp_dbg_fft,
                ! which the RHS now also sets around the seeding self-burgX OPR_Burgers_X call (the 2026-06-24
                ! origin). 'NWFR' = Node-Window Forward Real. a = input; node_recv_fptr_i = MY recv buffer =
                ! ims_npro_i source-slots of 'mas' each (slot s pushed by source rank s). Decision tree:
                !   NWFR:in clean on ALL I-comm peers but NWFR:wcpu/seg garbage -> the cross-rank GPU PUSH
                !       corrupts healthy data (node-window WRITE bug); seg index = which source slot.
                !   NWFR:wcpu/win/out all clean but self-burgX still garbage -> the FDM, not the transpose.
                if (trp_dbg_fft) DNS_PROBE('NWFR:in', a(1), size, -1)
                call MPI_Win_fence(0, node_win_i, ims_err)
#ifdef TRP_I_FUSEDFENCE
                ! V1: each peer's contiguous slot is WRITTEN AND fenced in ONE wavefront by the tested
                ! hip_write_with_fence kernel (write a -> peer's window slot, then __threadfence_system in the
                ! same wavefront). a(m*mas+1:..) is contiguous per peer; node_all_i is contiguous.
                ! 'a' is produced by the caller on the OpenMP offload stream; the HIP push runs on the HIP
                ! stream -> sync first so the push reads the final 'a', not an in-flight value.
                hip_sync_err = hipDeviceSynchronize()
                do m = 0, ims_npro_i - 1
                    off = int(node_lrank_i(m), 8)*int(apu_size_i, 8) + int(ims_pro_i, 8)*int(mas, 8)
                    call hip_write_with_fence(a(m*mas + 1:m*mas + mas), node_all_i(off + 1:off + mas), int(mas, c_int))
                end do
                ! CRITICAL: hip_write_with_fence is an ASYNC kernel launch. The in-kernel __threadfence_system
                ! orders the write but does NOT make the host wait. Without this sync, the MPI_Win_fence below
                ! (a CPU/MPI barrier that does not block on GPU kernels) closes the RMA epoch while the pushes
                ! are still in flight -> peers read stale/partial window data. Wait for completion here.
                hip_sync_err = hipDeviceSynchronize()
#else
                ! ONE fused GPU write over all I-peers (all intra-node) — collapse(2) over (m,i).
                !$omp target teams distribute parallel do collapse(2)
                do m = 0, ims_npro_i - 1
                    do i = 1, mas
                        node_all_i(int(node_lrank_i(m),8)*int(apu_size_i,8) + ims_pro_i*mas + i) = a(m*mas + i)
                    end do
                end do
                !$omp end target teams distribute parallel do
                hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#ifdef TRP_I_SYSFENCE
                call hip_system_fence()   ! FAST FIX: system-scope L2 write-back -> cross-rank push visible in MALL
#endif
#endif
                call MPI_Win_fence(0, node_win_i, ims_err)
                if (trp_dbg_fft) then
                    DNS_PROBE_CPU('NWFR:wcpu', node_recv_fptr_i(1), size, -1)   ! CPU read of HBM, first
                    DNS_PROBE('NWFR:win', node_recv_fptr_i(1), size, -1)
                    do m = 0, ims_npro_i - 1
                        DNS_PROBE('NWFR:seg', node_recv_fptr_i(m*mas + 1), mas, m)
                    end do
                end if
#ifdef TRP_I_READER_ACQUIRE
                call hip_invalidate_recv(node_recv_fptr_i, int(size, c_int))   ! reader acquire: invalidate L2 -> fresh MALL (pairs with SYSFENCE or FUSEDFENCE writer release)
#endif
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_i - 1
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            b(m*nmax_p + i*nmax_full + j + 1) = node_recv_fptr_i(m*mas + i*nmax_p + j + 1)
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                if (trp_dbg_fft) DNS_PROBE('NWFR:out', b(1), size, -1)
            else
                ! Fallback: all-MPI on fabric_mpi_comm_i (clean MPI_COMM_WORLD split). hip_write_with_fence
                ! flushes GPU L2 -> HBM before the CPU ISENDs.
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                l = 0
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(c_wrk_dp((nr - 1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ipr, ims_tag, fabric_mpi_comm_i, request(l), ims_err)
                end do
                call hip_write_with_fence(a(1:size), wrk_mpi_dp(size + 1), int(size, c_int))
                do m = 1, ims_npro_i
                    ns = maps_send_i(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(wrk_mpi_dp(size + (ns - 1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ips, ims_tag, fabric_mpi_comm_i, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1
                    flat_off = (nr - 1)*nmax_p*nlines_p
                    disp_nr  = trp_plan%disp_r(nr)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            b(disp_nr + i*nmax_full + j + 1) = c_wrk_dp(flat_off + i*nmax_p + j + 1)
                        end do
                    end do
                end do
                nullify (c_wrk_dp)
            end if

        else   ! CPU paths: ASYNCHRONOUS, SENDRECV, ALLTOALL
#endif
        ! ==================================================================== !
        ! CPU paths — MPI ISEND/IRECV (ASYNCHRONOUS), SENDRECV, ALLTOALL.     !
        ! ==================================================================== !
            if (trp_datatype_i == MPI_REAL4 .and. wp == dp) then
                ! Single-precision path: dp→sp before send, unpack merges sp→dp.
                size = trp_plan%size3d
                if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    a_wrk => wrk_mpi_fptr(1:size)           ! send: dp→sp copy of a
                    c_wrk => wrk_mpi_fptr(size + 1:2*size)  ! recv: flat staging
#ifdef USE_APU
                    !$omp target teams distribute parallel do
#endif
                    do i = 1, size
                        a_wrk(i) = real(a(i), sp)   ! dp→sp
                    end do
#ifdef USE_APU
                    !$omp end target teams distribute parallel do
#endif
                    ! Step 1: IRECVs
                    l = 0
                    do m = 1, ims_npro_i
                        nr = maps_recv_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_IRECV(c_wrk((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, MPI_REAL4, &
                                       ipr, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    ! Step 2: ISENDs from flat a_wrk
                    do m = 1, ims_npro_i
                        ns = maps_send_i(m) + 1; ips = ns - 1
                        l = l + 1
                        call MPI_ISEND(a_wrk(trp_plan%disp_s(ns) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                       ips, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    ! Step 3: WAITALL
                    call MPI_WAITALL(l, request, status, ims_err)
                    ! Step 4: GPU unpack flat c_wrk → strided b (sp→dp merged)
                    do m = 1, ims_npro_i
                        nr = maps_recv_i(m) + 1
                        flat_off = (nr - 1)*nmax_p*nlines_p
                        disp_nr  = trp_plan%disp_r(nr)
#ifdef USE_APU
                        !$omp target teams distribute parallel do collapse(2)
#endif
                        do i = 0, nlines_p - 1
                            do j = 0, nmax_p - 1
                                b(disp_nr + i*nmax_full + j + 1) = real(c_wrk(flat_off + i*nmax_p + j + 1), dp)
                            end do
                        end do
#ifdef USE_APU
                        !$omp end target teams distribute parallel do
#endif
                    end do
                    nullify (a_wrk, c_wrk)
                else
                    a_wrk => wrk_mpi_fptr(1:size)
                    b_wrk => wrk_mpi_fptr(size + 1:2*size)
#ifdef USE_APU
                    !$omp target teams distribute parallel do
#endif
                    do i = 1, size
                        a_wrk(i) = real(a(i), sp)   ! dp→sp
                    end do
#ifdef USE_APU
                    !$omp end target teams distribute parallel do
#endif
                    call Transpose_Kernel_Single(a_wrk, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                                 b_wrk, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                                 ims_comm_x, trp_sizBlock_i, trp_mode_i)
#ifdef USE_APU
                    !$omp target teams distribute parallel do
#endif
                    do i = 1, size
                        b(i) = real(b_wrk(i), dp)   ! sp→dp
                    end do
#ifdef USE_APU
                    !$omp end target teams distribute parallel do
#endif
                    nullify (a_wrk, b_wrk)
                end if

            else
                ! Double-precision path
                if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    ! a is flat; CPU pack into wrk_mpi_dp second half before ISEND.
                    ! Direct ISEND from GPU-written a would send stale HBM on cross-node paths;
                    ! CPU pack reads a via hardware cache coherency (correct on MI300A), then
                    ! MPI flushes CPU cache to HBM before NIC DMA.
                    size = trp_plan%size3d
                    call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                    ! Step 1: IRECVs into first half of wrk_mpi_dp
                    l = 0
                    do m = 1, ims_npro_i
                        nr = maps_recv_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_IRECV(c_wrk_dp((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    ! Step 2: CPU pack a → second half of wrk_mpi_dp
                    do m = 1, ims_npro_i
                        ns    = maps_send_i(m) + 1
                        flat_off = (ns - 1)*nmax_p*nlines_p
                        disp_nr  = trp_plan%disp_s(ns)
                        do i = 1, nmax_p*nlines_p
                            wrk_mpi_dp(size + flat_off + i) = a(disp_nr + i)
                        end do
                    end do
                    ! Step 3: ISENDs from second half (CPU-written, coherent with NIC DMA)
                    do m = 1, ims_npro_i
                        ns = maps_send_i(m) + 1; ips = ns - 1
                        l = l + 1
                        call MPI_ISEND(wrk_mpi_dp(size + (ns - 1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    ! Step 4: WAITALL
                    call MPI_WAITALL(l, request, status, ims_err)
                    ! Step 5: GPU scatter flat c_wrk_dp → strided b
                    do m = 1, ims_npro_i
                        nr = maps_recv_i(m) + 1
                        flat_off = (nr - 1)*nmax_p*nlines_p
                        disp_nr  = trp_plan%disp_r(nr)
#ifdef USE_APU
                        !$omp target teams distribute parallel do collapse(2)
#endif
                        do i = 0, nlines_p - 1
                            do j = 0, nmax_p - 1
                                b(disp_nr + i*nmax_full + j + 1) = c_wrk_dp(flat_off + i*nmax_p + j + 1)
                            end do
                        end do
#ifdef USE_APU
                        !$omp end target teams distribute parallel do
#endif
                    end do
                    nullify (c_wrk_dp)
                else
                    call Transpose_Kernel_Double(a(1:trp_plan%size3d), maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                                 b(1:trp_plan%size3d), maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                                 ims_comm_x, trp_sizBlock_i, trp_mode_i)
                end if
            end if
#ifdef USE_APU
        end if   ! end APU/CPU dispatch
#endif

        return
    end subroutine TLabMPI_Trp_ExecI_Forward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Forward_Complex(a, b, trp_plan)
        ! I-Forward (complex): flat I-space chunks in a → strided X-space layout in b.
        ! a(m*chunk + i) where chunk = nmax_p*nlines_p, m = peer index, i = local element.
        ! b(m*nmax_p + i*nmax_full + j + 1): strided with full-width stride nmax_full = nmax_p*npro.
        ! b must be passed full-size [nmax_full*nlines] by the caller (NOT a single column like
        ! c_out(:,1)); this routine indexes it linearly across all nlines. A single-column slice
        ! makes every write with i>=1 a formal out-of-bounds access (harmless at -O0,
        ! miscompiled/aborted at -O2). Kept assumed-shape (:) — Cray OpenMP target rejects (*).
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan
        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_nr
        integer :: send_to, recv_from, fbd_tag   ! FABRIC_DIRECT: global rank + distinct tag
#ifdef USE_APU
        ! apu_cx_all: complex view spanning all peers' shared windows (stride = apu_stride_i/2 complex units).
        complex(dp), pointer :: apu_cx_all(:) => null()
        integer(wi) :: mas
        integer(c_int) :: hip_sync_err   ! FABRIC_DIRECT fallback flush
#endif
        type(MPI_Comm) :: trp_comm_i   ! MPI_COMM_WORLD for FABRIC_DIRECT (see K-Forward_Complex comment).
#ifdef USE_APU
        if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
            trp_comm_i = fabric_mpi_comm_i   ! clean MPI_COMM_WORLD split; local rank = I-dir rank
            fbd_tag = ims_tag
        else
            trp_comm_i = ims_comm_x
            fbd_tag = ims_tag
        end if
#else
        trp_comm_i = ims_comm_x
        fbd_tag = ims_tag
#endif

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i

        ! ==================================================================== !
        ! APU paths — GPU direct writes between shared-memory windows.         !
        ! ==================================================================== !
#ifdef USE_APU
#ifndef TRP_CX_MPI
        if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT) then
            ! a is flat I-space chunks; fused GPU kernel pushes to all peers' recv buffers,
            ! then unpacks our own recv buffer (strided) into b.
            size = trp_plan%size3d
            mas  = nmax_p * nlines_p
            ! Build complex-typed view spanning all peers' contiguous window segments.
            call c_f_pointer(apu_peer_cptr_i(0), apu_cx_all, [apu_stride_i*ims_npro_i/2])
            ! Fence 1: open epoch — all ranks ready to receive direct writes.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! DEBUG (crash hunt): window state BEFORE our push (right after open fence). If this is already
            ! corrupt (>5) the bad value is a PRE-EXISTING clobber/leftover (a prior apu_win_i user, or a
            ! coverage gap), NOT produced by this epoch's push/fence. If clean, the corruption is born here.
            if (trp_dbg_fft) DNS_PROBE('X-FWC-pre', apu_recv_fptr_i(1), 2*size, -1)
            ! Push: write each peer's flat chunk into peer m's buffer at slot own_rank*chunk.
            !$omp target teams distribute parallel do collapse(2)
            do m = 0, ims_npro_i - 1
                do i = 1, nmax_p * nlines_p
                    apu_cx_all(m*(apu_stride_i/2) + ims_pro_i*nmax_p*nlines_p + i) = &
                        a(m * nmax_p * nlines_p + i)
                end do
            end do
            !$omp end target teams distribute parallel do
            ! ROOT FIX: flush the GPU push to HBM BEFORE the fence. MPI_Win_fence orders the RMA epoch but
            ! does NOT flush GPU writes on MI300A, so peers would read our cross-rank window writes stale
            ! -> the intermittent FFT blow-up. hipDeviceSynchronize makes it HBM-visible.
            hip_sync_err = hipDeviceSynchronize()
#ifdef TRP_I_FUSEDFENCE
            ! apudirect complex (the originally-documented apudirect seed, ExecI_Forward_Complex): device-scope
            ! hipDeviceSynchronize is NOT a system-scope L2 write-back, so commit L2->MALL before the close
            ! fence (writer release). The unpack below now reads the window on the GPU, paired with a reader-side
            ! hip_invalidate_recv acquire after the close fence (the all-GPU replacement for the old CPU read).
            call hip_system_fence()
#endif
            ! Fence 2: close epoch — all writes committed; recv buffers fully populated.
            call MPI_Win_fence(0, apu_win_i, ims_err)
#ifdef TRP_I_FUSEDFENCE
            ! reader acquire: invalidate this rank's GPU L2 before the GPU unpack reads the window, so it reloads
            ! the freshly-committed MALL data instead of a stale line from a previous substep. This is the proper
            ! GPU fix that supersedes the 2026-06-22 CPU-read workaround; mirrors ExecI_Backward_Complex. Real
            ! alias of the complex window, 2*size reals.
            call hip_invalidate_recv(apu_recv_fptr_i, int(2*size, c_int))
#endif
            ! DEBUG (X-FFT region): recv window AFTER cross-rank push+fence, BEFORE unpack. If this jumps the
            ! corruption is in the push/fence even WITH the flush (apu_recv_fptr_i = real alias; 2*size reals).
            if (trp_dbg_fft) DNS_PROBE('X-FWC-win', apu_recv_fptr_i(1), 2*size, -1)
            ! DEBUG (crash hunt): per-peer-segment window max, sub=segment index l. OUR window is npro_i
            ! contiguous complex segments of mas=nmax_p*nlines_p each; segment l holds the chunk pushed by
            ! I-rank l. All sources are <=5, so a >5 segment localizes the faulting push: l==ims_pro_i is our
            ! OWN local write (a local kernel/memory bug); l/=ims_pro_i is I-rank l's CROSS-rank push (a
            ! fence/visibility race). Real view: segment l spans reals [2*l*mas+1 .. 2*(l+1)*mas].
            if (trp_dbg_fft) then
                do l = 0, ims_npro_i - 1
                    DNS_PROBE('X-FWC-seg', apu_recv_fptr_i(2*l*mas + 1), 2*mas, l)
                end do
            end if
            ! Unpack on the GPU: scatter recv buffer (flat m*chunk+i layout) → b (strided m*nmax_p + i*nmax_full + j).
            ! HISTORY: the 2026-06-22 fix moved this READ to the CPU because the GPU shared-window cross-rank read
            ! returned STALE GPU L2 (rank0->rank2 seg0 read 204.7 from a <=5 source -> deterministic blow-up at
            ! it=241751). The writer push was SOUND; the bug was purely the reader's stale L2. The reader-side
            ! hip_invalidate_recv acquire added after the close fence (above) refreshes that L2, so the unpack now
            ! stays on the GPU (all-on-GPU, faster, exploits the single-node shared window) and feeds b to the CPU
            ! c2r FFTW after the flush below. Mirrors ExecI_Backward_Complex. (Revert to the CPU loop only if a
            ! long apudirect run past it=241751 re-blows -> the invalidate would be insufficient on this path.)
            !$omp target teams distribute parallel do collapse(3)
            do m = 0, ims_npro_i - 1
                do i = 0, nlines_p - 1
                    do j = 0, nmax_p - 1
                        b(m*nmax_p + i*nmax_full + j + 1) = &
                            apu_cx_recv_fptr_i(m*nmax_p*nlines_p + i*nmax_p + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            ! b is GPU-written, consumed by the CPU FFTW in OPR_Fourier_X_Forward -> flush for CPU coherence.
            hip_sync_err = hipDeviceSynchronize()
            nullify (apu_cx_all)

        else if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
#else
        if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT .or. trp_mode_i == TLAB_MPI_TRP_APU_DIRECT) then
#endif
            ! I-Forward complex FABRIC_DIRECT. For npro_i=6 all I-peers are intra-node (one XCD), so
            ! Option 1 (use_node_win_i .and. all intra) goes entirely through node-local complex
            ! shared-window GPU writes + MPI_Win_fence (apudirect-I pattern), zero MPI. b is GPU-unpacked
            ! and consumed by the LOCAL CPU FFTW (coherent, as the real-I node path already is). Else:
            ! all-MPI complex fallback on fabric_mpi_comm_i. Complex view: apu_size_i/2 cx per segment.
            size = trp_plan%size3d
            mas  = nmax_p * nlines_p
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
            if (use_node_win_i .and. all(is_intra_i(0:ims_npro_i - 1))) then
                call MPI_Win_fence(0, node_win_i, ims_err)
                ! ONE fused GPU write over all I-peers (all intra-node) — collapse(2) over (m,i).
                !$omp target teams distribute parallel do collapse(2)
                do m = 0, ims_npro_i - 1
                    do i = 1, mas
                        node_cx_all_i(int(node_lrank_i(m),8)*int(apu_size_i/2,8) + ims_pro_i*mas + i) = a(m*mas + i)
                    end do
                end do
                !$omp end target teams distribute parallel do
                hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#ifdef TRP_I_SYSFENCE
                call hip_system_fence()   ! FAST FIX: system-scope L2 write-back -> cross-rank push visible in MALL
#endif
                call MPI_Win_fence(0, node_win_i, ims_err)
#ifdef TRP_I_SYSFENCE
                call hip_invalidate_recv(node_recv_fptr_i, int(2*size, c_int))   ! reader-side (real view of cx): fresh MALL before unpack
#endif
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_i - 1
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            b(m*nmax_p + i*nmax_full + j + 1) = node_cx_recv_fptr_i(m*mas + i*nmax_p + j + 1)
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                ! b is GPU-unpacked but consumed by the CPU FFTW; flush so the GPU write is CPU-coherent.
                hip_sync_err = hipDeviceSynchronize()
            else
                ! Fallback: all-MPI complex on fabric_mpi_comm_i. Flush GPU-written a for the CPU ISEND
                ! (needed if the fallback ever carries an inter-node peer; harmless when all intra).
                hip_sync_err = hipDeviceSynchronize()
                do j = 1, ims_npro_i, trp_sizBlock_i
                    l = 0
                    do m = j, min(j + trp_sizBlock_i - 1, ims_npro_i)
                        ns = maps_send_i(m) + 1; ips = ns - 1
                        nr = maps_recv_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_ISEND(a(trp_plan%disp_s(ns) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, fabric_mpi_comm_i, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(c_wrk_cx((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, fabric_mpi_comm_i, request(l), ims_err)
                    end do
                    call MPI_WAITALL(l, request, status, ims_err)
                end do
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1
                    flat_off = (nr - 1)*nmax_p*nlines_p
                    disp_nr  = trp_plan%disp_r(nr)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            b(disp_nr + i*nmax_full + j + 1) = c_wrk_cx(flat_off + i*nmax_p + j + 1)
                        end do
                    end do
                end do
            end if
            nullify (c_wrk_cx)

        else   ! ASYNCHRONOUS complex; SENDRECV/ALLTOALL go to kernel.
#endif
        ! ==================================================================== !
        ! CPU paths — MPI ISEND/IRECV, SENDRECV, ALLTOALL                     !
        ! ==================================================================== !
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! Flat send of each peer's chunk from a; recv into flat staging c_wrk_cx;
                ! then scatter c_wrk_cx → strided b.
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
                ! Pipeline: batched ISEND+IRECV → WAITALL per trp_sizBlock_i batch.
                ! trp_comm_i is fabric_mpi_comm_i (local I-comm rank = dir rank directly).
                do j = 1, ims_npro_i, trp_sizBlock_i
                    l = 0
                    do m = j, min(j + trp_sizBlock_i - 1, ims_npro_i)
                        ns = maps_send_i(m) + 1; ips = ns - 1
                        nr = maps_recv_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_ISEND(a(trp_plan%disp_s(ns) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, fbd_tag, trp_comm_i, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(c_wrk_cx((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, fbd_tag, trp_comm_i, request(l), ims_err)
                    end do
                    call MPI_WAITALL(l, request, status, ims_err)
                end do
                ! Scatter flat c_wrk_cx → strided b.
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1
                    flat_off = (nr - 1)*nmax_p*nlines_p
                    disp_nr  = trp_plan%disp_r(nr)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            b(disp_nr + i*nmax_full + j + 1) = c_wrk_cx(flat_off + i*nmax_p + j + 1)
                        end do
                    end do
                end do
                nullify (c_wrk_cx)
            else
                call Transpose_Kernel_Complex(a, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                              b, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                              trp_comm_i, trp_sizBlock_i, trp_mode_i)
            end if
#ifdef USE_APU
        end if   ! end APU/CPU dispatch
#endif
        return
    end subroutine TLabMPI_Trp_ExecI_Forward_Complex

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Backward_Real(b, a, trp_plan)
        ! I-Backward: inverse of I-Forward. Reverses strided X-space (b) → flat I-space (a).
        ! b(m*nmax_p + i*nmax_full + j + 1): strided X-space layout, nmax_full = nmax_p*npro.
        ! a: flat output, layout mirrors the forward input (m*chunk + local elements).
        real(wp), intent(in) :: b(:)
        real(wp), intent(out) :: a(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan
        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_ns, mas
#ifdef USE_APU
        integer(c_int) :: hip_sync_err
        integer(8) :: off                     ! int64 window offset for the V1 fused write+fence push
#endif

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i
        mas       = nmax_p * nlines_p

        ! ==================================================================== !
        ! APU paths — GPU direct writes between shared-memory windows.         !
        ! ==================================================================== !
#ifdef USE_APU
        if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT) then
            ! b is strided X-space (disp_r stride). Each rank packs its slice from b into
            ! every peer m's recv buffer at slot own_rank*chunk. Single fused collapse(3)
            ! kernel covers all m in one HIP launch, eliminating per-peer launch overhead.
            size = trp_plan%size3d
            ! DEBUG (X-FFT region): input b BEFORE the transpose (if this jumps, the corruption is upstream).
            if (trp_dbg_fft) DNS_PROBE('X-BWR-in', b(1), size, -1)
            ! Fence 1: open epoch — all ranks ready to receive direct writes.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! Push: pack strided b[m] → peer m's recv buffer at slot own_rank*chunk (flat).
#ifdef TRP_I_FUSEDFENCE
            ! V1 (apudirect): strided source b -> GPU-gather into contiguous staging (wrk_mpi_dp), then fused
            ! hip_write_with_fence per peer (write + __threadfence_system in one wavefront). No hipDeviceSynchronize:
            ! the blocking !$omp end target below already completed the gather and hip_write_with_fence is itself
            ! synchronous, so the push reads committed staging and the epoch closes after it. (Re-add a sync ONLY
            ! if these target regions become nowait.)
            !$omp target teams distribute parallel do collapse(3)
            do m = 0, ims_npro_i - 1
                do i = 0, nlines_p - 1
                    do j = 0, nmax_p - 1
                        wrk_mpi_dp(m*mas + i*nmax_p + j + 1) = b(m*nmax_p + i*nmax_full + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            do m = 0, ims_npro_i - 1
                off = int(m, 8)*int(apu_stride_i, 8) + int(ims_pro_i, 8)*int(mas, 8)
                call hip_write_with_fence(wrk_mpi_dp(m*mas + 1:m*mas + mas), apu_all_i(off + 1:off + mas), int(mas, c_int))
            end do
#else
            !$omp target teams distribute parallel do collapse(3)
            do m = 0, ims_npro_i - 1
                do i = 0, nlines_p - 1
                    do j = 0, nmax_p - 1
                        apu_all_i(m*apu_stride_i + ims_pro_i*nmax_p*nlines_p + i*nmax_p + j + 1) = &
                            b(m*nmax_p + i*nmax_full + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            hip_sync_err = hipDeviceSynchronize()   ! ROOT FIX: flush GPU push to HBM before the fence (cross-rank visibility)
#endif
            ! Fence 2: close epoch — all writes committed; recv buffers fully populated.
            call MPI_Win_fence(0, apu_win_i, ims_err)
#ifdef TRP_I_FUSEDFENCE
            call hip_invalidate_recv(apu_recv_fptr_i, int(size, c_int))   ! reader acquire: invalidate L2 -> fresh MALL (pairs with the fused writer release)
#endif
            ! DEBUG (X-FFT region): recv window AFTER push+fence (if this jumps but X-BWR-in was clean -> push/fence).
            if (trp_dbg_fft) DNS_PROBE('X-BWR-win', apu_recv_fptr_i(1), size, -1)
            ! Flat copy: recv buffer layout is flat and matches a 1:1.
            !$omp target teams distribute parallel do
            do i = 1, size
                a(i) = apu_recv_fptr_i(i)
            end do
            !$omp end target teams distribute parallel do
            ! DEBUG (X-FFT region): output a AFTER unpack (if X-BWR-win clean and this jumps -> unpack kernel).
            if (trp_dbg_fft) DNS_PROBE('X-BWR-out', a(1), size, -1)
        else if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
            size = trp_plan%size3d
            if (use_node_win_i .and. all(is_intra_i(0:ims_npro_i - 1))) then
                ! Option 1: all I-peers intra-node -> node-window GPU push (apudirect-I backward), no MPI.
                ! ----- FINE LOCALIZATION of the 2026-06-24 fabricdirect seed (this node-window read-back).
                ! gated by trp_dbg_fft => fires ONLY in the Poisson OPR_Fourier_X_Backward (where it seeded),
                ! not in the per-substep OPR_Partial_X calls (keeps log volume sane). 'NWBR' = Node-Window
                ! Backward Real. b is the (flushed) c2r FFTW output = healthy input; node_recv_fptr_i is MY
                ! recv segment, = ims_npro_i source-slots of 'mas' each (slot s pushed by source rank s).
                if (trp_dbg_fft) DNS_PROBE('NWBR:in', b(1), size, -1)
                call MPI_Win_fence(0, node_win_i, ims_err)
#ifdef TRP_I_FUSEDFENCE
                ! V1: the backward source b is STRIDED, so gather it into a contiguous staging buffer
                ! (wrk_mpi_dp, per-peer contiguous), then WRITE+fence each peer slot in one wavefront via the
                ! tested hip_write_with_fence kernel. (Same data path as the !$omp push, but the cross-rank
                ! write and its __threadfence_system live in the same wavefront -> reliable system flush.)
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_i - 1
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            wrk_mpi_dp(m*mas + i*nmax_p + j + 1) = b(m*nmax_p + i*nmax_full + j + 1)
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                hip_sync_err = hipDeviceSynchronize()   ! gather (OMP stream) must complete before the HIP push reads wrk_mpi_dp
                do m = 0, ims_npro_i - 1
                    off = int(node_lrank_i(m), 8)*int(apu_size_i, 8) + int(ims_pro_i, 8)*int(mas, 8)
                    call hip_write_with_fence(wrk_mpi_dp(m*mas + 1:m*mas + mas), node_all_i(off + 1:off + mas), int(mas, c_int))
                end do
                ! CRITICAL: same async-launch race as ExecI_Forward_Real -- wait for the push+fence kernels to
                ! COMPLETE before the MPI_Win_fence closes the RMA epoch (CPU barrier does not block GPU kernels).
                hip_sync_err = hipDeviceSynchronize()
#else
                ! ONE fused GPU push over all I-peers (all intra-node) — collapse(3) over (m,i,j).
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_i - 1
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            node_all_i(int(node_lrank_i(m),8)*int(apu_size_i,8) + ims_pro_i*mas + i*nmax_p + j + 1) = &
                                b(m*nmax_p + i*nmax_full + j + 1)
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#ifdef TRP_I_SYSFENCE
                call hip_system_fence()   ! FAST FIX: system-scope L2 write-back -> cross-rank push visible in MALL
#endif
#endif
                call MPI_Win_fence(0, node_win_i, ims_err)
                ! WINDOW STATE after push+fence, BEFORE the GPU read-back. Decision tree on the next crash:
                !   NWBR:wcpu garbage              -> writer side: a peer pushed bad data / push mis-addressed.
                !   NWBR:wcpu clean but NWBR:win or NWBR:out garbage
                !                                  -> reader side: GPU read-back returns STALE L2 cache (HBM is
                !                                     correct); NWBR:seg index = which source slot is stale.
                ! NWBR:wcpu MUST come first (CPU read of HBM, uncontaminated by any GPU read of the window).
                if (trp_dbg_fft) then
                    DNS_PROBE_CPU('NWBR:wcpu', node_recv_fptr_i(1), size, -1)
                    DNS_PROBE('NWBR:win', node_recv_fptr_i(1), size, -1)
                    do m = 0, ims_npro_i - 1
                        DNS_PROBE('NWBR:seg', node_recv_fptr_i(m*mas + 1), mas, m)
                    end do
                end if
#ifdef TRP_I_READER_ACQUIRE
                call hip_invalidate_recv(node_recv_fptr_i, int(size, c_int))   ! reader acquire: invalidate L2 -> fresh MALL (pairs with SYSFENCE or FUSEDFENCE writer release)
#endif
                !$omp target teams distribute parallel do
                do i = 1, size
                    a(i) = node_recv_fptr_i(i)
                end do
                !$omp end target teams distribute parallel do
                if (trp_dbg_fft) DNS_PROBE('NWBR:out', a(1), size, -1)
            else
                ! Fallback: all-MPI on fabric_mpi_comm_i (clean MPI_COMM_WORLD split).
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                l = 0
                do m = 1, ims_npro_i
                    nr = maps_send_i(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(a(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ipr, ims_tag, fabric_mpi_comm_i, request(l), ims_err)
                end do
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_r(ns)
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            wrk_mpi_dp(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end do
                hip_sync_err = hipDeviceSynchronize()
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp((ns - 1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ips, ims_tag, fabric_mpi_comm_i, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
                nullify (c_wrk_dp)
            end if
        else   ! CPU paths
#endif
        ! ==================================================================== !
        ! CPU paths — MPI ISEND/IRECV, SENDRECV, ALLTOALL                     !
        ! ==================================================================== !
            if (trp_datatype_i == MPI_REAL4 .and. wp == dp) then
                ! Single-precision path: dp→sp pack before send, sp→dp unpack after recv.
                size = trp_plan%size3d
                if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    ! Pipeline: post IRECVs → GPU dp→sp+strided→flat pack → ISENDs → WAITALL → sp→dp.
                    a_wrk => wrk_mpi_fptr(1:size)           ! recv: flat sp staging → sp→dp to a
                    c_wrk => wrk_mpi_fptr(size + 1:2*size)  ! send: packed from b dp→sp+strided→flat
                    ! Step 1: post all IRECVs into flat a_wrk recv slots before pack starts.
                    l = 0
                    do m = 1, ims_npro_i
                        nr = maps_send_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_IRECV(a_wrk(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                       ipr, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    ! Step 2: GPU pack b→c_wrk (dp→sp+strided→flat) while recv buffers prepare.
                    do m = 1, ims_npro_i
                        ns = maps_recv_i(m) + 1
                        flat_off = (ns - 1)*nmax_p*nlines_p
                        disp_ns  = trp_plan%disp_r(ns)
#ifdef USE_APU
                        !$omp target teams distribute parallel do collapse(2)
#endif
                        do i = 0, nlines_p - 1
                            do j = 0, nmax_p - 1
                                c_wrk(flat_off + i*nmax_p + j + 1) = real(b(disp_ns + i*nmax_full + j + 1), sp)
                            end do
                        end do
#ifdef USE_APU
                        !$omp end target teams distribute parallel do
#endif
                    end do
                    ! Step 3: post all ISENDs from packed flat c_wrk.
                    do m = 1, ims_npro_i
                        ns = maps_recv_i(m) + 1; ips = ns - 1
                        l = l + 1
                        call MPI_ISEND(c_wrk((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, MPI_REAL4, &
                                       ips, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    ! Step 4: single WAITALL for all sends and receives.
                    call MPI_WAITALL(l, request, status, ims_err)
                    ! sp→dp: flat 1:1 conversion of recv buffer to a.
#ifdef USE_APU
                    !$omp target teams distribute parallel do
#endif
                    do i = 1, size
                        a(i) = real(a_wrk(i), dp)
                    end do
#ifdef USE_APU
                    !$omp end target teams distribute parallel do
#endif
                    nullify (a_wrk, c_wrk)
                else
                    ! SENDRECV/ALLTOALL sp path: dp→sp, kernel, sp→dp.
                    b_wrk => wrk_mpi_fptr(1:size)
                    a_wrk => wrk_mpi_fptr(size + 1:2*size)
                    ! dp→sp conversion of b into b_wrk.
#ifdef USE_APU
                    !$omp target teams distribute parallel do
#endif
                    do i = 1, size
                        b_wrk(i) = real(b(i), sp)
                    end do
#ifdef USE_APU
                    !$omp end target teams distribute parallel do
#endif
                    call Transpose_Kernel_Single(b_wrk, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                                 a_wrk, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                                 ims_comm_x, trp_sizBlock_i, trp_mode_i)
                    ! sp→dp conversion of a_wrk into a.
#ifdef USE_APU
                    !$omp target teams distribute parallel do
#endif
                    do i = 1, size
                        a(i) = real(a_wrk(i), dp)
                    end do
#ifdef USE_APU
                    !$omp end target teams distribute parallel do
#endif
                    nullify (a_wrk, b_wrk)
                end if
            else
                ! Double-precision path.
                if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                    ! Pipeline: post IRECVs → CPU strided→flat pack → ISENDs → WAITALL.
                    size = trp_plan%size3d
                    call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                    ! Step 1: post all IRECVs into flat a recv slots before pack starts.
                    l = 0
                    do m = 1, ims_npro_i
                        nr = maps_send_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_IRECV(a(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    ! Step 2: CPU pack b→c_wrk_dp (strided→flat); no !$omp target (Fix G pattern).
                    ! GPU-written b is readable by CPU via cache coherency on MI300A.
                    ! !$omp target pack followed by MPI_ISEND sends stale HBM on cross-node paths.
                    do m = 1, ims_npro_i
                        ns = maps_recv_i(m) + 1
                        flat_off = (ns - 1)*nmax_p*nlines_p
                        disp_ns  = trp_plan%disp_r(ns)
                        do i = 0, nlines_p - 1
                            do j = 0, nmax_p - 1
                                c_wrk_dp(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                            end do
                        end do
                    end do
                    ! Step 3: post all ISENDs from packed flat c_wrk_dp.
                    do m = 1, ims_npro_i
                        ns = maps_recv_i(m) + 1; ips = ns - 1
                        l = l + 1
                        call MPI_ISEND(c_wrk_dp((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    ! Step 4: single WAITALL for all sends and receives.
                    call MPI_WAITALL(l, request, status, ims_err)
                    nullify (c_wrk_dp)
                else
                    call Transpose_Kernel_Double(b(1:trp_plan%size3d), maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                                 a(1:trp_plan%size3d), maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                                 ims_comm_x, trp_sizBlock_i, trp_mode_i)
                end if
            end if
#ifdef USE_APU
        end if   ! end APU/CPU dispatch
#endif

        return
    end subroutine TLabMPI_Trp_ExecI_Backward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Backward_Complex(b, a, trp_plan)
        ! I-Backward (complex): inverse of I-Forward. Reverses strided X-space (b) → flat I-space (a).
        ! b(m*nmax_p + i*nmax_full + j + 1): strided X-space, nmax_full = nmax_p*npro.
        ! a: flat output, chunks of nmax_p*nlines_p per peer (mirrors forward input layout).
        ! b must be passed full-size [nmax_full*nlines] by the caller (NOT a single column like
        ! wrk1(:,1)); it is read linearly across all nlines here. A single-column slice makes
        ! every read with i>=1 a formal out-of-bounds access. Kept (:) — Cray target rejects (*).
        complex(wp), intent(in) :: b(:)
        complex(wp), intent(out) :: a(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan
        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_ns
        integer :: send_to, recv_from, fbd_tag   ! FABRIC_DIRECT: global rank + distinct tag
#ifdef USE_APU
        ! apu_cx_all: complex view spanning all peers' shared windows (stride = apu_stride_i/2 complex units).
        complex(dp), pointer :: apu_cx_all(:) => null()
        integer(wi) :: mas
        integer(c_int) :: hip_sync_err   ! flush GPU-assembled a before the CPU FFTW reads it
#endif
        type(MPI_Comm) :: trp_comm_i   ! MPI_COMM_WORLD for FABRIC_DIRECT (see K-Forward_Complex comment).
#ifdef USE_APU
        if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
            trp_comm_i = fabric_mpi_comm_i
            fbd_tag = ims_tag
        else
            trp_comm_i = ims_comm_x
            fbd_tag = ims_tag
        end if
#else
        trp_comm_i = ims_comm_x
        fbd_tag = ims_tag
#endif

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i

        ! ==================================================================== !
        ! APU paths — GPU direct writes between shared-memory windows.         !
        ! ==================================================================== !
#ifdef USE_APU
#ifndef TRP_CX_MPI
        if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT) then
            ! b is strided X-space. Fused collapse(3) kernel packs all peers in one HIP launch;
            ! each peer m's slot at own_rank*chunk receives this rank's strided b[m] data.
            size = trp_plan%size3d
            mas  = nmax_p * nlines_p
            ! Build complex-typed view spanning all peers' contiguous window segments.
            call c_f_pointer(apu_peer_cptr_i(0), apu_cx_all, [apu_stride_i*ims_npro_i/2])
            ! Fence 1: open epoch — all ranks ready to receive direct writes.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! Push: pack strided b[m] → peer m's recv buffer at slot own_rank*chunk.
            !$omp target teams distribute parallel do collapse(3)
            do m = 0, ims_npro_i - 1
                do i = 0, nlines_p - 1
                    do j = 0, nmax_p - 1
                        apu_cx_all(m*(apu_stride_i/2) + ims_pro_i*nmax_p*nlines_p + i*nmax_p + j + 1) = &
                            b(m*nmax_p + i*nmax_full + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            hip_sync_err = hipDeviceSynchronize()   ! ROOT FIX: flush GPU push to HBM before the fence (cross-rank visibility)
#ifdef TRP_I_FUSEDFENCE
            ! apudirect complex backward: writer-side system-scope commit (L2->MALL) before the close fence.
            call hip_system_fence()
#endif
            ! Fence 2: close epoch — all writes committed; recv buffers fully populated.
            call MPI_Win_fence(0, apu_win_i, ims_err)
#ifdef TRP_I_FUSEDFENCE
            ! ...and reader-side: this unpack reads the window on the GPU, so invalidate its L2 first
            ! (apu_recv_fptr_i = real alias of the complex window; 2*size reals).
            call hip_invalidate_recv(apu_recv_fptr_i, int(2*size, c_int))
#endif
            ! Flat copy: recv buffer is flat and matches a layout 1:1.
            !$omp target teams distribute parallel do
            do i = 1, size
                a(i) = apu_cx_recv_fptr_i(i)
            end do
            !$omp end target teams distribute parallel do
            ! a is GPU-written but consumed by the CPU FFTW in OPR_Fourier_X_Backward; flush so the GPU
            ! write is CPU-coherent before the FFTW reads it (same GPU->CPU hazard as the K-transpose).
            hip_sync_err = hipDeviceSynchronize()
            nullify (apu_cx_all)

        else if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
#else
        if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT .or. trp_mode_i == TLAB_MPI_TRP_APU_DIRECT) then
#endif
            ! I-Backward complex FABRIC_DIRECT. For npro_i=6 all I-peers are intra-node, so Option 1
            ! (use_node_win_i .and. all intra) goes entirely through node-local complex shared-window GPU
            ! push + MPI_Win_fence (apudirect-I backward pattern), zero MPI; flat-copy our segment → a.
            ! b is CPU-written (FFTW), GPU-read here (coherent on unified memory, as the apudirect path is);
            ! a is GPU-written (as apudirect/real-I backward). Else: all-MPI complex fallback (no flush
            ! needed — b and the CPU-packed c_wrk_cx are both CPU-resident). Complex view: apu_size_i/2 cx/seg.
            size = trp_plan%size3d
            mas  = nmax_p * nlines_p
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
            if (use_node_win_i .and. all(is_intra_i(0:ims_npro_i - 1))) then
                call MPI_Win_fence(0, node_win_i, ims_err)
                ! ONE fused GPU push over all I-peers (all intra-node) — collapse(3) over (m,i,j).
                !$omp target teams distribute parallel do collapse(3)
                do m = 0, ims_npro_i - 1
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            node_cx_all_i(int(node_lrank_i(m),8)*int(apu_size_i/2,8) + ims_pro_i*mas + i*nmax_p + j + 1) = &
                                b(m*nmax_p + i*nmax_full + j + 1)
                        end do
                    end do
                end do
                !$omp end target teams distribute parallel do
                hip_sync_err = hipDeviceSynchronize()   ! flush GPU push to HBM before the fence (cross-rank visibility)
#ifdef TRP_I_SYSFENCE
                call hip_system_fence()   ! FAST FIX: system-scope L2 write-back -> cross-rank push visible in MALL
#endif
                call MPI_Win_fence(0, node_win_i, ims_err)
#ifdef TRP_I_SYSFENCE
                call hip_invalidate_recv(node_recv_fptr_i, int(2*size, c_int))   ! reader-side (real view of cx): fresh MALL before unpack
#endif
                !$omp target teams distribute parallel do
                do i = 1, size
                    a(i) = node_cx_recv_fptr_i(i)
                end do
                !$omp end target teams distribute parallel do
                ! a is GPU-written but consumed by the CPU FFTW; flush so the GPU write is CPU-coherent.
                hip_sync_err = hipDeviceSynchronize()
            else
                ! Fallback: all-MPI complex on fabric_mpi_comm_i. Pack strided b → flat c_wrk_cx, ISEND;
                ! IRECV directly into flat a (disp_s(m+1) = m*mas on the fabric comm).
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_r(ns)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            c_wrk_cx(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                end do
                do j = 1, ims_npro_i, trp_sizBlock_i
                    l = 0
                    do m = j, min(j + trp_sizBlock_i - 1, ims_npro_i)
                        ns = maps_recv_i(m) + 1; ips = ns - 1
                        nr = maps_send_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_ISEND(c_wrk_cx((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, fabric_mpi_comm_i, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(a(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, fabric_mpi_comm_i, request(l), ims_err)
                    end do
                    call MPI_WAITALL(l, request, status, ims_err)
                end do
            end if
            nullify (c_wrk_cx)

        else   ! ASYNCHRONOUS complex; SENDRECV/ALLTOALL go to kernel.
#endif
        ! ==================================================================== !
        ! CPU paths — MPI ISEND/IRECV, SENDRECV, ALLTOALL                     !
        ! ==================================================================== !
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! Pack strided b → flat c_wrk_cx; batched ISEND+IRECV → WAITALL per block.
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
                ! Gather strided b → flat c_wrk_cx staging (per-peer).
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_r(ns)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            c_wrk_cx(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                end do
                ! Batched ISEND from flat c_wrk_cx + IRECV directly into a → WAITALL.
                ! For FABRIC_DIRECT: ips/ipr are I-ranks (= peer pro_i); global = ims_pro_k*npro_i + pro_i.
                do j = 1, ims_npro_i, trp_sizBlock_i
                    l = 0
                    do m = j, min(j + trp_sizBlock_i - 1, ims_npro_i)
                        ns = maps_recv_i(m) + 1; ips = ns - 1
                        nr = maps_send_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_ISEND(c_wrk_cx((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, fbd_tag, trp_comm_i, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(a(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, fbd_tag, trp_comm_i, request(l), ims_err)
                    end do
                    call MPI_WAITALL(l, request, status, ims_err)
                end do
                nullify (c_wrk_cx)
            else
                call Transpose_Kernel_Complex(b, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                              a, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                              trp_comm_i, trp_sizBlock_i, trp_mode_i)
            end if
#ifdef USE_APU
        end if   ! end APU/CPU dispatch
#endif

        return
    end subroutine TLabMPI_Trp_ExecI_Backward_Complex

    !########################################################################
    !########################################################################
    subroutine Transpose_Kernel_Double(a, msend, dsend, tsend, b, mrecv, drecv, trecv, comm, step, mode)
        ! Assumed-shape (not assumed-size) so Cray OpenMP target accepts these arrays.
        ! Callers must pass explicit array sections, e.g. a(1:size3d).
        real(wp), intent(in)  :: a(:)
        real(wp), intent(out) :: b(:)

        type(MPI_Comm), intent(in) :: comm                         ! communicator
        type(MPI_Datatype), intent(in) :: tsend, trecv             ! types send/receive
        integer(wi), intent(in) :: dsend(:), drecv(:)              ! displacements send/receive
        integer(wi), intent(in) :: msend(:), mrecv(:)              ! maps send/receive
        integer(wi), intent(in) :: step
        integer, intent(in) :: mode

        ! -----------------------------------------------------------------------
        integer(wi) npro
        integer(wi) j, l, m, ns, nr, ips, ipr

        ! #######################################################################
        npro = size(dsend(:))

        select case (mode)
        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, status(1), ims_err)
                end do
            end do

        case (TLAB_MPI_TRP_ALLTOALL)
            types_send(1:npro) = tsend
            types_recv(1:npro) = trecv
            call MPI_ALLTOALLW(a, counts, dsend*int(sizeof(1.0_wp)), types_send, &
                               b, counts, drecv*int(sizeof(1.0_wp)), types_recv, comm, ims_err)
        end select

        return
    end subroutine Transpose_Kernel_Double

    !########################################################################
    !########################################################################
    subroutine Transpose_Kernel_Single(a, msend, dsend, tsend, b, mrecv, drecv, trecv, comm, step, mode)
        ! Assumed-shape (not assumed-size) so Cray OpenMP target accepts these arrays.
        real(sp), intent(in) :: a(:)
        real(sp), intent(out) :: b(:)

        type(MPI_Comm), intent(in) :: comm                         ! communicator
        type(MPI_Datatype), intent(in) :: tsend, trecv             ! types send/receive
        integer(wi), intent(in) :: dsend(:), drecv(:)              ! displacements send/receive
        integer(wi), intent(in) :: msend(:), mrecv(:)              ! maps send/receive
        integer(wi), intent(in) :: step
        integer, intent(in) :: mode

        ! -----------------------------------------------------------------------
        integer(wi) npro
        integer(wi) j, l, m, ns, nr, ips, ipr

        ! #######################################################################
        npro = size(dsend(:))

        select case (mode)
        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, status(1), ims_err)
                end do
            end do

        case (TLAB_MPI_TRP_ALLTOALL)
            types_send(1:npro) = tsend
            types_recv(1:npro) = trecv
            call MPI_ALLTOALLW(a, counts, dsend*int(sizeof(1.0_sp)), types_send, &
                               b, counts, drecv*int(sizeof(1.0_sp)), types_recv, comm, ims_err)
        end select

        return
    end subroutine Transpose_Kernel_Single

    !########################################################################
    !########################################################################
    subroutine Transpose_Kernel_Complex(a, msend, dsend, tsend, b, mrecv, drecv, trecv, comm, step, mode)
        complex(wp), intent(in) :: a(*)
        complex(wp), intent(out) :: b(*)

        type(MPI_Comm), intent(in) :: comm                         ! communicator
        type(MPI_Datatype), intent(in) :: tsend, trecv             ! types send/receive
        integer(wi), intent(in) :: dsend(:), drecv(:)              ! displacements send/receive
        integer(wi), intent(in) :: msend(:), mrecv(:)              ! maps send/receive
        integer(wi), intent(in) :: step
        integer, intent(in) :: mode

        ! -----------------------------------------------------------------------
        integer(wi) npro
        integer(wi) j, l, m, ns, nr, ips, ipr

        ! #######################################################################
        npro = size(dsend(:))

        select case (mode)
        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, status(1), ims_err)
                end do
            end do

        case (TLAB_MPI_TRP_ALLTOALL)
            types_send(1:npro) = tsend
            types_recv(1:npro) = trecv
            call MPI_ALLTOALLW(a, counts, dsend*int(sizeof(1.0_sp)), types_send, &
                               b, counts, drecv*int(sizeof(1.0_sp)), types_recv, comm, ims_err)
        end select

        return
    end subroutine Transpose_Kernel_Complex

end module TLabMPI_Transpose
