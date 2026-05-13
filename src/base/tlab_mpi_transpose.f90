#include "dns_error.h"

! Circular transposition within directional communicators
module TLabMPI_Transpose
#ifdef USE_APU
    use omp_lib, only: omp_get_default_device, omp_target_alloc, omp_target_free
#endif
    use TLab_Constants, only: lfile, efile, wp, dp, sp, wi, sizeofreal
    use TLab_Memory, only: imax, jmax, kmax, isize_wrk3d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: TLab_Allocate_Real
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc, c_null_ptr, c_associated, c_size_t, c_sizeof
    ! c_ptr / c_intptr_t are accessible via mpi_f08 (which re-exports iso_c_binding); re-declaring causes ambiguity.
    ! For debug address dumps below we reuse MPI_ADDRESS_KIND (8-byte) instead of c_intptr_t.
    use TLabMPI_VARS
#ifdef USE_APU
    ! FABRIC_DIRECT debug taps (per-rank checksum logs in fort.500+rank or debug_thread_testing<rank>.log).
    ! Diff the per-rank files between an asynchronous run (working reference) and a fabricdirect run on
    ! Hunter; the first checkpoint where sums differ identifies where the data is going wrong.
    use Tlab_Debug, only: TLab_Debug_Print_1D
#endif
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
    integer, parameter :: TLAB_MPI_TRP_APU_DIRECT = 4    ! APU: direct GPU writes to peer device buffers, MPI_Barrier sync
    integer, parameter :: TLAB_MPI_TRP_APU_ASYNC  = 5    ! APU: intra-node direct writes + inter-node MPI ISEND/IRECV
    integer, parameter :: TLAB_MPI_TRP_FABRIC_DIRECT = 6 ! APU: intra-node direct writes + inter-node one-sided MPI_Put (RMA)

#ifdef USE_APU
    ! APU direct mode state: shared-memory MPI windows, one per direction.
    ! Each rank allocates its recv buffer as a shared segment so all ranks
    ! in the communicator can write directly to each other's buffers via
    ! MPI_Win_shared_query pointers (valid in every process's address space).
    integer(wi) :: apu_size_k = 0_wi, apu_size_i = 0_wi
    type(MPI_Win) :: apu_win_k, apu_win_i
    real(dp), pointer :: apu_recv_fptr_k(:) => null(), apu_recv_fptr_i(:) => null()
    type(c_ptr), allocatable :: apu_peer_cptr_k(:), apu_peer_cptr_i(:) ! per-rank pointers from Win_shared_query
    ! APU_ASYNC state: node-local shared windows + intra-node rank detection
    integer(wi) :: apu_async_size_k = 0_wi, apu_async_size_i = 0_wi
    type(MPI_Win) :: apu_async_win_k, apu_async_win_i
    ! contiguous: required so the buffer passed to MPI gets the real data address, not a
    ! copy-in/copy-out temporary (c_ptr address ≠ Fortran-pointer descriptor address).
    real(dp), pointer, contiguous :: apu_async_recv_k(:) => null(), apu_async_recv_i(:) => null()
    type(c_ptr), allocatable :: apu_async_peer_k(:), apu_async_peer_i(:)
    logical, allocatable :: apu_async_is_local_k(:), apu_async_is_local_i(:)
    ! FABRIC_DIRECT inter-node leg: plain two-sided MPI (ISEND/IRECV/WAITALL) into apu_async_recv_*,
    ! with all packing/unpacking done on the CPU. One-sided RMA into the shm/MPI_Win_allocate window
    ! was tried and does not deliver inter-node on Cray MPICH (remote puts land nowhere). Intra-node
    ! peers still use the GPU shared-window direct writes (the apudirect mechanism).
    ! Complex-typed aliases for the same shared windows; used by complex APU_DIRECT paths.
    complex(dp), pointer :: apu_cx_recv_fptr_k(:) => null(), apu_cx_recv_fptr_i(:) => null()
    ! Contiguous span over all peers' recv segments for fused single-kernel writes.
    ! apu_stride_k/i = apu_size_k/i (segment size in dp elements; peers are allocated contiguously).
    ! apu_all_k/i(m*stride + 1 : (m+1)*stride) = rank m's recv buffer.
    integer(wi) :: apu_stride_k = 0_wi, apu_stride_i = 0_wi
    real(dp), pointer :: apu_all_k(:) => null(), apu_all_i(:) => null()
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
        type(MPI_Comm)  :: apu_async_shmem_comm
        type(MPI_Group) :: apu_async_group_dir, apu_async_group_shmem
        integer, allocatable :: apu_async_src_ranks(:), apu_async_trans_ranks(:)
        integer(MPI_ADDRESS_KIND) :: dbg_addr   ! FABRIC_DIRECT debug: holds a transferred c_ptr address
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
        elseif (trim(adjustl(sRes)) == 'apuasync') then
#ifdef USE_APU
            trp_mode_i = TLAB_MPI_TRP_APU_ASYNC
#else
            call TLab_Write_ASCII(efile, __FILE__//'. TransposeModeI=apuasync requires USE_APU.')
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
        elseif (trim(adjustl(sRes)) == 'apuasync') then
#ifdef USE_APU
            trp_mode_k = TLAB_MPI_TRP_APU_ASYNC
#else
            call TLab_Write_ASCII(efile, __FILE__//'. TransposeModeK=apuasync requires USE_APU.')
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
            trp_mode_i == TLAB_MPI_TRP_APU_ASYNC    .or. trp_mode_k == TLAB_MPI_TRP_APU_ASYNC    .or. &
            trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT .or. trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
            allocate (wrk_mpi_dp(2*imax*jmax*kmax))
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated dp flat-MPI staging buffer (2 x size3d).')
        end if

#ifdef USE_APU
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
                                         int(c_sizeof(1.0_dp)), MPI_INFO_NULL, ims_comm_z, win_baseptr, apu_win_k, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU K recv buffer.')
                call TLab_Stop(DNS_ERROR_ALLOC)
            end if
            call c_f_pointer(win_baseptr, apu_recv_fptr_k, [apu_size_k])
            call c_f_pointer(win_baseptr, apu_cx_recv_fptr_k, [apu_size_k/2])
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
                                         int(c_sizeof(1.0_dp)), MPI_INFO_NULL, ims_comm_x, win_baseptr, apu_win_i, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU I recv buffer.')
                call TLab_Stop(DNS_ERROR_ALLOC)
            end if
            call c_f_pointer(win_baseptr, apu_recv_fptr_i, [apu_size_i])
            call c_f_pointer(win_baseptr, apu_cx_recv_fptr_i, [apu_size_i/2])
            allocate (apu_peer_cptr_i(0:ims_npro_i - 1))
            do ip = 0, ims_npro_i - 1
                call MPI_Win_shared_query(apu_win_i, ip, win_query_size, win_disp_unit, apu_peer_cptr_i(ip), ims_err)
            end do
            ! apu_all_i spans all peers' contiguous segments (MPI-3 shared windows are always contiguous).
            apu_stride_i = apu_size_i
            call c_f_pointer(apu_peer_cptr_i(0), apu_all_i, [apu_stride_i*ims_npro_i])
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated I APU direct recv buffer.')
        end if
        ! APU_ASYNC / FABRIC_DIRECT: node-local shared windows using MPI_COMM_TYPE_SHARED sub-communicators.
        ! Each rank splits ims_comm_z/x into a node-local communicator; the shared window
        ! is allocated only over that sub-communicator. Ranks absent from it are inter-node
        ! peers, handled via standard ISEND/IRECV (APU_ASYNC) or one-sided MPI_Put (FABRIC_DIRECT,
        ! which additionally exposes the recv buffer as an RMA window over the full transpose comm).
        if ((trp_mode_k == TLAB_MPI_TRP_APU_ASYNC .or. trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) .and. ims_npro_k > 1) then
            apu_async_size_k = int(imax, wi)*int(jmax, wi)*int(kmax, wi)
            allocate (apu_async_is_local_k(0:ims_npro_k - 1))
            allocate (apu_async_peer_k(0:ims_npro_k - 1))
            apu_async_peer_k = c_null_ptr
            ! Build node-local communicator and detect intra-node K peers
            call MPI_Comm_split_type(ims_comm_z, MPI_COMM_TYPE_SHARED, ims_pro_k, MPI_INFO_NULL, apu_async_shmem_comm, ims_err)
            call MPI_Comm_group(ims_comm_z, apu_async_group_dir, ims_err)
            call MPI_Comm_group(apu_async_shmem_comm, apu_async_group_shmem, ims_err)
            allocate (apu_async_src_ranks(0:ims_npro_k - 1), apu_async_trans_ranks(0:ims_npro_k - 1))
            do ip = 0, ims_npro_k - 1; apu_async_src_ranks(ip) = ip; end do
            call MPI_Group_translate_ranks(apu_async_group_dir, ims_npro_k, apu_async_src_ranks, &
                                           apu_async_group_shmem, apu_async_trans_ranks, ims_err)
            do ip = 0, ims_npro_k - 1
                apu_async_is_local_k(ip) = (apu_async_trans_ranks(ip) /= MPI_UNDEFINED)
            end do
            call MPI_Group_free(apu_async_group_dir, ims_err)
            call MPI_Group_free(apu_async_group_shmem, ims_err)
            ! Allocate shared window on the node-local sub-communicator
            call MPI_Win_allocate_shared(int(apu_async_size_k, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
                                         int(c_sizeof(1.0_dp)), MPI_INFO_NULL, apu_async_shmem_comm, win_baseptr, apu_async_win_k, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU_ASYNC K recv buffer.')
                call TLab_Stop(DNS_ERROR_ALLOC)
            end if
            call c_f_pointer(win_baseptr, apu_async_recv_k, [apu_async_size_k])
            ! The rank argument of MPI_Win_shared_query is a rank in the window's group
            ! (the node-local shmem comm), not in ims_comm_z — translate it.
            do ip = 0, ims_npro_k - 1
                if (apu_async_is_local_k(ip)) &
                    call MPI_Win_shared_query(apu_async_win_k, apu_async_trans_ranks(ip), &
                                              win_query_size, win_disp_unit, apu_async_peer_k(ip), ims_err)
            end do
            deallocate (apu_async_src_ranks, apu_async_trans_ranks)
            call MPI_Comm_free(apu_async_shmem_comm, ims_err)
            if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
                ! debug: shm recv-buffer addresses (c_ptr from MPI_Win_allocate_shared vs c_loc of the
                ! Fortran pointer over it — MUST match) and the intra/inter peer list.
                dbg_addr = transfer(win_baseptr, dbg_addr)
                write(500 + ims_pro, *) '[INIT_FBD_K] PE', ims_pro, ' shm_baseptr=', dbg_addr, &
                    ' size=', apu_async_size_k, ' npro_k=', ims_npro_k, ' pro_k=', ims_pro_k
                dbg_addr = transfer(c_loc(apu_async_recv_k(1)), dbg_addr)
                write(500 + ims_pro, *) '[INIT_FBD_K] PE', ims_pro, ' c_loc(apu_async_recv_k(1))=', dbg_addr, &
                    '  (must equal shm_baseptr above)'
                do ip = 0, ims_npro_k - 1
                    if (apu_async_is_local_k(ip)) then
                        dbg_addr = transfer(apu_async_peer_k(ip), dbg_addr)
                        write(500 + ims_pro, *) '[INIT_FBD_K] PE', ims_pro, ' INTRA-node peer k-rank', ip, &
                            ' shared_query_cptr=', dbg_addr
                    else
                        write(500 + ims_pro, *) '[INIT_FBD_K] PE', ims_pro, ' INTER-node peer k-rank', ip
                    end if
                end do
                flush(500 + ims_pro)
            end if
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated K APU_ASYNC/FABRIC_DIRECT recv buffer.')
        end if
        if ((trp_mode_i == TLAB_MPI_TRP_APU_ASYNC .or. trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) .and. ims_npro_i > 1) then
            apu_async_size_i = int(imax, wi)*int(jmax, wi)*int(kmax, wi)
            allocate (apu_async_is_local_i(0:ims_npro_i - 1))
            allocate (apu_async_peer_i(0:ims_npro_i - 1))
            apu_async_peer_i = c_null_ptr
            call MPI_Comm_split_type(ims_comm_x, MPI_COMM_TYPE_SHARED, ims_pro_i, MPI_INFO_NULL, apu_async_shmem_comm, ims_err)
            call MPI_Comm_group(ims_comm_x, apu_async_group_dir, ims_err)
            call MPI_Comm_group(apu_async_shmem_comm, apu_async_group_shmem, ims_err)
            allocate (apu_async_src_ranks(0:ims_npro_i - 1), apu_async_trans_ranks(0:ims_npro_i - 1))
            do ip = 0, ims_npro_i - 1; apu_async_src_ranks(ip) = ip; end do
            call MPI_Group_translate_ranks(apu_async_group_dir, ims_npro_i, apu_async_src_ranks, &
                                           apu_async_group_shmem, apu_async_trans_ranks, ims_err)
            do ip = 0, ims_npro_i - 1
                apu_async_is_local_i(ip) = (apu_async_trans_ranks(ip) /= MPI_UNDEFINED)
            end do
            call MPI_Group_free(apu_async_group_dir, ims_err)
            call MPI_Group_free(apu_async_group_shmem, ims_err)
            call MPI_Win_allocate_shared(int(apu_async_size_i, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
                                         int(c_sizeof(1.0_dp)), MPI_INFO_NULL, apu_async_shmem_comm, win_baseptr, apu_async_win_i, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU_ASYNC I recv buffer.')
                call TLab_Stop(DNS_ERROR_ALLOC)
            end if
            call c_f_pointer(win_baseptr, apu_async_recv_i, [apu_async_size_i])
            ! The rank argument of MPI_Win_shared_query is a rank in the window's group
            ! (the node-local shmem comm), not in ims_comm_x — translate it.
            do ip = 0, ims_npro_i - 1
                if (apu_async_is_local_i(ip)) &
                    call MPI_Win_shared_query(apu_async_win_i, apu_async_trans_ranks(ip), &
                                              win_query_size, win_disp_unit, apu_async_peer_i(ip), ims_err)
            end do
            deallocate (apu_async_src_ranks, apu_async_trans_ranks)
            call MPI_Comm_free(apu_async_shmem_comm, ims_err)
            if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
                dbg_addr = transfer(win_baseptr, dbg_addr)
                write(500 + ims_pro, *) '[INIT_FBD_I] PE', ims_pro, ' shm_baseptr=', dbg_addr, &
                    ' size=', apu_async_size_i, ' npro_i=', ims_npro_i, ' pro_i=', ims_pro_i
                dbg_addr = transfer(c_loc(apu_async_recv_i(1)), dbg_addr)
                write(500 + ims_pro, *) '[INIT_FBD_I] PE', ims_pro, ' c_loc(apu_async_recv_i(1))=', dbg_addr, &
                    '  (must equal shm_baseptr above)'
                do ip = 0, ims_npro_i - 1
                    if (apu_async_is_local_i(ip)) then
                        dbg_addr = transfer(apu_async_peer_i(ip), dbg_addr)
                        write(500 + ims_pro, *) '[INIT_FBD_I] PE', ims_pro, ' INTRA-node peer i-rank', ip, &
                            ' shared_query_cptr=', dbg_addr
                    else
                        write(500 + ims_pro, *) '[INIT_FBD_I] PE', ims_pro, ' INTER-node peer i-rank', ip
                    end if
                end do
                flush(500 + ims_pro)
            end if
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated I APU_ASYNC/FABRIC_DIRECT recv buffer.')
        end if
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
        real(dp), pointer :: apu_pfptr_k(:) => null()   ! scratch pointer for APU_ASYNC per-peer writes
        real(dp) :: dbg_intra, dbg_inter   ! FABRIC_DIRECT debug: recv-buffer split checksums
        integer(MPI_ADDRESS_KIND) :: dbg_addr
        integer(wi) :: fp_nerr               ! FINGERPRINT TEST: mismatch counter
        real(dp) :: fp_expected, fp_actual   ! FINGERPRINT TEST: expected and actual values
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

        ! Common-path debug tap: input checksum (runs for ALL modes, diff async vs fabricdirect).
#ifdef USE_APU
        call TLab_Debug_Print_1D('[KFR_pre] sum(a)=', a)
#endif

        ! ==================================================================== !
        ! APU paths — GPU direct writes between shared-memory windows.         !
        ! All ranks allocate their recv buffer via MPI_Win_allocate_shared so  !
        ! every peer can write directly into it without MPI send/recv.         !
        ! apu_all_k(m*stride + own_rank*chunk + 1:...) = our data for peer m. !
        ! MPI_Win_fence acts as the collective barrier between push and unpack. !
        ! ==================================================================== !
#ifdef USE_APU
        if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
            ! -- Push: one fused GPU kernel writes our chunk to ALL peers simultaneously.
            ! Eliminates per-peer kernel-launch overhead (npro separate launches → 1).
            size = trp_plan%size3d
            call MPI_Win_fence(0, apu_win_k, ims_err)
            !$omp target teams distribute parallel do collapse(3) if(mas * sizeofreal > 100000_wi)
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
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: our recv buffer is now fully populated
            ! -- Unpack: recv buffer is already in the flat K-space layout; one-to-one copy to b.
            !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
            do i = 1, size
                b(i) = apu_recv_fptr_k(i)
            end do
            !$omp end target teams distribute parallel do

        else if (trp_mode_k == TLAB_MPI_TRP_APU_ASYNC) then
            ! Hybrid: intra-node peers via direct shared-memory writes (as APU_DIRECT above);
            ! inter-node peers via MPI ISEND/IRECV. Both paths fill apu_async_recv_k,
            ! which is then flat-copied to b.
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
            l = 0
            ! Step 1: post IRECVs for inter-node peers before opening the shared-window epoch
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    l = l + 1
                    call MPI_IRECV(apu_async_recv_k(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_z, request(l), ims_err)
                end if
            end do
            ! Step 2: open shared-window epoch (collective over node-local communicator)
            call MPI_Win_fence(0, apu_async_win_k, ims_err)
            ! Step 3: intra-node — push our strided chunk directly to each peer's recv buffer slot
            do m = 0, ims_npro_k - 1
                if (apu_async_is_local_k(m)) then
                    call c_f_pointer(apu_async_peer_k(m), apu_pfptr_k, [apu_async_size_k])
                    flat_off = ims_pro_k * nmax_p * nlines_p   ! our write slot in peer's buffer
                    disp_ns  = m * nlines_p
                    !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            apu_pfptr_k(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end do
            ! Step 4: inter-node — pack strided chunk into flat staging buffer and ISEND
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    flat_off = m * nmax_p * nlines_p
                    disp_ns  = m * nlines_p
                    !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            c_wrk_dp(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp(flat_off + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_z, request(l), ims_err)
                end if
            end do
            ! Step 5: close epoch (all intra-node writes done) and wait for inter-node MPI
            call MPI_Win_fence(0, apu_async_win_k, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! Step 6: recv buffer is fully populated; flat copy to b
            !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
            do i = 1, size
                b(i) = apu_async_recv_k(i)
            end do
            !$omp end target teams distribute parallel do
            nullify (c_wrk_dp, apu_pfptr_k)

        else if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
            ! Two-level: intra-node peers via direct shared-window GPU writes (the apudirect mechanism);
            ! inter-node peers via plain two-sided MPI (ISEND/IRECV/WAITALL) — the same primitives the
            ! ASYNCHRONOUS mode uses, which work on Hunter (one-sided RMA into the shm/MPI_Win_allocate
            ! recv buffer does NOT deliver inter-node on Cray MPICH). Everything that touches the MPI
            ! buffers — the inter-node strided pack and the inter-node slots of the final copy — is done
            ! on the CPU to avoid reading DMA-written data through a stale GPU cache (the suspected APU
            ! corruption mode); the intra-node slots stay on the GPU.
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
            l = 0
            ! Step 1: post IRECVs for inter-node peers into their slots in the (shm) recv buffer
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    l = l + 1
                    call MPI_IRECV(apu_async_recv_k(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_z, request(l), ims_err)
                end if
            end do
            ! Step 2: open shared-window epoch (collective over the node-local communicator)
            call MPI_Win_fence(0, apu_async_win_k, ims_err)
            ! FINGERPRINT TEST step 3: intra-node — CPU writes fingerprint to peer's recv buffer at our slot.
            ! fingerprint(pos) = our_rank * 1e9 + pos  (pos is 1-based within the send chunk)
            do m = 0, ims_npro_k - 1
                if (apu_async_is_local_k(m)) then
                    call c_f_pointer(apu_async_peer_k(m), apu_pfptr_k, [apu_async_size_k])
                    flat_off = ims_pro_k * nmax_p * nlines_p
                    do i = 1, nmax_p * nlines_p
                        apu_pfptr_k(flat_off + i) = real(ims_pro_k, dp) * 1.0e9_dp + real(i, dp)
                    end do
                end if
            end do
            ! FINGERPRINT TEST step 4: inter-node — fill staging with fingerprint, then ISEND
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    flat_off = m * nmax_p * nlines_p
                    do i = 1, nmax_p * nlines_p
                        c_wrk_dp(flat_off + i) = real(ims_pro_k, dp) * 1.0e9_dp + real(i, dp)
                    end do
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp(flat_off + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_z, request(l), ims_err)
                end if
            end do
            ! Step 5: close shared-window epoch; wait for inter-node MPI
            call MPI_Win_fence(0, apu_async_win_k, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! FINGERPRINT VERIFICATION: slot m in apu_async_recv_k must hold data from rank m.
            ! expected(m, i) = m * 1e9 + i  for i = 1..nmax_p*nlines_p
            fp_nerr = 0
            do m = 0, ims_npro_k - 1
                flat_off = m * nmax_p * nlines_p
                do i = 1, nmax_p * nlines_p
                    fp_expected = real(m, dp) * 1.0e9_dp + real(i, dp)
                    fp_actual   = apu_async_recv_k(flat_off + i)
                    if (abs(fp_actual - fp_expected) > 0.5_dp) then
                        fp_nerr = fp_nerr + 1
                        if (fp_nerr <= 20) then
                            write(500 + ims_pro, *) '[KFR_FP_ERR] PE', ims_pro, &
                                ' from_rank=', m, ' pos=', i, &
                                ' expected=', fp_expected, ' got=', fp_actual, &
                                ' local=', apu_async_is_local_k(m)
                        end if
                    end if
                end do
            end do
            write(500 + ims_pro, *) '[KFR_FP] PE', ims_pro, &
                ' chunk=', nmax_p*nlines_p, ' total_errors=', fp_nerr
            flush(500 + ims_pro)
            nullify (c_wrk_dp, apu_pfptr_k)
            STOP 'KFR fingerprint test done -- check fort.500+ logs'

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
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                        !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
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
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                    ! Step 2: GPU gather: a (strided) → c_wrk_dp (flat)
                    do m = 1, ims_npro_k
                        ns = maps_send_k(m) + 1
                        flat_off = (ns - 1)*nmax_p*nlines_p
                        disp_ns  = trp_plan%disp_s(ns)
#ifdef USE_APU
                        !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
#endif
                        do i = 0, nmax_p - 1
                            do j = 0, nlines_p - 1
                                c_wrk_dp(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                            end do
                        end do
#ifdef USE_APU
                        !$omp end target teams distribute parallel do
#endif
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

        ! Common-path debug tap: output checksum (runs for ALL modes; first divergence vs the
        ! asynchronous-reference log identifies which K-forward call goes wrong in fabricdirect).
#ifdef USE_APU
        call TLab_Debug_Print_1D('[KFR_post] sum(b)=', b)
#endif
        return
    end subroutine TLabMPI_Trp_ExecK_Forward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Forward_Complex(a, b, trp_plan)
        ! K-Forward (complex): same strided-to-flat scatter as the real version but for complex(dp).
        ! apu_cx_all spans all peers' windows reinterpreted as complex; stride in complex units = apu_stride_k/2.
        ! APU_ASYNC has no dedicated complex path and falls back to ASYNCHRONOUS MPI.
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr
        integer(wi) :: nmax_p, nlines_p, npage, flat_off, disp_ns, mas
#ifdef USE_APU
        complex(dp), pointer :: apu_cx_all(:) => null()   ! complex view of apu_all_k across all peers
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
        if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
            size = trp_plan%size3d
            ! Reinterpret the contiguous real window as complex; size in complex units = apu_stride_k/2
            call c_f_pointer(apu_peer_cptr_k(0), apu_cx_all, [apu_stride_k*ims_npro_k/2])
            call MPI_Win_fence(0, apu_win_k, ims_err)
            ! Push: strided gather from a into all peers' recv buffers in one fused kernel
            !$omp target teams distribute parallel do collapse(3) if(mas * sizeofreal > 100000_wi)
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
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: recv buffer now fully populated
            ! Unpack: flat copy from complex-typed recv alias to b
            !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
            do i = 1, size
                b(i) = apu_cx_recv_fptr_k(i)
            end do
            !$omp end target teams distribute parallel do
            nullify (apu_cx_all)

        else   ! ASYNCHRONOUS, APU_ASYNC (no dedicated complex path), SENDRECV, ALLTOALL
#endif
        ! ==================================================================== !
        ! MPI path — pack into flat staging buffer c_wrk_cx, then ISEND/IRECV. !
        ! APU_ASYNC uses this same path for complex (hybrid not implemented).   !
        ! ==================================================================== !
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS .or. trp_mode_k == TLAB_MPI_TRP_APU_ASYNC .or. &
                trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
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
                ! ISEND/IRECV in batches of trp_sizBlock_k peers per WAITALL
                do j = 1, ims_npro_k, trp_sizBlock_k
                    l = 0
                    do m = j, min(j + trp_sizBlock_k - 1, ims_npro_k)
                        ns = maps_send_k(m) + 1; ips = ns - 1
                        nr = maps_recv_k(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_ISEND(c_wrk_cx((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, ims_comm_z, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(b(trp_plan%disp_r(nr) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, ims_comm_z, request(l), ims_err)
                    end do
                    call MPI_WAITALL(l, request, status, ims_err)
                end do
                nullify (c_wrk_cx)
            else
                call Transpose_Kernel_Complex(a, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                              b, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                              ims_comm_z, trp_sizBlock_k, trp_mode_k)
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
        real(dp), pointer :: apu_pfptr_k(:) => null()   ! scratch pointer for APU_ASYNC per-peer writes
        real(dp) :: dbg_intra, dbg_inter   ! FABRIC_DIRECT debug: recv-buffer split checksums
        integer(MPI_ADDRESS_KIND) :: dbg_addr
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

        ! Common-path debug tap: input checksum for K-backward.
#ifdef USE_APU
        call TLab_Debug_Print_1D('[KBR_pre] sum(b)=', b)
#endif

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
            !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
            do m = 0, ims_npro_k - 1
                do i = 1, nmax_p * nlines_p
                    ! Write our chunk (b[m*chunk]) into peer m's slot (own_rank*chunk) in their buffer
                    apu_all_k(m*apu_stride_k + ims_pro_k*nmax_p*nlines_p + i) = b(m * nmax_p * nlines_p + i)
                end do
            end do
            !$omp end target teams distribute parallel do
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: all peers have written to our buffer
            ! -- Unpack: recv buffer holds sorted chunks; scatter to strided a in one fused kernel.
            !$omp target teams distribute parallel do collapse(3) if(mas * sizeofreal > 100000_wi)
            do m = 0, ims_npro_k - 1
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        a(m*nlines_p + i*npage + j + 1) = &
                            apu_recv_fptr_k(m*nmax_p*nlines_p + i*nlines_p + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do

        else if (trp_mode_k == TLAB_MPI_TRP_APU_ASYNC) then
            ! Hybrid backward: each rank pushes b[m*chunk] to intra-node peer m directly,
            ! or ISENDs b[m*chunk] to inter-node peers. After sync, apu_async_recv_k is
            ! in the same layout as in APU_DIRECT and is scattered to a.
            size = trp_plan%size3d
            l = 0
            ! Step 1: post IRECVs for inter-node peers
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    l = l + 1
                    call MPI_IRECV(apu_async_recv_k(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_z, request(l), ims_err)
                end if
            end do
            ! Step 2: open shared-window epoch
            call MPI_Win_fence(0, apu_async_win_k, ims_err)
            ! Step 3: intra-node — push b[m*chunk] into peer m's recv buffer at our slot
            do m = 0, ims_npro_k - 1
                if (apu_async_is_local_k(m)) then
                    call c_f_pointer(apu_async_peer_k(m), apu_pfptr_k, [apu_async_size_k])
                    flat_off = ims_pro_k * nmax_p * nlines_p
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
                    do i = 1, nmax_p * nlines_p
                        apu_pfptr_k(flat_off + i) = b(m * nmax_p * nlines_p + i)
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end do
            ! Step 4: inter-node — ISEND b[m*chunk] directly (already flat, no packing needed)
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    l = l + 1
                    call MPI_ISEND(b(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_z, request(l), ims_err)
                end if
            end do
            ! Step 5: close epoch and wait for inter-node MPI
            call MPI_Win_fence(0, apu_async_win_k, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! Step 6: unpack recv buffer — recv_k[m*chunk] = data from rank m → scatter to a
            do m = 0, ims_npro_k - 1
                flat_off = m * nmax_p * nlines_p
                disp_nr  = m * nlines_p
                !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        a(disp_nr + i*npage + j + 1) = apu_async_recv_k(flat_off + i*nlines_p + j + 1)
                    end do
                end do
                !$omp end target teams distribute parallel do
            end do
            nullify (apu_pfptr_k)

        else if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
            ! Two-level backward: b is flat K-space (b[m*chunk] destined for peer m).
            ! Intra-node peers: GPU writes b[m*chunk] into peer m's recv buffer at our slot.
            ! Inter-node peers: plain two-sided MPI (ISEND b[m*chunk]; IRECV into our slot m).
            ! After sync, apu_async_recv_k[m*chunk] = data from rank m → scatter to strided a
            ! (GPU for intra-node source slots, CPU for inter-node ones).
            size = trp_plan%size3d
            l = 0
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    l = l + 1
                    call MPI_IRECV(apu_async_recv_k(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_z, request(l), ims_err)
                end if
            end do
            call MPI_Win_fence(0, apu_async_win_k, ims_err)
            do m = 0, ims_npro_k - 1
                if (apu_async_is_local_k(m)) then
                    call c_f_pointer(apu_async_peer_k(m), apu_pfptr_k, [apu_async_size_k])
                    flat_off = ims_pro_k * nmax_p * nlines_p
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
                    do i = 1, nmax_p * nlines_p
                        apu_pfptr_k(flat_off + i) = b(m * nmax_p * nlines_p + i)
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end do
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    l = l + 1
                    call MPI_ISEND(b(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_z, request(l), ims_err)
                end if
            end do
            call MPI_Win_fence(0, apu_async_win_k, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! debug: split checksum of the recv buffer before scattering to strided a.
            dbg_intra = 0.0_dp; dbg_inter = 0.0_dp
            do m = 0, ims_npro_k - 1
                if (apu_async_is_local_k(m)) then
                    dbg_intra = dbg_intra + sum(apu_async_recv_k(m*nmax_p*nlines_p + 1 : (m + 1)*nmax_p*nlines_p))
                else
                    dbg_inter = dbg_inter + sum(apu_async_recv_k(m*nmax_p*nlines_p + 1 : (m + 1)*nmax_p*nlines_p))
                end if
            end do
            write(500 + ims_pro, *) '[KBR_FBD] PE', ims_pro, ' recv intra=', dbg_intra, ' inter=', dbg_inter, &
                ' total=', sum(apu_async_recv_k(1:size))
            flush(500 + ims_pro)
            do m = 0, ims_npro_k - 1
                flat_off = m * nmax_p * nlines_p
                disp_nr  = m * nlines_p
                if (apu_async_is_local_k(m)) then
                    !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            a(disp_nr + i*npage + j + 1) = apu_async_recv_k(flat_off + i*nlines_p + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                else
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            a(disp_nr + i*npage + j + 1) = apu_async_recv_k(flat_off + i*nlines_p + j + 1)
                        end do
                    end do
                end if
            end do
            nullify (apu_pfptr_k)

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
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                        !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
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
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                        !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
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

        ! Common-path debug tap: K-backward output checksum.
#ifdef USE_APU
        call TLab_Debug_Print_1D('[KBR_post] sum(a)=', a)
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
#ifdef USE_APU
        complex(dp), pointer :: apu_cx_all(:) => null()
#endif

        nmax_p   = trp_plan%nmax
        nlines_p = trp_plan%nlines
        npage    = nlines_p * ims_npro_k
        mas      = nmax_p * nlines_p

        ! ==================================================================== !
        ! APU_DIRECT path — fused GPU kernels; inverse of K-Forward_Complex.   !
        ! ==================================================================== !
#ifdef USE_APU
        if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
            size = trp_plan%size3d
            call c_f_pointer(apu_peer_cptr_k(0), apu_cx_all, [apu_stride_k*ims_npro_k/2])
            call MPI_Win_fence(0, apu_win_k, ims_err)
            ! Push: b is flat K-space; write our chunk (b[m*chunk]) to every peer
            !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
            do m = 0, ims_npro_k - 1
                do i = 1, nmax_p * nlines_p
                    apu_cx_all(m*(apu_stride_k/2) + ims_pro_k*nmax_p*nlines_p + i) = &
                        b(m * nmax_p * nlines_p + i)
                end do
            end do
            !$omp end target teams distribute parallel do
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: recv buffer fully populated
            ! Unpack: complex recv buffer → strided Z-space a
            !$omp target teams distribute parallel do collapse(3) if(mas * sizeofreal > 100000_wi)
            do m = 0, ims_npro_k - 1
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        a(m*nlines_p + i*npage + j + 1) = &
                            apu_cx_recv_fptr_k(m*nmax_p*nlines_p + i*nlines_p + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            nullify (apu_cx_all)

        else   ! ASYNCHRONOUS, APU_ASYNC (fallback for complex), SENDRECV, ALLTOALL
#endif
        ! ==================================================================== !
        ! MPI path — ISEND/IRECV into flat c_wrk_cx, then scatter to a.       !
        ! ==================================================================== !
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS .or. trp_mode_k == TLAB_MPI_TRP_APU_ASYNC .or. &
                trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
                ! ISEND/IRECV in batches; recv into flat c_wrk_cx
                do j = 1, ims_npro_k, trp_sizBlock_k
                    l = 0
                    do m = j, min(j + trp_sizBlock_k - 1, ims_npro_k)
                        ns = maps_recv_k(m) + 1; ips = ns - 1   ! backward: send/recv maps swapped
                        nr = maps_send_k(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_ISEND(b(trp_plan%disp_r(ns) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, ims_comm_z, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(c_wrk_cx((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, ims_comm_z, request(l), ims_err)
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
                                              ims_comm_z, trp_sizBlock_k, trp_mode_k)
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
        real(dp), pointer :: apu_pfptr_i(:) => null()   ! scratch pointer for APU_ASYNC per-peer writes
        real(dp), pointer, contiguous :: c_recv_dp(:) => null()   ! FABRIC_DIRECT: non-shm IRECV staging
        real(dp) :: dbg_intra, dbg_inter   ! FABRIC_DIRECT debug: recv-buffer split checksums
        integer(MPI_ADDRESS_KIND) :: dbg_addr
        integer(wi) :: fp_nerr               ! FINGERPRINT TEST: mismatch counter
        real(dp) :: fp_expected, fp_actual   ! FINGERPRINT TEST: expected and actual values
#endif

        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i   ! total X-elements per line (stride in b)
        mas       = nmax_p * nlines_p

        ! Common-path debug tap: input checksum for I-forward.
#ifdef USE_APU
        call TLab_Debug_Print_1D('[IFR_pre] sum(a)=', a)
#endif

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
            !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
            do m = 0, ims_npro_i - 1
                do i = 1, nmax_p * nlines_p
                    ! a[m*chunk] is the flat chunk destined for peer m; write to their recv slot
                    apu_all_i(m*apu_stride_i + ims_pro_i*nmax_p*nlines_p + i) = a(m * nmax_p * nlines_p + i)
                end do
            end do
            !$omp end target teams distribute parallel do
            call MPI_Win_fence(0, apu_win_i, ims_err)   ! barrier: recv buffer fully populated
            ! -- Unpack: recv buffer holds sorted flat chunks; scatter to strided b in one fused kernel.
            !$omp target teams distribute parallel do collapse(3) if(mas * sizeofreal > 100000_wi)
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

        else if (trp_mode_i == TLAB_MPI_TRP_APU_ASYNC) then
            ! Hybrid I-forward: intra-node peers via direct writes; inter-node via MPI.
            ! a is flat (a[m*chunk] for peer m); both paths fill apu_async_recv_i, then scatter to b.
            size = trp_plan%size3d
            l = 0
            ! Step 1: post IRECVs for inter-node peers
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    l = l + 1
                    call MPI_IRECV(apu_async_recv_i(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                end if
            end do
            ! Step 2: open shared-window epoch
            call MPI_Win_fence(0, apu_async_win_i, ims_err)
            ! Step 3: intra-node — push a[m*chunk] directly to peer m's recv buffer at our slot
            do m = 0, ims_npro_i - 1
                if (apu_async_is_local_i(m)) then
                    call c_f_pointer(apu_async_peer_i(m), apu_pfptr_i, [apu_async_size_i])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
                    do i = 1, nmax_p * nlines_p
                        apu_pfptr_i(flat_off + i) = a(m * nmax_p * nlines_p + i)
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end do
            ! Step 4: inter-node — ISEND a[m*chunk] (already flat, no packing needed)
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    l = l + 1
                    call MPI_ISEND(a(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                end if
            end do
            ! Step 5: close epoch and wait for inter-node MPI
            call MPI_Win_fence(0, apu_async_win_i, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! Step 6: unpack recv buffer → strided b
            do m = 0, ims_npro_i - 1
                flat_off = m * nmax_p * nlines_p
                disp_nr  = m * nmax_p
                !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                do i = 0, nlines_p - 1
                    do j = 0, nmax_p - 1
                        b(disp_nr + i*nmax_full + j + 1) = apu_async_recv_i(flat_off + i*nmax_p + j + 1)
                    end do
                end do
                !$omp end target teams distribute parallel do
            end do
            nullify (apu_pfptr_i)

        else if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
            ! Two-level I-forward: fingerprint test.
            ! Sends fingerprint (our_rank * 1e9 + pos) instead of a(); verifies what arrives then STOPs.
            ! IRECV target is c_recv_dp (a non-shm Fortran-allocated buffer, upper half of wrk_mpi_dp),
            ! not apu_async_recv_i directly: previous attempt with IRECV into the MPI_Win_allocate_shared
            ! buffer delivered nothing inter-node on Cray MPICH (slot 6..11 came out as PE 0's own
            ! fingerprint pattern, not the sender's). Staging then copying matches what the working
            ! ASYNCHRONOUS path does and isolates the inter-node delivery from the shm window.
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)),        c_wrk_dp,  shape=[size])   ! ISEND staging
            call c_f_pointer(c_loc(wrk_mpi_dp(size + 1)), c_recv_dp, shape=[size])   ! IRECV staging
            l = 0
            ! Step 1: post IRECVs for inter-node peers INTO STAGING (not into apu_async_recv_i)
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    l = l + 1
                    call MPI_IRECV(c_recv_dp(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                end if
            end do
            call MPI_Win_fence(0, apu_async_win_i, ims_err)
            ! FINGERPRINT TEST step 3: intra-node — CPU writes fingerprint to peer's recv buffer at our slot
            do m = 0, ims_npro_i - 1
                if (apu_async_is_local_i(m)) then
                    call c_f_pointer(apu_async_peer_i(m), apu_pfptr_i, [apu_async_size_i])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    do i = 1, nmax_p * nlines_p
                        apu_pfptr_i(flat_off + i) = real(ims_pro_i, dp) * 1.0e9_dp + real(i, dp)
                    end do
                end if
            end do
            ! Sender-side probe: print ims_pro_i and the runtime value of the fingerprint expression.
            ! Goes to fort.500+rank, so PE 6's log (fort.506) will show what rank-prefix it used.
            write(500 + ims_pro, *) '[IFR_FP_SEND_RANK] PE', ims_pro, &
                ' ims_pro_i=', ims_pro_i, &
                ' real(ims_pro_i,dp)*1e9=', real(ims_pro_i, dp) * 1.0e9_dp
            flush(500 + ims_pro)
            ! FINGERPRINT TEST step 4: inter-node — fill staging with fingerprint, then ISEND
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    flat_off = m * nmax_p * nlines_p
                    do i = 1, nmax_p * nlines_p
                        c_wrk_dp(flat_off + i) = real(ims_pro_i, dp) * 1.0e9_dp + real(i, dp)
                    end do
                    ! Probe the values we just wrote, BEFORE the ISEND reads them.
                    ! If c_wrk_dp(flat_off+1..5) holds ims_pro_i*1e9+1..5, the fill is fine and
                    ! the bug is in transit. If it holds 1..5, the fill itself is wrong on this rank.
                    write(500 + ims_pro, *) '[IFR_FP_SEND] PE', ims_pro, &
                        ' to m=', m, ' c_wrk_dp(flat_off+1..5)=', &
                        c_wrk_dp(flat_off + 1), c_wrk_dp(flat_off + 2), &
                        c_wrk_dp(flat_off + 3), c_wrk_dp(flat_off + 4), &
                        c_wrk_dp(flat_off + 5)
                    flush(500 + ims_pro)
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp(flat_off + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                end if
            end do
            call MPI_Win_fence(0, apu_async_win_i, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! Address-aliasing probe: c_wrk_dp and c_recv_dp must point to different memory.
            ! If c_loc(wrk_mpi_dp(size+1)) is broken on this compiler, c_recv_addr == c_wrk_addr
            ! and the staging buffer is just the send buffer (which explains got=pos for inter-node).
            dbg_addr = transfer(c_loc(c_wrk_dp(1)),  dbg_addr)
            write(500 + ims_pro, *) '[IFR_FP_ADDR] PE', ims_pro, ' c_loc(c_wrk_dp(1))= ', dbg_addr
            dbg_addr = transfer(c_loc(c_recv_dp(1)), dbg_addr)
            write(500 + ims_pro, *) '[IFR_FP_ADDR] PE', ims_pro, ' c_loc(c_recv_dp(1))=', dbg_addr
            ! Probe: what IRECV actually delivered into c_recv_dp at slot 6 (the first inter-node peer).
            ! If staging works, c_recv_dp(6*chunk + 1..20) should hold 6e9+1..6e9+20 (from rank 6).
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    flat_off = m * nmax_p * nlines_p
                    write(500 + ims_pro, *) '[IFR_FP_PRECOPY] PE', ims_pro, &
                        ' m=', m, ' c_recv_dp(flat_off+1..5)=', &
                        c_recv_dp(flat_off + 1), c_recv_dp(flat_off + 2), &
                        c_recv_dp(flat_off + 3), c_recv_dp(flat_off + 4), &
                        c_recv_dp(flat_off + 5)
                    exit
                end if
            end do
            flush(500 + ims_pro)
            ! Copy inter-node slots from non-shm staging into apu_async_recv_i for verification.
            ! Intra-node slots are already in apu_async_recv_i via the shared-window writes in step 3.
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    flat_off = m * nmax_p * nlines_p
                    do i = 1, nmax_p * nlines_p
                        apu_async_recv_i(flat_off + i) = c_recv_dp(flat_off + i)
                    end do
                end if
            end do
            ! FINGERPRINT VERIFICATION: slot m in apu_async_recv_i must hold fingerprints from rank m.
            ! expected(m, i) = m * 1e9 + i  for i = 1..nmax_p*nlines_p
            fp_nerr = 0
            do m = 0, ims_npro_i - 1
                flat_off = m * nmax_p * nlines_p
                do i = 1, nmax_p * nlines_p
                    fp_expected = real(m, dp) * 1.0e9_dp + real(i, dp)
                    fp_actual   = apu_async_recv_i(flat_off + i)
                    if (abs(fp_actual - fp_expected) > 0.5_dp) then
                        fp_nerr = fp_nerr + 1
                        if (fp_nerr <= 20) then
                            write(500 + ims_pro, *) '[IFR_FP_ERR] PE', ims_pro, &
                                ' from_rank=', m, ' pos=', i, &
                                ' expected=', fp_expected, ' got=', fp_actual, &
                                ' local=', apu_async_is_local_i(m)
                        end if
                    end if
                end do
            end do
            write(500 + ims_pro, *) '[IFR_FP] PE', ims_pro, &
                ' chunk=', nmax_p*nlines_p, ' total_errors=', fp_nerr
            flush(500 + ims_pro)
            nullify (c_wrk_dp, c_recv_dp, apu_pfptr_i)
            STOP 'IFR fingerprint test done -- check fort.500+ logs'

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
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                        !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
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
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                    ! a is already flat; ISEND directly, recv into c_wrk_dp, scatter to strided b.
                    size = trp_plan%size3d
                    call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                    ! Step 1: IRECVs
                    l = 0
                    do m = 1, ims_npro_i
                        nr = maps_recv_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_IRECV(c_wrk_dp((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    ! Step 2: ISENDs from flat a
                    do m = 1, ims_npro_i
                        ns = maps_send_i(m) + 1; ips = ns - 1
                        l = l + 1
                        call MPI_ISEND(a(trp_plan%disp_s(ns) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    ! Step 3: WAITALL
                    call MPI_WAITALL(l, request, status, ims_err)
                    ! Step 4: GPU scatter flat c_wrk_dp → strided b
                    do m = 1, ims_npro_i
                        nr = maps_recv_i(m) + 1
                        flat_off = (nr - 1)*nmax_p*nlines_p
                        disp_nr  = trp_plan%disp_r(nr)
#ifdef USE_APU
                        !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
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

        ! Common-path debug tap: I-forward output checksum.
#ifdef USE_APU
        call TLab_Debug_Print_1D('[IFR_post] sum(b)=', b)
#endif
        return
    end subroutine TLabMPI_Trp_ExecI_Forward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Forward_Complex(a, b, trp_plan)
        ! I-Forward (complex): flat I-space chunks in a → strided X-space layout in b.
        ! a(m*chunk + i) where chunk = nmax_p*nlines_p, m = peer index, i = local element.
        ! b(m*nmax_p + i*nmax_full + j + 1): strided with full-width stride nmax_full = nmax_p*npro.
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan
        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_nr
#ifdef USE_APU
        ! apu_cx_all: complex view spanning all peers' shared windows (stride = apu_stride_i/2 complex units).
        complex(dp), pointer :: apu_cx_all(:) => null()
        integer(wi) :: mas
#endif

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i

        ! ==================================================================== !
        ! APU paths — GPU direct writes between shared-memory windows.         !
        ! ==================================================================== !
#ifdef USE_APU
        if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT) then
            ! a is flat I-space chunks; fused GPU kernel pushes to all peers' recv buffers,
            ! then unpacks our own recv buffer (strided) into b.
            size = trp_plan%size3d
            mas  = nmax_p * nlines_p
            ! Build complex-typed view spanning all peers' contiguous window segments.
            call c_f_pointer(apu_peer_cptr_i(0), apu_cx_all, [apu_stride_i*ims_npro_i/2])
            ! Fence 1: open epoch — all ranks ready to receive direct writes.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! Push: write each peer's flat chunk into peer m's buffer at slot own_rank*chunk.
            !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
            do m = 0, ims_npro_i - 1
                do i = 1, nmax_p * nlines_p
                    apu_cx_all(m*(apu_stride_i/2) + ims_pro_i*nmax_p*nlines_p + i) = &
                        a(m * nmax_p * nlines_p + i)
                end do
            end do
            !$omp end target teams distribute parallel do
            ! Fence 2: close epoch — all writes committed; recv buffers fully populated.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! Unpack: scatter recv buffer (flat m*chunk+i layout) → b (strided m*nmax_p + i*nmax_full + j).
            !$omp target teams distribute parallel do collapse(3) if(mas * sizeofreal > 100000_wi)
            do m = 0, ims_npro_i - 1
                do i = 0, nlines_p - 1
                    do j = 0, nmax_p - 1
                        b(m*nmax_p + i*nmax_full + j + 1) = &
                            apu_cx_recv_fptr_i(m*nmax_p*nlines_p + i*nmax_p + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            nullify (apu_cx_all)
        else   ! APU_ASYNC falls back to ASYNC for complex; SENDRECV/ALLTOALL go to kernel.
#endif
        ! ==================================================================== !
        ! CPU paths — MPI ISEND/IRECV, SENDRECV, ALLTOALL                     !
        ! ==================================================================== !
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS .or. trp_mode_i == TLAB_MPI_TRP_APU_ASYNC .or. &
                trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
                ! Flat send of each peer's chunk from a; recv into flat staging c_wrk_cx;
                ! then scatter c_wrk_cx → strided b.
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
                ! Pipeline: batched ISEND+IRECV → WAITALL per trp_sizBlock_i batch.
                do j = 1, ims_npro_i, trp_sizBlock_i
                    l = 0
                    do m = j, min(j + trp_sizBlock_i - 1, ims_npro_i)
                        ns = maps_send_i(m) + 1; ips = ns - 1
                        nr = maps_recv_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_ISEND(a(trp_plan%disp_s(ns) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(c_wrk_cx((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
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
                                              ims_comm_x, trp_sizBlock_i, trp_mode_i)
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
        ! apu_pfptr_i: per-peer pointer into apu_async_win_i for APU_ASYNC intra-node writes.
        real(dp), pointer :: apu_pfptr_i(:) => null()
        real(dp) :: dbg_intra, dbg_inter   ! FABRIC_DIRECT debug: recv-buffer split checksums
        integer(MPI_ADDRESS_KIND) :: dbg_addr
#endif

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i
        mas       = nmax_p * nlines_p

        ! Common-path debug tap: input checksum for I-backward.
#ifdef USE_APU
        call TLab_Debug_Print_1D('[IBR_pre] sum(b)=', b)
#endif

        ! ==================================================================== !
        ! APU paths — GPU direct writes between shared-memory windows.         !
        ! ==================================================================== !
#ifdef USE_APU
        if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT) then
            ! b is strided X-space (disp_r stride). Each rank packs its slice from b into
            ! every peer m's recv buffer at slot own_rank*chunk. Single fused collapse(3)
            ! kernel covers all m in one HIP launch, eliminating per-peer launch overhead.
            size = trp_plan%size3d
            ! Fence 1: open epoch — all ranks ready to receive direct writes.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! Push: pack strided b[m] → peer m's recv buffer at slot own_rank*chunk (flat).
            !$omp target teams distribute parallel do collapse(3) if(mas * sizeofreal > 100000_wi)
            do m = 0, ims_npro_i - 1
                do i = 0, nlines_p - 1
                    do j = 0, nmax_p - 1
                        apu_all_i(m*apu_stride_i + ims_pro_i*nmax_p*nlines_p + i*nmax_p + j + 1) = &
                            b(m*nmax_p + i*nmax_full + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            ! Fence 2: close epoch — all writes committed; recv buffers fully populated.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! Flat copy: recv buffer layout is flat and matches a 1:1.
            !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
            do i = 1, size
                a(i) = apu_recv_fptr_i(i)
            end do
            !$omp end target teams distribute parallel do
        else if (trp_mode_i == TLAB_MPI_TRP_APU_ASYNC) then
            ! Hybrid I-backward: b is strided X-space (disp_r stride = m*nmax_p).
            ! Intra-node peers: direct write strided b[m] into peer m's shared recv buffer at our slot.
            ! Inter-node peers: pack strided b[m] → flat c_wrk_dp → ISEND.
            ! After sync: recv buffer is flat → a.
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
            ! Step 1: post IRECVs for inter-node peers into their slots in the async recv window.
            l = 0
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    l = l + 1
                    call MPI_IRECV(apu_async_recv_i(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                end if
            end do
            ! Step 2: open shared-window epoch for intra-node writes.
            call MPI_Win_fence(0, apu_async_win_i, ims_err)
            ! Step 3: intra-node — pack strided b[m] into peer m's shared recv buffer at our slot.
            do m = 0, ims_npro_i - 1
                if (apu_async_is_local_i(m)) then
                    call c_f_pointer(apu_async_peer_i(m), apu_pfptr_i, [apu_async_size_i])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    disp_ns  = m * nmax_p
                    !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            apu_pfptr_i(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end do
            ! Step 4: inter-node — pack strided b[m] into flat c_wrk_dp then ISEND.
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    flat_off = m * nmax_p * nlines_p
                    disp_ns  = m * nmax_p
                    !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            c_wrk_dp(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp(flat_off + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                end if
            end do
            ! Step 5: close window epoch; wait for inter-node sends/recvs to complete.
            call MPI_Win_fence(0, apu_async_win_i, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! Step 6: flat copy from unified async recv buffer to a (1:1).
            !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
            do i = 1, size
                a(i) = apu_async_recv_i(i)
            end do
            !$omp end target teams distribute parallel do
            nullify (c_wrk_dp, apu_pfptr_i)
        else if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
            ! Two-level I-backward: b is strided X-space (disp_r stride = m*nmax_p).
            ! Intra-node peers: GPU writes strided b[m] into peer m's shared recv buffer at our slot.
            ! Inter-node peers: CPU-pack strided b[m] → flat c_wrk_dp slot → ISEND; IRECV into our slot m.
            ! After sync, apu_async_recv_i is flat → copy to a (GPU for intra-node slots, CPU for inter).
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
            l = 0
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    l = l + 1
                    call MPI_IRECV(apu_async_recv_i(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                end if
            end do
            call MPI_Win_fence(0, apu_async_win_i, ims_err)
            do m = 0, ims_npro_i - 1
                if (apu_async_is_local_i(m)) then
                    call c_f_pointer(apu_async_peer_i(m), apu_pfptr_i, [apu_async_size_i])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    disp_ns  = m * nmax_p
                    !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            apu_pfptr_i(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end do
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    flat_off = m * nmax_p * nlines_p
                    disp_ns  = m * nmax_p
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            c_wrk_dp(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp(flat_off + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                end if
            end do
            call MPI_Win_fence(0, apu_async_win_i, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! debug: split checksum of the assembled recv buffer; total must equal [IBR_post]'s sum(a).
            dbg_intra = 0.0_dp; dbg_inter = 0.0_dp
            do m = 0, ims_npro_i - 1
                if (apu_async_is_local_i(m)) then
                    dbg_intra = dbg_intra + sum(apu_async_recv_i(m*nmax_p*nlines_p + 1 : (m + 1)*nmax_p*nlines_p))
                else
                    dbg_inter = dbg_inter + sum(apu_async_recv_i(m*nmax_p*nlines_p + 1 : (m + 1)*nmax_p*nlines_p))
                end if
            end do
            write(500 + ims_pro, *) '[IBR_FBD] PE', ims_pro, ' recv intra=', dbg_intra, ' inter=', dbg_inter, &
                ' total=', sum(apu_async_recv_i(1:size))
            flush(500 + ims_pro)
            do m = 0, ims_npro_i - 1
                flat_off = m * nmax_p * nlines_p
                if (apu_async_is_local_i(m)) then
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
                    do i = 1, nmax_p*nlines_p
                        a(flat_off + i) = apu_async_recv_i(flat_off + i)
                    end do
                    !$omp end target teams distribute parallel do
                else
                    do i = 1, nmax_p*nlines_p
                        a(flat_off + i) = apu_async_recv_i(flat_off + i)
                    end do
                end if
            end do
            nullify (c_wrk_dp, apu_pfptr_i)
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
                        !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
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
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
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
                    ! Pipeline: post IRECVs → GPU strided→flat pack → ISENDs → WAITALL.
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
                    ! Step 2: GPU pack b→c_wrk_dp (strided→flat) while network prepares recv buffers.
                    do m = 1, ims_npro_i
                        ns = maps_recv_i(m) + 1
                        flat_off = (ns - 1)*nmax_p*nlines_p
                        disp_ns  = trp_plan%disp_r(ns)
#ifdef USE_APU
                        !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
#endif
                        do i = 0, nlines_p - 1
                            do j = 0, nmax_p - 1
                                c_wrk_dp(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                            end do
                        end do
#ifdef USE_APU
                        !$omp end target teams distribute parallel do
#endif
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
        ! Common-path debug tap: I-backward output checksum.
#ifdef USE_APU
        call TLab_Debug_Print_1D('[IBR_post] sum(a)=', a)
#endif
        return
    end subroutine TLabMPI_Trp_ExecI_Backward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Backward_Complex(b, a, trp_plan)
        ! I-Backward (complex): inverse of I-Forward. Reverses strided X-space (b) → flat I-space (a).
        ! b(m*nmax_p + i*nmax_full + j + 1): strided X-space, nmax_full = nmax_p*npro.
        ! a: flat output, chunks of nmax_p*nlines_p per peer (mirrors forward input layout).
        complex(wp), intent(in) :: b(:)
        complex(wp), intent(out) :: a(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan
        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_ns
#ifdef USE_APU
        ! apu_cx_all: complex view spanning all peers' shared windows (stride = apu_stride_i/2 complex units).
        complex(dp), pointer :: apu_cx_all(:) => null()
        integer(wi) :: mas
#endif

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i

        ! ==================================================================== !
        ! APU paths — GPU direct writes between shared-memory windows.         !
        ! ==================================================================== !
#ifdef USE_APU
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
            !$omp target teams distribute parallel do collapse(3) if(mas * sizeofreal > 100000_wi)
            do m = 0, ims_npro_i - 1
                do i = 0, nlines_p - 1
                    do j = 0, nmax_p - 1
                        apu_cx_all(m*(apu_stride_i/2) + ims_pro_i*nmax_p*nlines_p + i*nmax_p + j + 1) = &
                            b(m*nmax_p + i*nmax_full + j + 1)
                    end do
                end do
            end do
            !$omp end target teams distribute parallel do
            ! Fence 2: close epoch — all writes committed; recv buffers fully populated.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! Flat copy: recv buffer is flat and matches a layout 1:1.
            !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
            do i = 1, size
                a(i) = apu_cx_recv_fptr_i(i)
            end do
            !$omp end target teams distribute parallel do
            nullify (apu_cx_all)
        else   ! APU_ASYNC falls back to ASYNC for complex; SENDRECV/ALLTOALL go to kernel.
#endif
        ! ==================================================================== !
        ! CPU paths — MPI ISEND/IRECV, SENDRECV, ALLTOALL                     !
        ! ==================================================================== !
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS .or. trp_mode_i == TLAB_MPI_TRP_APU_ASYNC .or. &
                trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
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
                do j = 1, ims_npro_i, trp_sizBlock_i
                    l = 0
                    do m = j, min(j + trp_sizBlock_i - 1, ims_npro_i)
                        ns = maps_recv_i(m) + 1; ips = ns - 1
                        nr = maps_send_i(m) + 1; ipr = nr - 1
                        l = l + 1
                        call MPI_ISEND(c_wrk_cx((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
                        l = l + 1
                        call MPI_IRECV(a(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
                    end do
                    call MPI_WAITALL(l, request, status, ims_err)
                end do
                nullify (c_wrk_cx)
            else
                call Transpose_Kernel_Complex(b, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                              a, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                              ims_comm_x, trp_sizBlock_i, trp_mode_i)
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
