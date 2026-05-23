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
    ! APU_ASYNC / FABRIC_DIRECT state: node-local shmem comm + per-rank buffer pointers.
    integer(wi) :: apu_async_size_k = 0_wi, apu_async_size_i = 0_wi
    type(MPI_Comm) :: apu_async_node_comm_k, apu_async_node_comm_i   ! node-local shmem comm, kept alive for runtime barriers
    ! Recv buffers are backed by an MPI_Win_allocate_shared segment so node-local peers
    ! all see the SAME virtual addresses (contiguous shared mapping). Plain Fortran
    ! allocate() was previously tried under the assumption that APU unified memory plus
    ! XPMEM made private heap addresses cross-process-valid; the 2026-05-17 crashlog
    ! disproves that — addresses transmitted via MPI_Allgather pointed into private
    ! memory of each process, producing overlap/garbage in OpenACC's present-table
    ! ("Host region overlaps present region but is not contained for apu_pfptr_i(:)").
    ! Allocating a shared window on the node-local sub-comm of MPI_Comm_split_type DOES
    ! taint the parent ims_comm_z/x for two-sided traffic on Cray MPICH; the fix is to
    ! MPI_Comm_dup the parent BEFORE the split (the dup is outside the lineage that
    ! gets tainted) and route fabricdirect's inter-node MPI_ISEND/IRECV through the dup.
    type(MPI_Win) :: apu_async_win_k, apu_async_win_i              ! shared window backing recv buffers
    type(MPI_Comm) :: apu_async_mpi_comm_k, apu_async_mpi_comm_i   ! dup of ims_comm_z/x, used for two-sided traffic
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
    ! APU_ASYNC/FABRIC_DIRECT: contiguous all-local-shmem-peers view for fused single-kernel writes.
    ! Spans apu_async_shmem_size_i * apu_async_size_i elements starting at shmem rank 0's segment.
    ! apu_async_all_i(shmem_m*apu_async_size_i + 1 : (shmem_m+1)*apu_async_size_i) = shmem peer m's recv buffer.
    real(dp), pointer :: apu_async_all_i(:) => null()
    integer(wi) :: apu_async_shmem_size_i = 0_wi       ! number of local shmem peers in the I-window
    integer(wi) :: apu_async_shmem_base_pro_i = 0_wi   ! pro_i of shmem rank 0 in the local I-shmem group
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

#ifdef USE_APU
    interface
        subroutine hip_write_with_fence(src, dst, n) bind(C, name='hip_write_with_fence')
            use iso_c_binding
            real(c_double), intent(in)  :: src(*)
            real(c_double), intent(out) :: dst(*)
            integer(c_int), value       :: n
        end subroutine

        function hipHostRegister(ptr, sz, flags) bind(C, name='hipHostRegister') result(ierr)
            use iso_c_binding
            integer(c_int) :: ierr
            type(c_ptr), value       :: ptr
            integer(c_size_t), value :: sz
            integer(c_int), value    :: flags
        end function
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
        type(MPI_Comm)  :: apu_async_shmem_comm
        type(MPI_Comm)  :: ims_comm_x_dup2   ! throwaway dup of ims_comm_x — used only as split_type parent for I-window
        integer, allocatable :: apu_async_shmem_to_dir(:)  ! shmem-rank → dir-rank (Allgather result)
        integer :: apu_async_shmem_size
        integer(MPI_ADDRESS_KIND) :: win_segsize ! byte size of one peer's recv segment
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
        ! Take BOTH MPI_Comm_dup's FIRST, before ANY MPI_Win_allocate_shared.
        ! On Cray MPICH, allocating an MPI_Win_allocate_shared window on a MPI_Comm_split_type
        ! sub-comm taints the parent directional comm for two-sided traffic. A dup taken BEFORE
        ! any shared-window allocation lies outside the tainted lineage and stays clean.
        ! Both dups are taken together here so neither inherits the other direction's taint.
        if ((trp_mode_k == TLAB_MPI_TRP_APU_ASYNC .or. trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) .and. ims_npro_k > 1) then
            call MPI_Comm_dup(ims_comm_z, apu_async_mpi_comm_k, ims_err)
        end if
        if ((trp_mode_i == TLAB_MPI_TRP_APU_ASYNC .or. trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) .and. ims_npro_i > 1) then
            call MPI_Comm_dup(ims_comm_x, apu_async_mpi_comm_i, ims_err)
            ! Second dup: split_type parent for the I shared window. Only needed for APU_ASYNC;
            ! FABRIC_DIRECT uses no I shared window — allocating MPI_Win_allocate_shared on any
            ! I-direction sub-comm corrupts CXI cross-XCD intra-node routing on Hunter MI300A,
            ! causing WAITALL hangs regardless of comm ordering or timing.
            if (trp_mode_i == TLAB_MPI_TRP_APU_ASYNC) then
                call MPI_Comm_dup(ims_comm_x, ims_comm_x_dup2, ims_err)
            end if
        end if

        ! APU_ASYNC only: I-direction shared window — allocated BEFORE K-shmem.
        ! On Hunter MI300A, MPI_Win_allocate_shared on ANY communicator spanning more than
        ! one XCD creates per-XCD independent sub-windows. shared_query for cross-XCD ranks
        ! returns garbage VAs (confirmed: ranks 3-5 give random values when queried from
        ! ranks 0-2). Only intra-XCD peers (from MPI_Comm_split_type) give valid VAs.
        ! FABRIC_DIRECT does NOT allocate an I shared window: even a 3-member XCD sub-comm
        ! window corrupts CXI cross-XCD intra-node routing, causing WAITALL hangs (confirmed
        ! Hunter run 2026-05-23). FABRIC_DIRECT I-direction uses plain ISEND/IRECV on the
        ! untainted ims_comm_x for all 6 I-peers instead.
        if (trp_mode_i == TLAB_MPI_TRP_APU_ASYNC .and. ims_npro_i > 1) then
            apu_async_size_i = int(imax, wi)*int(jmax, wi)*int(kmax, wi)
            allocate (apu_async_is_local_i(0:ims_npro_i - 1))
            allocate (apu_async_peer_i(0:ims_npro_i - 1))
            apu_async_peer_i = c_null_ptr
            apu_async_is_local_i = .false.
            ! Identify intra-XCD peers via split_type on ims_comm_x_dup2 (throwaway dup).
            call MPI_Comm_split_type(ims_comm_x_dup2, MPI_COMM_TYPE_SHARED, ims_pro_i, MPI_INFO_NULL, &
                                     apu_async_shmem_comm, ims_err)
            call MPI_Comm_size(apu_async_shmem_comm, apu_async_shmem_size, ims_err)
            allocate (apu_async_shmem_to_dir(0:apu_async_shmem_size - 1))
            call MPI_Allgather(ims_pro_i, 1, MPI_INTEGER, apu_async_shmem_to_dir, 1, MPI_INTEGER, &
                               apu_async_shmem_comm, ims_err)
            ! Dup shmem comm BEFORE window allocation (stays clean for MPI_Barrier in APU_ASYNC).
            call MPI_Comm_dup(apu_async_shmem_comm, apu_async_node_comm_i, ims_err)
            ! Allocate window on the XCD-local sub-comm: valid VAs for all 3 intra-XCD peers.
            call MPI_Win_allocate_shared(int(apu_async_size_i, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
                                         int(c_sizeof(1.0_dp)), MPI_INFO_NULL, apu_async_shmem_comm, win_baseptr, apu_async_win_i, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU_ASYNC/FABRIC_DIRECT I recv buffer.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
            ! Query VAs for intra-XCD peers only; mark them local. Cross-XCD peers stay local=F.
            do ip = 0, apu_async_shmem_size - 1
                apu_async_is_local_i(apu_async_shmem_to_dir(ip)) = .true.
                call MPI_Win_shared_query(apu_async_win_i, ip, win_query_size, win_disp_unit, &
                                          apu_async_peer_i(apu_async_shmem_to_dir(ip)), ims_err)
            end do
            ! Bug A fix: bind recv_i to own-segment query result.
            call c_f_pointer(apu_async_peer_i(ims_pro_i), apu_async_recv_i, [apu_async_size_i])
            ! Fused kernel setup: intra-XCD peers only (apu_async_shmem_size members).
            apu_async_shmem_size_i = int(apu_async_shmem_size, wi)
            apu_async_shmem_base_pro_i = int(apu_async_shmem_to_dir(0), wi)
            call c_f_pointer(apu_async_peer_i(apu_async_shmem_to_dir(0)), apu_async_all_i, &
                             [apu_async_shmem_size_i * apu_async_size_i])
            deallocate (apu_async_shmem_to_dir)
            ! Debug: log locality + VAs + size sanity
            block
                integer(MPI_ADDRESS_KIND) :: dbg_va
                write(500+ims_pro,'(a,i4,a,i12,a,i4,a,i4)') &
                    '[INIT_FBD_I] PE', ims_pro, ' size=', apu_async_size_i, &
                    ' npro_i=', ims_npro_i, ' pro_i=', ims_pro_i
                dbg_va = transfer(c_loc(apu_async_recv_i(1)), dbg_va)
                write(500+ims_pro,'(a,i4,a,i22)') '[INIT_FBD_I] PE', ims_pro, ' own_VA=', dbg_va
                do ip = 0, ims_npro_i - 1
                    dbg_va = transfer(apu_async_peer_i(ip), dbg_va)
                    write(500+ims_pro,'(a,i4,a,i4,a,l1,a,i22)') &
                        '[INIT_FBD_I] PE', ims_pro, ' peer', ip, &
                        ' local=', apu_async_is_local_i(ip), ' VA=', dbg_va
                end do
                write(500+ims_pro,'(a,i4,a,i12,a,i12)') &
                    '[INIT_FBD_I] PE', ims_pro, ' max_flat_off=', &
                    (ims_npro_i-1)*apu_async_size_i/ims_npro_i, &
                    ' size=', apu_async_size_i
                flush(500+ims_pro)
            end block
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: I APU_ASYNC recv buffer setup complete.')
        end if
        if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT .and. ims_npro_i > 1) then
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: I FABRIC_DIRECT will use plain MPI (no shared window).')
        end if

        if ((trp_mode_k == TLAB_MPI_TRP_APU_ASYNC .or. trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) .and. ims_npro_k > 1) then
            apu_async_size_k = int(imax, wi)*int(jmax, wi)*int(kmax, wi)
            allocate (apu_async_is_local_k(0:ims_npro_k - 1))
            allocate (apu_async_peer_k(0:ims_npro_k - 1))
            apu_async_peer_k = c_null_ptr
            apu_async_is_local_k = .false.
            ! Step (b): identify node-local peers via MPI_Comm_split_type.
            call MPI_Comm_split_type(ims_comm_z, MPI_COMM_TYPE_SHARED, ims_pro_k, MPI_INFO_NULL, apu_async_shmem_comm, ims_err)
            call MPI_Comm_size(apu_async_shmem_comm, apu_async_shmem_size, ims_err)
            allocate (apu_async_shmem_to_dir(0:apu_async_shmem_size - 1))
            call MPI_Allgather(ims_pro_k, 1, MPI_INTEGER, apu_async_shmem_to_dir, 1, MPI_INTEGER, &
                               apu_async_shmem_comm, ims_err)
            ! Take a dup of the shmem comm BEFORE window allocation.
            ! MPI_Win_allocate_shared taints the shmem comm's own barrier internals on Cray MPICH
            ! (same mechanism that taints the parent two-sided comm). A dup taken here is outside
            ! that taint and provides a reliable MPI_Barrier for GPU write synchronization.
            call MPI_Comm_dup(apu_async_shmem_comm, apu_async_node_comm_k, ims_err)
            ! Step (c): allocate the recv buffer as ONE shared mapping across the node-local
            ! sub-comm. MPI_Win_shared_query then returns each peer's segment baseptr as a
            ! c_ptr valid in EVERY node-local process's address space — that is the only
            ! portable way to get cross-process pointers on Cray MPICH. (Plain allocate +
            ! MPI_Allgather of addresses was tried 2026-05-15..17 and fails: the addresses are
            ! private to each process; even when transmitted bit-correctly they point into the
            ! receiver's own heap, producing the OpenACC "overlaps present region" crash.)
            call MPI_Win_allocate_shared(int(apu_async_size_k, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
                                         int(c_sizeof(1.0_dp)), MPI_INFO_NULL, apu_async_shmem_comm, win_baseptr, apu_async_win_k, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU_ASYNC/FABRIC_DIRECT K recv buffer.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
            ! Peer segments: MPI_Win_shared_query with each shmem-rank, then route by K-rank.
            do ip = 0, apu_async_shmem_size - 1
                apu_async_is_local_k(apu_async_shmem_to_dir(ip)) = .true.
                call MPI_Win_shared_query(apu_async_win_k, ip, win_query_size, win_disp_unit, &
                                          apu_async_peer_k(apu_async_shmem_to_dir(ip)), ims_err)
            end do
            ! Bug A fix (2026-05-18): bind apu_async_recv_k to OUR query result, not win_baseptr.
            ! On Cray MPICH some shmem subcomms return win_baseptr that is 1 segsize BEFORE the
            ! caller's actual segment (segments turn out to be 2x segsize apart with a guard
            ! region between them). MPI_Win_shared_query(my_shmem_rank,...) returns the real
            ! own-segment baseptr — that's what was stored in apu_async_peer_k(ims_pro_k) by
            ! the loop above (since shmem_to_dir maps OUR shmem-rank → ims_pro_k).
            call c_f_pointer(apu_async_peer_k(ims_pro_k), apu_async_recv_k, [apu_async_size_k])
            deallocate (apu_async_shmem_to_dir)
            ! apu_async_node_comm_k already set (dup taken before window alloc above).
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated K APU_ASYNC/FABRIC_DIRECT recv buffer.')
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
        real(dp), pointer, contiguous :: apu_pfptr_k(:) => null()   ! scratch pointer for APU_ASYNC per-peer writes
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
        ! APU_ASYNC/FABRIC_DIRECT: recv buffers come from a node-local        !
        !   MPI_Win_allocate_shared on the shmem sub-comm — cross-process     !
        !   addresses (via MPI_Win_shared_query) are valid on every node-     !
        !   local process and live in ONE contiguous mapping (no overlap of   !
        !   the 47 MB per-peer views in OpenACC's present-table); inter-node  !
        !   peers go via two-sided MPI on the dup'd parent comm.              !
        ! ==================================================================== !
#ifdef USE_APU
        if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
            ! -- Push: one fused GPU kernel writes our chunk to ALL peers simultaneously.
            ! Eliminates per-peer kernel-launch overhead (npro separate launches → 1).
            size = trp_plan%size3d
            call MPI_Win_fence(0, apu_win_k, ims_err)
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
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: our recv buffer is now fully populated
            ! -- Unpack: recv buffer is already in the flat K-space layout; one-to-one copy to b.
            !$omp target teams distribute parallel do
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
                                   trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_k, request(l), ims_err)
                end if
            end do
            ! Step 2: barrier — ensure all IRECVs are posted before proceeding
            call MPI_Barrier(apu_async_node_comm_k, ims_err)
            ! Step 3: inter-node — GPU-pack strided chunk into flat staging, then ISEND.
            !   Done BEFORE intra-node GPU writes so network transfer overlaps with GPU work.
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    flat_off = m * nmax_p * nlines_p
                    disp_ns  = m * nlines_p
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            c_wrk_dp(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp(flat_off + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_k, request(l), ims_err)
                end if
            end do
            ! Step 4: intra-node — GPU direct writes while inter-node MPI is in-flight
            do m = 0, ims_npro_k - 1
                if (apu_async_is_local_k(m)) then
                    call c_f_pointer(apu_async_peer_k(m), apu_pfptr_k, [apu_async_size_k])
                    flat_off = ims_pro_k * nmax_p * nlines_p
                    disp_ns  = m * nlines_p
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            apu_pfptr_k(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end do
            ! Step 5: barrier (intra-node writes done) then WAITALL (inter-node likely complete)
            call MPI_Barrier(apu_async_node_comm_k, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! Step 6: recv buffer is fully populated; flat copy to b
            !$omp target teams distribute parallel do
            do i = 1, size
                b(i) = apu_async_recv_k(i)
            end do
            !$omp end target teams distribute parallel do
            nullify (c_wrk_dp, apu_pfptr_k)

        else if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
            ! K-Forward FABRIC_DIRECT: K-intra peers span different XCDs on MI300A
            ! (global ranks 0,6,12,18 → XCDs 0,1,2,3). Cross-XCD CPU writes to the shared
            ! window are NOT coherent — each XCD has independent CPU caches and MPI_Win_fence
            ! only synchronizes MPI RMA ops, not raw pointer writes. GPU writes in separate
            ! kernels on different XCDs are similarly incoherent. So ALL K-peers go through
            ! plain MPI ISEND/IRECV; apu_async_recv_k is just used as the IRECV target buffer.
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
            l = 0
            ! Post IRECVs for ALL K-peers into our recv buffer.
            ! Use local K-comm rank m on apu_async_mpi_comm_k (untainted dup of ims_comm_z).
            do m = 0, ims_npro_k - 1
                l = l + 1
                call MPI_IRECV(apu_async_recv_k(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                               trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_k, request(l), ims_err)
            end do
            ! Pack ALL K-peers (strided a → flat c_wrk_dp) and post ISENDs.
            do m = 0, ims_npro_k - 1
                flat_off = m * nmax_p * nlines_p
                disp_ns  = m * nlines_p
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        c_wrk_dp(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                    end do
                end do
                l = l + 1
                call MPI_ISEND(c_wrk_dp(flat_off + 1), nmax_p*nlines_p, &
                               trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_k, request(l), ims_err)
            end do
            call MPI_WAITALL(l, request, status, ims_err)
            ! CPU copy recv buffer → b.
            do i = 1, size
                b(i) = apu_async_recv_k(i)
            end do
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
                    ! Step 2: GPU gather: a (strided) → c_wrk_dp (flat)
                    do m = 1, ims_npro_k
                        ns = maps_send_k(m) + 1
                        flat_off = (ns - 1)*nmax_p*nlines_p
                        disp_ns  = trp_plan%disp_s(ns)
#ifdef USE_APU
                        !$omp target teams distribute parallel do collapse(2)
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
        integer :: send_to, recv_from, fbd_tag
#ifdef USE_APU
        complex(dp), pointer :: apu_cx_all(:) => null()   ! complex view of apu_all_k across all peers
#endif
        type(MPI_Comm) :: trp_comm_k   ! dup of ims_comm_z; K-comm local rank = K-dir rank
#ifdef USE_APU
        if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT .or. trp_mode_k == TLAB_MPI_TRP_APU_ASYNC) then
            trp_comm_k = apu_async_mpi_comm_k   ! untainted dup of ims_comm_z; local K-rank = dir rank
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
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: recv buffer now fully populated
            ! Unpack: flat copy from complex-typed recv alias to b
            !$omp target teams distribute parallel do
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
                ! ISEND/IRECV in batches of trp_sizBlock_k peers per WAITALL.
                ! trp_comm_k is apu_async_mpi_comm_k (dup of ims_comm_z) for FABRIC_DIRECT/APU_ASYNC;
                ! ips/ipr are local K-comm ranks (0..npro_k-1), valid in all K-direction comms.
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
        real(dp), pointer, contiguous :: apu_pfptr_k(:) => null()   ! scratch pointer for APU_ASYNC per-peer writes
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
            !$omp target teams distribute parallel do collapse(2)
            do m = 0, ims_npro_k - 1
                do i = 1, nmax_p * nlines_p
                    ! Write our chunk (b[m*chunk]) into peer m's slot (own_rank*chunk) in their buffer
                    apu_all_k(m*apu_stride_k + ims_pro_k*nmax_p*nlines_p + i) = b(m * nmax_p * nlines_p + i)
                end do
            end do
            !$omp end target teams distribute parallel do
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: all peers have written to our buffer
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

        else if (trp_mode_k == TLAB_MPI_TRP_APU_ASYNC) then
            ! Hybrid backward: each rank pushes b[m*chunk] to intra-node peer m directly,
            ! or ISENDs b[m*chunk] to inter-node peers. After sync, apu_async_recv_k is
            ! in the same layout as in APU_DIRECT and is scattered to a.
            size = trp_plan%size3d
            l = 0
            ! Step 1: post IRECVs for inter-node peers (on the dup'd comm — see init notes)
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    l = l + 1
                    call MPI_IRECV(apu_async_recv_k(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_k, request(l), ims_err)
                end if
            end do
            ! Step 2: barrier — ensure all IRECVs posted before proceeding
            call MPI_Barrier(apu_async_node_comm_k, ims_err)
            ! Step 3: inter-node — ISEND b[m*chunk] directly (flat, no packing needed).
            !   Done BEFORE GPU writes so network transfer overlaps with GPU work.
            do m = 0, ims_npro_k - 1
                if (.not. apu_async_is_local_k(m)) then
                    l = l + 1
                    call MPI_ISEND(b(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_k, request(l), ims_err)
                end if
            end do
            ! Step 4: intra-node — GPU direct writes while inter-node MPI is in-flight
            do m = 0, ims_npro_k - 1
                if (apu_async_is_local_k(m)) then
                    call c_f_pointer(apu_async_peer_k(m), apu_pfptr_k, [apu_async_size_k])
                    flat_off = ims_pro_k * nmax_p * nlines_p
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p * nlines_p
                        apu_pfptr_k(flat_off + i) = b(m * nmax_p * nlines_p + i)
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end do
            ! Step 5: barrier (intra-node writes done) then WAITALL (inter-node likely complete)
            call MPI_Barrier(apu_async_node_comm_k, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! Step 6: unpack recv buffer — recv_k[m*chunk] = data from rank m → scatter to a
            do m = 0, ims_npro_k - 1
                flat_off = m * nmax_p * nlines_p
                disp_nr  = m * nlines_p
                !$omp target teams distribute parallel do collapse(2)
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        a(disp_nr + i*npage + j + 1) = apu_async_recv_k(flat_off + i*nlines_p + j + 1)
                    end do
                end do
                !$omp end target teams distribute parallel do
            end do
            nullify (apu_pfptr_k)

        else if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT) then
            ! K-Backward FABRIC_DIRECT: same all-MPI approach as K-Forward (K-intra peers
            ! are cross-XCD; shared-window CPU/GPU writes are not cache-coherent across XCDs).
            ! b is already flat per-peer (b(m*chunk + i)), so ISEND directly from b.
            size = trp_plan%size3d
            l = 0
            ! Post IRECVs for ALL K-peers into our recv buffer.
            ! Use local K-comm rank m on apu_async_mpi_comm_k (untainted dup of ims_comm_z).
            do m = 0, ims_npro_k - 1
                l = l + 1
                call MPI_IRECV(apu_async_recv_k(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                               trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_k, request(l), ims_err)
            end do
            ! Post ISENDs for ALL K-peers from b (already flat).
            do m = 0, ims_npro_k - 1
                l = l + 1
                call MPI_ISEND(b(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                               trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_k, request(l), ims_err)
            end do
            call MPI_WAITALL(l, request, status, ims_err)
            ! Scatter recv buffer → strided a.
            do m = 0, ims_npro_k - 1
                flat_off = m * nmax_p * nlines_p
                disp_nr  = m * nlines_p
                do i = 0, nmax_p - 1
                    do j = 0, nlines_p - 1
                        a(disp_nr + i*npage + j + 1) = apu_async_recv_k(flat_off + i*nlines_p + j + 1)
                    end do
                end do
            end do

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
#endif
        type(MPI_Comm) :: trp_comm_k   ! dup of ims_comm_z; K-comm local rank = K-dir rank
#ifdef USE_APU
        if (trp_mode_k == TLAB_MPI_TRP_FABRIC_DIRECT .or. trp_mode_k == TLAB_MPI_TRP_APU_ASYNC) then
            trp_comm_k = apu_async_mpi_comm_k   ! untainted dup; local K-rank = dir rank m
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
            call MPI_Win_fence(0, apu_win_k, ims_err)   ! barrier: recv buffer fully populated
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
        real(dp), pointer, contiguous :: apu_pfptr_i(:) => null()   ! scratch pointer for APU_ASYNC/FABRIC_DIRECT per-peer writes
        real(dp), pointer, contiguous :: c_wrk_dp_recv(:) => null() ! clean recv buffer (NOT shared-window-aliased) for MPI IRECVs
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
            !$omp target teams distribute parallel do collapse(2)
            do m = 0, ims_npro_i - 1
                do i = 1, nmax_p * nlines_p
                    ! a[m*chunk] is the flat chunk destined for peer m; write to their recv slot
                    apu_all_i(m*apu_stride_i + ims_pro_i*nmax_p*nlines_p + i) = a(m * nmax_p * nlines_p + i)
                end do
            end do
            !$omp end target teams distribute parallel do
            call MPI_Win_fence(0, apu_win_i, ims_err)   ! barrier: recv buffer fully populated
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

        else if (trp_mode_i == TLAB_MPI_TRP_APU_ASYNC) then
            ! Hybrid I-forward: intra-node peers via direct writes; inter-node via MPI.
            ! a is flat (a[m*chunk] for peer m); both paths fill apu_async_recv_i, then scatter to b.
            size = trp_plan%size3d
            l = 0
            ! Step 1: post IRECVs for inter-node peers (on the dup'd, untainted comm)
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    l = l + 1
                    call MPI_IRECV(apu_async_recv_i(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_i, request(l), ims_err)
                end if
            end do
            ! Step 2: barrier — ensure all IRECVs posted before intra-node writes begin
            call MPI_Barrier(apu_async_node_comm_i, ims_err)
            ! Step 3: inter-node — ISEND a[m*chunk] (already flat, no packing needed).
            !   Done BEFORE GPU intra writes so network transfer overlaps with GPU work.
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    l = l + 1
                    call MPI_ISEND(a(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_i, request(l), ims_err)
                end if
            end do
            ! Step 4: intra-node — push a[m*chunk] directly to peer m's recv buffer at our slot.
            !   Happens while inter-node MPI is in-flight.
            do m = 0, ims_npro_i - 1
                if (apu_async_is_local_i(m)) then
                    call c_f_pointer(apu_async_peer_i(m), apu_pfptr_i, [apu_async_size_i])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p * nlines_p
                        apu_pfptr_i(flat_off + i) = a(m * nmax_p * nlines_p + i)
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end do
            ! Step 5: barrier (intra-node writes done) then wait for inter-node MPI
            call MPI_Barrier(apu_async_node_comm_i, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! Step 6: unpack recv buffer → strided b
            do m = 0, ims_npro_i - 1
                flat_off = m * nmax_p * nlines_p
                disp_nr  = m * nmax_p
                !$omp target teams distribute parallel do collapse(2)
                do i = 0, nlines_p - 1
                    do j = 0, nmax_p - 1
                        b(disp_nr + i*nmax_full + j + 1) = apu_async_recv_i(flat_off + i*nmax_p + j + 1)
                    end do
                end do
                !$omp end target teams distribute parallel do
            end do
            nullify (apu_pfptr_i)

        else if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
            ! No shared window for I-direction: MPI_Win_allocate_shared on any I-direction
            ! sub-comm (XCD-local or full I-comm) corrupts CXI cross-XCD intra-node routing,
            ! causing WAITALL hangs (confirmed Hunter runs 2026-05-22/23). Use plain
            ! ISEND/IRECV on ims_comm_x (never had a window; fully clean) for all 6 I-peers.
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
            write(500+ims_pro,'(a,i4,a,g20.6)') '[IFR_pre] PE', ims_pro, ' sum(a)=', sum(a)
            flush(500+ims_pro)
            ! IRECVs
            l = 0
            do m = 1, ims_npro_i
                nr = maps_recv_i(m) + 1; ipr = nr - 1
                l = l + 1
                call MPI_IRECV(c_wrk_dp((nr - 1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                               trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
            end do
            write(500+ims_pro,'(a)') '[IFR_S4]'
            flush(500+ims_pro)
            ! ISENDs
            do m = 1, ims_npro_i
                ns = maps_send_i(m) + 1; ips = ns - 1
                l = l + 1
                call MPI_ISEND(a(trp_plan%disp_s(ns) + 1), nmax_p*nlines_p, &
                               trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
            end do
            write(500+ims_pro,'(a)') '[IFR_S8]'
            flush(500+ims_pro)
            call MPI_WAITALL(l, request, status, ims_err)
            write(500+ims_pro,'(a)') '[IFR_S9]'
            flush(500+ims_pro)
            ! Scatter c_wrk_dp → strided b
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
            write(500+ims_pro,'(a,i4,a,g20.6)') '[IFR_post] PE', ims_pro, ' sum(b)=', sum(b)
            flush(500+ims_pro)
            nullify (c_wrk_dp)

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
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan
        integer(wi) :: size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_nr
        integer :: send_to, recv_from, fbd_tag   ! FABRIC_DIRECT: global rank + distinct tag
#ifdef USE_APU
        ! apu_cx_all: complex view spanning all peers' shared windows (stride = apu_stride_i/2 complex units).
        complex(dp), pointer :: apu_cx_all(:) => null()
        integer(wi) :: mas
#endif
        type(MPI_Comm) :: trp_comm_i   ! MPI_COMM_WORLD for FABRIC_DIRECT (see K-Forward_Complex comment).
#ifdef USE_APU
        if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT .or. trp_mode_i == TLAB_MPI_TRP_APU_ASYNC) then
            trp_comm_i = apu_async_mpi_comm_i   ! untainted dup of ims_comm_x; local rank = I-dir rank
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
            !$omp target teams distribute parallel do collapse(2)
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
                ! trp_comm_i is apu_async_mpi_comm_i (local I-comm rank = dir rank directly).
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
        real(dp), pointer, contiguous :: apu_pfptr_i(:) => null()   ! scratch pointer for APU_ASYNC/FABRIC_DIRECT per-peer writes
        real(dp), pointer, contiguous :: c_wrk_dp_recv(:) => null() ! clean recv buffer (NOT shared-window-aliased) for MPI IRECVs
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
            ! Fence 1: open epoch — all ranks ready to receive direct writes.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! Push: pack strided b[m] → peer m's recv buffer at slot own_rank*chunk (flat).
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
            ! Fence 2: close epoch — all writes committed; recv buffers fully populated.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! Flat copy: recv buffer layout is flat and matches a 1:1.
            !$omp target teams distribute parallel do
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
            ! Step 1: post IRECVs for inter-node peers into their slots in the async recv window
            ! (on the dup'd, untainted comm).
            l = 0
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    l = l + 1
                    call MPI_IRECV(apu_async_recv_i(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_i, request(l), ims_err)
                end if
            end do
            ! Step 2: barrier — ensure all IRECVs posted before intra-node writes begin.
            call MPI_Barrier(apu_async_node_comm_i, ims_err)
            ! Step 3: inter-node — GPU pack strided b[m] → flat c_wrk_dp then ISEND.
            !   Done BEFORE GPU intra writes so network transfer overlaps with GPU work.
            do m = 0, ims_npro_i - 1
                if (.not. apu_async_is_local_i(m)) then
                    flat_off = m * nmax_p * nlines_p
                    disp_ns  = m * nmax_p
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            c_wrk_dp(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp(flat_off + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, m, ims_tag, apu_async_mpi_comm_i, request(l), ims_err)
                end if
            end do
            ! Step 4: intra-node — GPU pack strided b[m] into peer m's shared recv buffer at our slot.
            !   Happens while inter-node MPI is in-flight.
            do m = 0, ims_npro_i - 1
                if (apu_async_is_local_i(m)) then
                    call c_f_pointer(apu_async_peer_i(m), apu_pfptr_i, [apu_async_size_i])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    disp_ns  = m * nmax_p
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            apu_pfptr_i(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end if
            end do
            ! Step 5: barrier (intra-node writes done) then wait for inter-node sends/recvs.
            call MPI_Barrier(apu_async_node_comm_i, ims_err)
            if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
            ! Step 6: flat copy from unified async recv buffer to a (1:1).
            !$omp target teams distribute parallel do
            do i = 1, size
                a(i) = apu_async_recv_i(i)
            end do
            !$omp end target teams distribute parallel do
            nullify (c_wrk_dp, apu_pfptr_i)
        else if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT) then
            ! No shared window for I-direction — same reason as IFR.
            ! Use plain ISEND/IRECV on ims_comm_x for all 6 I-peers.
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
            write(500+ims_pro,'(a,i4,a,g20.6)') '[IBR_pre] PE', ims_pro, ' sum(b)=', sum(b)
            flush(500+ims_pro)
            ! IRECVs into a directly (recv layout mirrors disp_s)
            l = 0
            do m = 1, ims_npro_i
                nr = maps_send_i(m) + 1; ipr = nr - 1
                l = l + 1
                call MPI_IRECV(a(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, &
                               trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
            end do
            write(500+ims_pro,'(a)') '[IBR_S4]'
            flush(500+ims_pro)
            ! Pack b → c_wrk_dp (strided → flat)
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
            ! ISENDs from packed flat c_wrk_dp
            do m = 1, ims_npro_i
                ns = maps_recv_i(m) + 1; ips = ns - 1
                l = l + 1
                call MPI_ISEND(c_wrk_dp((ns - 1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                               trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
            end do
            write(500+ims_pro,'(a)') '[IBR_S8]'
            flush(500+ims_pro)
            call MPI_WAITALL(l, request, status, ims_err)
            write(500+ims_pro,'(a)') '[IBR_S9]'
            flush(500+ims_pro)
            write(500+ims_pro,'(a,i4,a,g20.6)') '[IBR_post] PE', ims_pro, ' sum(a)=', sum(a)
            flush(500+ims_pro)
            nullify (c_wrk_dp)
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
                        !$omp target teams distribute parallel do collapse(2)
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
        integer :: send_to, recv_from, fbd_tag   ! FABRIC_DIRECT: global rank + distinct tag
#ifdef USE_APU
        ! apu_cx_all: complex view spanning all peers' shared windows (stride = apu_stride_i/2 complex units).
        complex(dp), pointer :: apu_cx_all(:) => null()
        integer(wi) :: mas
#endif
        type(MPI_Comm) :: trp_comm_i   ! MPI_COMM_WORLD for FABRIC_DIRECT (see K-Forward_Complex comment).
#ifdef USE_APU
        if (trp_mode_i == TLAB_MPI_TRP_FABRIC_DIRECT .or. trp_mode_i == TLAB_MPI_TRP_APU_ASYNC) then
            trp_comm_i = apu_async_mpi_comm_i
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
            ! Fence 2: close epoch — all writes committed; recv buffers fully populated.
            call MPI_Win_fence(0, apu_win_i, ims_err)
            ! Flat copy: recv buffer is flat and matches a layout 1:1.
            !$omp target teams distribute parallel do
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
