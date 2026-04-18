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
    ! c_ptr is accessible via mpi_f08 (which re-exports iso_c_binding); declaring it here again causes ambiguity
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
    real(dp), pointer :: apu_async_recv_k(:) => null(), apu_async_recv_i(:) => null()
    type(c_ptr), allocatable :: apu_async_peer_k(:), apu_async_peer_i(:)
    logical, allocatable :: apu_async_is_local_k(:), apu_async_is_local_i(:)
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
    real(dp), pointer    :: c_wrk_dp(:) => null()
    complex(dp), pointer :: c_wrk_cx(:) => null()
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
        ! local PE mappings for explicit send/recv
        allocate (maps_send_i(ims_npro_i))
        allocate (maps_recv_i(ims_npro_i))
        do ip = 0, ims_npro_i - 1
            maps_send_i(ip + 1) = ip
            maps_recv_i(ip + 1) = mod(ims_npro_i - ip, ims_npro_i)
        end do
        maps_send_i = cshift(maps_send_i, ims_pro_i)
        maps_recv_i = cshift(maps_recv_i, -ims_pro_i)

        allocate (maps_send_k(ims_npro_k))
        allocate (maps_recv_k(ims_npro_k))
        do ip = 0, ims_npro_k - 1
            maps_send_k(ip + 1) = ip
            maps_recv_k(ip + 1) = mod(ims_npro_k - ip, ims_npro_k)
        end do
        maps_send_k = cshift(maps_send_k, ims_pro_k)
        maps_recv_k = cshift(maps_recv_k, -ims_pro_k)

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
            trp_mode_i == TLAB_MPI_TRP_APU_ASYNC    .or. trp_mode_k == TLAB_MPI_TRP_APU_ASYNC) then
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
            apu_size_k = int(imax, wi)*int(jmax, wi)*int(kmax, wi)
            call MPI_Win_allocate_shared(int(apu_size_k, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
                                         int(c_sizeof(1.0_dp)), MPI_INFO_NULL, ims_comm_z, win_baseptr, apu_win_k, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU K recv buffer.')
                call TLab_Stop(DNS_ERROR_ALLOC)
            end if
            call c_f_pointer(win_baseptr, apu_recv_fptr_k, [apu_size_k])
            allocate (apu_peer_cptr_k(0:ims_npro_k - 1))
            do ip = 0, ims_npro_k - 1
                call MPI_Win_shared_query(apu_win_k, ip, win_query_size, win_disp_unit, apu_peer_cptr_k(ip), ims_err)
            end do
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated K APU direct recv buffer.')
        end if
        if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT .and. ims_npro_i > 1) then
            apu_size_i = int(imax, wi)*int(jmax, wi)*int(kmax, wi)
            call MPI_Win_allocate_shared(int(apu_size_i, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
                                         int(c_sizeof(1.0_dp)), MPI_INFO_NULL, ims_comm_x, win_baseptr, apu_win_i, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU I recv buffer.')
                call TLab_Stop(DNS_ERROR_ALLOC)
            end if
            call c_f_pointer(win_baseptr, apu_recv_fptr_i, [apu_size_i])
            allocate (apu_peer_cptr_i(0:ims_npro_i - 1))
            do ip = 0, ims_npro_i - 1
                call MPI_Win_shared_query(apu_win_i, ip, win_query_size, win_disp_unit, apu_peer_cptr_i(ip), ims_err)
            end do
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated I APU direct recv buffer.')
        end if
        ! APU_ASYNC: node-local shared windows using MPI_COMM_TYPE_SHARED sub-communicators.
        ! Each rank splits ims_comm_z/x into a node-local communicator; the shared window
        ! is allocated only over that sub-communicator. Ranks absent from it are inter-node
        ! peers handled via standard ISEND/IRECV.
        if (trp_mode_k == TLAB_MPI_TRP_APU_ASYNC .and. ims_npro_k > 1) then
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
            deallocate (apu_async_src_ranks, apu_async_trans_ranks)
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
            do ip = 0, ims_npro_k - 1
                if (apu_async_is_local_k(ip)) &
                    call MPI_Win_shared_query(apu_async_win_k, ip, win_query_size, win_disp_unit, apu_async_peer_k(ip), ims_err)
            end do
            call MPI_Comm_free(apu_async_shmem_comm, ims_err)
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated K APU_ASYNC recv buffer.')
        end if
        if (trp_mode_i == TLAB_MPI_TRP_APU_ASYNC .and. ims_npro_i > 1) then
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
            deallocate (apu_async_src_ranks, apu_async_trans_ranks)
            call MPI_Group_free(apu_async_group_dir, ims_err)
            call MPI_Group_free(apu_async_group_shmem, ims_err)
            call MPI_Win_allocate_shared(int(apu_async_size_i, MPI_ADDRESS_KIND)*int(c_sizeof(1.0_dp), MPI_ADDRESS_KIND), &
                                         int(c_sizeof(1.0_dp)), MPI_INFO_NULL, apu_async_shmem_comm, win_baseptr, apu_async_win_i, ims_err)
            if (ims_err /= MPI_SUCCESS) then
                call TLab_Write_ASCII(efile, __FILE__//'. MPI_Win_allocate_shared failed for APU_ASYNC I recv buffer.')
                call TLab_Stop(DNS_ERROR_ALLOC)
            end if
            call c_f_pointer(win_baseptr, apu_async_recv_i, [apu_async_size_i])
            do ip = 0, ims_npro_i - 1
                if (apu_async_is_local_i(ip)) &
                    call MPI_Win_shared_query(apu_async_win_i, ip, win_query_size, win_disp_unit, apu_async_peer_i(ip), ims_err)
            end do
            call MPI_Comm_free(apu_async_shmem_comm, ims_err)
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated I APU_ASYNC recv buffer.')
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
        real(wp), intent(in) :: a(:)
        real(wp), intent(out) :: b(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        ! -----------------------------------------------------------------------
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, npage, flat_off, disp_ns, mas
#ifdef USE_APU
        real(dp), pointer :: apu_pfptr_k(:) => null()
#endif

#ifdef PROFILE_ON
        real(wp) time_loc_1, time_loc_2
#endif

        ! #######################################################################
#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif

        nmax_p   = trp_plan%nmax
        nlines_p = trp_plan%nlines
        npage    = nlines_p * ims_npro_k   ! total lines spanning all K ranks
        mas      = nmax_p * nlines_p       ! elements per peer message

        if (trp_datatype_k == MPI_REAL4 .and. wp == dp .and. &
            trp_mode_k /= TLAB_MPI_TRP_APU_DIRECT .and. trp_mode_k /= TLAB_MPI_TRP_APU_ASYNC) then
            size = trp_plan%size3d
            a_wrk => wrk_mpi_fptr(1:size)
            b_wrk => wrk_mpi_fptr(size + 1:2*size)
            c_wrk => wrk_mpi_fptr(2*size + 1:3*size)
            ! dp→sp conversion
#ifdef USE_APU
            !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
#endif
            do i = 1, size
                a_wrk(i) = real(a(i), sp)
            end do
#ifdef USE_APU
            !$omp end target teams distribute parallel do
#endif
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! Pipeline: IRECVs posted before GPU pack, ISENDs posted after pack, single WAITALL.
                ! Recv buffers (b_wrk) are independent of the pack target (c_wrk), so IRECVs can
                ! be registered with MPI while the GPU gathers the strided send data into c_wrk.
                ! Step 1: post all IRECVs — recv buffers are ready before pack starts
                l = 0
                do m = 1, ims_npro_k
                    nr = maps_recv_k(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(b_wrk(trp_plan%disp_r(nr) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ipr, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                ! Step 2: GPU pack a_wrk→c_wrk (strided→flat) while network prepares recv buffers
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
                ! Step 3: post all ISENDs from packed flat c_wrk
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
            ! sp→dp conversion
#ifdef USE_APU
            !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
#endif
            do i = 1, size
                b(i) = real(b_wrk(i), dp)
            end do
#ifdef USE_APU
            !$omp end target teams distribute parallel do
#endif
            nullify (a_wrk, b_wrk, c_wrk)
        else
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! dp path: pipeline IRECVs → GPU pack a→c_wrk_dp → ISENDs → WAITALL
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                ! Step 1: post all IRECVs before pack starts
                l = 0
                do m = 1, ims_npro_k
                    nr = maps_recv_k(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(b(trp_plan%disp_r(nr) + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ipr, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                ! Step 2: GPU pack a→c_wrk_dp (strided→flat) while network prepares recv buffers
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
                ! Step 3: post all ISENDs from packed flat c_wrk_dp
                do m = 1, ims_npro_k
                    ns = maps_send_k(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ips, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                ! Step 4: single WAITALL for all sends and receives
                call MPI_WAITALL(l, request, status, ims_err)
                nullify (c_wrk_dp)
#ifdef USE_APU
            else if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
                ! Each rank writes its own chunk into every peer's recv buffer at slot own_rank*chunk.
                ! After Win_fence, apu_recv_fptr_k[r*chunk] = data from rank r (all peers complete).
                size = trp_plan%size3d
                call MPI_Win_fence(0, apu_win_k, ims_err)
                do m = 0, ims_npro_k - 1
                    call c_f_pointer(apu_peer_cptr_k(m), apu_pfptr_k, [apu_size_k])
                    flat_off = ims_pro_k * nmax_p * nlines_p
                    disp_ns  = m * nlines_p
                    !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            apu_pfptr_k(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end do
                call MPI_Win_fence(0, apu_win_k, ims_err)
                ! recv buffer is now fully populated: flat copy to b
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
                do i = 1, size
                    b(i) = apu_recv_fptr_k(i)
                end do
                !$omp end target teams distribute parallel do
                nullify (apu_pfptr_k)
            else if (trp_mode_k == TLAB_MPI_TRP_APU_ASYNC) then
                ! Hybrid: intra-node peers receive via direct shared-memory writes (APU_DIRECT style);
                ! inter-node peers receive via MPI ISEND/IRECV. Unified recv buffer apu_async_recv_k
                ! holds all data after both synchronisation barriers complete.
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                ! Step 1: post IRECVs for inter-node peers before Win_fence (non-blocking, independent)
                l = 0
                do m = 0, ims_npro_k - 1
                    if (.not. apu_async_is_local_k(m)) then
                        l = l + 1
                        call MPI_IRECV(apu_async_recv_k(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, m, ims_tag, ims_comm_z, request(l), ims_err)
                    end if
                end do
                ! Step 2: open shared-window epoch (collective over node-local comm)
                call MPI_Win_fence(0, apu_async_win_k, ims_err)
                ! Step 3: intra-node — push our chunk to each intra-node peer's recv buffer slot
                do m = 0, ims_npro_k - 1
                    if (apu_async_is_local_k(m)) then
                        call c_f_pointer(apu_async_peer_k(m), apu_pfptr_k, [apu_async_size_k])
                        flat_off = ims_pro_k * nmax_p * nlines_p
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
                ! Step 4: inter-node — pack strided chunk into flat buffer and ISEND
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
                ! Step 5: close epoch (all intra-node writes complete) and wait for inter-node
                call MPI_Win_fence(0, apu_async_win_k, ims_err)
                if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
                ! Step 6: flat copy from unified recv buffer to b
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
                do i = 1, size
                    b(i) = apu_async_recv_k(i)
                end do
                !$omp end target teams distribute parallel do
                nullify (c_wrk_dp, apu_pfptr_k)
#endif
            else
                call Transpose_Kernel_Double(a(1:trp_plan%size3d), maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                             b(1:trp_plan%size3d), maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                             ims_comm_z, trp_sizBlock_k, trp_mode_k)
            end if
        end if

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine TLabMPI_Trp_ExecK_Forward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Forward_Complex(a, b, trp_plan)
        complex(wp), intent(in) :: a(*)
        complex(wp), intent(out) :: b(*)
        type(tmpi_transpose_dt), intent(in) :: trp_plan
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, npage, flat_off, disp_ns

        ! #######################################################################
        nmax_p   = trp_plan%nmax
        nlines_p = trp_plan%nlines
        npage    = nlines_p * ims_npro_k
        ! APU_DIRECT and APU_ASYNC fall back to ASYNC for complex (dp path, no sp conversion needed)
        if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS .or. trp_mode_k == TLAB_MPI_TRP_APU_DIRECT .or. &
            trp_mode_k == TLAB_MPI_TRP_APU_ASYNC) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
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
        return
    end subroutine TLabMPI_Trp_ExecK_Forward_Complex

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Backward_Real(b, a, trp_plan)
        real(wp), intent(in) :: b(:)
        real(wp), intent(out) :: a(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        ! -----------------------------------------------------------------------
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, npage, flat_off, disp_nr, mas
#ifdef USE_APU
        real(dp), pointer :: apu_pfptr_k(:) => null()
#endif
#ifdef PROFILE_ON
        real(wp) time_loc_1, time_loc_2
#endif

        ! #######################################################################
#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif

        nmax_p   = trp_plan%nmax
        nlines_p = trp_plan%nlines
        npage    = nlines_p * ims_npro_k
        mas      = nmax_p * nlines_p

        if (trp_datatype_k == MPI_REAL4 .and. wp == dp .and. &
            trp_mode_k /= TLAB_MPI_TRP_APU_DIRECT .and. trp_mode_k /= TLAB_MPI_TRP_APU_ASYNC) then
            size = trp_plan%size3d
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! Q1: merge unpack + sp→dp into one pass; only 2 sp slots needed.
                b_wrk => wrk_mpi_fptr(1:size)           ! send: dp→sp copy of b
                c_wrk => wrk_mpi_fptr(size + 1:2*size)  ! recv: flat staging
                ! dp→sp conversion
#ifdef USE_APU
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
#endif
                do i = 1, size
                    b_wrk(i) = real(b(i), sp)
                end do
#ifdef USE_APU
                !$omp end target teams distribute parallel do
#endif
                ! K-backward: send flat from b_wrk, recv into c_wrk, then unpack+sp→dp directly to a.
                ! Step 1: post all IRECVs into flat c_wrk recv slots
                l = 0
                do m = 1, ims_npro_k
                    nr = maps_send_k(m) + 1; ipr = nr - 1   ! backward: mrecv=maps_send_k
                    l = l + 1
                    call MPI_IRECV(c_wrk((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ipr, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                ! Step 2: post all ISENDs from flat b_wrk send slots
                do m = 1, ims_npro_k
                    ns = maps_recv_k(m) + 1; ips = ns - 1   ! backward: msend=maps_recv_k
                    l = l + 1
                    call MPI_ISEND(b_wrk(trp_plan%disp_r(ns) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ips, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                ! Step 3: single WAITALL for all sends and receives
                call MPI_WAITALL(l, request, status, ims_err)
                ! GPU unpack: flat c_wrk → a (dp), unpack+sp→dp merged in one pass
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
                ! dp→sp conversion
#ifdef USE_APU
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
#endif
                do i = 1, size
                    b_wrk(i) = real(b(i), sp)
                end do
#ifdef USE_APU
                !$omp end target teams distribute parallel do
#endif
                call Transpose_Kernel_Single(b_wrk, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                             a_wrk, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                             ims_comm_z, trp_sizBlock_k, trp_mode_k)
                ! sp→dp conversion
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
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! dp path: flat send from b, flat recv into c_wrk_dp, GPU unpack c_wrk_dp→a
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                ! Step 1: post all IRECVs into flat c_wrk_dp recv slots
                l = 0
                do m = 1, ims_npro_k
                    nr = maps_send_k(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(c_wrk_dp((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ipr, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                ! Step 2: post all ISENDs from flat b send slots
                do m = 1, ims_npro_k
                    ns = maps_recv_k(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(b(trp_plan%disp_r(ns) + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ips, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                ! Step 3: single WAITALL for all sends and receives
                call MPI_WAITALL(l, request, status, ims_err)
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
#ifdef USE_APU
            else if (trp_mode_k == TLAB_MPI_TRP_APU_DIRECT) then
                ! Reverse of K-Forward: b is in the fully-assembled K-space layout.
                ! Each rank writes its own flat chunk into every peer's recv buffer at slot own_rank*chunk.
                size = trp_plan%size3d
                call MPI_Win_fence(0, apu_win_k, ims_err)
                do m = 0, ims_npro_k - 1
                    call c_f_pointer(apu_peer_cptr_k(m), apu_pfptr_k, [apu_size_k])
                    flat_off = ims_pro_k * nmax_p * nlines_p
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
                    do i = 1, nmax_p * nlines_p
                        apu_pfptr_k(flat_off + i) = b(m * nmax_p * nlines_p + i)
                    end do
                    !$omp end target teams distribute parallel do
                end do
                call MPI_Win_fence(0, apu_win_k, ims_err)
                do m = 0, ims_npro_k - 1
                    flat_off = m * nmax_p * nlines_p
                    disp_nr  = m * nlines_p
                    !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            a(disp_nr + i*npage + j + 1) = apu_recv_fptr_k(flat_off + i*nlines_p + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end do
                nullify (apu_pfptr_k)
            else if (trp_mode_k == TLAB_MPI_TRP_APU_ASYNC) then
                ! Hybrid backward: each rank pushes b[m*chunk] to intra-node peer m's recv buffer at
                ! our slot (ims_pro_k*chunk), and ISENDs b[m*chunk] to inter-node peers.
                ! After sync, apu_async_recv_k[m*chunk] = b[ims_pro_k*chunk] from rank m → unpack to a.
                size = trp_plan%size3d
                ! Step 1: post IRECVs from inter-node peers into recv buffer (they send b[ims_pro_k*chunk])
                l = 0
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
                ! Step 4: inter-node — ISEND b[m*chunk] to peer m
                do m = 0, ims_npro_k - 1
                    if (.not. apu_async_is_local_k(m)) then
                        l = l + 1
                        call MPI_ISEND(b(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, m, ims_tag, ims_comm_z, request(l), ims_err)
                    end if
                end do
                ! Step 5: close epoch and wait for inter-node
                call MPI_Win_fence(0, apu_async_win_k, ims_err)
                if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
                ! Step 6: unpack recv buffer — recv_k[m*chunk] contains data from rank m
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
#endif
            else
                call Transpose_Kernel_Double(b(1:trp_plan%size3d), maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                             a(1:trp_plan%size3d), maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                             ims_comm_z, trp_sizBlock_k, trp_mode_k)
            end if
        end if

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine TLabMPI_Trp_ExecK_Backward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Backward_Complex(b, a, trp_plan)
        complex(wp), intent(in) :: b(*)
        complex(wp), intent(out) :: a(*)
        type(tmpi_transpose_dt), intent(in) :: trp_plan
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, npage, flat_off, disp_nr

        ! #######################################################################
        nmax_p   = trp_plan%nmax
        nlines_p = trp_plan%nlines
        npage    = nlines_p * ims_npro_k
        ! APU_DIRECT and APU_ASYNC fall back to ASYNC for complex
        if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS .or. trp_mode_k == TLAB_MPI_TRP_APU_DIRECT .or. &
            trp_mode_k == TLAB_MPI_TRP_APU_ASYNC) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
            do j = 1, ims_npro_k, trp_sizBlock_k
                l = 0
                do m = j, min(j + trp_sizBlock_k - 1, ims_npro_k)
                    ns = maps_recv_k(m) + 1; ips = ns - 1
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
        return
    end subroutine TLabMPI_Trp_ExecK_Backward_Complex

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Forward_Real(a, b, trp_plan)
        real(wp), dimension(:), intent(in) :: a
        real(wp), dimension(:), intent(out) :: b
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        ! -----------------------------------------------------------------------
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_nr, mas
#ifdef USE_APU
        real(dp), pointer :: apu_pfptr_i(:) => null()
#endif

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i   ! total elements per line across all I ranks (stride in recv)
        mas       = nmax_p * nlines_p

        if (trp_datatype_i == MPI_REAL4 .and. wp == dp .and. &
            trp_mode_i /= TLAB_MPI_TRP_APU_DIRECT .and. trp_mode_i /= TLAB_MPI_TRP_APU_ASYNC) then
            size = trp_plan%size3d
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! Q1: merge unpack + sp→dp into one pass; only 2 sp slots needed.
                a_wrk => wrk_mpi_fptr(1:size)           ! send: dp→sp copy of a
                c_wrk => wrk_mpi_fptr(size + 1:2*size)  ! recv: flat staging
                ! dp→sp
#ifdef USE_APU
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
#endif
                do i = 1, size
                    a_wrk(i) = real(a(i), sp)
                end do
#ifdef USE_APU
                !$omp end target teams distribute parallel do
#endif
                ! I-forward: send flat from a_wrk, recv into c_wrk, then unpack+sp→dp directly to b.
                ! Step 1: post all IRECVs into flat c_wrk recv slots
                l = 0
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(c_wrk((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ipr, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                ! Step 2: post all ISENDs from flat a_wrk send slots
                do m = 1, ims_npro_i
                    ns = maps_send_i(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(a_wrk(trp_plan%disp_s(ns) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ips, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                ! Step 3: single WAITALL for all sends and receives
                call MPI_WAITALL(l, request, status, ims_err)
                ! GPU unpack: flat c_wrk → b (dp), unpack+sp→dp merged in one pass
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
                ! dp→sp
#ifdef USE_APU
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
#endif
                do i = 1, size
                    a_wrk(i) = real(a(i), sp)
                end do
#ifdef USE_APU
                !$omp end target teams distribute parallel do
#endif
                call Transpose_Kernel_Single(a_wrk, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                             b_wrk, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                             ims_comm_x, trp_sizBlock_i, trp_mode_i)
                ! sp→dp
#ifdef USE_APU
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
#endif
                do i = 1, size
                    b(i) = real(b_wrk(i), dp)
                end do
#ifdef USE_APU
                !$omp end target teams distribute parallel do
#endif
                nullify (a_wrk, b_wrk)
            end if
        else
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! dp path: send flat from a, recv into c_wrk_dp, GPU unpack c_wrk_dp→b.
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                ! Step 1: post all IRECVs into flat c_wrk_dp recv slots
                l = 0
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(c_wrk_dp((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                ! Step 2: post all ISENDs from flat a send slots
                do m = 1, ims_npro_i
                    ns = maps_send_i(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(a(trp_plan%disp_s(ns) + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                ! Step 3: single WAITALL for all sends and receives
                call MPI_WAITALL(l, request, status, ims_err)
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
#ifdef USE_APU
            else if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT) then
                ! Each rank writes its own flat chunk into every peer's recv buffer at slot own_rank*chunk.
                ! After Win_fence, apu_recv_fptr_i[r*chunk] = data from rank r; unpack into strided b.
                size = trp_plan%size3d
                call MPI_Win_fence(0, apu_win_i, ims_err)
                do m = 0, ims_npro_i - 1
                    call c_f_pointer(apu_peer_cptr_i(m), apu_pfptr_i, [apu_size_i])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
                    do i = 1, nmax_p * nlines_p
                        apu_pfptr_i(flat_off + i) = a(m * nmax_p * nlines_p + i)
                    end do
                    !$omp end target teams distribute parallel do
                end do
                call MPI_Win_fence(0, apu_win_i, ims_err)
                do m = 0, ims_npro_i - 1
                    flat_off = m * nmax_p * nlines_p
                    disp_nr  = m * nmax_p
                    !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            b(disp_nr + i*nmax_full + j + 1) = apu_recv_fptr_i(flat_off + i*nmax_p + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end do
                nullify (apu_pfptr_i)
            else if (trp_mode_i == TLAB_MPI_TRP_APU_ASYNC) then
                ! Hybrid I-forward: a is flat (rank m's chunk at a[m*chunk]).
                ! Intra-node: push a[m*chunk] directly to peer m's recv buffer at our slot.
                ! Inter-node: IRECV then ISEND a[m*chunk] via MPI.
                ! After sync, recv_i[m*chunk] = a[m*chunk] from rank m → strided unpack to b.
                size = trp_plan%size3d
                ! Step 1: post IRECVs for inter-node peers
                l = 0
                do m = 0, ims_npro_i - 1
                    if (.not. apu_async_is_local_i(m)) then
                        l = l + 1
                        call MPI_IRECV(apu_async_recv_i(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                    end if
                end do
                ! Step 2: open shared-window epoch
                call MPI_Win_fence(0, apu_async_win_i, ims_err)
                ! Step 3: intra-node — push a[m*chunk] to peer m's recv buffer at our slot
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
                ! Step 4: inter-node — ISEND a[m*chunk] to peer m (contiguous, no packing needed)
                do m = 0, ims_npro_i - 1
                    if (.not. apu_async_is_local_i(m)) then
                        l = l + 1
                        call MPI_ISEND(a(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                    end if
                end do
                ! Step 5: close epoch and wait for inter-node
                call MPI_Win_fence(0, apu_async_win_i, ims_err)
                if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
                ! Step 6: strided unpack from unified recv buffer to b
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
#endif
            else
                call Transpose_Kernel_Double(a(1:trp_plan%size3d), maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                             b(1:trp_plan%size3d), maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                             ims_comm_x, trp_sizBlock_i, trp_mode_i)
            end if
        end if

        return
    end subroutine TLabMPI_Trp_ExecI_Forward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Forward_Complex(a, b, trp_plan)
        complex(wp), intent(in) :: a(*)
        complex(wp), intent(out) :: b(*)
        type(tmpi_transpose_dt), intent(in) :: trp_plan
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_nr

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i
        ! APU_DIRECT and APU_ASYNC fall back to ASYNC for complex
        if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS .or. trp_mode_i == TLAB_MPI_TRP_APU_DIRECT .or. &
            trp_mode_i == TLAB_MPI_TRP_APU_ASYNC) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
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
        return
    end subroutine TLabMPI_Trp_ExecI_Forward_Complex

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Backward_Real(b, a, trp_plan)
        real(wp), intent(in) :: b(:)
        real(wp), intent(out) :: a(:)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        ! -----------------------------------------------------------------------
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_ns, mas
#ifdef USE_APU
        real(dp), pointer :: apu_pfptr_i(:) => null()
#endif

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i
        mas       = nmax_p * nlines_p

        if (trp_datatype_i == MPI_REAL4 .and. wp == dp .and. &
            trp_mode_i /= TLAB_MPI_TRP_APU_DIRECT .and. trp_mode_i /= TLAB_MPI_TRP_APU_ASYNC) then
            size = trp_plan%size3d
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! Q1: pack directly from b (dp→sp+strided→flat in one pass); only 2 sp slots needed.
                a_wrk => wrk_mpi_fptr(1:size)           ! recv: flat staging → sp→dp to a
                c_wrk => wrk_mpi_fptr(size + 1:2*size)  ! send: packed from b dp→sp+strided→flat
                ! Pipeline: IRECVs posted before GPU pack, ISENDs posted after pack, single WAITALL.
                ! Step 1: post all IRECVs into flat a_wrk recv slots before pack starts
                l = 0
                do m = 1, ims_npro_i
                    nr = maps_send_i(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(a_wrk(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ipr, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                ! Step 2: GPU pack b→c_wrk (dp→sp+strided→flat in one pass) while recv buffers prepare
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
                ! Step 3: post all ISENDs from packed flat c_wrk
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(c_wrk((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ips, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                ! Step 4: single WAITALL for all sends and receives
                call MPI_WAITALL(l, request, status, ims_err)
                ! sp→dp: 1:1 conversion of flat recv buffer to a
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
                b_wrk => wrk_mpi_fptr(1:size)
                a_wrk => wrk_mpi_fptr(size + 1:2*size)
                ! dp→sp
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
                ! sp→dp
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
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! dp path: pipeline IRECVs → GPU pack b→c_wrk_dp → ISENDs → WAITALL
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                ! Step 1: post all IRECVs into flat a recv slots before pack starts
                l = 0
                do m = 1, ims_npro_i
                    nr = maps_send_i(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(a(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                ! Step 2: GPU pack b→c_wrk_dp (strided→flat) while network prepares recv buffers
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
                ! Step 3: post all ISENDs from packed flat c_wrk_dp
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                ! Step 4: single WAITALL for all sends and receives
                call MPI_WAITALL(l, request, status, ims_err)
                nullify (c_wrk_dp)
#ifdef USE_APU
            else if (trp_mode_i == TLAB_MPI_TRP_APU_DIRECT) then
                ! Reverse of I-Forward: b is in fully-assembled X-space layout (strided).
                ! Each rank packs its strided chunk from b into every peer's recv buffer at slot own_rank*chunk.
                size = trp_plan%size3d
                call MPI_Win_fence(0, apu_win_i, ims_err)
                do m = 0, ims_npro_i - 1
                    call c_f_pointer(apu_peer_cptr_i(m), apu_pfptr_i, [apu_size_i])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    disp_ns  = m * nmax_p
                    !$omp target teams distribute parallel do collapse(2) if(mas * sizeofreal > 100000_wi)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            apu_pfptr_i(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end do
                call MPI_Win_fence(0, apu_win_i, ims_err)
                ! recv buffer layout matches a (flat, 1:1 copy)
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
                do i = 1, size
                    a(i) = apu_recv_fptr_i(i)
                end do
                !$omp end target teams distribute parallel do
                nullify (apu_pfptr_i)
            else if (trp_mode_i == TLAB_MPI_TRP_APU_ASYNC) then
                ! Hybrid I-backward: b is in X-space layout (strided, disp_r stride).
                ! Each rank packs strided b[m] into peer m's recv buffer at our slot (intra)
                ! or into c_wrk_dp for ISEND (inter). After sync, recv_i is flat → a.
                size = trp_plan%size3d
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                ! Step 1: post IRECVs for inter-node peers (they send b[ims_pro_i*chunk] strided→flat)
                l = 0
                do m = 0, ims_npro_i - 1
                    if (.not. apu_async_is_local_i(m)) then
                        l = l + 1
                        call MPI_IRECV(apu_async_recv_i(m*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                       trp_plan%base_type, m, ims_tag, ims_comm_x, request(l), ims_err)
                    end if
                end do
                ! Step 2: open shared-window epoch
                call MPI_Win_fence(0, apu_async_win_i, ims_err)
                ! Step 3: intra-node — pack strided b[m] into peer m's recv buffer at our slot
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
                ! Step 4: inter-node — pack strided b[m] into flat buffer and ISEND
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
                ! Step 5: close epoch and wait for inter-node
                call MPI_Win_fence(0, apu_async_win_i, ims_err)
                if (l > 0) call MPI_WAITALL(l, request, status, ims_err)
                ! Step 6: flat copy from unified recv buffer to a
                !$omp target teams distribute parallel do if(mas * sizeofreal > 100000_wi)
                do i = 1, size
                    a(i) = apu_async_recv_i(i)
                end do
                !$omp end target teams distribute parallel do
                nullify (c_wrk_dp, apu_pfptr_i)
#endif
            else
                call Transpose_Kernel_Double(b(1:trp_plan%size3d), maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                             a(1:trp_plan%size3d), maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                             ims_comm_x, trp_sizBlock_i, trp_mode_i)
            end if
        end if

        return
    end subroutine TLabMPI_Trp_ExecI_Backward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Backward_Complex(b, a, trp_plan)
        complex(wp), intent(in) :: b(*)
        complex(wp), intent(out) :: a(*)
        type(tmpi_transpose_dt), intent(in) :: trp_plan
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_ns

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i
        ! APU_DIRECT and APU_ASYNC fall back to ASYNC for complex
        if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS .or. trp_mode_i == TLAB_MPI_TRP_APU_DIRECT .or. &
            trp_mode_i == TLAB_MPI_TRP_APU_ASYNC) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
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
