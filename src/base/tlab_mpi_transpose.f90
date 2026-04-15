#include "dns_error.h"

! Circular transposition within directional communicators
module TLabMPI_Transpose
    use TLab_Constants, only: lfile, efile, wp, dp, sp, wi, sizeofreal
    use TLab_Memory, only: imax, jmax, kmax, isize_wrk3d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: TLab_Allocate_Real
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc, c_ptr
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

    integer(wi) :: trp_sizBlock_i, trp_sizBlock_k                   ! explicit sed/recv: group sizes of rend/recv messages
    integer(wi), allocatable :: maps_send_i(:), maps_recv_i(:)      ! PE maps to use explicit sed/recv
    integer(wi), allocatable :: maps_send_k(:), maps_recv_k(:)
    type(MPI_Datatype), allocatable :: types_send(:), types_recv(:) ! alltoallw
    integer, allocatable :: counts(:)

    type(MPI_Datatype) :: trp_datatype_i, trp_datatype_k            ! Transposition in double or single precision

    ! sp work buffer: three sections of size3d sp elements each.
    !   [1..size]       = a_wrk: sp send/pack buffer
    !   [size+1..2*size]= b_wrk: sp recv buffer
    !   [2*size+1..3*size] = c_wrk: flat staging buffer for strided send/recv (sp)
    ! Declared real(sp) allocatable target so OpenMP target can map it directly.
    real(sp), allocatable, target :: wrk_mpi(:)
    real(sp), pointer :: a_wrk(:) => null(), b_wrk(:) => null(), c_wrk(:) => null()

    ! dp/complex staging buffer for flat MPI send/recv (replaces strided MPI_TYPE_VECTOR).
    ! Size = 2*imax*jmax*kmax: covers real(dp) (size3d elements) and complex(dp) (size3d/2 complex = size3d reals).
    real(dp), allocatable, target :: wrk_mpi_dp(:)
    real(dp), pointer    :: c_wrk_dp(:) => null()
    complex(dp), pointer :: c_wrk_cx(:) => null()
    type(MPI_Status) status(128)
    type(MPI_Request) request(128)

#ifdef USE_APU
    ! Shared-memory MPI windows for zero-copy ASYNCHRONOUS transpose on APU.
    ! Each rank allocates one window slot of imax*jmax*kmax elements; peers read directly.
    ! K-direction windows (ims_comm_z), I-direction windows (ims_comm_x).
    type(MPI_Win) :: win_k_sp, win_k_dp   ! K-direction windows
    type(MPI_Win) :: win_i_sp, win_i_dp   ! I-direction windows
    real(sp), pointer :: win_k_sp_own(:) => null()  ! own K sp window slice
    real(dp), pointer :: win_k_dp_own(:) => null()  ! own K dp window slice
    real(sp), pointer :: win_i_sp_own(:) => null()  ! own I sp window slice
    real(dp), pointer :: win_i_dp_own(:) => null()  ! own I dp window slice
    ! Per-peer base pointers from MPI_Win_shared_query (0-based indexing = MPI rank)
    type(c_ptr), allocatable :: peer_base_k_sp(:)   ! size 0:ims_npro_k-1
    type(c_ptr), allocatable :: peer_base_k_dp(:)
    type(c_ptr), allocatable :: peer_base_i_sp(:)   ! size 0:ims_npro_i-1
    type(c_ptr), allocatable :: peer_base_i_dp(:)
#endif

    integer(wi), allocatable :: test(:)
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

        ! -----------------------------------------------------------------------TLabMPI_Trp_Initialize
        integer(wi) ip, npage, dummy
        character(len=32) bakfile, block
        character(len=512) sRes, line
        character*64 lstr
#ifdef USE_APU
        integer(8) :: win_size_bytes
        type(c_ptr) :: own_cptr
        integer(8) :: sz_out
        integer :: du_out
#endif

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
        ! On hawk, we tested that 192 yields optimum performace;
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
        ! wrk_mpi (sp, 3 sections): a_wrk/b_wrk dp<->sp conversion, c_wrk flat staging.
        !   Only needed when ASYNCHRONOUS transport uses single precision (MPI_REAL4).
        if ((trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS .and. trp_datatype_i == MPI_REAL4) .or. &
            (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS .and. trp_datatype_k == MPI_REAL4)) then
            allocate (wrk_mpi(3*imax*jmax*kmax))
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated sp flat-MPI staging buffer (3 x size3d).')
        end if

        ! wrk_mpi_dp (dp): staging buffer for dp-real and complex ASYNCHRONOUS paths.
        !   Not needed for SENDRECV or ALLTOALL which use MPI_TYPE_VECTOR via the kernels.
        if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS .or. trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
            allocate (wrk_mpi_dp(2*imax*jmax*kmax))
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated dp flat-MPI staging buffer (2 x size3d).')
        end if

#ifdef USE_APU
        ! -----------------------------------------------------------------------
        ! Shared-memory MPI windows for ASYNCHRONOUS transpose (APU only).
        ! Window size = imax*jmax*kmax elements per rank (covers all transpose plans).
        ! After allocation, MPI_Win_shared_query retrieves each peer's base pointer.
        win_size_bytes = int(imax, 8) * int(jmax, 8) * int(kmax, 8)

        if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
            if (trp_datatype_k == MPI_REAL4) then
                call MPI_Win_allocate_shared(win_size_bytes * 4_8, 4, MPI_INFO_NULL, &
                                             ims_comm_z, own_cptr, win_k_sp, ims_err)
                call c_f_pointer(own_cptr, win_k_sp_own, [imax*jmax*kmax])
                allocate(peer_base_k_sp(0:ims_npro_k - 1))
                do ip = 0, ims_npro_k - 1
                    call MPI_Win_shared_query(win_k_sp, ip, sz_out, du_out, peer_base_k_sp(ip), ims_err)
                end do
            end if
            call MPI_Win_allocate_shared(win_size_bytes * 8_8, 8, MPI_INFO_NULL, &
                                         ims_comm_z, own_cptr, win_k_dp, ims_err)
            call c_f_pointer(own_cptr, win_k_dp_own, [imax*jmax*kmax])
            allocate(peer_base_k_dp(0:ims_npro_k - 1))
            do ip = 0, ims_npro_k - 1
                call MPI_Win_shared_query(win_k_dp, ip, sz_out, du_out, peer_base_k_dp(ip), ims_err)
            end do
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated K shared MPI windows.')
        end if

        if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
            if (trp_datatype_i == MPI_REAL4) then
                call MPI_Win_allocate_shared(win_size_bytes * 4_8, 4, MPI_INFO_NULL, &
                                             ims_comm_x, own_cptr, win_i_sp, ims_err)
                call c_f_pointer(own_cptr, win_i_sp_own, [imax*jmax*kmax])
                allocate(peer_base_i_sp(0:ims_npro_i - 1))
                do ip = 0, ims_npro_i - 1
                    call MPI_Win_shared_query(win_i_sp, ip, sz_out, du_out, peer_base_i_sp(ip), ims_err)
                end do
            end if
            call MPI_Win_allocate_shared(win_size_bytes * 8_8, 8, MPI_INFO_NULL, &
                                         ims_comm_x, own_cptr, win_i_dp, ims_err)
            call c_f_pointer(own_cptr, win_i_dp_own, [imax*jmax*kmax])
            allocate(peer_base_i_dp(0:ims_npro_i - 1))
            do ip = 0, ims_npro_i - 1
                call MPI_Win_shared_query(win_i_dp, ip, sz_out, du_out, peer_base_i_dp(ip), ims_err)
            end do
            call TLab_Write_ASCII(lfile, 'TLabMPI_Trp_Initialize: allocated I shared MPI windows.')
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
            datatype = trp_datatype_i
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
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, npage, flat_off, disp_ns
        real(sp), pointer :: peer_k_sp_ptr(:)
        real(dp), pointer :: peer_k_dp_ptr(:)

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

        if (trp_datatype_k == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi(1)),             a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi(size + 1)),      b_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi(2*size + 1)),    c_wrk, shape=[size])
            ! dp→sp conversion
#ifdef USE_APU
            !$omp target teams distribute parallel do
#endif
            do i = 1, size
                a_wrk(i) = real(a(i), sp)
            end do
#ifdef USE_APU
            !$omp end target teams distribute parallel do
#endif
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
#ifdef USE_APU
                ! APU path: GPU pack into own shared window, Barrier, GPU read from peer windows.
                ! Replaces ISEND/IRECV/WAITALL with two MPI_Barriers (microseconds).
                ! Step 1: GPU pack a_wrk (strided) → own K sp window slot (flat, slot ns-1 = dest rank)
                !$omp target data use_device_addr(a_wrk, win_k_sp_own)
                do m = 1, ims_npro_k
                    ns = maps_send_k(m) + 1
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_s(ns)
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            win_k_sp_own(flat_off + i*nlines_p + j + 1) = a_wrk(disp_ns + i*npage + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end do
                !$omp end target data
                ! Step 2: memory fence + barrier — all ranks done packing their window slots
                call MPI_Win_sync(win_k_sp, ims_err)
                call MPI_Barrier(ims_comm_z, ims_err)
                ! Step 3: GPU read from each peer's window slot for current rank (ims_pro_k)
                !         into b_wrk at the appropriate displacement
                do m = 1, ims_npro_k
                    nr = maps_recv_k(m) + 1; ipr = nr - 1
                    call c_f_pointer(peer_base_k_sp(ipr), peer_k_sp_ptr, [trp_plan%size3d])
                    flat_off = ims_pro_k * nmax_p * nlines_p   ! slot for current rank in peer's window
                    disp_ns  = trp_plan%disp_r(nr)
                    !$omp target data use_device_addr(peer_k_sp_ptr, b_wrk)
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p*nlines_p
                        b_wrk(disp_ns + i) = peer_k_sp_ptr(flat_off + i)
                    end do
                    !$omp end target teams distribute parallel do
                    !$omp end target data
                end do
                ! Step 4: barrier — all ranks done reading; safe to reuse window on next call
                call MPI_Win_sync(win_k_sp, ims_err)
                call MPI_Barrier(ims_comm_z, ims_err)
#else
                ! CPU path: pipeline IRECVs → GPU pack → ISENDs → WAITALL
                ! Step 1: post all IRECVs — recv buffers are ready before pack starts
                l = 0
                do m = 1, ims_npro_k
                    nr = maps_recv_k(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(b_wrk(trp_plan%disp_r(nr) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ipr, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                ! Step 2: pack a_wrk→c_wrk (strided→flat)
                do m = 1, ims_npro_k
                    ns = maps_send_k(m) + 1
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_s(ns)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            c_wrk(flat_off + i*nlines_p + j + 1) = a_wrk(disp_ns + i*npage + j + 1)
                        end do
                    end do
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
#endif
            else
                call Transpose_Kernel_Single(a_wrk, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                             b_wrk, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                             ims_comm_z, trp_sizBlock_k, trp_mode_k)
            end if
            ! sp→dp conversion
#ifdef USE_APU
            !$omp target teams distribute parallel do
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
                ! dp path
                size = trp_plan%size3d
#ifdef USE_APU
                ! APU: GPU pack into own K dp window, Barrier, GPU read from peer windows.
                ! Step 1: GPU pack a (strided) → own K dp window slot (flat)
                !$omp target data use_device_addr(a, win_k_dp_own)
                do m = 1, ims_npro_k
                    ns = maps_send_k(m) + 1
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_s(ns)
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            win_k_dp_own(flat_off + i*nlines_p + j + 1) = a(disp_ns + i*npage + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end do
                !$omp end target data
                call MPI_Win_sync(win_k_dp, ims_err)
                call MPI_Barrier(ims_comm_z, ims_err)
                ! Step 2: GPU read from each peer's window slot for current rank into b
                do m = 1, ims_npro_k
                    nr = maps_recv_k(m) + 1; ipr = nr - 1
                    call c_f_pointer(peer_base_k_dp(ipr), peer_k_dp_ptr, [size])
                    flat_off = ims_pro_k * nmax_p * nlines_p
                    disp_ns  = trp_plan%disp_r(nr)
                    !$omp target data use_device_addr(peer_k_dp_ptr, b)
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p*nlines_p
                        b(disp_ns + i) = peer_k_dp_ptr(flat_off + i)
                    end do
                    !$omp end target teams distribute parallel do
                    !$omp end target data
                end do
                call MPI_Win_sync(win_k_dp, ims_err)
                call MPI_Barrier(ims_comm_z, ims_err)
#else
                ! CPU path: pipeline IRECVs → pack → ISENDs → WAITALL
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                l = 0
                do m = 1, ims_npro_k
                    nr = maps_recv_k(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(b(trp_plan%disp_r(nr) + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ipr, ims_tag, ims_comm_z, request(l), ims_err)
                end do
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
                do m = 1, ims_npro_k
                    ns = maps_send_k(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ips, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
                nullify (c_wrk_dp)
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
        if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
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
#ifdef USE_APU
            !$omp target data use_device_addr(c_wrk_cx)
#endif
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
#ifdef USE_APU
            !$omp end target data
#endif
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
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, npage, flat_off, disp_nr
        real(sp), pointer :: peer_k_sp_ptr(:)
        real(dp), pointer :: peer_k_dp_ptr(:)
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

        if (trp_datatype_k == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi(1)),          b_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi(size + 1)),   a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi(2*size + 1)), c_wrk, shape=[size])
            ! dp→sp conversion
#ifdef USE_APU
            !$omp target teams distribute parallel do
#endif
            do i = 1, size
                b_wrk(i) = real(b(i), sp)
            end do
#ifdef USE_APU
            !$omp end target teams distribute parallel do
#endif
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! K-Backward: b_wrk (flat K-layout) → communicate → a_wrk (strided original layout)
                ! Send: flat chunks from b_wrk (slot ns = maps_recv_k destination)
                ! Recv: flat chunks for self from peer windows (slot ims_pro_k in peer's window)
#ifdef USE_APU
                ! APU: pack own send data into K sp window, Barrier, GPU read from peer windows.
                ! Step 1: GPU pack flat b_wrk chunks → own window (slot ns-1 = dest rank)
                !$omp target data use_device_addr(b_wrk, win_k_sp_own)
                do m = 1, ims_npro_k
                    ns = maps_recv_k(m) + 1   ! backward msend = maps_recv_k
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_nr  = trp_plan%disp_r(ns)
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p*nlines_p
                        win_k_sp_own(flat_off + i) = b_wrk(disp_nr + i)
                    end do
                    !$omp end target teams distribute parallel do
                end do
                !$omp end target data
                call MPI_Win_sync(win_k_sp, ims_err)
                call MPI_Barrier(ims_comm_z, ims_err)
                ! Step 2: GPU read from each peer's window slot for current rank → c_wrk
                !$omp target data use_device_addr(c_wrk)
                do m = 1, ims_npro_k
                    nr = maps_send_k(m) + 1; ipr = nr - 1   ! backward mrecv = maps_send_k
                    call c_f_pointer(peer_base_k_sp(ipr), peer_k_sp_ptr, [trp_plan%size3d])
                    flat_off = ims_pro_k * nmax_p * nlines_p
                    disp_nr  = (nr - 1)*nmax_p*nlines_p
                    !$omp target data use_device_addr(peer_k_sp_ptr)
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p*nlines_p
                        c_wrk(disp_nr + i) = peer_k_sp_ptr(flat_off + i)
                    end do
                    !$omp end target teams distribute parallel do
                    !$omp end target data
                end do
                !$omp end target data
                call MPI_Win_sync(win_k_sp, ims_err)
                call MPI_Barrier(ims_comm_z, ims_err)
#else
                ! CPU path: IRECVs → ISENDs → WAITALL
                l = 0
                do m = 1, ims_npro_k
                    nr = maps_send_k(m) + 1; ipr = nr - 1   ! backward: mrecv=maps_send_k
                    l = l + 1
                    call MPI_IRECV(c_wrk((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ipr, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                do m = 1, ims_npro_k
                    ns = maps_recv_k(m) + 1; ips = ns - 1   ! backward: msend=maps_recv_k
                    l = l + 1
                    call MPI_ISEND(b_wrk(trp_plan%disp_r(ns) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ips, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
#endif
                ! GPU unpack: flat c_wrk → strided a_wrk (K layout, stride=npage)
                do m = 1, ims_npro_k
                    nr = maps_send_k(m) + 1   ! backward mrecv=maps_send_k
                    flat_off = (nr - 1)*nmax_p*nlines_p
                    disp_nr  = trp_plan%disp_s(nr)
#ifdef USE_APU
                    !$omp target teams distribute parallel do collapse(2)
#endif
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            a_wrk(disp_nr + i*npage + j + 1) = c_wrk(flat_off + i*nlines_p + j + 1)
                        end do
                    end do
#ifdef USE_APU
                    !$omp end target teams distribute parallel do
#endif
                end do
            else
                call Transpose_Kernel_Single(b_wrk, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                             a_wrk, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                             ims_comm_z, trp_sizBlock_k, trp_mode_k)
            end if
            ! sp→dp conversion
#ifdef USE_APU
            !$omp target teams distribute parallel do
#endif
            do i = 1, size
                a(i) = real(a_wrk(i), dp)
            end do
#ifdef USE_APU
            !$omp end target teams distribute parallel do
#endif
            nullify (a_wrk, b_wrk, c_wrk)
        else
            if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! dp path: K-Backward — send flat from b, recv into c_wrk_dp, unpack to a
                size = trp_plan%size3d
#ifdef USE_APU
                ! APU: pack own send data into K dp window, Barrier, GPU read from peer windows,
                !      then GPU unpack into a.
                ! Step 1: GPU pack flat b chunks → own window (slot ns-1 = dest rank)
                !$omp target data use_device_addr(b, win_k_dp_own)
                do m = 1, ims_npro_k
                    ns = maps_recv_k(m) + 1   ! backward msend = maps_recv_k
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_nr  = trp_plan%disp_r(ns)
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p*nlines_p
                        win_k_dp_own(flat_off + i) = b(disp_nr + i)
                    end do
                    !$omp end target teams distribute parallel do
                end do
                !$omp end target data
                call MPI_Win_sync(win_k_dp, ims_err)
                call MPI_Barrier(ims_comm_z, ims_err)
                ! Step 2: GPU read from each peer's window slot for current rank → a (direct unpack)
                do m = 1, ims_npro_k
                    nr = maps_send_k(m) + 1; ipr = nr - 1   ! backward mrecv = maps_send_k
                    call c_f_pointer(peer_base_k_dp(ipr), peer_k_dp_ptr, [size])
                    flat_off = ims_pro_k * nmax_p * nlines_p
                    disp_nr  = trp_plan%disp_s(nr)
                    !$omp target data use_device_addr(peer_k_dp_ptr, a)
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            a(disp_nr + i*npage + j + 1) = peer_k_dp_ptr(flat_off + i*nlines_p + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                    !$omp end target data
                end do
                call MPI_Win_sync(win_k_dp, ims_err)
                call MPI_Barrier(ims_comm_z, ims_err)
#else
                ! CPU path: IRECVs → ISENDs → WAITALL → unpack
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                l = 0
                do m = 1, ims_npro_k
                    nr = maps_send_k(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(c_wrk_dp((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ipr, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                do m = 1, ims_npro_k
                    ns = maps_recv_k(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(b(trp_plan%disp_r(ns) + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ips, ims_tag, ims_comm_z, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
                do m = 1, ims_npro_k
                    nr = maps_send_k(m) + 1
                    flat_off = (nr - 1)*nmax_p*nlines_p
                    disp_nr  = trp_plan%disp_s(nr)
                    do i = 0, nmax_p - 1
                        do j = 0, nlines_p - 1
                            a(disp_nr + i*npage + j + 1) = c_wrk_dp(flat_off + i*nlines_p + j + 1)
                        end do
                    end do
                end do
                nullify (c_wrk_dp)
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
        if (trp_mode_k == TLAB_MPI_TRP_ASYNCHRONOUS) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
#ifdef USE_APU
            !$omp target data use_device_addr(c_wrk_cx)
#endif
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
#ifdef USE_APU
            !$omp end target data
#endif
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
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_nr
        real(sp), pointer :: peer_i_sp_ptr(:)
        real(dp), pointer :: peer_i_dp_ptr(:)

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i   ! total elements per line across all I ranks (stride in recv)

        if (trp_datatype_i == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi(1)),          a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi(size + 1)),   b_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi(2*size + 1)), c_wrk, shape=[size])
            ! dp→sp
#ifdef USE_APU
            !$omp target teams distribute parallel do
#endif
            do i = 1, size
                a_wrk(i) = real(a(i), sp)
            end do
#ifdef USE_APU
            !$omp end target teams distribute parallel do
#endif
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
#ifdef USE_APU
                ! APU: GPU pack a_wrk slots into own I sp window, Barrier,
                !      GPU read from peer windows directly into b_wrk (strided unpack).
                ! Step 1: GPU copy a_wrk flat slots → own I sp window (slot ns-1 = dest rank)
                !$omp target data use_device_addr(a_wrk, win_i_sp_own)
                do m = 1, ims_npro_i
                    ns = maps_send_i(m) + 1
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p*nlines_p
                        win_i_sp_own(flat_off + i) = a_wrk(trp_plan%disp_s(ns) + i)
                    end do
                    !$omp end target teams distribute parallel do
                end do
                !$omp end target data
                call MPI_Win_sync(win_i_sp, ims_err)
                call MPI_Barrier(ims_comm_x, ims_err)
                ! Step 2: GPU read from each peer's window slot for current rank → b_wrk (strided)
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1; ipr = nr - 1
                    call c_f_pointer(peer_base_i_sp(ipr), peer_i_sp_ptr, [trp_plan%size3d])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    disp_nr  = trp_plan%disp_r(nr)
                    !$omp target data use_device_addr(peer_i_sp_ptr, b_wrk)
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            b_wrk(disp_nr + i*nmax_full + j + 1) = peer_i_sp_ptr(flat_off + i*nmax_p + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                    !$omp end target data
                end do
                call MPI_Win_sync(win_i_sp, ims_err)
                call MPI_Barrier(ims_comm_x, ims_err)
#else
                ! CPU path: IRECVs → ISENDs → WAITALL → unpack
                l = 0
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(c_wrk((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ipr, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                do m = 1, ims_npro_i
                    ns = maps_send_i(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(a_wrk(trp_plan%disp_s(ns) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ips, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1
                    flat_off = (nr - 1)*nmax_p*nlines_p
                    disp_nr  = trp_plan%disp_r(nr)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            b_wrk(disp_nr + i*nmax_full + j + 1) = c_wrk(flat_off + i*nmax_p + j + 1)
                        end do
                    end do
                end do
#endif
            else
                call Transpose_Kernel_Single(a_wrk, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                             b_wrk, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                             ims_comm_x, trp_sizBlock_i, trp_mode_i)
            end if
            ! sp→dp
#ifdef USE_APU
            !$omp target teams distribute parallel do
#endif
            do i = 1, size
                b(i) = real(b_wrk(i), dp)
            end do
#ifdef USE_APU
            !$omp end target teams distribute parallel do
#endif
            nullify (a_wrk, b_wrk, c_wrk)
        else
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! dp path
                size = trp_plan%size3d
#ifdef USE_APU
                ! APU: GPU pack a slots → own I dp window, Barrier, GPU read peers → b (strided).
                ! Step 1: GPU copy a flat slots → own I dp window (slot ns-1 = dest rank)
                !$omp target data use_device_addr(a, win_i_dp_own)
                do m = 1, ims_npro_i
                    ns = maps_send_i(m) + 1
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p*nlines_p
                        win_i_dp_own(flat_off + i) = a(trp_plan%disp_s(ns) + i)
                    end do
                    !$omp end target teams distribute parallel do
                end do
                !$omp end target data
                call MPI_Win_sync(win_i_dp, ims_err)
                call MPI_Barrier(ims_comm_x, ims_err)
                ! Step 2: GPU read from each peer's window slot for current rank → b (strided)
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1; ipr = nr - 1
                    call c_f_pointer(peer_base_i_dp(ipr), peer_i_dp_ptr, [size])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    disp_nr  = trp_plan%disp_r(nr)
                    !$omp target data use_device_addr(peer_i_dp_ptr, b)
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            b(disp_nr + i*nmax_full + j + 1) = peer_i_dp_ptr(flat_off + i*nmax_p + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                    !$omp end target data
                end do
                call MPI_Win_sync(win_i_dp, ims_err)
                call MPI_Barrier(ims_comm_x, ims_err)
#else
                ! CPU path: IRECVs → ISENDs → WAITALL → unpack
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                l = 0
                do m = 1, ims_npro_i
                    nr = maps_recv_i(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(c_wrk_dp((nr-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                do m = 1, ims_npro_i
                    ns = maps_send_i(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(a(trp_plan%disp_s(ns) + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
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
        if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_cx, shape=[size])
#ifdef USE_APU
            !$omp target data use_device_addr(c_wrk_cx)
#endif
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
#ifdef USE_APU
            !$omp end target data
#endif
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
        integer(wi) size, i, j, l, m, ns, nr, ips, ipr, nmax_p, nlines_p, nmax_full, flat_off, disp_ns
        real(sp), pointer :: peer_i_sp_ptr(:)
        real(dp), pointer :: peer_i_dp_ptr(:)

        ! #######################################################################
        nmax_p    = trp_plan%nmax
        nlines_p  = trp_plan%nlines
        nmax_full = nmax_p * ims_npro_i

        if (trp_datatype_i == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(wrk_mpi(1)),          b_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi(size + 1)),   a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi(2*size + 1)), c_wrk, shape=[size])
            ! dp→sp
#ifdef USE_APU
            !$omp target teams distribute parallel do
#endif
            do i = 1, size
                b_wrk(i) = real(b(i), sp)
            end do
#ifdef USE_APU
            !$omp end target teams distribute parallel do
#endif
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
#ifdef USE_APU
                ! APU: GPU pack b_wrk (strided) → own I sp window (flat), Barrier,
                !      GPU read from peer windows → a_wrk (flat).
                ! Step 1: GPU pack b_wrk (strided) → own I sp window (slot ns-1 = dest rank)
                !$omp target data use_device_addr(b_wrk, win_i_sp_own)
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1   ! backward msend = maps_recv_i
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_r(ns)
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            win_i_sp_own(flat_off + i*nmax_p + j + 1) = b_wrk(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end do
                !$omp end target data
                call MPI_Win_sync(win_i_sp, ims_err)
                call MPI_Barrier(ims_comm_x, ims_err)
                ! Step 2: GPU read from each peer's window slot for current rank → a_wrk (flat)
                do m = 1, ims_npro_i
                    nr = maps_send_i(m) + 1; ipr = nr - 1   ! backward mrecv = maps_send_i
                    call c_f_pointer(peer_base_i_sp(ipr), peer_i_sp_ptr, [trp_plan%size3d])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    disp_ns  = trp_plan%disp_s(nr)
                    !$omp target data use_device_addr(peer_i_sp_ptr, a_wrk)
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p*nlines_p
                        a_wrk(disp_ns + i) = peer_i_sp_ptr(flat_off + i)
                    end do
                    !$omp end target teams distribute parallel do
                    !$omp end target data
                end do
                call MPI_Win_sync(win_i_sp, ims_err)
                call MPI_Barrier(ims_comm_x, ims_err)
#else
                ! CPU path: IRECVs → pack → ISENDs → WAITALL
                l = 0
                do m = 1, ims_npro_i
                    nr = maps_send_i(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(a_wrk(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ipr, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1   ! backward msend=maps_recv_i
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_r(ns)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            c_wrk(flat_off + i*nmax_p + j + 1) = b_wrk(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                end do
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(c_wrk((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, MPI_REAL4, &
                                   ips, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
#endif
            else
                call Transpose_Kernel_Single(b_wrk, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                             a_wrk, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                             ims_comm_x, trp_sizBlock_i, trp_mode_i)
            end if
            ! sp→dp
#ifdef USE_APU
            !$omp target teams distribute parallel do
#endif
            do i = 1, size
                a(i) = real(a_wrk(i), dp)
            end do
#ifdef USE_APU
            !$omp end target teams distribute parallel do
#endif
            nullify (a_wrk, b_wrk, c_wrk)
        else
            if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
                ! dp path
                size = trp_plan%size3d
#ifdef USE_APU
                ! APU: GPU pack b (strided) → own I dp window (flat), Barrier,
                !      GPU read from peer windows → a (flat).
                ! Step 1: GPU pack b (strided) → own I dp window (slot ns-1 = dest rank)
                !$omp target data use_device_addr(b, win_i_dp_own)
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1   ! backward msend = maps_recv_i
                    flat_off = (ns - 1)*nmax_p*nlines_p
                    disp_ns  = trp_plan%disp_r(ns)
                    !$omp target teams distribute parallel do collapse(2)
                    do i = 0, nlines_p - 1
                        do j = 0, nmax_p - 1
                            win_i_dp_own(flat_off + i*nmax_p + j + 1) = b(disp_ns + i*nmax_full + j + 1)
                        end do
                    end do
                    !$omp end target teams distribute parallel do
                end do
                !$omp end target data
                call MPI_Win_sync(win_i_dp, ims_err)
                call MPI_Barrier(ims_comm_x, ims_err)
                ! Step 2: GPU read from each peer's window slot for current rank → a (flat)
                do m = 1, ims_npro_i
                    nr = maps_send_i(m) + 1; ipr = nr - 1   ! backward mrecv = maps_send_i
                    call c_f_pointer(peer_base_i_dp(ipr), peer_i_dp_ptr, [size])
                    flat_off = ims_pro_i * nmax_p * nlines_p
                    disp_ns  = trp_plan%disp_s(nr)
                    !$omp target data use_device_addr(peer_i_dp_ptr, a)
                    !$omp target teams distribute parallel do
                    do i = 1, nmax_p*nlines_p
                        a(disp_ns + i) = peer_i_dp_ptr(flat_off + i)
                    end do
                    !$omp end target teams distribute parallel do
                    !$omp end target data
                end do
                call MPI_Win_sync(win_i_dp, ims_err)
                call MPI_Barrier(ims_comm_x, ims_err)
#else
                ! CPU path: IRECVs → pack → ISENDs → WAITALL
                call c_f_pointer(c_loc(wrk_mpi_dp(1)), c_wrk_dp, shape=[size])
                l = 0
                do m = 1, ims_npro_i
                    nr = maps_send_i(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_IRECV(a(trp_plan%disp_s(nr) + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ipr, ims_tag, ims_comm_x, request(l), ims_err)
                end do
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
                do m = 1, ims_npro_i
                    ns = maps_recv_i(m) + 1; ips = ns - 1
                    l = l + 1
                    call MPI_ISEND(c_wrk_dp((ns-1)*nmax_p*nlines_p + 1), nmax_p*nlines_p, &
                                   trp_plan%base_type, ips, ims_tag, ims_comm_x, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
                nullify (c_wrk_dp)
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
        if (trp_mode_i == TLAB_MPI_TRP_ASYNCHRONOUS) then
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
#ifdef USE_APU
            !$omp target data use_device_addr(c_wrk_cx)
#endif
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
#ifdef USE_APU
            !$omp end target data
#endif
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
        type(MPI_Datatype), intent(in) :: tsend, trecv                 ! types send/receive
        integer(wi), intent(in) :: dsend(:), drecv(:)       ! displacements send/receive
        integer(wi), intent(in) :: msend(:), mrecv(:)       ! maps send/receive
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
        type(MPI_Datatype), intent(in) :: tsend, trecv                 ! types send/receive
        integer(wi), intent(in) :: dsend(:), drecv(:)       ! displacements send/receive
        integer(wi), intent(in) :: msend(:), mrecv(:)       ! maps send/receive
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
            ! call MPI_ALLTOALLW(a, spread(1, 1, npro), dsend*int(sizeof(1.0_wp)), spread(tsend, 1, npro), &
            !                    b, spread(1, 1, npro), drecv*int(sizeof(1.0_wp)), spread(trecv, 1, npro), comm, ims_err)

        end select

        return
    end subroutine Transpose_Kernel_Single

    !########################################################################
    !########################################################################
    subroutine Transpose_Kernel_Complex(a, msend, dsend, tsend, b, mrecv, drecv, trecv, comm, step, mode)
        complex(wp), intent(in) :: a(*)
        complex(wp), intent(out) :: b(*)

        type(MPI_Comm), intent(in) :: comm                         ! communicator
        type(MPI_Datatype), intent(in) :: tsend, trecv                 ! types send/receive
        integer(wi), intent(in) :: dsend(:), drecv(:)       ! displacements send/receive
        integer(wi), intent(in) :: msend(:), mrecv(:)       ! maps send/receive
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
            ! call MPI_ALLTOALLW(a, spread(1, 1, npro), dsend*int(sizeof(1.0_wp)), spread(tsend, 1, npro), &
            !                    b, spread(1, 1, npro), drecv*int(sizeof(1.0_wp)), spread(trecv, 1, npro), comm, ims_err)

        end select

        return
    end subroutine Transpose_Kernel_Complex

end module TLabMPI_Transpose