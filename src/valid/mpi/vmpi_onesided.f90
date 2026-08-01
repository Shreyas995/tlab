! vmpi_onesided.f90
! ============================================================================
! Standalone harness to RE-ADJUDICATE the project's oldest un-rooted verdict:
!
!     "One-sided MPI (MPI_Put / MPI_Win_create) does NOT deliver inter-node on
!      Cray MPICH (returns MPI_SUCCESS, data lands nowhere)."  [claimed 2026-05-12/13]
!
! That entry (CLAUDE.md "Hunter-specific MPI/OS behaviors" + CLAUDE.log rejected list)
! is a ONE-LINE symptom with NO harness and NO diagnostic trail behind it -- unlike every
! other transpose finding, each of which has a standalone src/valid/mpi/vmpi_* test
! (commtest_p1-p4, nodewin, gpuaware, xcd_write, hip_shmwrite). Two *completely different*
! mechanisms produce that identical symptom, and they were never told apart:
!
!   (a) USAGE BUG, not an MPI defect.  MPI_Put is NOT required to have updated the target
!       when the call returns -- RMA lets an implementation buffer/defer the transfer until
!       the ACCESS EPOCH CLOSES. The only things that prove delivery are MPI_Win_fence
!       (active target), MPI_Win_flush / MPI_Win_unlock (passive target), or
!       MPI_Win_complete + MPI_Win_wait (PSCW). If the 2026-05 test read the target buffer
!       before a correct close -- or fenced the wrong window -- then "MPI_SUCCESS, data
!       nowhere" is EXACTLY what a perfectly conforming MPI looks like.
!       => one-sided is viable, and the rejected-approaches entry is wrong.
!
!   (b) GENUINE Cray MPICH GPU-RMA GAP.  Device-memory RMA (Put/Get targeting a window
!       backed by GPU memory) is a far less mature path than device-memory TWO-SIDED
!       Isend/Irecv (which IS proven working here -- see vmpi_gpuaware, MPICH_GPU_SUPPORT_ENABLED).
!       If the library mis-detects the target window's memory type it can fail to program a
!       real RDMA transfer while still returning MPI_SUCCESS locally (nothing was wrong on
!       THIS rank); the failure is invisible until you inspect the target buffer remotely.
!       => one-sided is dead for GPU buffers; next scope is ROC_SHMEM / libfabric.
!
!   (c) [third candidate this harness also covers] SHARED-WINDOW TAINT.  This project has a
!       PROVEN Cray behavior: allocating any MPI_Win_allocate_shared retroactively breaks
!       two-sided traffic on ims_comm_x/ims_comm_z and any dup of them (vmpi_commtest_p1
!       HANGS, p3/p4 PASS). The 2026-05 one-sided test may have run inside the production
!       binary, where node windows exist. If Put works standalone but fails with a shared
!       window alive, the verdict was real but MISATTRIBUTED.
!
! ---------------------------------------------------------------------------
! DECISION TABLE  (this is the whole point of the harness -- read the [VERDICT] lines)
! ---------------------------------------------------------------------------
!   A1 FAIL  + A2..A5 PASS                 -> (a) usage bug. The old verdict is WRONG.
!                                             One-sided becomes a legitimate lower-overhead
!                                             alternative to the two-sided inter-node leg.
!   A2..A5 PASS + B1/B2 PASS + B3/B4 FAIL  -> (b) real GPU-RMA gap, confined to *device*
!                                             allocations. USM buffers still usable.
!   A* PASS + B* FAIL                      -> (b) full GPU-RMA gap. Verdict stands; go
!                                             ROC_SHMEM / libfabric or stay two-sided.
!   A2 PASS but C1 FAIL                    -> (c) shared-window taint, NOT an RMA defect.
!   A2 FAIL                                -> the harness itself is suspect (host-memory Put
!                                             with a correct fence must work). Do NOT ship a
!                                             verdict from this run; debug the harness first.
!
! Every case is independently timeout-guarded by run_onesided.pbs, so a HANG in one mode
! cannot hide the results of the others (a hang is itself a datum -- exit=124).
!
! ---------------------------------------------------------------------------
! WHAT A "PASS" MEANS HERE -- the verification is remote, not local
! ---------------------------------------------------------------------------
! Origin rank with K-coordinate pk_me Puts `mas` doubles into EVERY K-peer's window at
! target displacement pk_me*mas. Encoded value:  v(j) = world_rank_of_origin*ENC + j.
! After the epoch close, each rank reads its OWN window and checks slot m against
! peer_world(m)*ENC + j. Three distinguishable outcomes per element -- this is what makes
! the harness diagnostic rather than pass/fail:
!     n_untouched : still the FILL sentinel        -> nothing was ever written (the claimed symptom)
!     n_wrong     : a valid encoding, WRONG slot   -> delivered but mis-placed (cf. the known
!                   "shared-window memory is not a valid IRECV target" pathology)
!     garbage     : decodes to no valid rank       -> genuine corruption
! The first mismatch is decoded and logged (src=, want=, got=) per rank to fort.<400+rank>.
!
! ---------------------------------------------------------------------------
! GOAL CONTEXT (why we care -- speed, at equal accuracy)
! ---------------------------------------------------------------------------
! The inter-node K leg of the transpose is two-sided Isend/Irecv on fabric_mpi_comm_k.
! If (a) holds, one-sided Put removes tag matching, receive posting, and the unexpected-message
! path. HONEST EXPECTATION, set in advance so we do not fool ourselves: vmpi_gpuaware's M5
! sweep already showed this leg is FABRIC-BANDWIDTH-BOUND above ~16 MB/peer. One-sided will
! NOT raise the bandwidth ceiling. The only wins available are latency/overhead and freed CPU
! cores. Mode 20-23 measure exactly that against the production two-sided baseline, at the
! production message size, so the GO/NO-GO is a measurement and not a hope.
!
! ---------------------------------------------------------------------------
! MODES (argv1)
! ---------------------------------------------------------------------------
!   A-series  host memory, isolate the SYNCHRONIZATION axis  -> discriminates (a)
!     1  A1  Win_create,   NO closing sync (read right after Put)  [CONTROL: expected FAIL]
!     2  A2  Win_create,   MPI_Win_fence
!     3  A3  Win_create,   Win_lock/Win_unlock        (passive target)
!     4  A4  Win_create,   Win_lock_all/Win_flush_all (passive target)
!     5  A5  Win_create,   PSCW post/start/complete/wait
!     6  A6  Win_allocate, MPI_Win_fence              (library-allocated == pre-registered?)
!   B-series  same correct sync, escalate the MEMORY axis     -> discriminates (b)
!     7  B1  USM buffer, GPU-written origin, Win_create,        fence
!     8  B2  USM buffer, GPU-written origin, Win_allocate,      fence
!     9  B3  device (omp_target_alloc) origin AND window, Win_create,  fence
!    10  B4  device origin AND window, Win_create_dynamic/attach,      lock_all+flush_all
!   C-series  with an MPI_Win_allocate_shared alive           -> discriminates (c)
!    11  C1  = A2 with taint window
!    12  C2  = B1 with taint window
!   PERF
!    20  P1  two-sided Isend/Irecv baseline (production shape), inter-node peers
!    21  P2  one-sided Put + fence
!    22  P3  one-sided Put + lock_all/flush_all
!    23  P-SWEEP: P1/P2/P3 at 1, 4, 16, 64 MB per peer
!     0  run the full correctness matrix 1..12 in one job (convenience; PBS prefers separate)
!
! ---------------------------------------------------------------------------
! HUNTER BUILD/RUN GOTCHAS OBEYED HERE (hard-won, see CLAUDE.md)
! ---------------------------------------------------------------------------
!   - `use mpi`, NOT `use mpi_f08`: under -fopenmp offload Cray miscompiles mpi_f08's
!     predefined type(MPI_Datatype) constants -> "Invalid datatype".
!   - NO `module load craype-accel-amd-gfx942` in the PBS script: it reloads cce and breaks
!     MPI predefined datatypes. Inherit the .bashrc HLRS/APU env.
!   - NO block-local array + array-constructor under -fopenmp + unified_shared_memory
!     (Cray zeroed them at runtime). All sweep sizes are passed as literal arguments.
!   - Test comm is MPI_Comm_split(MPI_COMM_WORLD, pro_i, pro_k) -- the SAME clean split
!     production uses for fabric_mpi_comm_k. Never a Cartesian sub-comm or a dup of one:
!     that would inject the known taint and we would misread it as an RMA failure.
!   - Per-rank output is a FLUSHED unit (fort.<400+rank>); write(*,*) is buffered and lost
!     if a rank dies. Read the LAST line of a dead rank's file.
!
!   NOTE on B4 / MPI_Get_address: exchanging MPI_Aint addresses via Allgather for a DYNAMIC
!   window is the standard, CORRECT mechanism (MPI translates them at the target). Do not
!   confuse it with the project's rejected "Allgather of c_ptr for cross-process pointers"
!   -- that failed because it used the raw addresses for direct load/store, which are
!   process-private. Here MPI, not our code, dereferences them.
!
! Build:  make -f Makefile.onesided        Run: qsub run_onesided.pbs
! ============================================================================
module onesided_mod
    use mpi
    use iso_c_binding
#ifdef USE_APU
    use omp_lib
#endif
    implicit none
#ifdef USE_APU
    !$omp requires unified_shared_memory
#endif

    integer, parameter :: dp = kind(1.0d0)

    ! Production decomposition (1152x432x1152, npro_i=6 x npro_k=8 = 48 ranks, 2 nodes x 24).
    integer, parameter :: NPRO_I = 6, NPRO_K = 8, RANKS_PER_NODE = 24
    integer, parameter :: FUNIT_BASE = 400          ! fort.400 .. fort.447

    real(dp), parameter :: ENC  = 1.0d9             ! v = world*ENC + idx  (exact in dp: < 2^53)
    real(dp), parameter :: FILL = -999.0_dp         ! "nothing arrived" sentinel
    real(dp), parameter :: TOL  = 0.25_dp           ! all encodings are integers; 0.25 is generous

    ! --- window construction ---
    integer, parameter :: WIN_CREATE  = 1           ! MPI_Win_create over memory WE allocated
    integer, parameter :: WIN_ALLOC   = 2           ! MPI_Win_allocate (library-allocated/registered)
    integer, parameter :: WIN_DYNAMIC = 3           ! MPI_Win_create_dynamic + MPI_Win_attach

    ! --- epoch synchronization ---
    integer, parameter :: SYNC_NONE  = 0            ! open fence, Put, read WITHOUT closing (the control)
    integer, parameter :: SYNC_FENCE = 1
    integer, parameter :: SYNC_LOCK  = 2            ! per-target lock/unlock + barrier
    integer, parameter :: SYNC_FLUSH = 3            ! lock_all / flush_all / barrier
    integer, parameter :: SYNC_PSCW  = 4

    ! --- memory type of origin buffer AND window ---
    integer, parameter :: MEM_HOST = 1              ! plain allocate, CPU-written  (host-resident)
    integer, parameter :: MEM_USM  = 2              ! plain allocate, GPU-written  (unified, GPU-dirty)
    integer, parameter :: MEM_DEV  = 3              ! omp_target_alloc             (true device pointer)

    ! --- perf methods ---
    integer, parameter :: PERF_TWOSIDED = 1, PERF_PUT_FENCE = 2, PERF_PUT_FLUSH = 3

    integer :: ims_rank, ims_nprocs, ims_pro_i, ims_pro_k, my_node, funit
    integer :: fabric_comm_k

    ! shared-window taint source (hypothesis (c)); kept alive across a case when requested
    integer :: taint_win = MPI_WIN_NULL, taint_comm = MPI_COMM_NULL
    logical :: taint_live = .false.

#ifdef USE_APU
    interface
        function hipDeviceSynchronize() bind(C, name='hipDeviceSynchronize') result(r)
            use iso_c_binding
            integer(c_int) :: r
        end function
    end interface
#endif

contains

    !-----------------------------------------------------------------------
    subroutine setup_topology()
        ! Cartesian ONLY to reproduce production's rank->(pro_k,pro_i) mapping
        ! (world = pro_k*NPRO_I + pro_i). All test traffic then runs on a clean
        ! MPI_COMM_WORLD split -- never on the Cartesian comm or a dup of it.
        integer :: comm_xz, dims(2), coord(2), ierr
        logical :: period(2)
        dims(1) = NPRO_K; dims(2) = NPRO_I; period = .true.
        call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, .false., comm_xz, ierr)
        call MPI_Cart_coords(comm_xz, ims_rank, 2, coord, ierr)
        ims_pro_k = coord(1); ims_pro_i = coord(2)
        call MPI_Comm_split(MPI_COMM_WORLD, ims_pro_i, ims_pro_k, fabric_comm_k, ierr)
        my_node = ims_rank/RANKS_PER_NODE
        call MPI_Comm_free(comm_xz, ierr)
    end subroutine setup_topology

    integer function peer_world(pk)             ! local K-rank -> world rank
        integer, intent(in) :: pk
        peer_world = pk*NPRO_I + ims_pro_i
    end function peer_world

    integer function peer_node(pk)
        integer, intent(in) :: pk
        peer_node = (pk*NPRO_I + ims_pro_i)/RANKS_PER_NODE
    end function peer_node

    !-----------------------------------------------------------------------
    subroutine report_env()
        ! Prove the !$omp target regions really run on the GPU, so a B-series PASS
        ! cannot be a false positive from a silent host fallback.
        logical :: on_host
        on_host = .true.
#ifdef USE_APU
        !$omp target map(tofrom: on_host)
        on_host = omp_is_initial_device()
        !$omp end target
        if (ims_rank == 0) &
            write (*, '(a,i0,a,l1,a)') '[ENV] omp_get_num_devices=', omp_get_num_devices(), &
            '   target_on_host=', on_host, '   (target_on_host=F means real GPU offload)'
#else
        if (ims_rank == 0) write (*, '(a)') '[ENV] built WITHOUT USE_APU: B-series (GPU memory) will be SKIPPED'
#endif
    end subroutine report_env

    !-----------------------------------------------------------------------
    subroutine taint_on()
        ! Reproduce the production environment: a live MPI_Win_allocate_shared.
        ! In production these are node_win_i/node_win_k. Proven to retroactively break
        ! two-sided traffic on Cartesian comms; whether it also perturbs RMA on a clean
        ! split comm is exactly what C1/C2 ask.
        integer :: ierr
        type(c_ptr) :: bp
        integer(MPI_ADDRESS_KIND) :: nb
        if (taint_live) return
        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, ims_rank, &
                                 MPI_INFO_NULL, taint_comm, ierr)
        nb = 4096_MPI_ADDRESS_KIND*8_MPI_ADDRESS_KIND
        call MPI_Win_allocate_shared(nb, 8, MPI_INFO_NULL, taint_comm, bp, taint_win, ierr)
        taint_live = .true.
        write (funit, '(a,i0)') '[TAINT] MPI_Win_allocate_shared active, ierr=', ierr
        flush (funit)
    end subroutine taint_on

    subroutine taint_off()
        integer :: ierr
        if (.not. taint_live) return
        call MPI_Win_free(taint_win, ierr)
        call MPI_Comm_free(taint_comm, ierr)
        taint_live = .false.
    end subroutine taint_off

    !=======================================================================
    ! CORRECTNESS CASE
    !=======================================================================
    subroutine run_case(label, win_kind, sync_kind, mem_kind, do_taint, nmax_p, nlines_p)
        character(*), intent(in) :: label
        integer, intent(in) :: win_kind, sync_kind, mem_kind, nmax_p, nlines_p
        logical, intent(in) :: do_taint

        integer :: mas, wsize, pk, j, ierr, win, grp, disp_unit
        integer :: put_ierr, put_ierr_max, first_src, g_first_src
        integer :: n_untouched, n_wrong, n_inter_bad, n_intra_bad
        integer :: g_untouched, g_wrong, g_inter_bad, g_intra_bad, g_put_ierr
        integer(MPI_ADDRESS_KIND) :: nbytes
        integer(MPI_ADDRESS_KIND), allocatable :: peer_base(:)
        integer(MPI_ADDRESS_KIND) :: my_base_aint
        real(dp), pointer :: sbuf(:) => null(), wbuf(:) => null()
        real(dp), allocatable, target :: host_stage(:)
        type(c_ptr) :: sptr, wptr
        real(dp) :: first_want, first_got

        mas = nmax_p*nlines_p                 ! doubles pushed to ONE peer
        wsize = NPRO_K*mas                    ! window = one slot per K-peer
        nbytes = int(wsize, MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
        disp_unit = 8
        sptr = c_null_ptr; wptr = c_null_ptr
        win = MPI_WIN_NULL

#ifndef USE_APU
        if (mem_kind /= MEM_HOST) then
            if (ims_rank == 0) write (*, '(3a)') '[SKIP] ', label, ' -- needs USE_APU build'
            return
        end if
#endif

        if (do_taint) call taint_on()

        write (funit, '(a)') '================================================================'
        write (funit, '(3a,i0,a,i0,a,i0)') '[CASE] ', label, '  win=', win_kind, ' sync=', sync_kind, ' mem=', mem_kind
        write (funit, '(a,i0,a,i0,a,i0,a,i0)') '[CASE] pro_k=', ims_pro_k, ' pro_i=', ims_pro_i, &
            ' node=', my_node, ' mas=', mas
        flush (funit)

        !--- allocate origin buffer + window memory ------------------------
        call alloc_buf(sptr, sbuf, wsize, mem_kind)
        if (win_kind /= WIN_ALLOC) call alloc_buf(wptr, wbuf, wsize, mem_kind)

        !--- fill origin buffer: every target slot carries OUR encoding ----
        allocate (host_stage(wsize))
        do pk = 0, NPRO_K - 1
            do j = 1, mas
                host_stage(pk*mas + j) = real(ims_rank, dp)*ENC + real(j, dp)
            end do
        end do
        call put_host_to_buf(host_stage, sptr, sbuf, wsize, mem_kind)

        !--- prime the window with the "nothing arrived" sentinel ----------
        host_stage(:) = FILL

        !--- create the window --------------------------------------------
        select case (win_kind)
        case (WIN_CREATE)
            call put_host_to_buf(host_stage, wptr, wbuf, wsize, mem_kind)
            call MPI_Win_create(wbuf, nbytes, disp_unit, MPI_INFO_NULL, fabric_comm_k, win, ierr)
            write (funit, '(a,i0)') '[WIN] MPI_Win_create ierr=', ierr
        case (WIN_ALLOC)
            call MPI_Win_allocate(nbytes, disp_unit, MPI_INFO_NULL, fabric_comm_k, wptr, win, ierr)
            write (funit, '(a,i0)') '[WIN] MPI_Win_allocate ierr=', ierr
            call c_f_pointer(wptr, wbuf, [wsize])
            wbuf(:) = FILL                     ! library-allocated memory is always host-addressable
        case (WIN_DYNAMIC)
            call put_host_to_buf(host_stage, wptr, wbuf, wsize, mem_kind)
            call MPI_Win_create_dynamic(MPI_INFO_NULL, fabric_comm_k, win, ierr)
            write (funit, '(a,i0)') '[WIN] MPI_Win_create_dynamic ierr=', ierr
            call MPI_Win_attach(win, wbuf, nbytes, ierr)
            write (funit, '(a,i0)') '[WIN] MPI_Win_attach ierr=', ierr
            ! Standard dynamic-window address exchange (MPI dereferences these, we never do).
            allocate (peer_base(0:NPRO_K - 1))
            call MPI_Get_address(wbuf, my_base_aint, ierr)
            call MPI_Allgather(my_base_aint, 1, MPI_AINT, peer_base, 1, MPI_AINT, fabric_comm_k, ierr)
        end select
        flush (funit)

        call MPI_Barrier(fabric_comm_k, ierr)

        !--- THE EXPERIMENT: one epoch, NPRO_K Puts, then close ------------
        put_ierr_max = 0
        select case (sync_kind)

        case (SYNC_NONE)
            ! CONTROL. Open an epoch, Put, and read the target WITHOUT closing it.
            ! A conforming MPI is free to have delivered nothing yet. If this is the only
            ! failing case, the 2026-05 verdict was hypothesis (a) -- a usage bug.
            call MPI_Win_fence(0, win, ierr)
            call do_puts(win, sbuf, mas, WIN_CREATE, peer_base, put_ierr_max)
            call MPI_Barrier(fabric_comm_k, ierr)     ! processes arrive; RMA is NOT completed
            call verify(wptr, wbuf, mas, mem_kind, n_untouched, n_wrong, n_inter_bad, n_intra_bad, &
                        first_src, first_want, first_got)
            call MPI_Win_fence(0, win, ierr)          ! close properly so the window can be freed

        case (SYNC_FENCE)
            call MPI_Win_fence(0, win, ierr)
            call do_puts(win, sbuf, mas, win_kind, peer_base, put_ierr_max)
            call MPI_Win_fence(0, win, ierr)          ! <-- the completion guarantee
            call verify(wptr, wbuf, mas, mem_kind, n_untouched, n_wrong, n_inter_bad, n_intra_bad, &
                        first_src, first_want, first_got)

        case (SYNC_LOCK)
            do pk = 0, NPRO_K - 1
                call MPI_Win_lock(MPI_LOCK_SHARED, pk, 0, win, ierr)
                call one_put(win, sbuf, pk, mas, win_kind, peer_base, put_ierr)
                put_ierr_max = max(put_ierr_max, put_ierr)
                call MPI_Win_unlock(pk, win, ierr)    ! <-- the completion guarantee
            end do
            call MPI_Barrier(fabric_comm_k, ierr)     ! targets must know all origins finished
            call verify(wptr, wbuf, mas, mem_kind, n_untouched, n_wrong, n_inter_bad, n_intra_bad, &
                        first_src, first_want, first_got)

        case (SYNC_FLUSH)
            call MPI_Win_lock_all(0, win, ierr)
            call do_puts(win, sbuf, mas, win_kind, peer_base, put_ierr_max)
            call MPI_Win_flush_all(win, ierr)         ! <-- the completion guarantee
            call MPI_Barrier(fabric_comm_k, ierr)
            call verify(wptr, wbuf, mas, mem_kind, n_untouched, n_wrong, n_inter_bad, n_intra_bad, &
                        first_src, first_want, first_got)
            call MPI_Win_unlock_all(win, ierr)

        case (SYNC_PSCW)
            call MPI_Comm_group(fabric_comm_k, grp, ierr)
            call MPI_Win_post(grp, 0, win, ierr)
            call MPI_Win_start(grp, 0, win, ierr)
            call do_puts(win, sbuf, mas, win_kind, peer_base, put_ierr_max)
            call MPI_Win_complete(win, ierr)          ! origin side done
            call MPI_Win_wait(win, ierr)              ! <-- target side completion guarantee
            call verify(wptr, wbuf, mas, mem_kind, n_untouched, n_wrong, n_inter_bad, n_intra_bad, &
                        first_src, first_want, first_got)
            call MPI_Group_free(grp, ierr)
        end select

        !--- reduce + verdict ---------------------------------------------
        call MPI_Reduce(n_untouched, g_untouched, 1, MPI_INTEGER, MPI_SUM, 0, fabric_comm_k, ierr)
        call MPI_Reduce(n_wrong, g_wrong, 1, MPI_INTEGER, MPI_SUM, 0, fabric_comm_k, ierr)
        call MPI_Reduce(n_inter_bad, g_inter_bad, 1, MPI_INTEGER, MPI_SUM, 0, fabric_comm_k, ierr)
        call MPI_Reduce(n_intra_bad, g_intra_bad, 1, MPI_INTEGER, MPI_SUM, 0, fabric_comm_k, ierr)
        call MPI_Reduce(put_ierr_max, g_put_ierr, 1, MPI_INTEGER, MPI_MAX, 0, fabric_comm_k, ierr)
        call MPI_Reduce(first_src, g_first_src, 1, MPI_INTEGER, MPI_MAX, 0, fabric_comm_k, ierr)

        write (funit, '(a,i0,a,i0,a,i0,a,i0,a,i0)') '[LOCAL] untouched=', n_untouched, &
            ' wrong=', n_wrong, ' inter_bad=', n_inter_bad, ' intra_bad=', n_intra_bad, &
            ' put_ierr=', put_ierr_max
        if (n_untouched + n_wrong > 0) &
            write (funit, '(a,i0,a,es24.16,a,es24.16)') '[FIRSTBAD] decoded_src=', first_src, &
            ' want=', first_want, ' got=', first_got
        flush (funit)

        if (ims_pro_k == 0) then
            write (*, '(a)') '----------------------------------------------------------------'
            write (*, '(2a)') '[CASE] ', label
            write (*, '(a,i0,a,i0,a,i0,a,i0,a,i0)') '   put_ierr(max)=', g_put_ierr, &
                '   untouched=', g_untouched, '   wrong-slot=', g_wrong, &
                '   inter_bad=', g_inter_bad, '   intra_bad=', g_intra_bad
            if (g_untouched + g_wrong == 0) then
                write (*, '(2a)') '   [VERDICT] PASS -- every Put was visible at the target after the epoch close'
            else if (g_untouched > 0 .and. g_wrong == 0) then
                write (*, '(2a)') '   [VERDICT] FAIL: DATA NEVER ARRIVED (the historical symptom).', &
                    '  put_ierr above tells you whether MPI also claimed success.'
            else
                write (*, '(2a)') '   [VERDICT] FAIL: data DELIVERED BUT MIS-PLACED', &
                    '  -- an offset/displacement problem, NOT a delivery problem.'
            end if
            if (g_inter_bad > 0 .and. g_intra_bad == 0) &
                write (*, '(a)') '   [HINT] failures are INTER-NODE ONLY -> fabric/RDMA registration, not RMA semantics.'
            if (g_intra_bad > 0 .and. g_inter_bad == 0) &
                write (*, '(a)') '   [HINT] failures are INTRA-NODE ONLY -> shared-memory RMA path, not the NIC.'
        end if

        !--- teardown ------------------------------------------------------
        if (win_kind == WIN_DYNAMIC) then
            call MPI_Win_detach(win, wbuf, ierr)
            deallocate (peer_base)
        end if
        call MPI_Win_free(win, ierr)
        call free_buf(sptr, sbuf, mem_kind)
        if (win_kind /= WIN_ALLOC) call free_buf(wptr, wbuf, mem_kind)
        deallocate (host_stage)
        if (do_taint) call taint_off()
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine run_case

    !-----------------------------------------------------------------------
    subroutine do_puts(win, sbuf, mas, win_kind, peer_base, ierr_max)
        integer, intent(in) :: win, mas, win_kind
        real(dp), intent(in) :: sbuf(*)
        integer(MPI_ADDRESS_KIND), intent(in), allocatable :: peer_base(:)
        integer, intent(inout) :: ierr_max
        integer :: pk, ie
        do pk = 0, NPRO_K - 1
            call one_put(win, sbuf, pk, mas, win_kind, peer_base, ie)
            ierr_max = max(ierr_max, ie)
        end do
    end subroutine do_puts

    subroutine one_put(win, sbuf, pk, mas, win_kind, peer_base, ierr)
        integer, intent(in) :: win, pk, mas, win_kind
        real(dp), intent(in) :: sbuf(*)
        integer(MPI_ADDRESS_KIND), intent(in), allocatable :: peer_base(:)
        integer, intent(out) :: ierr
        integer(MPI_ADDRESS_KIND) :: tdisp
        ! Target displacement: our K-coordinate selects the slot we own in EVERY peer's window.
        ! For a dynamic window the displacement is an absolute address at the target.
        if (win_kind == WIN_DYNAMIC) then
            tdisp = peer_base(pk) + int(ims_pro_k, MPI_ADDRESS_KIND)*int(mas, MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
        else
            tdisp = int(ims_pro_k, MPI_ADDRESS_KIND)*int(mas, MPI_ADDRESS_KIND)
        end if
        call MPI_Put(sbuf(pk*mas + 1), mas, MPI_DOUBLE_PRECISION, &
                     pk, tdisp, mas, MPI_DOUBLE_PRECISION, win, ierr)
    end subroutine one_put

    !-----------------------------------------------------------------------
    subroutine verify(wptr, wbuf, mas, mem_kind, n_untouched, n_wrong, n_inter_bad, n_intra_bad, &
                      first_src, first_want, first_got)
        ! Read our OWN window and check slot m against source K-rank m. The read happens on the
        ! HOST; for a device-allocated window we stage it back with omp_target_memcpy first, so
        ! the check never dereferences a device pointer from the CPU.
        type(c_ptr), intent(in) :: wptr
        real(dp), pointer, intent(in) :: wbuf(:)
        integer, intent(in) :: mas, mem_kind
        integer, intent(out) :: n_untouched, n_wrong, n_inter_bad, n_intra_bad, first_src
        real(dp), intent(out) :: first_want, first_got

        real(dp), allocatable, target :: chk(:)
        integer :: pk, j, wsize
        real(dp) :: want, got
        logical :: bad

        wsize = NPRO_K*mas
        allocate (chk(wsize))
        call get_buf_to_host(wptr, wbuf, chk, wsize, mem_kind)

        n_untouched = 0; n_wrong = 0; n_inter_bad = 0; n_intra_bad = 0
        first_src = -1; first_want = 0.0_dp; first_got = 0.0_dp

        do pk = 0, NPRO_K - 1
            do j = 1, mas
                want = real(peer_world(pk), dp)*ENC + real(j, dp)
                got = chk(pk*mas + j)
                bad = .false.
                if (abs(got - FILL) < TOL) then
                    n_untouched = n_untouched + 1; bad = .true.
                else if (abs(got - want) > TOL) then
                    n_wrong = n_wrong + 1; bad = .true.
                end if
                if (bad) then
                    if (peer_node(pk) /= my_node) then
                        n_inter_bad = n_inter_bad + 1
                    else
                        n_intra_bad = n_intra_bad + 1
                    end if
                    if (first_src < 0) then
                        ! Decode what DID land: (got - j)/ENC recovers the origin world rank if
                        ! the payload is a valid encoding that simply arrived at the wrong slot.
                        first_src = nint((got - real(j, dp))/ENC)
                        first_want = want; first_got = got
                    end if
                end if
            end do
        end do
        deallocate (chk)
    end subroutine verify

    !=======================================================================
    ! memory helpers -- one place where MEM_HOST / MEM_USM / MEM_DEV differ
    !=======================================================================
    subroutine alloc_buf(cp, fp, n, mem_kind)
        type(c_ptr), intent(out) :: cp
        real(dp), pointer, intent(out) :: fp(:)
        integer, intent(in) :: n, mem_kind
        if (mem_kind == MEM_DEV) then
#ifdef USE_APU
            cp = omp_target_alloc(int(n, c_size_t)*8_c_size_t, omp_get_default_device())
            call c_f_pointer(cp, fp, [n])          ! address only -- never dereferenced on the host
#endif
        else
            allocate (fp(n))                        ! plain allocate == unified memory on MI300A
            cp = c_loc(fp(1))
        end if
    end subroutine alloc_buf

    subroutine free_buf(cp, fp, mem_kind)
        type(c_ptr), intent(inout) :: cp
        real(dp), pointer, intent(inout) :: fp(:)
        integer, intent(in) :: mem_kind
        if (mem_kind == MEM_DEV) then
#ifdef USE_APU
            call omp_target_free(cp, omp_get_default_device())
#endif
        else
            if (associated(fp)) deallocate (fp)
        end if
        fp => null(); cp = c_null_ptr
    end subroutine free_buf

    subroutine put_host_to_buf(src, cp, fp, n, mem_kind)
        ! host array -> the buffer, using the memory-appropriate route.
        ! MEM_USM deliberately writes through a GPU kernel: that is the production shape
        ! (plain allocate, GPU-dirty) and is what makes B1/B2 a real test of GPU-RMA.
        real(dp), intent(in), target :: src(:)
        type(c_ptr), intent(in) :: cp
        real(dp), pointer, intent(inout) :: fp(:)
        integer, intent(in) :: n, mem_kind
        integer :: i
#ifdef USE_APU
        integer(c_int) :: hip_err, ie
#endif
        select case (mem_kind)
        case (MEM_HOST)
            do i = 1, n
                fp(i) = src(i)
            end do
        case (MEM_USM)
#ifdef USE_APU
            do i = 1, n
                fp(i) = src(i)
            end do
            !$omp target teams distribute parallel do
            do i = 1, n
                fp(i) = fp(i)                        ! touch on device: pages/L2 become GPU-dirty
            end do
            hip_err = hipDeviceSynchronize()
#endif
        case (MEM_DEV)
#ifdef USE_APU
            ie = omp_target_memcpy(cp, c_loc(src(1)), int(n, c_size_t)*8_c_size_t, &
                                   0_c_size_t, 0_c_size_t, &
                                   omp_get_default_device(), omp_get_initial_device())
#endif
        end select
    end subroutine put_host_to_buf

    subroutine get_buf_to_host(cp, fp, dst, n, mem_kind)
        type(c_ptr), intent(in) :: cp
        real(dp), pointer, intent(in) :: fp(:)
        real(dp), intent(out), target :: dst(:)
        integer, intent(in) :: n, mem_kind
        integer :: i
#ifdef USE_APU
        integer(c_int) :: hip_err, ie
#endif
        if (mem_kind == MEM_DEV) then
#ifdef USE_APU
            ie = omp_target_memcpy(c_loc(dst(1)), cp, int(n, c_size_t)*8_c_size_t, &
                                   0_c_size_t, 0_c_size_t, &
                                   omp_get_initial_device(), omp_get_default_device())
#endif
        else
#ifdef USE_APU
            hip_err = hipDeviceSynchronize()         ! make any GPU-side state visible to this CPU read
#endif
            do i = 1, n
                dst(i) = fp(i)
            end do
        end if
    end subroutine get_buf_to_host

    !=======================================================================
    ! PERFORMANCE: one-sided vs the production two-sided inter-node leg
    !=======================================================================
    subroutine run_perf(label, method, nmax_p, nlines_p, niter)
        ! Inter-node K-peers only -- that is the leg one-sided would replace. Same message
        ! size and same comm as production. Reports per-iteration wall time and per-rank
        ! effective bandwidth so the GO/NO-GO is a number, not an impression.
        character(*), intent(in) :: label
        integer, intent(in) :: method, nmax_p, nlines_p, niter

        integer :: mas, wsize, pk, it, ierr, win, l, np_inter
        integer, allocatable :: req(:)
        real(dp), allocatable, target :: sbuf(:), rbuf(:)
        real(dp) :: t0, t1, tloc, tmin, tmax, tsum, mb
        integer(MPI_ADDRESS_KIND) :: nbytes, tdisp

        mas = nmax_p*nlines_p
        wsize = NPRO_K*mas
        nbytes = int(wsize, MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
        allocate (sbuf(wsize), rbuf(wsize), req(2*NPRO_K))
        sbuf(:) = real(ims_rank, dp); rbuf(:) = FILL
        win = MPI_WIN_NULL

        np_inter = 0
        do pk = 0, NPRO_K - 1
            if (peer_node(pk) /= my_node) np_inter = np_inter + 1
        end do

        if (method /= PERF_TWOSIDED) &
            call MPI_Win_create(rbuf, nbytes, 8, MPI_INFO_NULL, fabric_comm_k, win, ierr)

        call MPI_Barrier(fabric_comm_k, ierr)
        t0 = MPI_Wtime()
        do it = 1, niter
            select case (method)
            case (PERF_TWOSIDED)
                l = 0
                do pk = 0, NPRO_K - 1
                    if (peer_node(pk) == my_node) cycle
                    l = l + 1
                    call MPI_Irecv(rbuf(pk*mas + 1), mas, MPI_DOUBLE_PRECISION, pk, 99, &
                                   fabric_comm_k, req(l), ierr)
                end do
                do pk = 0, NPRO_K - 1
                    if (peer_node(pk) == my_node) cycle
                    l = l + 1
                    call MPI_Isend(sbuf(pk*mas + 1), mas, MPI_DOUBLE_PRECISION, pk, 99, &
                                   fabric_comm_k, req(l), ierr)
                end do
                call MPI_Waitall(l, req, MPI_STATUSES_IGNORE, ierr)
            case (PERF_PUT_FENCE)
                call MPI_Win_fence(0, win, ierr)
                do pk = 0, NPRO_K - 1
                    if (peer_node(pk) == my_node) cycle
                    tdisp = int(ims_pro_k, MPI_ADDRESS_KIND)*int(mas, MPI_ADDRESS_KIND)
                    call MPI_Put(sbuf(pk*mas + 1), mas, MPI_DOUBLE_PRECISION, pk, tdisp, &
                                 mas, MPI_DOUBLE_PRECISION, win, ierr)
                end do
                call MPI_Win_fence(0, win, ierr)
            case (PERF_PUT_FLUSH)
                call MPI_Win_lock_all(0, win, ierr)
                do pk = 0, NPRO_K - 1
                    if (peer_node(pk) == my_node) cycle
                    tdisp = int(ims_pro_k, MPI_ADDRESS_KIND)*int(mas, MPI_ADDRESS_KIND)
                    call MPI_Put(sbuf(pk*mas + 1), mas, MPI_DOUBLE_PRECISION, pk, tdisp, &
                                 mas, MPI_DOUBLE_PRECISION, win, ierr)
                end do
                call MPI_Win_flush_all(win, ierr)
                call MPI_Win_unlock_all(win, ierr)
                call MPI_Barrier(fabric_comm_k, ierr)   ! the two-sided path is also target-synchronizing
            end select
        end do
        t1 = MPI_Wtime()
        tloc = (t1 - t0)/real(max(niter, 1), dp)

        call MPI_Reduce(tloc, tmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, fabric_comm_k, ierr)
        call MPI_Reduce(tloc, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, fabric_comm_k, ierr)
        call MPI_Reduce(tloc, tsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, fabric_comm_k, ierr)

        if (ims_pro_k == 0) then
            mb = real(np_inter, dp)*real(mas, dp)*8.0_dp/1024.0_dp/1024.0_dp
            write (*, '(2a,f10.3,a,f10.3,a,f10.3,a,f8.2,a,f8.3,a)') '[PERF] ', label, &
                tsum/real(NPRO_K, dp)*1.0d3, ' ms avg  (min ', tmin*1.0d3, ' max ', tmax*1.0d3, &
                ' ms)   ', mb, ' MB/rank   ', mb/(tsum/real(NPRO_K, dp))/1024.0_dp, ' GB/s/rank'
        end if

        if (method /= PERF_TWOSIDED) call MPI_Win_free(win, ierr)
        deallocate (sbuf, rbuf, req)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end subroutine run_perf

end module onesided_mod

!===========================================================================
program vmpi_onesided
    use onesided_mod
    implicit none
    integer :: ierr, mode, nargs, nmax_p, nlines_p, niter
    character(len=32) :: arg

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, ims_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, ims_nprocs, ierr)
    funit = FUNIT_BASE + ims_rank

    if (ims_nprocs /= NPRO_I*NPRO_K) then
        if (ims_rank == 0) write (*, *) 'ERROR: need ', NPRO_I*NPRO_K, ' ranks, got ', ims_nprocs
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if

    call setup_topology()
    call report_env()

    ! Correctness default: 8 x 1024 = 8192 doubles = 64 KB/peer. Small enough that the
    ! element-wise verification is cheap, large enough to cross the eager/rendezvous
    ! threshold (so a rendezvous-only bug cannot hide).
    mode = 0; nmax_p = 8; nlines_p = 1024; niter = 20
    nargs = command_argument_count()
    if (nargs >= 1) then; call get_command_argument(1, arg); read (arg, *) mode; end if
    if (nargs >= 2) then; call get_command_argument(2, arg); read (arg, *) nmax_p; end if
    if (nargs >= 3) then; call get_command_argument(3, arg); read (arg, *) nlines_p; end if
    if (nargs >= 4) then; call get_command_argument(4, arg); read (arg, *) niter; end if

    select case (mode)
        ! ---- A: isolate the SYNCHRONIZATION axis on host memory -> tests hypothesis (a)
    case (1); call run_case('A1 host  Win_create   NO-CLOSE   [CONTROL: expect FAIL]', &
                            WIN_CREATE, SYNC_NONE, MEM_HOST, .false., nmax_p, nlines_p)
    case (2); call run_case('A2 host  Win_create   fence', &
                            WIN_CREATE, SYNC_FENCE, MEM_HOST, .false., nmax_p, nlines_p)
    case (3); call run_case('A3 host  Win_create   lock/unlock', &
                            WIN_CREATE, SYNC_LOCK, MEM_HOST, .false., nmax_p, nlines_p)
    case (4); call run_case('A4 host  Win_create   lock_all/flush_all', &
                            WIN_CREATE, SYNC_FLUSH, MEM_HOST, .false., nmax_p, nlines_p)
    case (5); call run_case('A5 host  Win_create   PSCW', &
                            WIN_CREATE, SYNC_PSCW, MEM_HOST, .false., nmax_p, nlines_p)
    case (6); call run_case('A6 host  Win_allocate fence', &
                            WIN_ALLOC, SYNC_FENCE, MEM_HOST, .false., nmax_p, nlines_p)
        ! ---- B: correct sync, escalate the MEMORY axis -> tests hypothesis (b)
    case (7); call run_case('B1 USM(GPU-written) Win_create   fence', &
                            WIN_CREATE, SYNC_FENCE, MEM_USM, .false., nmax_p, nlines_p)
    case (8); call run_case('B2 USM(GPU-written) Win_allocate fence', &
                            WIN_ALLOC, SYNC_FENCE, MEM_USM, .false., nmax_p, nlines_p)
    case (9); call run_case('B3 device ptr       Win_create   fence', &
                            WIN_CREATE, SYNC_FENCE, MEM_DEV, .false., nmax_p, nlines_p)
    case (10); call run_case('B4 device ptr       Win_dynamic  lock_all/flush_all', &
                             WIN_DYNAMIC, SYNC_FLUSH, MEM_DEV, .false., nmax_p, nlines_p)
        ! ---- C: shared-window taint alive -> tests hypothesis (c)
    case (11); call run_case('C1 = A2 WITH MPI_Win_allocate_shared alive', &
                             WIN_CREATE, SYNC_FENCE, MEM_HOST, .true., nmax_p, nlines_p)
    case (12); call run_case('C2 = B1 WITH MPI_Win_allocate_shared alive', &
                             WIN_CREATE, SYNC_FENCE, MEM_USM, .true., nmax_p, nlines_p)
        ! ---- PERF: is one-sided actually faster than the production two-sided leg?
    case (20); call run_perf('P1 two-sided Isend/Irecv  ', PERF_TWOSIDED, nmax_p, nlines_p, niter)
    case (21); call run_perf('P2 Put + fence            ', PERF_PUT_FENCE, nmax_p, nlines_p, niter)
    case (22); call run_perf('P3 Put + lock_all/flush   ', PERF_PUT_FLUSH, nmax_p, nlines_p, niter)
    case (23)
        ! Literal args only: a block-local array + constructor comes out ALL ZEROS under
        ! Cray -fopenmp + requires unified_shared_memory (see vmpi_gpuaware M5).
        if (ims_rank == 0) write (*, '(a)') '### P-SWEEP: per-peer 1 -> 64 MB, two-sided vs Put ###'
        call run_perf('P1 two-sided   1MB', PERF_TWOSIDED, 128, 1024, 20)
        call run_perf('P2 Put+fence   1MB', PERF_PUT_FENCE, 128, 1024, 20)
        call run_perf('P3 Put+flush   1MB', PERF_PUT_FLUSH, 128, 1024, 20)
        call run_perf('P1 two-sided   4MB', PERF_TWOSIDED, 128, 4096, 20)
        call run_perf('P2 Put+fence   4MB', PERF_PUT_FENCE, 128, 4096, 20)
        call run_perf('P3 Put+flush   4MB', PERF_PUT_FLUSH, 128, 4096, 20)
        call run_perf('P1 two-sided  16MB', PERF_TWOSIDED, 128, 16384, 20)
        call run_perf('P2 Put+fence  16MB', PERF_PUT_FENCE, 128, 16384, 20)
        call run_perf('P3 Put+flush  16MB', PERF_PUT_FLUSH, 128, 16384, 20)
        call run_perf('P1 two-sided  64MB', PERF_TWOSIDED, 128, 65536, 10)
        call run_perf('P2 Put+fence  64MB', PERF_PUT_FENCE, 128, 65536, 10)
        call run_perf('P3 Put+flush  64MB', PERF_PUT_FLUSH, 128, 65536, 10)
    case (0)
        ! Full correctness matrix in one job. Convenient, but a HANG in an early case hides
        ! the rest -- the PBS script runs each mode as its own timeout-guarded mpirun.
        call run_case('A1 host  Win_create   NO-CLOSE   [CONTROL: expect FAIL]', &
                      WIN_CREATE, SYNC_NONE, MEM_HOST, .false., nmax_p, nlines_p)
        call run_case('A2 host  Win_create   fence', &
                      WIN_CREATE, SYNC_FENCE, MEM_HOST, .false., nmax_p, nlines_p)
        call run_case('A3 host  Win_create   lock/unlock', &
                      WIN_CREATE, SYNC_LOCK, MEM_HOST, .false., nmax_p, nlines_p)
        call run_case('A4 host  Win_create   lock_all/flush_all', &
                      WIN_CREATE, SYNC_FLUSH, MEM_HOST, .false., nmax_p, nlines_p)
        call run_case('A5 host  Win_create   PSCW', &
                      WIN_CREATE, SYNC_PSCW, MEM_HOST, .false., nmax_p, nlines_p)
        call run_case('A6 host  Win_allocate fence', &
                      WIN_ALLOC, SYNC_FENCE, MEM_HOST, .false., nmax_p, nlines_p)
        call run_case('B1 USM(GPU-written) Win_create   fence', &
                      WIN_CREATE, SYNC_FENCE, MEM_USM, .false., nmax_p, nlines_p)
        call run_case('B2 USM(GPU-written) Win_allocate fence', &
                      WIN_ALLOC, SYNC_FENCE, MEM_USM, .false., nmax_p, nlines_p)
        call run_case('B3 device ptr       Win_create   fence', &
                      WIN_CREATE, SYNC_FENCE, MEM_DEV, .false., nmax_p, nlines_p)
        call run_case('B4 device ptr       Win_dynamic  lock_all/flush_all', &
                      WIN_DYNAMIC, SYNC_FLUSH, MEM_DEV, .false., nmax_p, nlines_p)
        call run_case('C1 = A2 WITH MPI_Win_allocate_shared alive', &
                      WIN_CREATE, SYNC_FENCE, MEM_HOST, .true., nmax_p, nlines_p)
        call run_case('C2 = B1 WITH MPI_Win_allocate_shared alive', &
                      WIN_CREATE, SYNC_FENCE, MEM_USM, .true., nmax_p, nlines_p)
    case default
        if (ims_rank == 0) write (*, *) 'unknown mode ', mode
    end select

    call MPI_Finalize(ierr)
end program vmpi_onesided
