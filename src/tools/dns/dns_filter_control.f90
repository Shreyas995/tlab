!########################################################################
!#
!# Adaptive in-run controller for the compact [Filter] strength alpha_f.
!#
!# Ramps alpha_f from ParameterStart (strong filtering, survives the
!# polluted-IC transient) toward ParameterEnd (mild, minimal statistics
!# bias) and optionally disables the filter, gated by the DILATATION:
!#   * relax  (alpha += Step) only when the dilatation moving-average is
!#     flat/decaying over the window  (trend ratio <= 1 + RiseTol);
!#   * revert (alpha -= Step) when it is growing on average
!#     (trend ratio >= 1 + RevertTol) -- the milder filter was too weak.
!# ALL decisions use the WINDOWED AVERAGE of the dilatation, never the
!# instantaneous value. A checkpoint is written before each relax step so
!# the change is revertible. Mirrors the viscosity-ramp pattern in dns_main.
!#
!# Default OFF ([FilterControl] Adaptive=no) => no change to behaviour.
!#
!########################################################################
module DNS_FILTER_CONTROL_M
    use TLab_Constants, only: wp, wi, lfile, wfile, tag_flow, tag_scal
    use TLab_Memory, only: imax, jmax, kmax, inb_flow, inb_scal
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use FDM, only: g
    use OPR_FILTERS, only: FilterDomain, OPR_FILTER_REINIT, DNS_FILTER_COMPACT, DNS_FILTER_NONE
    use IO_Fields, only: io_header_q, io_header_s, IO_Write_Fields
    use DNS_LOCAL, only: logs_data
    use TLab_Time, only: itime, rtime
    use NavierStokes, only: nse_eqns, DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC
    implicit none
    save
    private

    public :: DNS_FILTER_CONTROL_INITIALIZE, DNS_FILTER_CONTROL, DNS_FILTER_CONTROL_ALPHA

    logical :: ctrl_active = .false.
    real(wp) :: alpha_cur = -1.0_wp   ! current filter strength; <0 => adaptive control never active
    real(wp) :: alpha_start, alpha_end, alpha_step, rise_tol, revert_tol
    integer :: win, minhold
    integer :: hold_count, buf_pos, buf_count
    real(wp), allocatable :: dil_buf(:)
    logical :: disable_at_end, save_before_change

contains

    !########################################################################
    subroutine DNS_FILTER_CONTROL_INITIALIZE(inifile)
        character(len=*), intent(in) :: inifile

        character(len=32) :: bakfile
        character(len=512) :: sRes
        integer :: ig
        logical :: have_compact

        bakfile = trim(adjustl(inifile))//'.bak'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#[FilterControl]')
        call TLab_Write_ASCII(bakfile, '#Adaptive=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#ParameterStart=<alpha_f at start>')
        call TLab_Write_ASCII(bakfile, '#ParameterEnd=<alpha_f target>')
        call TLab_Write_ASCII(bakfile, '#Step=<alpha_f increment>')
        call TLab_Write_ASCII(bakfile, '#Window=<iters in the moving average>')
        call TLab_Write_ASCII(bakfile, '#MinHold=<min iters between changes>')
        call TLab_Write_ASCII(bakfile, '#RiseTol=<relax if avg-trend <= 1+RiseTol>')
        call TLab_Write_ASCII(bakfile, '#RevertTol=<revert if avg-trend >= 1+RevertTol>')
        call TLab_Write_ASCII(bakfile, '#DisableAtEnd=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#SaveBeforeChange=<yes/no>')

        ctrl_active = .false.
        call ScanFile_Char(bakfile, inifile, 'FilterControl', 'Adaptive', 'no', sRes)
        if (trim(adjustl(sRes)) == 'yes') ctrl_active = .true.
        if (.not. ctrl_active) return

        ! Dilatation-gated control only meaningful for incompressible/anelastic,
        ! where logs_data(10:11) hold DilMin/DilMax.
        if (nse_eqns /= DNS_EQNS_INCOMPRESSIBLE .and. nse_eqns /= DNS_EQNS_ANELASTIC) then
            call TLab_Write_ASCII(wfile, 'DNS_FILTER_CONTROL. Needs incompressible/anelastic; adaptive control disabled.')
            ctrl_active = .false.
            return
        end if

        ! Needs at least one active compact [Filter] direction to ramp.
        have_compact = .false.
        do ig = 1, 3
            if (FilterDomain(ig)%type == DNS_FILTER_COMPACT) have_compact = .true.
        end do
        if (.not. have_compact) then
            call TLab_Write_ASCII(wfile, 'DNS_FILTER_CONTROL. No active compact [Filter]; adaptive control disabled.')
            ctrl_active = .false.
            return
        end if

        call ScanFile_Real(bakfile, inifile, 'FilterControl', 'ParameterStart', '0.40', alpha_start)
        call ScanFile_Real(bakfile, inifile, 'FilterControl', 'ParameterEnd', '0.49', alpha_end)
        call ScanFile_Real(bakfile, inifile, 'FilterControl', 'Step', '0.03', alpha_step)
        call ScanFile_Int(bakfile, inifile, 'FilterControl', 'Window', '150', win)
        call ScanFile_Int(bakfile, inifile, 'FilterControl', 'MinHold', '250', minhold)
        call ScanFile_Real(bakfile, inifile, 'FilterControl', 'RiseTol', '0.02', rise_tol)
        call ScanFile_Real(bakfile, inifile, 'FilterControl', 'RevertTol', '0.10', revert_tol)
        call ScanFile_Char(bakfile, inifile, 'FilterControl', 'DisableAtEnd', 'yes', sRes)
        disable_at_end = (trim(adjustl(sRes)) == 'yes')
        call ScanFile_Char(bakfile, inifile, 'FilterControl', 'SaveBeforeChange', 'yes', sRes)
        save_before_change = (trim(adjustl(sRes)) == 'yes')

        ! ---- sanitise ----
        if (win < 4) win = 4
        if (mod(win, 2) /= 0) win = win + 1          ! even, so the two halves are equal
        if (minhold < win) minhold = win             ! averages must reflect the current alpha
        if (alpha_step <= 0.0_wp) alpha_step = 0.03_wp
        if (alpha_end < alpha_start) alpha_end = alpha_start
        if (revert_tol <= rise_tol) revert_tol = rise_tol + 0.05_wp

        allocate (dil_buf(win))
        dil_buf(:) = 0.0_wp
        buf_pos = 0; buf_count = 0; hold_count = 0

        alpha_cur = alpha_start
        call apply_alpha()                           ! start the run at ParameterStart

        write (sRes, 1000) alpha_start, alpha_end, alpha_step, win, minhold
        call TLab_Write_ASCII(lfile, 'DNS_FILTER_CONTROL active: '//trim(adjustl(sRes)))

        return
1000    format('start=', F5.3, ' end=', F5.3, ' step=', F5.3, ' window=', I0, ' minhold=', I0)
    end subroutine DNS_FILTER_CONTROL_INITIALIZE

    !########################################################################
    ! Run ONCE per completed iteration, after the dilatation is computed
    ! (dns_main, just after DNS_BOUNDS_CONTROL). q,s are passed so a
    ! checkpoint can be written before a relax step.
    !########################################################################
    subroutine DNS_FILTER_CONTROL(q, s)
        real(wp), dimension(:, :), intent(in) :: q, s

        real(wp) :: dil, mean_recent, mean_older, ratio
        integer :: nh, i, idx
        character(len=128) :: line

        if (.not. ctrl_active) return

        ! ---- push current dilatation magnitude into the ring buffer ----
        dil = max(abs(logs_data(10)), abs(logs_data(11)))
        buf_pos = buf_pos + 1
        if (buf_pos > win) buf_pos = 1
        dil_buf(buf_pos) = dil
        if (buf_count < win) buf_count = buf_count + 1
        hold_count = hold_count + 1

        ! ---- gate: need a full window AND a minimum hold since the last change ----
        if (buf_count < win) return
        if (hold_count < minhold) return

        ! ---- two half-window averages: newest nh vs the nh before them ----
        nh = win/2
        mean_recent = 0.0_wp
        do i = 0, nh - 1
            idx = buf_pos - i; if (idx < 1) idx = idx + win
            mean_recent = mean_recent + dil_buf(idx)
        end do
        mean_older = 0.0_wp
        do i = nh, 2*nh - 1
            idx = buf_pos - i; if (idx < 1) idx = idx + win
            mean_older = mean_older + dil_buf(idx)
        end do
        mean_recent = mean_recent/real(nh, wp)
        mean_older = mean_older/real(nh, wp)
        ratio = mean_recent/max(mean_older, 1.0e-30_wp)

        ! ---- decide, all on averaged quantities ----
        if (ratio >= 1.0_wp + revert_tol .and. alpha_cur > alpha_start + 1.0e-6_wp) then
            ! growing on average -> milder filter too weak: re-strengthen in place
            alpha_cur = max(alpha_start, alpha_cur - alpha_step)
            call apply_alpha()
            call reset_window()
            write (line, 2000) 'REVERT ', itime, alpha_cur, mean_recent, ratio
            call TLab_Write_ASCII(lfile, trim(adjustl(line)))

        else if (ratio <= 1.0_wp + rise_tol .and. alpha_cur < alpha_end - 1.0e-6_wp) then
            ! flat/decaying on average -> safe to relax; checkpoint first (revertible)
            if (save_before_change) call save_checkpoint(q, s)
            alpha_cur = min(alpha_end, alpha_cur + alpha_step)
            call apply_alpha()
            call reset_window()
            write (line, 2000) 'RELAX  ', itime, alpha_cur, mean_recent, ratio
            call TLab_Write_ASCII(lfile, trim(adjustl(line)))

        else if (ratio <= 1.0_wp + rise_tol .and. alpha_cur >= alpha_end - 1.0e-6_wp .and. disable_at_end) then
            ! at the mild end and calm -> disable the domain filter, stop controlling
            call disable_filter()
            ctrl_active = .false.
            write (line, 2000) 'DISABLE', itime, alpha_cur, mean_recent, ratio
            call TLab_Write_ASCII(lfile, trim(adjustl(line)))
        end if

        return
2000    format('DNS_FILTER_CONTROL ', A7, ' it=', I0, ' alpha=', F5.3, ' dil_avg=', ES10.3, ' trend=', F6.3)
    end subroutine DNS_FILTER_CONTROL

    !########################################################################
    ! Current compact-filter strength alpha_f (for logging). Holds the last
    ! value even after DISABLE; -1.0 if adaptive control was never active.
    !########################################################################
    function DNS_FILTER_CONTROL_ALPHA() result(a)
        real(wp) :: a
        a = alpha_cur
    end function DNS_FILTER_CONTROL_ALPHA

    !########################################################################
    ! Push alpha_cur into every active compact [Filter] direction and rebuild
    ! its coefficients in place (no re-allocation). On MI300A unified memory
    ! the host coeff update is visible to the GPU filter without an explicit map.
    !########################################################################
    subroutine apply_alpha()
        integer :: ig
        do ig = 1, 3
            if (FilterDomain(ig)%type == DNS_FILTER_COMPACT) then
                FilterDomain(ig)%parameters(1) = alpha_cur
                call OPR_FILTER_REINIT(g(ig), FilterDomain(ig))
            end if
        end do
        return
    end subroutine apply_alpha

    !########################################################################
    subroutine disable_filter()
        integer :: ig
        do ig = 1, 3
            if (FilterDomain(ig)%type == DNS_FILTER_COMPACT) FilterDomain(ig)%type = DNS_FILTER_NONE
        end do
        return
    end subroutine disable_filter

    !########################################################################
    subroutine reset_window()
        ! Discard pre-change history so the next decision uses only data under
        ! the new alpha (refills over `win` iters; minhold >= win guarantees it).
        buf_count = 0; buf_pos = 0; hold_count = 0
        return
    end subroutine reset_window

    !########################################################################
    ! Standard restart save (flow.<itime>/scal.<itime>), mirroring the
    ! check-pointing block in dns_main; the saved itime is a valid restart point.
    !########################################################################
    subroutine save_checkpoint(q, s)
        real(wp), dimension(:, :), intent(in) :: q, s
        character(len=32) :: fname

        if (inb_flow > 0) then
            write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
            io_header_q(1)%params(1) = rtime
            call IO_Write_Fields(fname, imax, jmax, kmax, itime, inb_flow, q, io_header_q(1:1))
        end if
        if (inb_scal > 0) then
            write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
            io_header_s(:)%params(1) = rtime
            call IO_Write_Fields(fname, imax, jmax, kmax, itime, inb_scal, s, io_header_s(1:inb_scal))
        end if
        return
    end subroutine save_checkpoint

end module DNS_FILTER_CONTROL_M
