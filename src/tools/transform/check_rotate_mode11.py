#!/usr/bin/env python3
"""
check_rotate_mode11.py -- verify that transfields.f90 ParamTransform=11
(rotation of the horizontal velocity about Oy) did what it was supposed to do.

It compares the ORIGINAL flow field against the rotated '.trn' field and certifies,
independently of the Fortran, the properties the transform claims:

  C1  v (component 2) and any extra components (4,5,...) are BIT-IDENTICAL
      -> the transform touched only u and w.
  C2  the plane means turned correctly:
          ubar_new(j) =  ubar(j)*cos(beta) + wbar(j)*sin(beta)
          wbar_new(j) = -ubar(j)*sin(beta) + wbar(j)*cos(beta)
      with beta = alpha_new - alpha_old.
  C3  the rotation is orthogonal: |(ubar,wbar)|(j) is unchanged.
  C4  the mean turns by -beta:  atan2(wbar_new,ubar_new) - atan2(wbar,ubar) = -beta.
  C5  no turbulence is created or lost: <u'^2 + w'^2>(j) is unchanged.
  C6  the requested MODE was actually applied:
          mode 0 -> fluctuations untouched:  u_new - u  is CONSTANT in each y-plane
                    (and equal to ubar_new - ubar); same for w.
          mode 1 -> pointwise rotation:      u_new = u*cos(b) + w*sin(b) everywhere.
      Both residuals are always measured, so the script also REPORTS which mode the
      data is consistent with -- it detects a mode mix-up even if you pass the wrong one.

TLab field-file format (little-endian stream, ONE FILE PER COMPONENT):
    int32   header_offset          # bytes; field data starts here
    int32   nx, ny, nz, nt
    float64 params[(header_offset-20)//8]
    float64 data[nx*ny*nz]         # Fortran order: i fastest, then j (=y), then k

Usage:
    python3 check_rotate_mode11.py ORIG_BASE TRN_BASE ALPHA_OLD ALPHA_NEW [MODE]

      ORIG_BASE   e.g. flow.100000       (reads flow.100000.1 / .2 / .3 ...)
      TRN_BASE    e.g. flow.100000.trn   (reads flow.100000.trn.1 / .2 / .3 ...)
      ALPHA_OLD   the 2nd number of ParamTransform, in radians
      ALPHA_NEW   the 3rd number, in radians
      MODE        optional 4th number: 0 = mean only (default), 1 = mean+fluctuations

    Options:
      --degrees        ALPHA_OLD/ALPHA_NEW are given in degrees
      --scal A B       additionally check two scalar bases are bit-identical
      --sample N       max points to sample per component without numpy (default 400000)
      --full           force a full-field scan even without numpy (slow on big grids)

Stdlib-only (struct + array); numpy is used automatically if importable.
Exit code 0 if every check passes, 1 otherwise.
"""
import sys
import os
import math
import struct
import array

try:
    import numpy as np
    HAVE_NUMPY = True
except Exception:
    HAVE_NUMPY = False

# Tolerances are RELATIVE to the field scale (max|u|,max|w|). All the identities above
# are exact in exact arithmetic; what is left is a couple of rounded operations per
# point, i.e. a few ULP. 1e-12 sits ~3 decades above that and ~orders below any real bug.
TOL_REL = 1.0e-12
# A residual above this (relative) means the identity plainly does not hold -- used to
# decide which mode the data is consistent with.
TOL_MODE = 1.0e-6


# ----------------------------------------------------------------------------- io
def read_field(path):
    """Return (nx, ny, nz, nt, data) with data as array('d') or numpy array."""
    with open(path, 'rb') as f:
        head = f.read(20)
        if len(head) < 20:
            raise ValueError("%s: too small to hold a 20-byte header" % path)
        header_offset, nx, ny, nz, nt = struct.unpack('<5i', head)
        f.seek(header_offset)
        raw = f.read()
    npts = nx * ny * nz
    if len(raw) != npts * 8:
        raise ValueError("%s: expected %d doubles (%dx%dx%d), file holds %.3f" %
                         (path, npts, nx, ny, nz, len(raw) / 8.0))
    if HAVE_NUMPY:
        data = np.frombuffer(raw, dtype='<f8').reshape(nz, ny, nx)   # [k][j][i]
    else:
        data = array.array('d')
        data.frombytes(raw)
        if sys.byteorder != 'little':
            data.byteswap()
    return nx, ny, nz, nt, data


def component_files(base, ncomp_max=8):
    """Existing 'base.N' component files, in order."""
    out = []
    for n in range(1, ncomp_max + 1):
        p = "%s.%d" % (base, n)
        if os.path.exists(p):
            out.append((n, p))
    return out


def files_identical(pa, pb):
    """Byte-compare the field payloads (headers may legitimately differ in nt)."""
    with open(pa, 'rb') as fa, open(pb, 'rb') as fb:
        ha = struct.unpack('<5i', fa.read(20))
        hb = struct.unpack('<5i', fb.read(20))
        if ha[1:4] != hb[1:4]:
            return False, "grid mismatch %s vs %s" % (ha[1:4], hb[1:4])
        fa.seek(ha[0])
        fb.seek(hb[0])
        while True:
            ba = fa.read(1 << 20)
            bb = fb.read(1 << 20)
            if ba != bb:
                return False, "field data differs"
            if not ba:
                return True, "bit-identical"


# ------------------------------------------------------------------- measurements
class Profiles(object):
    """Per-y-plane sums needed by every check, plus the mode residuals."""

    def __init__(self, ny):
        self.ny = ny
        self.n = [0] * ny            # points per plane actually sampled
        self.su = [0.0] * ny         # sum u_old
        self.sw = [0.0] * ny
        self.su2 = [0.0] * ny        # sum u_old^2 + w_old^2
        self.sun = [0.0] * ny        # sum u_new
        self.swn = [0.0] * ny
        self.su2n = [0.0] * ny
        self.scale = 0.0             # max |u|,|w| over old+new
        self.res_rot = 0.0           # max |mode-1 pointwise residual|
        self.res_shift = 0.0         # max |mode-0 residual| (increment not constant)


def measure(u_o, w_o, u_n, w_n, nx, ny, nz, cb, sb, stride_k, stride_i):
    p = Profiles(ny)
    if HAVE_NUMPY:
        # full field, no sampling: numpy makes it cheap
        p.n = [nx * nz] * ny
        p.su = u_o.sum(axis=(0, 2))
        p.sw = w_o.sum(axis=(0, 2))
        p.su2 = (u_o * u_o + w_o * w_o).sum(axis=(0, 2))
        p.sun = u_n.sum(axis=(0, 2))
        p.swn = w_n.sum(axis=(0, 2))
        p.su2n = (u_n * u_n + w_n * w_n).sum(axis=(0, 2))
        p.scale = max(float(np.abs(u_o).max()), float(np.abs(w_o).max()),
                      float(np.abs(u_n).max()), float(np.abs(w_n).max()))
        # mode 1: pointwise rotation residual
        p.res_rot = max(float(np.abs(u_n - (u_o * cb + w_o * sb)).max()),
                        float(np.abs(w_n - (-u_o * sb + w_o * cb)).max()))
        # mode 0: the increment must be constant inside each plane
        du = u_n - u_o
        dw = w_n - w_o
        p.res_shift = max(
            float(np.abs(du - du.mean(axis=(0, 2))[None, :, None]).max()),
            float(np.abs(dw - dw.mean(axis=(0, 2))[None, :, None]).max()))
        return p

    # ---- stdlib path: sample rows, accumulate per plane, two passes over the rows
    du_mean = [0.0] * ny
    dw_mean = [0.0] * ny
    rows = []                                     # (j, base) of every sampled row
    for k in range(0, nz, stride_k):
        for j in range(ny):
            rows.append((j, (k * ny + j) * nx))

    for j, base in rows:
        ro_u = u_o[base:base + nx:stride_i]
        ro_w = w_o[base:base + nx:stride_i]
        rn_u = u_n[base:base + nx:stride_i]
        rn_w = w_n[base:base + nx:stride_i]
        p.n[j] += len(ro_u)
        p.su[j] += sum(ro_u)
        p.sw[j] += sum(ro_w)
        p.sun[j] += sum(rn_u)
        p.swn[j] += sum(rn_w)
        p.su2[j] += sum(a * a + b * b for a, b in zip(ro_u, ro_w))
        p.su2n[j] += sum(a * a + b * b for a, b in zip(rn_u, rn_w))
        du_mean[j] += sum(rn_u) - sum(ro_u)
        dw_mean[j] += sum(rn_w) - sum(ro_w)
        m = max(max(map(abs, ro_u)), max(map(abs, ro_w)),
                max(map(abs, rn_u)), max(map(abs, rn_w)))
        if m > p.scale:
            p.scale = m
        for a, b, an, bn in zip(ro_u, ro_w, rn_u, rn_w):
            r = abs(an - (a * cb + b * sb))
            if r > p.res_rot:
                p.res_rot = r
            r = abs(bn - (-a * sb + b * cb))
            if r > p.res_rot:
                p.res_rot = r
    for j in range(ny):
        if p.n[j]:
            du_mean[j] /= p.n[j]
            dw_mean[j] /= p.n[j]
    for j, base in rows:
        ro_u = u_o[base:base + nx:stride_i]
        ro_w = w_o[base:base + nx:stride_i]
        rn_u = u_n[base:base + nx:stride_i]
        rn_w = w_n[base:base + nx:stride_i]
        for a, an in zip(ro_u, rn_u):
            r = abs((an - a) - du_mean[j])
            if r > p.res_shift:
                p.res_shift = r
        for b, bn in zip(ro_w, rn_w):
            r = abs((bn - b) - dw_mean[j])
            if r > p.res_shift:
                p.res_shift = r
    return p


# ----------------------------------------------------------------------- checking
def worst(values):
    """(max |value|, index of it)."""
    m, im = 0.0, 0
    for i, v in enumerate(values):
        if abs(v) > m:
            m, im = abs(v), i
    return m, im


def report(name, err, tol, extra=""):
    ok = err <= tol
    print("  [%s] %-46s  max err %10.3e  (tol %7.1e)%s"
          % ("PASS" if ok else "FAIL", name, err, tol, extra))
    return ok


def main(argv):
    args = [a for a in argv[1:] if not a.startswith('--')]
    flags = [a for a in argv[1:] if a.startswith('--')]
    degrees = '--degrees' in flags
    full = '--full' in flags
    sample = 400000
    scal_pair = None
    for i, a in enumerate(argv[1:]):
        if a == '--sample' and i + 2 < len(argv):
            sample = int(argv[i + 2])
        if a == '--scal' and i + 3 < len(argv):
            scal_pair = (argv[i + 2], argv[i + 3])
    # strip option values out of the positional list
    drop = set()
    for i, a in enumerate(argv):
        if a == '--sample':
            drop.add(argv[i + 1])
        if a == '--scal':
            drop.add(argv[i + 1])
            drop.add(argv[i + 2])
    args = [a for a in args if a not in drop]

    if len(args) < 4:
        print(__doc__)
        return 2

    orig_base, trn_base = args[0], args[1]
    alpha_old, alpha_new = float(args[2]), float(args[3])
    mode = int(float(args[4])) if len(args) > 4 else 0
    if degrees:
        alpha_old = math.radians(alpha_old)
        alpha_new = math.radians(alpha_new)
    beta = alpha_new - alpha_old
    cb, sb = math.cos(beta), math.sin(beta)

    print("transfields ParamTransform=11 verification")
    print("  original : %s.{1..}" % orig_base)
    print("  rotated  : %s.{1..}" % trn_base)
    print("  alpha_old=%.10g  alpha_new=%.10g  beta=%.10g rad (%.6g deg)  mode=%d"
          % (alpha_old, alpha_new, beta, math.degrees(beta), mode))
    print("  backend  : %s" % ("numpy (full field)" if HAVE_NUMPY else "stdlib (sampled)"))

    comps_o = component_files(orig_base)
    comps_t = component_files(trn_base)
    have_o = set(n for n, _ in comps_o)
    have_t = set(n for n, _ in comps_t)
    if not ({1, 3} <= have_o and {1, 3} <= have_t):
        print("ERROR: need components 1 (u) and 3 (w) in both bases; "
              "found %s and %s" % (sorted(have_o), sorted(have_t)))
        return 2
    if have_o != have_t:
        print("WARNING: component sets differ: original %s, rotated %s"
              % (sorted(have_o), sorted(have_t)))

    ok = True

    # --- C1: untouched components ------------------------------------------------
    print("\nC1  components other than u,w are untouched")
    others = sorted((have_o & have_t) - {1, 3})
    if not others:
        print("  [ -- ] no components besides u,w to check")
    for n in others:
        same, why = files_identical("%s.%d" % (orig_base, n), "%s.%d" % (trn_base, n))
        print("  [%s] component %d %s" % ("PASS" if same else "FAIL", n, why))
        ok = ok and same
    if scal_pair:
        for n, p in component_files(scal_pair[0]):
            q = "%s.%d" % (scal_pair[1], n)
            if not os.path.exists(q):
                print("  [FAIL] scalar %d missing in %s" % (n, scal_pair[1]))
                ok = False
                continue
            same, why = files_identical(p, q)
            print("  [%s] scalar %d %s" % ("PASS" if same else "FAIL", n, why))
            ok = ok and same

    # --- load u,w ----------------------------------------------------------------
    nx, ny, nz, nt_o, u_o = read_field("%s.1" % orig_base)
    nxt, nyt, nzt, nt_t, u_n = read_field("%s.1" % trn_base)
    if (nx, ny, nz) != (nxt, nyt, nzt):
        print("ERROR: grid mismatch %dx%dx%d vs %dx%dx%d" % (nx, ny, nz, nxt, nyt, nzt))
        return 2
    _, _, _, _, w_o = read_field("%s.3" % orig_base)
    _, _, _, _, w_n = read_field("%s.3" % trn_base)
    print("\n  grid %dx%dx%d, itime %d -> %d" % (nx, ny, nz, nt_o, nt_t))

    stride_k = stride_i = 1
    if not HAVE_NUMPY and not full:
        # keep the sampled point count near `sample`, striding k first then i
        total = nx * ny * nz
        need = max(1, total // max(1, sample))
        stride_k = min(nz, int(math.sqrt(need)) or 1)
        stride_i = max(1, need // stride_k)
        stride_i = min(nx, stride_i)
        if stride_k > 1 or stride_i > 1:
            print("  sampling every %d-th k and %d-th i (~%d points/plane pair); "
                  "use --full for a complete scan"
                  % (stride_k, stride_i, (nz // stride_k) * (nx // stride_i)))

    p = measure(u_o, w_o, u_n, w_n, nx, ny, nz, cb, sb, stride_k, stride_i)
    scale = p.scale if p.scale > 0 else 1.0
    tol = TOL_REL * scale

    ubar = [p.su[j] / p.n[j] for j in range(ny)]
    wbar = [p.sw[j] / p.n[j] for j in range(ny)]
    ubn = [p.sun[j] / p.n[j] for j in range(ny)]
    wbn = [p.swn[j] / p.n[j] for j in range(ny)]
    hvar = [p.su2[j] / p.n[j] - ubar[j] ** 2 - wbar[j] ** 2 for j in range(ny)]
    hvarn = [p.su2n[j] / p.n[j] - ubn[j] ** 2 - wbn[j] ** 2 for j in range(ny)]

    # --- C2: the means turned by the rotation matrix ------------------------------
    print("\nC2  plane means follow the rotation matrix (beta=%.6g rad)" % beta)
    eu = [ubn[j] - (ubar[j] * cb + wbar[j] * sb) for j in range(ny)]
    ew = [wbn[j] - (-ubar[j] * sb + wbar[j] * cb) for j in range(ny)]
    m1, j1 = worst(eu)
    m2, j2 = worst(ew)
    ok &= report("ubar_new == ubar*cos+wbar*sin", m1, tol, "  worst plane j=%d" % (j1 + 1))
    ok &= report("wbar_new == -ubar*sin+wbar*cos", m2, tol, "  worst plane j=%d" % (j2 + 1))

    # --- C3: orthogonality --------------------------------------------------------
    print("\nC3  rotation is orthogonal: |mean| unchanged")
    em = [math.hypot(ubn[j], wbn[j]) - math.hypot(ubar[j], wbar[j]) for j in range(ny)]
    m, jm = worst(em)
    ok &= report("|(ubar,wbar)| unchanged", m, tol, "  worst plane j=%d" % (jm + 1))

    # --- C4: the mean turned by -beta ---------------------------------------------
    print("\nC4  the mean vector turns by -beta")
    ea = []
    small = 1.0e-8 * scale
    for j in range(ny):
        if math.hypot(ubar[j], wbar[j]) < small:
            ea.append(0.0)                     # angle undefined for a null mean
            continue
        d = math.atan2(wbn[j], ubn[j]) - math.atan2(wbar[j], ubar[j])
        d = (d + math.pi) % (2 * math.pi) - math.pi
        ea.append(d + beta)                    # expected d = -beta
    m, jm = worst(ea)
    ok &= report("angle_after - angle_before == -beta", m, 1.0e-9,
                 "  worst plane j=%d" % (jm + 1))
    jtop = ny - 1
    print("       top plane j=%d: angle %.10g -> %.10g rad, expected change %.10g"
          % (ny, math.atan2(wbar[jtop], ubar[jtop]),
             math.atan2(wbn[jtop], ubn[jtop]), -beta))
    print("       mean at top plane: (u,w) = (%.10g, %.10g) -> (%.10g, %.10g)"
          % (ubar[jtop], wbar[jtop], ubn[jtop], wbn[jtop]))

    # --- C5: no turbulence gained or lost -----------------------------------------
    print("\nC5  horizontal fluctuation energy <u'^2+w'^2>(y) is unchanged")
    peak = max(abs(v) for v in hvar) if ny else 0.0
    eh = [hvarn[j] - hvar[j] for j in range(ny)]
    m, jm = worst(eh)
    ok &= report("<u'^2+w'^2> unchanged", m, max(TOL_REL * peak, 1e-14 * scale * scale),
                 "  worst plane j=%d, peak %.4e" % (jm + 1, peak))

    # --- C6: the requested mode was applied ---------------------------------------
    print("\nC6  the requested mode was actually applied")
    print("       mode-0 residual (increment constant in plane) : %.3e" % p.res_shift)
    print("       mode-1 residual (pointwise rotation)          : %.3e" % p.res_rot)
    if mode == 0:
        this_ok = report("mode 0: fluctuations untouched", p.res_shift, tol)
        if this_ok and p.res_rot <= TOL_MODE * scale and abs(beta) > 1e-12:
            print("  [WARN] the data ALSO satisfies the pointwise rotation -- the "
                  "fluctuations are too weak here to tell the two modes apart")
        elif not this_ok and p.res_rot <= TOL_MODE * scale:
            print("  [WARN] the data IS a pointwise rotation -- mode 1 was applied, "
                  "not the mode 0 you asked about")
    else:
        this_ok = report("mode 1: pointwise rotation everywhere", p.res_rot, tol)
        if this_ok and p.res_shift <= TOL_MODE * scale and abs(beta) > 1e-12:
            print("  [WARN] the data ALSO looks like a pure plane shift -- the "
                  "fluctuations are too weak here to tell the two modes apart")
        elif not this_ok and p.res_shift <= TOL_MODE * scale:
            print("  [WARN] the data IS a pure plane shift -- mode 0 was applied, "
                  "not the mode 1 you asked about")
    ok &= this_ok
    consistent = []
    if p.res_shift <= TOL_MODE * scale:
        consistent.append(0)
    if p.res_rot <= TOL_MODE * scale:
        consistent.append(1)
    print("       data is consistent with mode(s): %s"
          % (consistent if consistent else "NEITHER -- the field is not a rotation of the original"))

    # --- profile excerpt -----------------------------------------------------------
    print("\n  profile excerpt (j, ubar, wbar -> ubar_new, wbar_new, <u'^2+w'^2>)")
    step = max(1, ny // 10)
    for j in list(range(0, ny, step)) + [ny - 1]:
        print("    j=%5d  %12.6e %12.6e  ->  %12.6e %12.6e   %12.6e"
              % (j + 1, ubar[j], wbar[j], ubn[j], wbn[j], hvar[j]))

    print("\n%s" % ("ALL CHECKS PASSED -- mode 11 did what it was supposed to do."
                    if ok else "SOME CHECKS FAILED -- see the [FAIL] lines above."))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv))
