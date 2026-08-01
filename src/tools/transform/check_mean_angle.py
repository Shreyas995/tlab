#!/usr/bin/env python3
"""
check_mean_angle.py -- measure the angle of the MEAN horizontal wind in a TLab flow
field, plane by plane, and check the free stream against the geostrophic direction.

Use this when you only have ONE field (e.g. the '.trn' output of transfields
ParamTransform=11) and want to confirm the mean wind ended up at the intended angle.
No reference/original field is needed.

CONVENTION (from src/physics/rotation.f90, Rotation_Coriolis, EQNS_COR_NORMALIZED):
        geo_u =  cos(alpha)*G
        geo_w = -sin(alpha)*G
    with alpha = [Rotation] Parameters(1) [rad] and G = Parameters(2). So the
    geostrophic wind points at atan2(geo_w, geo_u) = -alpha with respect to Ox:
    a POSITIVE alpha means the free-stream wind is rotated CLOCKWISE (towards -z).
    transfields mode 11 with (alpha_old, alpha_new) turns the mean by
    beta = -(alpha_new - alpha_old), i.e. it re-aims the field at -alpha_new.

The script reports, for every y-plane j:
    ubar(j), wbar(j), |mean|(j), the angle atan2(wbar,ubar), and the implied
    alpha = -angle; then it checks the free stream (the top planes) against the
    expected -alpha and prints the Ekman veering across the profile.

TLab field-file format (little-endian stream, ONE FILE PER COMPONENT):
    int32 header_offset; int32 nx, ny, nz, nt; float64 params[]; float64 data[nx*ny*nz]
    data is Fortran order: i fastest, then j (=y), then k.  base.1 = u, base.3 = w.

Usage:
    python3 check_mean_angle.py BASE [--alpha A | --ini tlab.ini] [options]

      BASE            e.g. flow.100000.trn   (reads flow.100000.trn.1 and .3)

    Options:
      --alpha A       expected alpha in radians (the alpha_new you gave transfields)
      --ini FILE      read alpha and G from the [Rotation] Parameters line of FILE
      --degrees       --alpha is given in degrees (output always shows both)
      --tol T         tolerance on the free-stream angle, in degrees (default 0.5)
      --top N         average the free stream over the top N planes (default 5)
      --profile       print the full per-plane table instead of an excerpt
      --sample N      max points sampled per component without numpy (default 400000)
      --full          force a complete scan even without numpy (slow on big grids)

Stdlib-only (struct + array); numpy is used automatically if importable.
Exit code 0 if the free-stream angle matches the expectation (or if none was given
and the field could be read), 1 otherwise.
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


# ----------------------------------------------------------------------------- io
def read_field(path):
    """Return (nx, ny, nz, nt, data) with data as numpy [k][j][i] or array('d')."""
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
        return nx, ny, nz, nt, np.frombuffer(raw, dtype='<f8').reshape(nz, ny, nx)
    data = array.array('d')
    data.frombytes(raw)
    if sys.byteorder != 'little':
        data.byteswap()
    return nx, ny, nz, nt, data


def plane_means(u, w, nx, ny, nz, stride_k, stride_i):
    """Per-y-plane means of u and w."""
    if HAVE_NUMPY:
        return list(u.mean(axis=(0, 2))), list(w.mean(axis=(0, 2)))
    su = [0.0] * ny
    sw = [0.0] * ny
    n = [0] * ny
    for k in range(0, nz, stride_k):
        for j in range(ny):
            base = (k * ny + j) * nx
            ru = u[base:base + nx:stride_i]
            rw = w[base:base + nx:stride_i]
            su[j] += sum(ru)
            sw[j] += sum(rw)
            n[j] += len(ru)
    return [su[j] / n[j] for j in range(ny)], [sw[j] / n[j] for j in range(ny)]


def read_ini_rotation(path):
    """Return (alpha, G) from the [Rotation] Parameters line, or (None, None)."""
    alpha = geo = None
    section = None
    with open(path, 'r', errors='replace') as f:
        for line in f:
            s = line.split('#')[0].split('!')[0].strip()
            if not s:
                continue
            if s.startswith('[') and s.endswith(']'):
                section = s[1:-1].strip().lower()
                continue
            if section == 'rotation' and '=' in s:
                key, val = s.split('=', 1)
                if key.strip().lower() == 'parameters':
                    vals = [v for v in val.replace(',', ' ').split() if v]
                    if vals:
                        alpha = float(vals[0])
                    if len(vals) > 1:
                        geo = float(vals[1])
    return alpha, geo


def wrap_pi(a):
    return (a + math.pi) % (2 * math.pi) - math.pi


# --------------------------------------------------------------------------- main
def main(argv):
    args = []
    opts = {}
    i = 1
    valued = {'--alpha', '--ini', '--tol', '--top', '--sample'}
    while i < len(argv):
        a = argv[i]
        if a in valued:
            opts[a] = argv[i + 1]
            i += 2
        elif a.startswith('--'):
            opts[a] = True
            i += 1
        else:
            args.append(a)
            i += 1

    if not args:
        print(__doc__)
        return 2

    base = args[0]
    degrees = '--degrees' in opts
    tol_deg = float(opts.get('--tol', 0.5))
    ntop = int(opts.get('--top', 5))
    sample = int(opts.get('--sample', 400000))

    alpha = None
    geo = None
    if '--ini' in opts:
        alpha, geo = read_ini_rotation(opts['--ini'])
        if alpha is None:
            print("WARNING: no [Rotation] Parameters found in %s" % opts['--ini'])
    if '--alpha' in opts:
        alpha = float(opts['--alpha'])
        if degrees:
            alpha = math.radians(alpha)
    if len(args) > 1 and alpha is None:          # allow the alpha as a bare 2nd argument
        alpha = float(args[1])
        if degrees:
            alpha = math.radians(alpha)

    fu = "%s.1" % base
    fw = "%s.3" % base
    for p in (fu, fw):
        if not os.path.exists(p):
            print("ERROR: %s not found (expected the u and w component files "
                  "'%s.1' and '%s.3')" % (p, base, base))
            return 2

    nx, ny, nz, nt, u = read_field(fu)
    nxw, nyw, nzw, _, w = read_field(fw)
    if (nx, ny, nz) != (nxw, nyw, nzw):
        print("ERROR: u and w grids differ: %dx%dx%d vs %dx%dx%d"
              % (nx, ny, nz, nxw, nyw, nzw))
        return 2

    print("mean horizontal wind angle")
    print("  field   : %s.{1,3}   grid %dx%dx%d   itime %d" % (base, nx, ny, nz, nt))
    print("  backend : %s" % ("numpy (full field)" if HAVE_NUMPY else "stdlib (sampled)"))
    print("  sign    : geo_u=cos(alpha)*G, geo_w=-sin(alpha)*G  =>  wind angle = -alpha")

    stride_k = stride_i = 1
    if not HAVE_NUMPY and '--full' not in opts:
        need = max(1, (nx * ny * nz) // max(1, sample))
        stride_k = min(nz, max(1, int(math.sqrt(need))))
        stride_i = min(nx, max(1, need // stride_k))
        if stride_k > 1 or stride_i > 1:
            print("  sampling: every %d-th k, every %d-th i (~%d points/plane); "
                  "--full for a complete scan"
                  % (stride_k, stride_i, (nz // stride_k) * (nx // stride_i)))

    ubar, wbar = plane_means(u, w, nx, ny, nz, stride_k, stride_i)
    mag = [math.hypot(ubar[j], wbar[j]) for j in range(ny)]
    ang = [math.atan2(wbar[j], ubar[j]) for j in range(ny)]

    # ---- profile table -----------------------------------------------------------
    print("\n  j        ubar          wbar         |mean|      angle[deg]   implied alpha[deg]")
    idxs = range(ny) if '--profile' in opts else \
        sorted(set(list(range(0, ny, max(1, ny // 12))) + list(range(max(0, ny - 3), ny))))
    for j in idxs:
        print("  %5d %13.6e %13.6e %12.6e %12.4f %14.4f"
              % (j + 1, ubar[j], wbar[j], mag[j],
                 math.degrees(ang[j]), math.degrees(-ang[j])))

    # ---- free stream ---------------------------------------------------------------
    ntop = max(1, min(ntop, ny))
    su = sum(ubar[ny - ntop:]) / ntop
    sw = sum(wbar[ny - ntop:]) / ntop
    ang_top = math.atan2(sw, su)
    mag_top = math.hypot(su, sw)
    jmax = max(range(ny), key=lambda j: mag[j])

    print("\n  free stream (mean over the top %d planes, j=%d..%d):" % (ntop, ny - ntop + 1, ny))
    print("    (ubar,wbar) = (%.10g, %.10g)   |mean| = %.10g" % (su, sw, mag_top))
    print("    angle       = %.6f deg  (%.10g rad)" % (math.degrees(ang_top), ang_top))
    print("    implied alpha = -angle = %.6f deg  (%.10g rad)"
          % (math.degrees(-ang_top), -ang_top))
    print("    top plane j=%d alone: angle %.6f deg, |mean| %.6g"
          % (ny, math.degrees(ang[ny - 1]), mag[ny - 1]))
    print("    strongest mean at j=%d: |mean| %.6g, angle %.6f deg"
          % (jmax + 1, mag[jmax], math.degrees(ang[jmax])))

    # ---- Ekman veering --------------------------------------------------------------
    jlow = next((j for j in range(ny) if mag[j] > 0.05 * mag_top), 0)
    print("\n  veering across the profile:")
    print("    lowest plane with |mean| > 5%% of the free stream: j=%d, angle %.4f deg"
          % (jlow + 1, math.degrees(ang[jlow])))
    print("    turning from there to the free stream: %.4f deg"
          % math.degrees(wrap_pi(ang_top - ang[jlow])))

    # ---- verdict ----------------------------------------------------------------------
    ok = True
    if alpha is None:
        print("\n  no expected alpha given (--alpha A / --ini tlab.ini): "
              "reporting the measured angle only.")
    else:
        expected = wrap_pi(-alpha)
        err = wrap_pi(ang_top - expected)
        ok = abs(math.degrees(err)) <= tol_deg
        print("\n  expected: alpha = %.10g rad (%.6f deg)  =>  wind angle %.6f deg"
              % (alpha, math.degrees(alpha), math.degrees(expected)))
        print("  measured: %.6f deg      difference: %.6f deg   (tol %.3f deg)"
              % (math.degrees(ang_top), math.degrees(err), tol_deg))
        print("  [%s] free-stream mean wind %s aligned with the geostrophic direction"
              % ("PASS" if ok else "FAIL", "is" if ok else "is NOT"))
        # Diagnosis only (not a gate), so it uses a deliberately loose window: the
        # free-stream mean of a real field carries a degree or so of scatter.
        if not ok and abs(math.degrees(wrap_pi(ang_top - alpha))) <= max(5 * tol_deg, 2.0):
            print("  [WARN] the measured angle is ~ +alpha, not -alpha -- "
                  "this is exactly what a SIGN ERROR in the rotation looks like")
        if geo is not None:
            rel = abs(mag_top - geo) / geo if geo else 0.0
            print("  geostrophic magnitude G=%.10g from the ini, measured |mean|=%.10g "
                  "(%.2f%% off)" % (geo, mag_top, 100 * rel))

    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv))
