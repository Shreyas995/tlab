# HLRS HUNTER-GPU (AMD MI300A APU) — AMD compiler (flang / amdflang-new) variant.
#
# This is the AMD-flang counterpart of config/hunter-cray.cmake (HPE Cray CCE, the verified default).
# It produces the SAME -D macros and the SAME build toggles; only the compiler driver and the
# compiler-specific flags differ (Cray `-eZ`/`-h*` -> flang `-cpp`/`-O*`).
#
# ⚠️ STATUS [2026-06-30] — AMD-compiler build is BLOCKED on Hunter (system provisioning gap, not code):
#   cray-mpich ships `mpi_f08.mod` ONLY for the GNU frontend (`.../ofi/gnu/11.2/include`); the AMD flavor
#   `.../ofi/amd/6.0/include` has only the F90 `mpi.mod` (+ mpi_base/constants/sizeofs) — NO `mpi_f08.mod`
#   in ANY mpich version (8.1.30-9.0.1, confirmed by `find /opt/cray/pe/mpich -name mpi_f08.mod`). The code
#   requires `use mpi_f08` (tlab_mpi_vars.f90). So:
#     * `ftn` under PrgEnv-amd + amd/6.4.1 drives `amdflang` = CLASSIC flang (flang-legacy, ROCm 6.4.1) — it
#       both lacks mpi_f08 AND can't parse OpenMP-5 `!$omp declare mapper` (tlab_type.f90). Dead end.
#     * The LLVM `amdflang-new` exists ONLY in un-modulized /opt/rocm-7.0.2 — it parses OpenMP-5 fine but
#       has NO MPI module it can read (no amd mpi_f08; the amd `mpi.mod` is classic-format, and LLVM-flang
#       and classic-flang .mod files are mutually unreadable). Also dead end.
#   FIX = HLRS must modulize rocm-7.0.2 / amdflang-new AND build cray-mpich F2008 modules for it. Until then
#   use CCE (config/hunter-cray.cmake, verified). Closest non-CCE alternative on the AMD GPU today is
#   PrgEnv-gnu-amd (gfortran frontend — HAS gnu mpi_f08.mod — + AMD GPU offload), untested for this code.
#
# HOW TO SELECT THE AMD COMPILER ON HUNTER (Cray/HPE system):
#   Hunter's `ftn` is the Cray compiler *wrapper*; it dispatches to whichever PrgEnv is loaded and
#   provides Cray-MPI + module integration regardless. To make `ftn` drive AMD flang instead of CCE:
#       module swap PrgEnv-cray PrgEnv-amd     # (or: module load PrgEnv-amd)
#   then `ftn` wraps the AMD compiler and you still get MPI for free. NOTE [2026-06-30]: under PrgEnv-amd +
#   amd/6.4.1, `ftn`->`amdflang` = CLASSIC flang, NOT amdflang-new (see the STATUS block above). This is why
#   we keep the compiler as `ftn` below (NOT a bare `amdflang-new`, which has no MPI). Do NOT additionally
#   `module load craype-accel-amd-gfx942` — per the project notes that reloads cce and breaks MPI
#   predefined datatypes; the .bashrc HLRS/APU env already provides the offload toolchain.
#   (If you ever build OUTSIDE a Cray PrgEnv, e.g. the RAC-Plano cluster, set FC to an MPI wrapper
#    around amdflang-new instead — see the commented block at the bottom.)
#
# Build, exactly like the Cray config:
#   cd build_parallel && cmake ../src -DSYST=hunter-flang -DBUILD_TYPE=PARALLEL -DACCELERATE=TRUE && make dns.x

if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE PARALLEL)
endif()

if ( NOT HYBRID )
   set(HYBRID FALSE)
else()
   message(WARNING "Compiling for hybrid openMP/MPI usage (Only in combination! Not tested!)")
endif()

if ( NOT PROFILE )
   set(PROFILE FALSE)
endif()

if (${PROFILE} STREQUAL "TRUE" )
   set(USER_profile_FLAGS "-g")
endif()

# ACCELERATE=TRUE -> GPU offload (ALL transposes on the GPU). flang needs the offload arch passed
# EXPLICITLY (--offload-arch=gfx942); the Cray config gets it from the loaded craype-accel module.
# A plain -DACCELERATE=TRUE gives the COMPLETE all-on-GPU build: the per-workgroup write+fence
# coherence fix is folded in automatically (it is REQUIRED for a correct GPU run, not a choice).
# See the matching note in config/hunter-cray.cmake: the original apudirect was fast because it had
# NO fences/syncs, NOT because the complex transposes ran on the CPU -- a CPU/MPI transpose is always
# slower on this unified APU, so the all-GPU fenced path is the production path.
if ( NOT ACCELERATE )
  set(ACCELERATE "FALSE")
elseif( ${ACCELERATE} STREQUAL "TRUE" )
  set(USER_APU_FLAGS "-fopenmp --offload-arch=gfx942")
  add_definitions(-DUSE_APU)

  # Default GPU-transpose coherence fix = per-workgroup fused write+fence + reader L2 invalidate.
  set(_TRP_FUSED_DEFAULT TRUE)

  # ----- Advanced / experimental transpose overrides (default OFF; A/B diagnosis only; all slower) -----
  if (TRP_I_SYSFENCE)         # LEGACY separate fence kernel; superseded by, and exclusive with, the fused fence
    add_definitions(-DTRP_I_SYSFENCE)
    set(_TRP_FUSED_DEFAULT FALSE)
  endif ()
  if (_TRP_FUSED_DEFAULT)
    add_definitions(-DTRP_I_FUSEDFENCE)
  endif ()
  if (TRP_I_MEMCPY)           # K pushes via blocking hipMemcpy instead of the fused fence
    add_definitions(-DTRP_I_MEMCPY)
  endif ()
  if (TRP_I_FORCE_MPI)        # route the intra-node I-transpose onto MPI (legacy stabilizer, slower)
    add_definitions(-DTRP_I_FORCE_MPI)
  endif ()
  if (TRP_CX_MPI)             # route the complex (Poisson) transposes onto MPI (slower, A/B only)
    add_definitions(-DTRP_CX_MPI)
  endif ()
endif()

# Debug probe LEVEL (-DDNS_DEBUG=<n>; default 0 = OFF, zero cost). Replaces -DDNS_DEBUG_PROBES + -DPROBESYNC.
#   1 = full probes (NWFR/NWBR/MAXVAL/POIS-trace -> fort.5xx, "print all").
#   2 = Heisenbug isolation (probes do ONLY hipDeviceSynchronize, no read/log).
if (NOT DNS_DEBUG)
  set(DNS_DEBUG 0)
endif ()
if (DNS_DEBUG GREATER_EQUAL 1)
  add_definitions(-DDNS_DEBUG_PROBES)
endif ()
if (DNS_DEBUG GREATER_EQUAL 2)
  add_definitions(-DPROBE_SYNC_ONLY)
endif ()

# compiler for parallel build
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" OR ${BUILD_TYPE} STREQUAL "NONBLOCKING" )
   set(ENV{FC} ftn) # Cray wrapper -> AMD flang when PrgEnv-amd is loaded; provides Cray-MPI.
   add_definitions(-DUSE_MPI -DUSE_MPI_IO -DUSE_ALLTOALL)

# OpenMP (hybrid) flags
   if ( ${HYBRID} STREQUAL "TRUE" )
     set(USER_OMP_FLAGS " -fopenmp")
     add_definitions(-DUSE_OPENMP)
   endif()

# compiler for serial build
else( ${BUILD_TYPE} STREQUAL "SERIAL" )
  set(ENV{FC} ftn)
endif()

# Fortran flags. flang: -cpp for the C preprocessor (the Cray config uses -eZ instead).
set(USER_Fortran_FLAGS         "-cpp ${USER_OMP_FLAGS} ${USER_APU_FLAGS} ${USER_profile_FLAGS}")

# Release flags. The Cray config's -hipa2/-hfp2/-hunroll2/-hfusion2/-hscalar1/-m4 are CCE-only and have
# NO flang equivalent; use a portable optimization level. Start conservative (-O2) and bump to -O3 only
# after a bit-match check vs the apudirect reference (compact FDM is sensitive to aggressive reassoc).
set(USER_Fortran_FLAGS_RELEASE "-O2")
set(USER_Fortran_FLAGS_DEBUG   "-O0 -g -ffpe-trap=invalid,zero,overflow")

if ( NOT CMAKE_BUILD_TYPE )
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

# FFT backend: FFTW (default) or hipFFT (-DHIPFFT=TRUE). Same switch as the Cray config. Uses the
# Hunter module-provided FFTW/NetCDF (cray-fftw + cray module env), NOT a hardcoded cluster path.
if (HIPFFT)
  add_definitions(-DNO_ASSUMED_RANKS -DUSE_HIPFFT -DUSE_NETCDF)
  set(FFTW_LIB "-lhipfft -lamdhip64 -lfftw3")
else ()
  add_definitions(-DNO_ASSUMED_RANKS -DUSE_FFTW -DUSE_NETCDF)
  set(FFTW_LIB "-lfftw3")
endif ()
set(NCDF_LIB "-lnetcdff")
set(LIBS     "${NCDF_LIB} ${FFTW_LIB} -lm")

set(GNU_SED  "gsed")

# --- ALTERNATIVE: standalone amdflang-new (no Cray PrgEnv, e.g. AMD RAC-Plano) -----------------------
# If you are NOT on a Cray PrgEnv, `ftn` is unavailable. Point FC at an MPI wrapper that uses
# amdflang-new (you need MPI — a bare amdflang-new compile of the PARALLEL build will fail to link MPI),
# and supply the local FFTW path. Example (adjust paths to your site):
#   set(ENV{FC} "/path/to/mpi/bin/mpifort")        # MPICH/OpenMPI wrapper built against amdflang-new
#   set(FFTW_PATH "/shared/midgard/home/jkostele_nld/")
#   set(FFTW_LIB  "-L${FFTW_PATH}/lib -lfftw3")
#   include_directories("${FFTW_PATH}/include/")
# (Set OMPI_FC / MPICH_FC=amdflang-new in the environment so the wrapper drives the AMD compiler.)
