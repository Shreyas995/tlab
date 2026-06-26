# HLRS HUNTER-GPU (AMD MI300A APU) — AMD compiler (flang / amdflang-new) variant.
#
# This is the AMD-flang counterpart of config/hunter-cray.cmake (HPE Cray CCE, the verified default).
# It produces the SAME -D macros and the SAME build toggles; only the compiler driver and the
# compiler-specific flags differ (Cray `-eZ`/`-h*` -> flang `-cpp`/`-O*`).
#
# HOW TO SELECT THE AMD COMPILER ON HUNTER (Cray/HPE system):
#   Hunter's `ftn` is the Cray compiler *wrapper*; it dispatches to whichever PrgEnv is loaded and
#   provides Cray-MPI + module integration regardless. To make `ftn` drive AMD flang instead of CCE:
#       module swap PrgEnv-cray PrgEnv-amd     # (or: module load PrgEnv-amd)
#   then `ftn` wraps amdflang/amdflang-new and you still get MPI for free. This is why we keep the
#   compiler as `ftn` below (NOT a bare `amdflang-new`, which has no MPI). Do NOT additionally
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

# ACCELERATE=TRUE -> GPU offload. flang needs the offload arch passed EXPLICITLY (--offload-arch=gfx942);
# the Cray config gets it from the loaded craype-accel module instead, hence the difference.
if ( NOT ACCELERATE )
  set(ACCELERATE "FALSE")
elseif( ${ACCELERATE} STREQUAL "TRUE" )
  set(USER_APU_FLAGS "-fopenmp --offload-arch=gfx942")
  add_definitions(-DUSE_APU)
endif()

# HEISENBUG ISOLATION: -DPROBESYNC=TRUE makes DNS_PRINT_MAXVAL do ONLY hipDeviceSynchronize() (no array
# read, no log). Default off. (Identical semantics to the Cray config.)
if (PROBESYNC)
  add_definitions(-DPROBE_SYNC_ONLY)
endif ()

# I-transpose node-window coherence fixes (the documented cross-rank GPU-push race). Mutually exclusive;
# default off leaves the current node-window path unchanged. For a production multinode run you almost
# certainly want ONE of these (see CLAUDE.md "Production blow-up"):
#   -DTRP_I_FORCE_MPI=TRUE  : always-correct STABILIZER (CPU/MPI fallback, ~1.5-2x slower on the I leg).
#   -DTRP_I_SYSFENCE=TRUE   : separate hip_system_fence after each push + reader invalidate (~baseline).
#   -DTRP_I_FUSEDFENCE=TRUE : V1 fused write+__threadfence_system + reader invalidate (current candidate).
if (TRP_I_FORCE_MPI)
  add_definitions(-DTRP_I_FORCE_MPI)
endif ()
if (TRP_I_SYSFENCE)
  add_definitions(-DTRP_I_SYSFENCE)
endif ()
if (TRP_I_FUSEDFENCE)
  add_definitions(-DTRP_I_FUSEDFENCE)
endif ()

# Compile-time debug probes (NWFR/NWBR/MAXVAL/POIS-trace). Default OFF = zero cost.
if (DNS_DEBUG_PROBES)
  add_definitions(-DDNS_DEBUG_PROBES)
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
