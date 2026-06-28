# HPE Cray Compiler (AAC7 Plano cluster with APUs: AMD Instinct MI300A Accelerator) 

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
   set(USER_profile_FLAGS "-g ")#-h profile_generate")
endif()

# ============================================================================
#  GPU / APU acceleration  (-DACCELERATE=TRUE)  --  ALL TRANSPOSES RUN ON THE GPU
# ============================================================================
# A plain -DACCELERATE=TRUE gives the COMPLETE, coherent, all-on-GPU build: every
# transpose (real I/K and complex I/K) runs on the GPU shared window, with the
# per-workgroup write+fence coherence fix folded in automatically -- it is REQUIRED
# for a correct GPU run, so there is no separate flag to remember (was -DTRP_I_FUSEDFENCE).
#
# DESIGN NOTE -- why everything stays on the GPU (do NOT "optimize" by moving transposes
# to the CPU/MPI): the original apudirect was fast NOT because the complex/Poisson
# transposes ran on the CPU -- they did not, every transpose ran on the GPU shared
# window -- but because it carried NO per-push system fences or synchronizations. The
# fences added later are a CORRECTNESS tax: the cross-rank node-window push needs a
# system-scope __threadfence_system + a reader-side L2 invalidate, or it delivers stale
# data and the run blows up. So the optimization target is a CHEAPER fence (the fused
# per-workgroup kernel), NEVER routing transposes through the CPU. On this unified APU a
# CPU/MPI transpose is always slower; the all-GPU fenced path IS the production path.
if ( NOT ACCELERATE )
  set(ACCELERATE "FALSE")
elseif( ${ACCELERATE} STREQUAL "TRUE" )
  set(USER_APU_FLAGS "-fopenmp ")#-L/opt/rh/gcc-toolset-12/root/usr/lib/gcc/x86_64-redhat-linux/12")
  add_definitions(-DUSE_APU)

  # Default GPU-transpose coherence fix = the per-workgroup fused write+fence (real I/K +
  # complex-K) + reader L2 invalidate. Always on for an APU build (see the note above).
  set(_TRP_FUSED_DEFAULT TRUE)

  # ----- Advanced / experimental transpose overrides (default OFF; A/B diagnosis only) -----
  # Every one of these is SLOWER than, or equal to, the all-GPU default -- keep them OFF for
  # production. Kept here only so a failing GPU run can be bisected against a CPU/MPI leg.
  if (TRP_I_SYSFENCE)         # LEGACY separate hip_system_fence kernel; superseded by the
    add_definitions(-DTRP_I_SYSFENCE)    # fused fence and mutually exclusive with it.
    set(_TRP_FUSED_DEFAULT FALSE)
  endif ()
  if (_TRP_FUSED_DEFAULT)
    add_definitions(-DTRP_I_FUSEDFENCE)
  endif ()
  if (TRP_I_MEMCPY)           # K pushes via a blocking system-coherent hipMemcpy instead of
    add_definitions(-DTRP_I_MEMCPY)      # the fused fence (combines with the fused I fence).
  endif ()
  if (TRP_I_FORCE_MPI)        # route the single-XCD intra-node I-transpose onto MPI; legacy
    add_definitions(-DTRP_I_FORCE_MPI)   # stabilizer, slower than the GPU fence.
  endif ()
  if (TRP_CX_MPI)             # route the complex (Poisson) transposes onto MPI; slower, A/B
    add_definitions(-DTRP_CX_MPI)        # only (real transposes stay on the GPU).
  endif ()
endif()

# ============================================================================
#  Debug probe LEVEL  (-DDNS_DEBUG=<n>;  default 0 = OFF, zero production cost)
# ============================================================================
# Crash-localization probes (NWFR/NWBR/MAXVAL/POIS-trace -> per-rank fort.5xx). When you
# debug you PRINT ALL: one level switch turns the whole probe family on (this replaces the
# two separate -DDNS_DEBUG_PROBES and -DPROBESYNC flags).
#   0 (default) : OFF. Probes compile out to `continue` -- no GPU reduction, no I/O.
#   1           : full probes ("print all") -- GPU/CPU max+nbad reductions to fort.5xx.
#   2           : Heisenbug isolation -- probes do ONLY hipDeviceSynchronize (no read/log);
#                 the exception used to test whether the full device read masks the bug.
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
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" )
   set(ENV{FC} ftn) # instead of running "export FC=ftn" in the terminal
   add_definitions(-DUSE_MPI -DUSE_MPI_IO -DUSE_ALLTOALL) # -DUSE_NETCDF (already later defined)
   # NOTE: debug probes are toggled ONLY by the `-DDNS_DEBUG=<n>` level block above (1=on, 2=Heisenbug).
   # Do NOT hard-add -DDNS_DEBUG_PROBES here — that would force probes always-on and defeat the toggle.
  
# OpenMP (hybrid) flags
   if ( ${HYBRID} STREQUAL "TRUE" )
     set(USER_OMP_FLAGS " -fopenmp")
     add_definitions(-DUSE_OPENMP) 
   endif()

# compiler for serial build
else( ${BUILD_TYPE} STREQUAL "SERIAL" )
  set(ENV{FC} ftn )
endif()     

# set(DRAGONEGG_FLAGS "-finline-aggressive -fslp-vectorize  -fmerge-all-constants") #  -mmadd4 -mfp64 -enable-strided-vectorization")
set(USER_Fortran_FLAGS         "-eZ ${USER_OMP_FLAGS} ${USER_APU_FLAGS} ${USER_profile_FLAGS} ") #-fallow-argument-mismatch from gnu-version10

# Production optimization flags for Hunter (Cray CCE 20.0.0, MI300A).
# Debug alternative: "-g -O0 -R abc -m4"
set(USER_Fortran_FLAGS_RELEASE "-hipa2 -hfp2 -hunroll2 -hfusion2 -hscalar1 -m4")

if ( NOT CMAKE_BUILD_TYPE ) 
  set(CMAKE_BUILD_TYPE RELEASE)  
endif() 

# FFT backend: FFTW (default; CPU + the GPU bit-reference) or hipFFT (GPU X/Z FFTs). Mutually exclusive
# for the OPR_Fourier module (opr_fourier.f90 vs opr_fourier_hipfft.f90). hipFFT still links FFTW because
# the Y direction (used only by OPR_Fourier_F/B) stays on CPU FFTW. Select hipFFT with -DHIPFFT=TRUE.
if (HIPFFT)
  add_definitions(-DNO_ASSUMED_RANKS -DUSE_HIPFFT -DUSE_NETCDF)
  set(FFTW_LIB "-lhipfft -lamdhip64 -lfftw3")
else ()
  add_definitions(-DNO_ASSUMED_RANKS -DUSE_FFTW -DUSE_NETCDF) # -DHLRS_HAWK -DUSE_BLAS -DUSE_MKL)
  set(FFTW_LIB "-lfftw3")
endif ()
set(NCDF_LIB "-lnetcdff")
set(LIBS     "${NCDF_LIB} ${FFTW_LIB} -lm")

set(GNU_SED  "gsed")
