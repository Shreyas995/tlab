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

if ( NOT ACCELERATE )
  set(ACCELERATE "FALSE")
elseif( ${ACCELERATE} STREQUAL "TRUE" )
  set(USER_APU_FLAGS "-fopenmp ")#-L/opt/rh/gcc-toolset-12/root/usr/lib/gcc/x86_64-redhat-linux/12")
  add_definitions(-DUSE_APU)
endif()

# HEISENBUG ISOLATION: -DPROBESYNC=TRUE makes DNS_PRINT_MAXVAL do ONLY hipDeviceSynchronize() (no array
# read, no log) instead of the full GPU reduction. Used to test whether the sentinels mask the crash via
# the stream sync (=> missing flush) or via the full device read (=> host->GPU coherency / UB). Default off.
if (PROBESYNC)
  add_definitions(-DPROBE_SYNC_ONLY)
endif ()

# I-transpose node-window corruption fixes (fabricdirect; the cross-rank GPU push delivers stale data because
# hipDeviceSynchronize is device-scope, not a system-scope L2 write-back). Mutually exclusive; default off
# leaves the current node-window path unchanged.
#   -DTRP_I_FORCE_MPI=TRUE : TEMP STABILIZER. Forces the 4 I-transposes onto the all-MPI fallback (correct,
#                            ~1.5-2x slower). Keeps the sim running.
#   -DTRP_I_SYSFENCE=TRUE  : FAST FIX. Adds a system-scope fence (hip_system_fence) after each node-window I
#                            push, committing the write to MALL before the close fence. ~baseline speed.
if (TRP_I_FORCE_MPI)
  add_definitions(-DTRP_I_FORCE_MPI)
endif ()
if (TRP_I_SYSFENCE)
  add_definitions(-DTRP_I_SYSFENCE)
endif ()

# compiler for parallel build	  
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" )
   set(ENV{FC} ftn) # instead of running "export FC=ftn" in the terminal
   add_definitions(-DUSE_MPI -DUSE_MPI_IO -DUSE_ALLTOALL) # -DUSE_NETCDF (already later defined)
  
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
