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

# ============================================================================
# Optimization flag tests (the LAST uncommented set() wins; earlier ones are overridden)
# ============================================================================
# Test 1: Baseline. Works correctly. Slow.
# set(USER_Fortran_FLAGS_RELEASE "-g -O0 -R abc -m4")

# Test 2: Light scalar opt. Explodes (production turbulent case).
# set(USER_Fortran_FLAGS_RELEASE "-O2 -hscalar1 -m4")

# Test 3: + loop unroll/fusion. Explodes.
# set(USER_Fortran_FLAGS_RELEASE "-O2 -hscalar1 -hunroll2 -hfusion2 -m4")

# Test 4: + IPA. Explodes.
# set(USER_Fortran_FLAGS_RELEASE "-O2 -hscalar1 -hunroll2 -hfusion2 -hipa2 -m4")

# Test 5: + IEEE-equivalent FP. Explodes (rules out FP-precision as the cause).
# set(USER_Fortran_FLAGS_RELEASE "-O2 -hscalar1 -hunroll2 -hfusion2 -hipa2 -hfp1 -m4")

# Test 6: Aggressive (original target). Explodes.
# set(USER_Fortran_FLAGS_RELEASE "-hipa2 -hfp2 -hunroll2 -hfusion2 -hscalar1 -m4")

# ============================================================================
# Test 7 (ACTIVE): -O1 only.
# -O1 does constant folding, dead-code elim, light inlining. NO vectorization,
# NO loop transformations, NO IPA across the host/device boundary.
# This is the cleanest "is the bug -O2-specific?" test using only standard flags.
#   - If this WORKS  -> the bug is in some -O2 transformation. Vectorization
#                       is the prime suspect for an iterative Poisson solver.
#                       Next step: -O2 -h vector0 to confirm.
#   - If this CRASHES-> bug is in even basic optimization. Then we move to
#                       per-subroutine !dir$ optimize(0) directives on the
#                       APU kernels to localize the culprit.
# ============================================================================
set(USER_Fortran_FLAGS_RELEASE "-O1 -m4")

# ============================================================================
# Test 8 (FALLBACK if Test 7 crashes): -O2 with vectorization disabled.
# If -O1 crashes too, skip to Test 9 instead.
# ============================================================================
# set(USER_Fortran_FLAGS_RELEASE "-O2 -h vector0 -m4")

# ============================================================================
# Test 9 (FALLBACK if -O1 also crashes): per-subroutine optimization control.
# Use Cray directive !dir$ optimize(0) inside specific .f90 files to disable
# optimization on individual APU kernels. Build with -O2 globally.
# ============================================================================
# set(USER_Fortran_FLAGS_RELEASE "-O2 -m4")

if ( NOT CMAKE_BUILD_TYPE ) 
  set(CMAKE_BUILD_TYPE RELEASE)  
endif() 

add_definitions(-DNO_ASSUMED_RANKS -DUSE_FFTW -DUSE_NETCDF) # -DHLRS_HAWK -DUSE_BLAS -DUSE_MKL)
set(FFTW_LIB "-lfftw3")
set(NCDF_LIB "-lnetcdff") 
set(LIBS     "${NCDF_LIB} ${FFTW_LIB} -lm")

set(GNU_SED  "gsed")
