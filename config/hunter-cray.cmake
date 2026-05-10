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
# Test 7: -O1 only. RESULT: Crashed (NaN at first timestep after restart).
# ============================================================================
# set(USER_Fortran_FLAGS_RELEASE "-O1 -m4")

# ============================================================================
# Test 9: APU kernels forced to -O0 via per-file CMakeLists.txt overrides.
# RESULT: Crashed (NaN at first timestep). Bug is NOT in just the four kernels.
#
# Files at -O0 in Test 9 (kept):
#     src/utils/LinearDss.f90       (PENTADSS_APU + friends)
#     src/utils/tlab_transpose.f90  (TLab_Transpose_COMPLEX_APU)
#     src/fdm/fdm_matmul.f90        (MatMul_3d_APU, MatMul_5d_APU)
#     src/fdm/fdm_integral.f90      (FDM_Int2_Solve_APU)
# ============================================================================

# ============================================================================
# Test 10: 7 files at -O0 (kernels + Poisson driver + MPI transpose).
# RESULT: Crashed. Bug is OUTSIDE the Poisson solver path.
# ============================================================================

# ============================================================================
# Test 11: 17 files (all with !$omp target). RESULT: Crashed but PROGRESS!
#   - dt, D#, visc are now correct (no longer NaN).
#   - CFL = 0 (suspicious: implies velocity field is zero).
#   - DilMin/DilMax = NaN (only the dilatation calc is broken now).
#   The bug is narrowed. Most likely culprit: opr_partial.f90 (computes ALL
#   derivatives every step; not in the -O0 list because no target regions).
#
# ============================================================================
# Test 12: All src/operators + src/utils at -O0 (43 files total).
# RESULT: dt+D# correct, CFL=0, dilatation=±Infinity (was NaN before).
# Velocity field zeroed except for isolated huge spikes. Bug is in code
# that runs every RK substage but isn't covered yet.
#
# Key insight: only the FIRST RK substage files (_1.f90) were at -O0.
# Substages 2 and 3 (rhs_*_2.f90, rhs_*_3.f90) were still at -O2.
# ============================================================================

# ============================================================================
# Test 13: 88 files (Test 12 + all of src/tools/dns + src/base).
# RESULT: SAME as Test 12. CFL=0, dilatation=±Inf. Bug is NOT in any of
# the per-timestep RK / boundary / halo / pointer-setup files.
# ============================================================================

# ============================================================================
# Test 14 (ACTIVE): Adds src/mappings + src/filters + src/statistics at -O0.
# Filters are the prime suspect: a buggy dealiasing filter applied every step
# would zero out velocity (matching CFL=0) and produce Inf at certain
# wavenumbers (matching DilMax=±Inf).
#
#   - If this WORKS  -> bug is in mappings, filters, or statistics.
#                       Bisect within these directories.
#   - If this CRASHES-> bug is in src/particles or src/thermodynamics
#                       (or some module file we haven't located yet).
# ============================================================================
set(USER_Fortran_FLAGS_RELEASE "-O2 -m4")

if ( NOT CMAKE_BUILD_TYPE ) 
  set(CMAKE_BUILD_TYPE RELEASE)  
endif() 

add_definitions(-DNO_ASSUMED_RANKS -DUSE_FFTW -DUSE_NETCDF) # -DHLRS_HAWK -DUSE_BLAS -DUSE_MKL)
set(FFTW_LIB "-lfftw3")
set(NCDF_LIB "-lnetcdff") 
set(LIBS     "${NCDF_LIB} ${FFTW_LIB} -lm")

set(GNU_SED  "gsed")
