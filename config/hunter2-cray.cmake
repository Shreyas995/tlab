#----------------------------
# Default build type
#----------------------------
if(NOT BUILD_TYPE)
    message(WARNING "No BUILD_TYPE provided. Defaulting to PARALLEL.")
    set(BUILD_TYPE "PARALLEL")
endif()

#----------------------------
# Optional features
#----------------------------
if(NOT ACCELERATE)
    set(ACCELERATE "FALSE")
endif()

if(NOT PROFILE)
    set(PROFILE "FALSE")
endif()

if(NOT HYBRID)
    set(HYBRID "FALSE")
endif()

#----------------------------
# Initialize flags
#----------------------------
set(USER_Fortran_FLAGS "")
set(USER_Fortran_FLAGS_DEBUG "")
set(USER_Fortran_FLAGS_RELEASE "")
set(USER_OMP_FLAGS "")
set(USER_APU_FLAGS "")
set(USER_profile_FLAGS "")

# Profile flags
if(PROFILE STREQUAL "TRUE")
    set(USER_profile_FLAGS "-g -h profile_generate")
endif()

# APU/OpenMP flags
if(ACCELERATE STREQUAL "TRUE")
    set(USER_APU_FLAGS "-fopenmp")
    add_definitions(-DUSE_APU)
endif()

# Debug flags (no optimization)
set(USER_Fortran_FLAGS_DEBUG "-O0 -g -Rb -Rbc -Rsv -h bounds -h msgs -h stop_on_fpe -h fp_trap=invalid,zero,overflow")

# Release optimization flags
set(USER_Fortran_FLAGS_RELEASE "-hipa2 -hfp2 -hunroll2 -hfusion2 -hscalar1 -m4")

# OpenMP for hybrid builds
if(HYBRID STREQUAL "TRUE")
    set(USER_OMP_FLAGS "-fopenmp")
    add_definitions(-DUSE_OPENMP)
endif()

#----------------------------
# Compiler selection and flags
#----------------------------
set(CMAKE_Fortran_COMPILER ftn CACHE FILEPATH "Cray Fortran wrapper")

if(BUILD_TYPE STREQUAL "PARALLEL")
    message(WARNING "Compiling in PARALLEL mode (MPI enabled)")
    add_definitions(-DUSE_MPI -DUSE_MPI_IO -DUSE_ALLTOALL)

    set(USER_Fortran_FLAGS "${USER_Fortran_FLAGS_RELEASE} ${USER_OMP_FLAGS} ${USER_APU_FLAGS} ${USER_profile_FLAGS}")

elseif(BUILD_TYPE STREQUAL "SERIAL")
    message(WARNING "Compiling in SERIAL mode (no MPI)")
    set(USER_Fortran_FLAGS "${USER_Fortran_FLAGS_RELEASE} ${USER_OMP_FLAGS} ${USER_APU_FLAGS} ${USER_profile_FLAGS}")

elseif(BUILD_TYPE STREQUAL "DEBUG")
    message(WARNING "Compiling in DEBUG mode (no optimization)")
    set(USER_Fortran_FLAGS "${USER_Fortran_FLAGS_DEBUG} ${USER_OMP_FLAGS} ${USER_APU_FLAGS} ${USER_profile_FLAGS}")

else()
    message(FATAL_ERROR "Unknown BUILD_TYPE: ${BUILD_TYPE}. Allowed: PARALLEL, SERIAL, DEBUG")
endif()

#----------------------------
# Libraries
#----------------------------
add_definitions(-DNO_ASSUMED_RANKS -DUSE_FFTW -DUSE_NETCDF)
set(FFTW_LIB "-lfftw3")
set(NCDF_LIB "-lnetcdff") 
set(LIBS "${NCDF_LIB} ${FFTW_LIB} -lm")

# For debugging
message(STATUS "Final Fortran flags: ${USER_Fortran_FLAGS}")
