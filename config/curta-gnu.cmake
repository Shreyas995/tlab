# ZEDAT CURTA - AMD (FU Berlin) 

if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE PARALLEL)
endif()

if ( NOT HYBRID ) 
   set(HYBRID FALSE) 
else()
   message(WARNING "Compiling for hybrid openMP/MPI usage") 
endif() 

if ( NOT PROFILE ) 
   set(PROFILE FALSE) 
endif() 

if ( ${PROFILE} STREQUAL "TRUE" )  
   set(USER_PROFILE_FLAGS "-pg")
endif() 

# compiler for parallel build and hybrid flags	  
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" OR ${BUILD_TYPE} STREQUAL "NONBLOCKING" )
   set(ENV{FC} mpif90) 

   add_definitions(-DUSE_MPI -DUSE_MPI_IO -DUSE_ALLTOALL)

   if ( ${HYBRID} STREQUAL "TRUE" )
     set(USER_omp_FLAGS " -fopenmp " )
     add_definitions(-DUSE_OPENMP) 
   endif()

# compiler for serial build
else( ${BUILD_TYPE} STREQUAL "SERIAL" )  
  set(ENV{FC} gfortran)
  set(CMAKE_Fortran_COMPILER gfortran)
endif()     

set(NCDF_INC "-I/trinity/shared/easybuild/arch/x86_64/amd/zen3/software/netCDF-Fortran/4.6.0-gompi-2022a/include" )
set(FFTW_INC "-I/trinity/shared/easybuild/arch/x86_64/amd/zen3/software/FFTW/3.3.10-GCC-11.3.0/include" ) 

set(USER_Fortran_FLAGS         "${USER_PROFILE_FLAGS} -cpp -ffree-form -ffree-line-length-2048 -fno-automatic -fallow-argument-mismatch ${NCDF_INC} ${FFTW_INC}")
set(USER_Fortran_FLAGS_RELEASE "-O3 -march=znver3 -mtune=znver3 -ffinite-math-only -fprefetch-loop-arrays")
# uncommented -fprefetch-loop-arrays --param prefetch-latency=300
# Do not use -funroll-all-loops / -funroll-loops

set(USER_Fortran_FLAGS_DEBUG   "-g -traceback -debug all -ffpe-trap=all") 

if ( NOT CMAKE_BUILD_TYPE ) 
  set(CMAKE_BUILD_TYPE RELEASE)  
endif() 

add_definitions(-DUSE_FFTW  -DZEDAT_AMD) # -DUSE_BLAS -DUSE_MKL -DNO_ASSUMED_RANKS -DUSE_MKL)
set(FFTW_LIB "-lfftw3")
set(NCDF_LIB "-lnetcdff") 
set(LIBS     "${FFTW_LIB} -lm")
#set(BLAS_LIB "-L/trinity/shared/easybuild/arch/x86_64/amd/zen3/software/FlexiBLAS/3.2.0-GCC-11.3.0/lib64 -lflexiblas")
#set(MKL_LIB  "-lmkl")

set(GNU_SED "gsed")

if ( ${BUILD_TYPE} STREQUAL "NONBLOCKING" ) # use NB3DFFT library
  add_definitions(-DUSE_PSFFT)

  set(PSFFT_INCLUDE_DIR   "$ENV{HOME}/nb3dFFT/src")
  set(PSFFT_INCLUDE_DIR   ${PSFFT_INCLUDE_DIR} "$ENV{HOME}/nb3dFFT/include")
  set(PSFFT_LIB           "$ENV{HOME}/nb3dFFT/lib/libnb3dfft.a")
  set(PSFFT_LINK_FLAGS    "-qopenmp")
  set(PSFFT_COMPILE_FLAGS "-qopenmp")

endif()
