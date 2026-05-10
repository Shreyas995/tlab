FILE(GLOB sourcefiles "*.f90")

if ( ${BUILD_TYPE} STREQUAL "PARALLEL" OR ${BUILD_TYPE} STREQUAL "NONBLOCKING" )
else()
FILE(GLOB to_remove "particle_mpi_sort.f90")
list(REMOVE_ITEM sourcefiles ${to_remove})
FILE(GLOB to_remove "particle_mpi_send_recv.f90")
list(REMOVE_ITEM sourcefiles ${to_remove})
endif()

add_library(particles STATIC ${sourcefiles})
add_dependencies(particles thermodynamics base)

# Force -O0 on all particles files.
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)
