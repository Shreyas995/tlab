FILE(GLOB sourcefiles "*.f90")
get_directory_property(defs COMPILE_DEFINITIONS)
if ( "USE_MPI" IN_LIST defs )
else()
  FILE(GLOB to_remove "tlab_mpi*.f90") # get the complete path
  list(REMOVE_ITEM sourcefiles ${to_remove})
endif()
add_library(base STATIC ${sourcefiles})

# Force -O0 on all base files. tlab_memory, tlab_arrays, tlab_pointers etc.
# set up shared array layout — corruption here would silently feed bad data
# everywhere downstream.
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)
