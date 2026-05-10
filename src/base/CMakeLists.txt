FILE(GLOB sourcefiles "*.f90")
get_directory_property(defs COMPILE_DEFINITIONS)
if ( "USE_MPI" IN_LIST defs )
else()
  FILE(GLOB to_remove "tlab_mpi*.f90") # get the complete path
  list(REMOVE_ITEM sourcefiles ${to_remove})
endif()
add_library(base STATIC ${sourcefiles})

# APU MPI transpose: force -O0. Has many !$omp target regions and
# c_f_pointer aliases in the active path.
if ( "USE_MPI" IN_LIST defs )
  set_source_files_properties(
      ${CMAKE_CURRENT_SOURCE_DIR}/tlab_mpi_transpose.f90
      PROPERTIES COMPILE_FLAGS "-O0"
  )
endif()
