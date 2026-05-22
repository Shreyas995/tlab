FILE(GLOB sourcefiles "*.f90")
get_directory_property(defs COMPILE_DEFINITIONS)
if ( "USE_MPI" IN_LIST defs )
else()
  FILE(GLOB to_remove "tlab_mpi*.f90") # get the complete path
  list(REMOVE_ITEM sourcefiles ${to_remove})
endif()
add_library(base STATIC ${sourcefiles})

if ( "USE_APU" IN_LIST defs )
  set(HIP_FENCE_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../valid/mpi/hip_write_fence.hip")
  set(HIP_FENCE_OBJ "${CMAKE_CURRENT_BINARY_DIR}/hip_write_fence.o")
  add_custom_command(
    OUTPUT  ${HIP_FENCE_OBJ}
    COMMAND hipcc --offload-arch=gfx942 -O2 -fPIC -c ${HIP_FENCE_SRC} -o ${HIP_FENCE_OBJ}
    DEPENDS ${HIP_FENCE_SRC}
    COMMENT "Compiling hip_write_fence.hip"
  )
  add_custom_target(hip_fence_obj DEPENDS ${HIP_FENCE_OBJ})
  add_dependencies(base hip_fence_obj)
  # Insert the compiled HIP object into the static archive so all executables
  # that link base (dns.x, inirand.x, ...) resolve hip_write_with_fence.
  add_custom_command(TARGET base POST_BUILD
    COMMAND ${CMAKE_AR} r $<TARGET_FILE:base> ${HIP_FENCE_OBJ}
    COMMENT "Adding hip_write_fence.o to libbase.a"
  )
  # Propagate -lamdhip64 transitively to every consumer of base.
  target_link_libraries(base INTERFACE -lamdhip64)
endif()
