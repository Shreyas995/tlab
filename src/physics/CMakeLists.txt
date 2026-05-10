FILE(GLOB sourcefiles "*.f90")
add_library(physics STATIC ${sourcefiles})
add_dependencies(physics mappings operators thermodynamics utils base)

# APU target regions in physics: force -O0 for the optimization-bug bisection.
set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/opr_burgers.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/rotation.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tlab_sources.f90
    PROPERTIES COMPILE_FLAGS "-O0"
)