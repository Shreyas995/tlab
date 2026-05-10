FILE(GLOB sourcefiles "../ibm/*.f90")
add_library(ibm STATIC ${sourcefiles})
add_dependencies(ibm fdm utils base)

# APU target region in ibm_bcs.f90: force -O0 for the optimization-bug bisection.
set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/../ibm/ibm_bcs.f90
    PROPERTIES COMPILE_FLAGS "-O0"
)