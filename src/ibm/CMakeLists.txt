FILE(GLOB sourcefiles "../ibm/*.f90")
add_library(ibm STATIC ${sourcefiles})
add_dependencies(ibm fdm utils base)

# Force -O0 on all ibm files (immersed boundary).
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)