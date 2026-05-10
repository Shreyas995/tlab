FILE(GLOB sourcefiles "*.f90")
add_library(thermodynamics STATIC ${sourcefiles})
add_dependencies(thermodynamics fdm base)

# Force -O0 on all thermodynamics files.
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)
