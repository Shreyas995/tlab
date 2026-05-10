FILE(GLOB sourcefiles "*.f90")
add_library(physics STATIC ${sourcefiles})
add_dependencies(physics mappings operators thermodynamics utils base)

# Force -O0 on all physics files (gravity, chemistry, largescale, etc.).
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)