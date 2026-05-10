FILE(GLOB sourcefiles "*.f90")
add_library(statistics STATIC ${sourcefiles})
add_dependencies(statistics physics mappings operators thermodynamics utils ibm base)

# Force -O0 on all statistics files. Some routines (e.g., avg_phase) may be
# called every step; even though statistics shouldn't modify the flow field,
# we put them at -O0 to fully exclude this category as a culprit.
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)