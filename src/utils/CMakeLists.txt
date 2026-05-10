FILE(GLOB sourcefiles "*.f90")
add_library(utils STATIC ${sourcefiles})
add_dependencies(utils base)

# Force -O0 on ALL utils files. linear3/5/7/9, averages, etc. are all
# called every timestep via various paths.
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)
