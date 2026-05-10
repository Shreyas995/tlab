FILE(GLOB sourcefiles "*.f90")
add_library(utils STATIC ${sourcefiles})
add_dependencies(utils base)

# APU-kernel files: force -O0 to bypass optimizer bug that produces NaN at first
# timestep in the production turbulent Poisson solve. Toggle by commenting out.
set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/LinearDss.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tlab_transpose.f90
    PROPERTIES COMPILE_FLAGS "-O0"
)
