FILE(GLOB sourcefiles "*.f90")
add_library(operators STATIC ${sourcefiles})
add_dependencies(operators ibm filters fdm base)

# APU driver-layer files: force -O0. These coordinate the APU kernels and
# compute the inputs (lambda, transposes, FFTs) that feed into them.
set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/opr_elliptic.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/opr_fourier.f90
    PROPERTIES COMPILE_FLAGS "-O0"
)