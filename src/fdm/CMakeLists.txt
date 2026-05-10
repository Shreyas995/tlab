FILE(GLOB sourcefiles "*.f90")
add_library(fdm STATIC ${sourcefiles})
add_dependencies(fdm base)

# Force -O0 on ALL fdm files. fdm_derivative.f90 contains the actual derivative
# kernels (called by opr_partial every step). fdm_base, fdm_*_jacobian etc.
# are also called from the per-step compute path.
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)