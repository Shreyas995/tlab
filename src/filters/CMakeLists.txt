FILE(GLOB sourcefiles "*.f90")
add_library(filters STATIC ${sourcefiles})
add_dependencies(filters base)

# Force -O0 on all filter files. Dealiasing filters are applied every timestep
# to the nonlinear term — a buggy filter could zero out velocity (CFL=0) and
# produce Inf at specific wavenumbers (DilMax=±Inf), matching our symptoms.
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)
