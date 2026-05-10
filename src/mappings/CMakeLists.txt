FILE(GLOB sourcefiles "*.f90")
add_library(mappings STATIC ${sourcefiles})
add_dependencies(mappings operators thermodynamics ibm base)

# Force -O0 on all mappings files. fi_gradient, fi_vectorcalculus etc. compute
# derived quantities used in time integration; corruption here propagates.
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)