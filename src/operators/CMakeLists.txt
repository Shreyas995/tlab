FILE(GLOB sourcefiles "*.f90")
add_library(operators STATIC ${sourcefiles})
add_dependencies(operators ibm filters fdm base)

# Force -O0 on ALL operator files. opr_partial in particular is called every
# timestep for ∂u/∂x etc.; if miscompiled at -O2, all derivatives are wrong.
set_source_files_properties(
    ${sourcefiles}
    PROPERTIES COMPILE_FLAGS "-O0"
)