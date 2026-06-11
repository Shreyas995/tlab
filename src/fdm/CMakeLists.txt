FILE(GLOB sourcefiles "*.f90")
add_library(fdm STATIC ${sourcefiles})
add_dependencies(fdm base utils)   # fdm_integral.f90 `use LinearDss` (utils) -> utils must build first for lineardss.mod