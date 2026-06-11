FILE(GLOB sourcefiles "*.f90")
add_library(fdm STATIC ${sourcefiles})
add_dependencies(fdm base)
# fdm_integral.f90 'use LinearDss' (added with the GPU elliptic Y-solve) pulls a
# module from the utils library, so utils must build first. Without this, a clean
# build can compile fdm before utils and fail with "No rule to make target
# modules/lineardss.mod". (Hunter rebuilds via cmake each time and masks this;
# the local CPU `make`-only flow does not.)
add_dependencies(fdm utils)