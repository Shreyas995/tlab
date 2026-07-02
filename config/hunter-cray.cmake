# HPE Cray Compiler (AAC7 Plano cluster with APUs: AMD Instinct MI300A Accelerator) 

if ( NOT BUILD_TYPE )
   message( WARNING "Setting CMAKE_BUILD_TYPE to default value." )
   set(BUILD_TYPE PARALLEL)
endif()

if ( NOT HYBRID ) 
   set(HYBRID FALSE) 
else()
   message(WARNING "Compiling for hybrid openMP/MPI usage (Only in combination! Not tested!)") 
endif() 

if ( NOT PROFILE ) 
   set(PROFILE FALSE) 
endif()

if (${PROFILE} STREQUAL "TRUE" )  
   set(USER_profile_FLAGS "-g ")#-h profile_generate")
endif()

# ============================================================================
#  GPU / APU acceleration  (-DACCELERATE=TRUE)  --  ALL TRANSPOSES RUN ON THE GPU
# ============================================================================
# A plain -DACCELERATE=TRUE gives the COMPLETE, coherent, all-on-GPU build: every
# transpose (real I/K and complex I/K) runs on the GPU shared window, with the node-window
# cross-XCD push coherence fix folded in automatically -- it is REQUIRED for a correct GPU
# run, so there is no separate flag to remember (was -DTRP_I_FUSEDFENCE, then -DTRP_I_MEMCPY).
# The fix lives entirely in #ifdef USE_APU code, which compiles ONLY under ACCELERATE, so it
# belongs here -- nothing to merge or pass on the command line.
#
# Two mechanisms, both folded in:
#   * TRP_I_MEMCPY  (writer commit, DEFAULT) -- a BLOCKING, system-coherent hipMemcpy for every
#       node-window cross-XCD push (real I/K + complex I/K). Committed to MALL on return, so it
#       does NOT rely on a manual fence. The per-workgroup fence only CUT the writer-side push
#       race ~14x, never closed it (the _750 it=236601 real-I blow-up: held ~1600 iters then a
#       1-iter detonation at CFL 0.516). hipMemcpy closes it. ~0.01-0.05 s/iter, still on-GPU.
#   * TRP_I_FUSEDFENCE (per-workgroup write+fence + reader L2 invalidate) -- kept ON: it provides
#       the staging/reader-acquire scaffolding the memcpy reuses, AND is the commit for the
#       apudirect-only paths that have no memcpy branch (hip_pushseg_fence).
#
# DESIGN NOTE -- why everything stays on the GPU (do NOT "optimize" by moving transposes to the
# CPU/MPI): the original apudirect was fast NOT because the complex/Poisson transposes ran on the
# CPU -- they did not, every transpose ran on the GPU shared window -- but because it carried NO
# per-push system fences or synchronizations. The coherence machinery added later is a CORRECTNESS
# tax: the cross-rank node-window push needs a system-scope commit (the blocking hipMemcpy, or a
# __threadfence_system) + a reader-side L2 invalidate, or it delivers stale data and the run blows
# up. The optimization target is a CHEAPER on-GPU commit, NEVER routing transposes through the CPU.
# On this unified APU a CPU/MPI transpose is always slower; the all-GPU path IS the production path.
if ( NOT ACCELERATE )
  set(ACCELERATE "FALSE")
elseif( ${ACCELERATE} STREQUAL "TRUE" )
  set(USER_APU_FLAGS "-fopenmp ")#-L/opt/rh/gcc-toolset-12/root/usr/lib/gcc/x86_64-redhat-linux/12")
  add_definitions(-DUSE_APU)

  # Default GPU-transpose coherence fix, BOTH folded into -DACCELERATE (see the note above):
  #   _TRP_MEMCPY_DEFAULT -> blocking system-coherent hipMemcpy on every node-window push (the fix).
  #   _TRP_FUSED_DEFAULT  -> per-workgroup write+fence + reader L2 invalidate (scaffolding + apudirect).
  set(_TRP_FUSED_DEFAULT TRUE)
  set(_TRP_MEMCPY_DEFAULT TRUE)

  # ----- Advanced / experimental transpose overrides (default OFF; A/B diagnosis only) -----
  # Every one of these is SLOWER than, or equal to, the all-GPU default -- keep them OFF for
  # production. Kept here only so a failing GPU run can be bisected against a fence/CPU/MPI leg.
  if (TRP_I_FENCE_ONLY)       # A/B: drop the hipMemcpy, keep the per-wg fence (to measure the memcpy
    set(_TRP_MEMCPY_DEFAULT FALSE)        # cost / show it closes the residual the fence leaks).
  endif ()
  if (TRP_I_SYSFENCE)         # LEGACY weak separate hip_system_fence kernel (A/B only); disables BOTH
    add_definitions(-DTRP_I_SYSFENCE)    # GPU defaults so the weak raw-push path is actually exercised.
    set(_TRP_FUSED_DEFAULT FALSE)
    set(_TRP_MEMCPY_DEFAULT FALSE)
  endif ()
  if (_TRP_MEMCPY_DEFAULT)
    add_definitions(-DTRP_I_MEMCPY)      # writer commit = blocking hipMemcpy (real I/K + complex I/K).
  endif ()
  if (_TRP_FUSED_DEFAULT)
    add_definitions(-DTRP_I_FUSEDFENCE)  # staging/reader-acquire scaffolding + apudirect fence commit.
  endif ()
  if (TRP_I_MEMCPY_ASYNC)     # A/B (F1): pipeline the node-window pushes -- N hipMemcpyAsync on ONE stream + 1
    add_definitions(-DTRP_I_MEMCPY -DTRP_I_MEMCPY_ASYNC)   # stream-sync, vs N blocking hipMemcpy. Implies MEMCPY.
  endif ()                                                 # Default OFF; do NOT combine with FENCE_ONLY/SYSFENCE.
  if (TRP_I_FORCE_MPI)        # route the single-XCD intra-node I-transpose onto MPI; legacy stabilizer,
    add_definitions(-DTRP_I_FORCE_MPI)   # now unnecessary -- hipMemcpy fixes the real-I push on the GPU.
  endif ()
  if (TRP_CX_MPI)             # route the complex (Poisson) transposes onto MPI; slower, A/B
    add_definitions(-DTRP_CX_MPI)        # only (real transposes stay on the GPU).
  endif ()

  # apudirect push optimization (2026-07-01, DEFAULT ON): replace the per-peer writer-release loop
  # (N blocking hipMemcpy) AND the per-workgroup __threadfence_system push (hip_pushseg_fence /
  # hip_write_with_fence, the profiler's 46%-of-GPU hot spot) with ONE strided coherent hipMemcpy2D
  # per apudirect transpose (real + complex). Same MALL commit (blocking, system-coherent on return),
  # single host round-trip. apudirect-only (#ifdef USE_APU + trp_mode==APU_DIRECT); fabricdirect paths
  # untouched. A/B: -DTRP_APU_NO_MEMCPY2D=TRUE reverts to the per-peer memcpy/fence push.
  set(_TRP_APU_MEMCPY2D_DEFAULT TRUE)
  if (NOT _TRP_MEMCPY_DEFAULT)   # FENCE_ONLY / SYSFENCE want the per-wg fence on ALL apudirect routines,
    set(_TRP_APU_MEMCPY2D_DEFAULT FALSE)   # so the coherent-copy 2D path follows the memcpy default off.
  endif ()
  if (TRP_APU_NO_MEMCPY2D)
    set(_TRP_APU_MEMCPY2D_DEFAULT FALSE)
  endif ()
  if (_TRP_APU_MEMCPY2D_DEFAULT)
    add_definitions(-DTRP_APU_MEMCPY2D)
  endif ()

  # apudirect coherence-trim A/B switches (2026-07-01, default OFF = keep the op). Single-node apudirect is
  # coherent via MPI_Win_fence + the blocking hipMemcpy2D writer release; these two extra ops MAY be redundant
  # on one node (the old bit-reference apudirect had neither). Each is its own flag so a LONG Hunter run can
  # attribute any silent corruption. Only affects apudirect (apu_win_*/apu_recv_fptr_*); fabricdirect untouched.
  if (TRP_APU_NO_INVALIDATE)  # T5: drop the reader-side hip_invalidate_recv (L2 acquire) after the close fence.
    add_definitions(-DTRP_APU_NO_INVALIDATE)
  endif ()
  if (TRP_APU_NO_PRESYNC)     # T6: drop the pre-push hipDeviceSynchronize that orders the OMP gather/caller
    add_definitions(-DTRP_APU_NO_PRESYNC)   # data before the hipMemcpy2D (OMP target regions are host-synchronous).
  endif ()

  # T7 (2026-07-01, DEFAULT ON): port the apudirect N->1 single-DMA push to the fabricdirect all-intra I
  # node-window. When node_lrank_i is a constant-stride grid (runtime check node_win_i_linear), the per-peer
  # hip_memcpy_push loop over the intra I-peers collapses to ONE strided hipMemcpy2D. Reader invalidate +
  # writer sync + inter-node MPI + K node-window are UNCHANGED (the _750 coherence fix stays). A/B off-switch:
  # -DTRP_FBD_NO_MEMCPY2D=TRUE reverts to the per-peer loop.
  set(_TRP_FBD_MEMCPY2D_DEFAULT TRUE)
  if (TRP_FBD_NO_MEMCPY2D)
    set(_TRP_FBD_MEMCPY2D_DEFAULT FALSE)
  endif ()
  if (_TRP_FBD_MEMCPY2D_DEFAULT)
    add_definitions(-DTRP_FBD_MEMCPY2D)
  endif ()
endif()

# ============================================================================
#  Debug probe LEVEL  (-DDNS_DEBUG=<n>;  default 0 = OFF, zero production cost)
# ============================================================================
# Crash-localization probes (NWFR/NWBR/MAXVAL/POIS-trace -> per-rank fort.5xx). When you
# debug you PRINT ALL: one level switch turns the whole probe family on (this replaces the
# two separate -DDNS_DEBUG_PROBES and -DPROBESYNC flags).
#   0 (default) : OFF. Probes compile out to `continue` -- no GPU reduction, no I/O.
#   1           : full probes ("print all") -- GPU/CPU max+nbad reductions to fort.5xx.
#   2           : Heisenbug isolation -- probes do ONLY hipDeviceSynchronize (no read/log);
#                 the exception used to test whether the full device read masks the bug.
if (NOT DNS_DEBUG)
  set(DNS_DEBUG 0)
endif ()
if (DNS_DEBUG GREATER_EQUAL 1)
  add_definitions(-DDNS_DEBUG_PROBES)
endif ()
if (DNS_DEBUG GREATER_EQUAL 2)
  add_definitions(-DPROBE_SYNC_ONLY)
endif ()

# -------- Transpose leak-hunt build  (-DTRP_LEAK=TRUE) --------
# Single knob for the decisive long debug run that hunts a T1-T7 push-coherence seed: it turns on
#   * DNS_DEBUG_PROBES  -> the modest per-iteration DNS_PROBE family (RHS/POIS/TIME/LOOP) + makes
#                          DNS_PRINT_MAXVAL actually reduce (it early-returns without this).
#   * TRP_LEAK_PROBE    -> the AGGRESSIVE TRP_LEAK probes at in/win/out of EVERY apudirect AND
#                          fabricdirect real+complex transpose push (fires on every call).
# Same binary works for BOTH the apudirect and the fabricdirect decisive runs (only the active branch
# logs). Heavy fort.5xx I/O by design; grep the files for `nbad>0` to find the first NaN in exec order.
# Independent of the DNS_DEBUG level so you don't have to also raise that.
if (TRP_LEAK)
  add_definitions(-DDNS_DEBUG_PROBES -DTRP_LEAK_PROBE)
endif ()

# compiler for parallel build	  
if ( ${BUILD_TYPE} STREQUAL "PARALLEL" )
   set(ENV{FC} ftn) # instead of running "export FC=ftn" in the terminal
   add_definitions(-DUSE_MPI -DUSE_MPI_IO -DUSE_ALLTOALL) # -DUSE_NETCDF (already later defined)
   # NOTE: debug probes are toggled ONLY by the `-DDNS_DEBUG=<n>` level block above (1=on, 2=Heisenbug).
   # Do NOT hard-add -DDNS_DEBUG_PROBES here — that would force probes always-on and defeat the toggle.
  
# OpenMP (hybrid) flags
   if ( ${HYBRID} STREQUAL "TRUE" )
     set(USER_OMP_FLAGS " -fopenmp")
     add_definitions(-DUSE_OPENMP) 
   endif()

# compiler for serial build
else( ${BUILD_TYPE} STREQUAL "SERIAL" )
  set(ENV{FC} ftn )
endif()     

# set(DRAGONEGG_FLAGS "-finline-aggressive -fslp-vectorize  -fmerge-all-constants") #  -mmadd4 -mfp64 -enable-strided-vectorization")
set(USER_Fortran_FLAGS         "-eZ ${USER_OMP_FLAGS} ${USER_APU_FLAGS} ${USER_profile_FLAGS} ") #-fallow-argument-mismatch from gnu-version10

# Production optimization flags for Hunter (Cray CCE 20.0.0, MI300A).
# Debug alternative: "-g -O0 -R abc -m4"
set(USER_Fortran_FLAGS_RELEASE "-hipa2 -hfp2 -hunroll2 -hfusion2 -hscalar1 -m4")

if ( NOT CMAKE_BUILD_TYPE ) 
  set(CMAKE_BUILD_TYPE RELEASE)  
endif() 

# FFT backend: FFTW (default; CPU + the GPU bit-reference) or hipFFT (GPU X/Z FFTs). Mutually exclusive
# for the OPR_Fourier module (opr_fourier.f90 vs opr_fourier_hipfft.f90). hipFFT still links FFTW because
# the Y direction (used only by OPR_Fourier_F/B) stays on CPU FFTW. Select hipFFT with -DHIPFFT=TRUE.
if (HIPFFT)
  add_definitions(-DNO_ASSUMED_RANKS -DUSE_HIPFFT -DUSE_NETCDF)
  set(FFTW_LIB "-lhipfft -lamdhip64 -lfftw3")
else ()
  add_definitions(-DNO_ASSUMED_RANKS -DUSE_FFTW -DUSE_NETCDF) # -DHLRS_HAWK -DUSE_BLAS -DUSE_MKL)
  set(FFTW_LIB "-lfftw3")
endif ()
set(NCDF_LIB "-lnetcdff")
set(LIBS     "${NCDF_LIB} ${FFTW_LIB} -lm")

set(GNU_SED  "gsed")
