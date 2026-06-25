// hip_write_fence.hip
//
// HIP C wrapper providing system-scope fenced writes for the fabricdirect validation.
// Callable from Fortran via iso_c_binding.
//
// Build on Hunter (Cray/ROCm):
//   hipcc -c hip_write_fence.hip -o hip_write_fence.o
//   ftn vmpi_hip_shmwrite.f90 hip_write_fence.o -lamdhip64 -o vmpi_hip_shmwrite.x
//
// What it does vs the current !$omp target approach
// --------------------------------------------------
// The current fabricdirect code uses:
//     !$omp target teams distribute parallel do
//     do i = 1, n
//         dst(i) = src(i)      ! direct write into peer's shared-window segment
//     end do
//     !$omp end target teams distribute parallel do
// followed by MPI_Win_fence on apu_async_win_i (collective on the tainted shmem subcomm).
//
// This file replaces BOTH the !$omp target writes AND the MPI_Win_fence:
//   1. HIP kernel writes src[0:n) to dst[0:n) on the GPU.
//   2. __threadfence_system() inside the kernel flushes the writer's L2 cache to
//      MALL/DRAM (AMD GFX94: buffer_wbl2 sc0:1 sc1:1 + s_waitcnt vmcnt(0)).
//      After this, all agents (other XCDs, CPUs) see the fresh data from MALL.
//   3. hipDeviceSynchronize() on the host ensures the kernel + fence complete before
//      returning to Fortran, so no MPI_Win_fence on any MPI comm is needed for
//      coherency.  Process-level ordering is handled by MPI_Barrier(MPI_COMM_WORLD)
//      in the Fortran caller.
//
// Assumption being tested by vmpi_hip_shmwrite.f90:
//   __threadfence_system() is SUFFICIENT for cross-XCD visibility on MI300A when
//   MPI_Win_allocate_shared maps memory as SVM_COHERENT (fine-grained).
//   If the test shows stale reads, add a reader-side buffer_gl2_inv (see comment
//   in hip_read_invalidate below) and call it before the unpack step.

#include <hip/hip_runtime.h>
#include <stdio.h>

// ------------------------------------------------------------------
// Writer kernel: copy src[0:n) to dst[0:n) then flush to system memory.
// The caller passes a pre-offset dst pointer so this kernel always
// starts at dst[0] — no offset parameter needed here.
// ------------------------------------------------------------------
__global__ void hip_write_fence_kernel(const double* __restrict__ src,
                                        double* __restrict__ dst,
                                        int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dst[idx] = src[idx];
    }
    // System-scope fence: after this all agents in the system see our writes.
    // On AMD GFX94 (MI300A) this compiles to:
    //   s_waitcnt vmcnt(0)        <- drain all in-flight vector stores
    //   buffer_wbl2 sc0:1 sc1:1   <- write-back L2 at system scope
    //   s_waitcnt vmcnt(0)        <- ensure wbl2 completes
    __threadfence_system();
}

// ------------------------------------------------------------------
// Optional reader-side invalidation kernel.
// Call on the RECEIVING rank just before reading its recv window if
// stale L2 data is suspected (i.e., vmpi_hip_shmwrite reports wrong
// values despite the writer flushing correctly).
// In that case, call this on the recv buffer before the verify loop.
// ------------------------------------------------------------------
__global__ void hip_invalidate_kernel(volatile double* __restrict__ buf, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    // Touch every cache-line boundary to force a fresh load from MALL.
    // On AMD GFX94 a simpler way is: __builtin_amdgcn_buffer_gl2_inv()
    // but that requires raw AMDGCN intrinsics.  The volatile read achieves
    // the same effect for our purposes.
    if (idx < n) {
        (void)buf[idx];
    }
    __threadfence_system();  // re-acquire fence to drain all loads
}

// ------------------------------------------------------------------
// Fortran-callable entry points
// ------------------------------------------------------------------

// hip_write_with_fence(src, dst, n)
//   src  — pointer to first element of the flat send-chunk (Fortran array)
//   dst  — pointer to first element of the TARGET slot in peer's recv window
//          (caller must offset: dst => peer_win_base(my_dir_rank * chunk + 1))
//   n    — number of double elements to write
//
// Fortran interface:
//   subroutine hip_write_with_fence(src, dst, n) bind(C, name='hip_write_with_fence')
//     use iso_c_binding
//     real(c_double), intent(in)  :: src(*)
//     real(c_double), intent(out) :: dst(*)
//     integer(c_int), value       :: n
//   end subroutine
extern "C" void hip_write_with_fence(const double* src, double* dst, int n)
{
    if (n <= 0) return;
    const int threads = 256;
    const int blocks  = (n + threads - 1) / threads;
    hipLaunchKernelGGL(hip_write_fence_kernel,
                       dim3(blocks), dim3(threads), 0, 0,
                       src, dst, n);
    // Block until kernel + fence complete.  After this returns,
    // MPI_Barrier(MPI_COMM_WORLD) in Fortran is sufficient for ordering.
    hipError_t err = hipDeviceSynchronize();
    if (err != hipSuccess) {
        fprintf(stderr, "[hip_write_with_fence] hipDeviceSynchronize failed: %s\n",
                hipGetErrorString(err));
    }
}

// hip_invalidate_recv(buf, n)
//   Call on the receiving rank just before reading the shared recv window
//   if vmpi_hip_shmwrite.f90 reports mismatches despite correct writer flushing.
//
// Fortran interface:
//   subroutine hip_invalidate_recv(buf, n) bind(C, name='hip_invalidate_recv')
//     use iso_c_binding
//     real(c_double), intent(inout) :: buf(*)
//     integer(c_int), value         :: n
//   end subroutine
extern "C" void hip_invalidate_recv(double* buf, int n)
{
    if (n <= 0) return;
    const int threads = 256;
    const int blocks  = (n + threads - 1) / threads;
    hipLaunchKernelGGL(hip_invalidate_kernel,
                       dim3(blocks), dim3(threads), 0, 0,
                       buf, n);
    (void)hipDeviceSynchronize();
}

// ------------------------------------------------------------------
// hip_system_fence()
//   Standalone system-scope L2 write-back. Call AFTER a fused !$omp target
//   push (whose writes are dirty in L2 once the preceding hipDeviceSynchronize
//   has returned): __threadfence_system() emits buffer_wbl2 sc0:1 sc1:1, which
//   writes back the L2's dirty lines to MALL so all agents (other ranks/XCDs,
//   CPU) see the fresh data. Used by the production node-window I-transpose
//   (Fix 2 / TRP_I_SYSFENCE), where the strided push cannot be expressed as a
//   contiguous hip_write_with_fence.
//
// Fortran interface:
//   subroutine hip_system_fence() bind(C, name='hip_system_fence')
//   end subroutine
// ------------------------------------------------------------------
__global__ void hip_system_fence_kernel()
{
    __threadfence_system();
}

extern "C" void hip_system_fence()
{
    hipLaunchKernelGGL(hip_system_fence_kernel, dim3(1), dim3(1), 0, 0);
    hipError_t err = hipDeviceSynchronize();
    if (err != hipSuccess) {
        fprintf(stderr, "[hip_system_fence] hipDeviceSynchronize failed: %s\n",
                hipGetErrorString(err));
    }
}

// ------------------------------------------------------------------
// DIAGNOSTIC wrapper: writes a detailed trace to hip_trace_<rank>.log
// so we can see exactly where inside the HIP layer the wedge happens.
// Use this in place of hip_write_with_fence to debug.
// ------------------------------------------------------------------
extern "C" void hip_write_with_fence_diag(int rank,
                                           const double* src,
                                           double* dst,
                                           int n)
{
    char fname[64];
    snprintf(fname, sizeof(fname), "hip_trace_%04d.log", rank);
    FILE* f = fopen(fname, "a");
    if (!f) return;

    fprintf(f, "[ENTRY] n=%d src=%p dst=%p\n", n, (const void*)src, (void*)dst); fflush(f);

    if (n <= 0) { fclose(f); return; }

    int dev_count = -1;
    hipError_t e_count = hipGetDeviceCount(&dev_count);
    fprintf(f, "[DEV] hipGetDeviceCount=%d err=%d (%s)\n",
            dev_count, (int)e_count, hipGetErrorString(e_count)); fflush(f);

    int dev = -1;
    hipError_t e_getdev = hipGetDevice(&dev);
    fprintf(f, "[DEV] current device=%d err=%d (%s)\n",
            dev, (int)e_getdev, hipGetErrorString(e_getdev)); fflush(f);

    hipDeviceProp_t props;
    hipError_t e_prop = hipGetDeviceProperties(&props, dev);
    if (e_prop == hipSuccess) {
        fprintf(f, "[DEV] name=%s\n", props.name); fflush(f);
    }

    // Probe: can the GPU actually see this 'dst' pointer?
    hipPointerAttribute_t attr;
    hipError_t e_attr = hipPointerGetAttributes(&attr, (const void*)dst);
    fprintf(f, "[PTR] hipPointerGetAttributes err=%d (%s)\n",
            (int)e_attr, hipGetErrorString(e_attr)); fflush(f);

    const int threads = 256;
    const int blocks  = (n + threads - 1) / threads;
    fprintf(f, "[KERNEL] launching blocks=%d threads=%d\n", blocks, threads); fflush(f);

    hipLaunchKernelGGL(hip_write_fence_kernel,
                       dim3(blocks), dim3(threads), 0, 0,
                       src, dst, n);

    hipError_t launch_err = hipGetLastError();
    fprintf(f, "[KERNEL] post-launch hipGetLastError=%d (%s)\n",
            (int)launch_err, hipGetErrorString(launch_err)); fflush(f);

    fprintf(f, "[SYNC] calling hipDeviceSynchronize ...\n"); fflush(f);
    hipError_t sync_err = hipDeviceSynchronize();
    fprintf(f, "[SYNC] hipDeviceSynchronize returned, err=%d (%s)\n",
            (int)sync_err, hipGetErrorString(sync_err)); fflush(f);

    fprintf(f, "[EXIT] OK\n"); fflush(f);
    fclose(f);
}
