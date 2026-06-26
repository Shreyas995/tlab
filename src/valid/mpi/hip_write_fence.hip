// hip_write_fence.hip
//
// HIP C wrapper providing system-scope fenced writes + reader-side invalidation for the
// node-window / apudirect I-transpose coherence fixes (TRP_I_FUSEDFENCE / TRP_I_SYSFENCE).
// Callable from Fortran via iso_c_binding.
//
// Build on Hunter (Cray/ROCm):  hipcc -c hip_write_fence.hip -o hip_write_fence.o
//
// PERFORMANCE NOTE (fixed 2026-06-26)
// -----------------------------------
// The previous hip_write_fence_kernel ran __threadfence_system() in EVERY thread. On
// MI300A (GFX94) that compiles to  buffer_wbl2 sc0:1 sc1:1 + s_waitcnt  — a *full*
// system-scope L2 write-back — so a million-element copy issued a million full-L2
// write-backs (~71 ms/call, 42% of all GPU time in a pat_report). A single
// __threadfence_system() already writes back the WHOLE L2, so it only needs to run ONCE,
// after the copy completes. hip_write_with_fence now does: copy kernel (no fence) -> a
// 1-thread fence kernel on the SAME in-order stream (runs after the copy, flushes the
// copied lines) -> one hipDeviceSynchronize. This is the cheap SYSFENCE mechanism.

#include <hip/hip_runtime.h>
#include <stdio.h>

// ------------------------------------------------------------------
// Kernels
// ------------------------------------------------------------------

// Copy src[0:n) -> dst[0:n).  NO per-thread fence (see performance note).
__global__ void hip_copy_kernel(const double* __restrict__ src,
                                 double* __restrict__ dst,
                                 int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dst[idx] = src[idx];
    }
}

// Single system-scope L2 write-back: __threadfence_system() emits buffer_wbl2 sc0:1 sc1:1,
// writing the WHOLE L2's dirty lines back to MALL so all agents (other ranks/XCDs, CPU) see
// the fresh data. Launched with <<<1,1>>>; on the default (in-order) stream it runs AFTER
// the preceding copy kernel completes, so it covers all of that copy's writes.
__global__ void hip_system_fence_kernel()
{
    __threadfence_system();
}

// Reader-side invalidation: touch every element (volatile) so the subsequent GPU unpack
// reloads fresh data instead of a stale L2 line. Fence ONCE (one thread), not per-thread.
__global__ void hip_invalidate_kernel(volatile double* __restrict__ buf, int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        (void)buf[idx];
    }
    if (idx == 0) {
        __threadfence_system();
    }
}

// LEGACY per-thread fenced write — kept ONLY for the diagnostic wrapper below. Do NOT use
// on a hot path: it fences in every thread (catastrophically slow). See performance note.
__global__ void hip_write_fence_kernel(const double* __restrict__ src,
                                        double* __restrict__ dst,
                                        int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dst[idx] = src[idx];
    }
    __threadfence_system();
}

// ------------------------------------------------------------------
// Fortran-callable entry points
// ------------------------------------------------------------------

// hip_write_with_fence(src, dst, n) — copy src -> dst then make it system-visible.
//   subroutine hip_write_with_fence(src, dst, n) bind(C, name='hip_write_with_fence')
extern "C" void hip_write_with_fence(const double* src, double* dst, int n)
{
    if (n <= 0) return;
    const int threads = 256;
    const int blocks  = (n + threads - 1) / threads;
    // 1. copy (no per-thread fence)
    hipLaunchKernelGGL(hip_copy_kernel, dim3(blocks), dim3(threads), 0, 0, src, dst, n);
    // 2. ONE system-scope L2 write-back, same default stream -> runs after the copy,
    //    flushes the copied lines to MALL. O(1) instead of O(n) fences.
    hipLaunchKernelGGL(hip_system_fence_kernel, dim3(1), dim3(1), 0, 0);
    // 3. block until the copy + fence complete (callers — node-window push, all-MPI
    //    fallback flush — rely on this being synchronous before MPI_Win_fence / ISEND).
    hipError_t err = hipDeviceSynchronize();
    if (err != hipSuccess) {
        fprintf(stderr, "[hip_write_with_fence] hipDeviceSynchronize failed: %s\n",
                hipGetErrorString(err));
    }
}

// hip_invalidate_recv(buf, n) — reader acquire before a GPU unpack of a shared window.
//   subroutine hip_invalidate_recv(buf, n) bind(C, name='hip_invalidate_recv')
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

// hip_system_fence() — standalone single system-scope L2 write-back (1 thread).
//   subroutine hip_system_fence() bind(C, name='hip_system_fence')
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
// DIAGNOSTIC wrapper (debug only): per-rank trace of the HIP layer. Uses the legacy
// per-thread fenced kernel; not for production use.
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

    int dev = -1;
    hipError_t e_getdev = hipGetDevice(&dev);
    fprintf(f, "[DEV] current device=%d err=%d (%s)\n",
            dev, (int)e_getdev, hipGetErrorString(e_getdev)); fflush(f);

    const int threads = 256;
    const int blocks  = (n + threads - 1) / threads;
    fprintf(f, "[KERNEL] launching blocks=%d threads=%d\n", blocks, threads); fflush(f);

    hipLaunchKernelGGL(hip_write_fence_kernel,
                       dim3(blocks), dim3(threads), 0, 0,
                       src, dst, n);

    hipError_t launch_err = hipGetLastError();
    fprintf(f, "[KERNEL] post-launch hipGetLastError=%d (%s)\n",
            (int)launch_err, hipGetErrorString(launch_err)); fflush(f);

    hipError_t sync_err = hipDeviceSynchronize();
    fprintf(f, "[SYNC] hipDeviceSynchronize returned, err=%d (%s)\n",
            (int)sync_err, hipGetErrorString(sync_err)); fflush(f);

    fprintf(f, "[EXIT] OK\n"); fflush(f);
    fclose(f);
}
