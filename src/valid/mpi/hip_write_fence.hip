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
// THREE variants of the fence, only the last is both correct AND fast:
//   (1) per-THREAD fence (original hip_write_fence_kernel): __threadfence_system() in every
//       thread -> a million-element copy issued a million full-L2 write-backs (buffer_wbl2
//       sc0:1 sc1:1), ~71 ms/call = 42% of GPU time. Correct but catastrophically slow.
//   (2) single <<<1,1>>> fence in a SEPARATE kernel after the copy (the 2026-06-26 perf
//       "fix"): O(1) and fast, but this IS the SYSFENCE mechanism — a separate single-thread
//       fence that only CUT (~20x) and never eliminated the cross-rank push race (the writer-
//       side window corruption returned: Re125000 fabricdirect crashed again at ~31-53 iters).
//   (3) FUSED copy + PER-WORKGROUP fence in ONE kernel (hip_copy_fence_kernel, current):
//       each workgroup copies then ONE thread fences after __syncthreads. ~n/256 fences (fast),
//       fused with the write (no inter-kernel L1->L2 gap), and spread across CUs (covers every
//       L2 slice). This is the correct-and-cheap middle ground. hip_write_with_fence uses it.

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

// FUSED copy + PER-WORKGROUP system-scope fence (the correct+fast middle ground).
// Each workgroup copies its chunk, __syncthreads() so the whole workgroup's stores are
// issued, then ONE thread (threadIdx.x==0) does __threadfence_system() (buffer_wbl2 sc0:1
// sc1:1 = whole-L2 write-back to MALL). So every store is followed, within its own
// workgroup, by a system fence -> all copied lines are made cross-rank/CPU visible BEFORE
// the host's hipDeviceSynchronize returns. Cost = (n/256) fences, not n (per-thread, the
// 71 ms/call catastrophe) and not 1 (a single <<<1,1>>> fence in a SEPARATE kernel = the
// SYSFENCE mechanism that only CUT, never eliminated, the cross-rank push race). Being in
// the SAME kernel as the copy removes the inter-kernel L1->L2 timing gap, and running on
// many CUs covers every L2 slice that took a write.
__global__ void hip_copy_fence_kernel(const double* __restrict__ src,
                                      double* __restrict__ dst,
                                      int n)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        dst[idx] = src[idx];
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        __threadfence_system();
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
    // FUSED copy + per-workgroup system fence in ONE kernel: correct (every store is
    // system-fenced within its workgroup, no inter-kernel L1->L2 gap, all L2 slices covered)
    // AND fast (~n/256 fences, not n and not 1). Replaces the copy + separate <<<1,1>>> fence
    // (= SYSFENCE) that only cut, never eliminated, the cross-rank push race.
    hipLaunchKernelGGL(hip_copy_fence_kernel, dim3(blocks), dim3(threads), 0, 0, src, dst, n);
    // block until the copy+fence kernel completes (callers — node-window/apudirect push and
    // the all-MPI fallback flush — rely on this being synchronous before MPI_Win_fence/ISEND).
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

// hip_memcpy_push(dst, src, n) — V2 (TRP_I_MEMCPY): blocking, system-coherent push of n doubles
// src -> dst via hipMemcpy(...hipMemcpyDefault). A hipMemcpy is a runtime-managed transfer that is
// SYNCHRONOUS on the host and SYSTEM-COHERENT on completion: the destination lines are committed to
// MALL before the call returns, so we need NO in-kernel __threadfence_system on the writer side (this
// is the alternative to the per-workgroup fused-fence kernel hip_write_with_fence, for testing whether
// a runtime-managed coherent copy fully closes the cross-XCD node-window push race that the GPU fence
// only cut, never eliminated). The READER still acquires via hip_invalidate_recv before its GPU unpack
// because another XCD's L2 may still hold a stale line for the destination window address (hipMemcpy
// guarantees the writer's release, not the remote reader's acquire). hipMemcpyDefault dispatches on the
// pointer attributes (both ends are MI300A unified/coarse-grained shared-window memory). Because the
// copy is synchronous, no extra hipDeviceSynchronize is required after it.
//   subroutine hip_memcpy_push(dst, src, n) bind(C, name='hip_memcpy_push')
extern "C" void hip_memcpy_push(double* dst, const double* src, int n)
{
    if (n <= 0) return;
    hipError_t err = hipMemcpy(dst, src, (size_t)n * sizeof(double), hipMemcpyDefault);
    if (err != hipSuccess) {
        fprintf(stderr, "[hip_memcpy_push] hipMemcpy failed: %s\n",
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
