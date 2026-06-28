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
#include <stdint.h>

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

// FUSED multi-segment strided SCATTER + PER-WORKGROUP system fence in ONE kernel. This collapses
// the per-peer hip_write_with_fence loop (and its strided staging gather) used by the apudirect /
// node-window REAL transpose pushes: instead of N contiguous copies, each its own kernel launch
// AND its own hipDeviceSynchronize, a single launch scatters the (possibly strided) source into all
// N peer window-slots and fences per workgroup.
//
//   element (m,i,j), m in [0,n_peer), i in [0,n_rows), j in [0,n_cols):
//     dst[m*dst_stride + dst_base + i*n_cols + j] = src[m*src_seg_stride + i*src_row_stride + j]
//
// Each peer's destination slot is contiguous (i*n_cols+j); the SOURCE is strided when
// src_row_stride != n_cols (the gather-on-the-fly forward/backward transposes) or contiguous when
// src_seg_stride = n_rows*n_cols, src_row_stride = 1, n_cols = 1 (the flat per-peer transposes).
// Coherency is BYTE-IDENTICAL to hip_copy_fence_kernel: an unconditional __syncthreads() then a
// thread-0 __threadfence_system() (whole-L2 write-back to MALL) per workgroup. __threadfence_system
// is NOT address-scoped, so a scattered set of stores is committed exactly as a contiguous one is —
// the per-peer loop / staging buffer were never the source of coherency, the per-workgroup whole-L2
// write-back is. dst offsets are int64 (peer windows are many GB apart on the big grids).
__global__ void hip_pushseg_fence_kernel(const double* __restrict__ src,
                                         double* __restrict__ dst,
                                         int n_peer, int n_rows, int n_cols,
                                         int src_seg_stride, int src_row_stride,
                                         int64_t dst_stride, int64_t dst_base)
{
    const int64_t seg_len = (int64_t)n_rows * n_cols;
    const int64_t total   = (int64_t)n_peer * seg_len;
    const int64_t tid     = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < total) {
        const int64_t m = tid / seg_len;
        const int64_t r = tid - m * seg_len;                 // contiguous position within the peer slot
        const int64_t i = r / n_cols;
        const int64_t j = r - i * n_cols;
        const int64_t s = m * src_seg_stride + i * src_row_stride + j;
        const int64_t d = m * dst_stride + dst_base + i * (int64_t)n_cols + j;
        dst[d] = src[s];
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

// hip_pushseg_fence(...) — fused multi-segment scatter + per-workgroup system fence (ONE launch,
// ONE sync). Drop-in replacement for the per-peer hip_write_with_fence loop on the real transposes.
//   subroutine hip_pushseg_fence(src, dst, n_peer, n_rows, n_cols, src_seg_stride, src_row_stride,
//                                dst_stride, dst_base) bind(C, name='hip_pushseg_fence')
extern "C" void hip_pushseg_fence(const double* src, double* dst,
                                  int n_peer, int n_rows, int n_cols,
                                  int src_seg_stride, int src_row_stride,
                                  int64_t dst_stride, int64_t dst_base)
{
    const int64_t total = (int64_t)n_peer * n_rows * n_cols;
    if (total <= 0) return;
    const int threads = 256;
    const int64_t blocks = (total + threads - 1) / threads;
    hipLaunchKernelGGL(hip_pushseg_fence_kernel, dim3((unsigned int)blocks), dim3(threads), 0, 0,
                       src, dst, n_peer, n_rows, n_cols,
                       src_seg_stride, src_row_stride, dst_stride, dst_base);
    // block until the scatter+fence kernel completes (caller relies on this being synchronous before
    // MPI_Win_fence closes the RMA epoch — the V1 iter-53 lesson: the fence does not wait on GPU work).
    hipError_t err = hipDeviceSynchronize();
    if (err != hipSuccess) {
        fprintf(stderr, "[hip_pushseg_fence] hipDeviceSynchronize failed: %s\n",
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
