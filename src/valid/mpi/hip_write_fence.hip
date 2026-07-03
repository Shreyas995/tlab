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

// FUSED multi-segment scatter + per-workgroup system fence (ONE launch, ONE sync). Drop-in
// replacement for the per-peer hip_write_with_fence LOOP on the real transposes: it scatters the
// STRIDED source directly into every peer's window slot, so it also subsumes the staging gather
// kernel. For peer m, element (i,j): dst[m*dst_stride + dst_base + i*n_cols + j] = src[m*src_seg_stride
// + i*src_row_stride + j]. Coherency is BYTE-IDENTICAL to hip_copy_fence_kernel: one __syncthreads()
// then a thread-0 __threadfence_system() (whole-L2 write-back to MALL) per workgroup — the fence is
// NOT address-scoped, so a scattered set of stores commits exactly as a contiguous one does (the
// per-peer loop / staging buffer were never the source of coherency, the per-workgroup write-back is).
// dst offsets are int64 (peer windows are many GB apart on the big grids).
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
// ONE sync). Drop-in replacement for the per-peer hip_write_with_fence LOOP on the real apudirect
// transposes (it reads the strided source directly, so it also replaces the staging gather kernel).
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

// F1 (TRP_I_MEMCPY_ASYNC): pipelined variant of hip_memcpy_push. The per-peer node-window push loop issues N
// non-blocking hipMemcpyAsync onto ONE persistent stream, then drains them with a SINGLE hip_memcpy_push_sync()
// after the loop (before MPI_Win_fence). Identical runtime DMA / MALL commit as the blocking hip_memcpy_push --
// async vs blocking changes only when the HOST is notified, not the memory scope -- so the cross-XCD window data
// is committed on stream-sync exactly as the blocking copy commits on return. Benefit: the host queues all N
// copies instead of round-tripping per copy, freeing the CPU during the DMAs (and letting the K paths post the
// inter-node MPI_ISEND while the intra copies run). One persistent stream per rank (OMP_NUM_THREADS=1 -> no
// concurrency), lazily created, never destroyed (1 stream/run = negligible; matches the no-teardown style).
static hipStream_t g_push_stream = nullptr;

// hip_memcpy_push_async(dst, src, n) -- queue one coherent copy on g_push_stream WITHOUT blocking the host.
// The caller MUST call hip_memcpy_push_sync() once after the per-peer loop to complete the batch.
extern "C" void hip_memcpy_push_async(double* dst, const double* src, int n)
{
    if (n <= 0) return;
    if (g_push_stream == nullptr) {
        hipError_t cerr = hipStreamCreate(&g_push_stream);
        if (cerr != hipSuccess) {
            fprintf(stderr, "[hip_memcpy_push_async] hipStreamCreate failed: %s\n",
                    hipGetErrorString(cerr));
            return;
        }
    }
    hipError_t err = hipMemcpyAsync(dst, src, (size_t)n * sizeof(double),
                                    hipMemcpyDefault, g_push_stream);
    if (err != hipSuccess) {
        fprintf(stderr, "[hip_memcpy_push_async] hipMemcpyAsync failed: %s\n",
                hipGetErrorString(err));
    }
}

// hip_memcpy_push_sync() -- block until ALL queued async pushes on g_push_stream have completed and committed
// to MALL. Pairs with hip_memcpy_push_async; safe no-op if the stream was never created.
extern "C" void hip_memcpy_push_sync()
{
    if (g_push_stream == nullptr) return;
    hipError_t err = hipStreamSynchronize(g_push_stream);
    if (err != hipSuccess) {
        fprintf(stderr, "[hip_memcpy_push_sync] hipStreamSynchronize failed: %s\n",
                hipGetErrorString(err));
    }
}

// hip_memcpy2d_push(dst, dpitch, src, spitch, width, height) — apudirect optimization (2026-07-01).
// ONE blocking, system-coherent strided DMA (hipMemcpy2D, hipMemcpyDefault) that pushes ALL peers at once,
// replacing (a) the per-peer hip_memcpy_push LOOP (N blocking copies -> 1) and (b) the per-workgroup
// __threadfence_system push kernels (hip_pushseg_fence / hip_write_with_fence). The N peer recv-slots are a
// clean 2D block: consecutive peers are dpitch doubles apart in dst, each row is 'width' contiguous doubles,
// source rows are spitch doubles apart. hipMemcpy2D is a runtime-managed transfer: SYNCHRONOUS on the host and
// SYSTEM-COHERENT on completion (destination committed to MALL before return) — identical writer-release
// semantics to the 1D hip_memcpy_push, so NO manual writer-side fence is required. All pitches/width/height
// arrive in DOUBLES and are converted to bytes here; they are int64 because peer window segments are many GB
// apart on big grids.
//   subroutine hip_memcpy2d_push(dst, dpitch, src, spitch, width, height) bind(C, name='hip_memcpy2d_push')
// dst/src arrive as whole-array base pointers with 0-based element offsets (the hip_pushseg_fence idiom).
extern "C" void hip_memcpy2d_push(double* dst, int64_t dst_off, int64_t dpitch,
                                  const double* src, int64_t src_off, int64_t spitch,
                                  int64_t width, int64_t height)
{
    if (width <= 0 || height <= 0) return;
    hipError_t err = hipMemcpy2D(dst + dst_off, (size_t)dpitch * sizeof(double),
                                 src + src_off, (size_t)spitch * sizeof(double),
                                 (size_t)width * sizeof(double),
                                 (size_t)height, hipMemcpyDefault);
    if (err != hipSuccess) {
        fprintf(stderr, "[hip_memcpy2d_push] hipMemcpy2D failed: %s\n",
                hipGetErrorString(err));
    }
}

// T30 (apudirect forward pilot, 2026-07-03): gather + 2D peer push on ONE persistent stream.
// The forward real transpose gathers a STRIDED source into a contiguous per-peer staging buffer, then pushes
// all peers with a strided 2D DMA. The default path runs the gather as an OMP target region (OMP offload
// stream) and the push as hipMemcpy2D (default HIP stream), so a full-device hipDeviceSynchronize is needed
// between them to order that cross-stream dependency. Here BOTH the gather kernel and the 2D copy run on the
// SAME persistent stream (g_push_stream, shared with the F1 async path): HIP stream ordering makes the gather's
// writes visible to the copy with NO explicit sync, and a single hipStreamSynchronize at the end commits the
// coherent async DMA to MALL before the close MPI_Win_fence (identical writer-release to the blocking copy --
// async vs blocking changes only host-return timing, not memory scope; see the F1 note above).
// gather index (0-based, mirrors the Fortran OMP loop exactly):
//   stage[m*mas + i*nlines_p + j] = src[m*nlines_p + i*npage + j],  mas = nmax_p*nlines_p,
//   over m in [0,npro), i in [0,nmax_p), j in [0,nlines_p).
__global__ void hip_gather3d_kernel(double* __restrict__ stage, const double* __restrict__ src,
                                    int npro, int nmax_p, int nlines_p, int npage)
{
    const long mas = (long)nmax_p * (long)nlines_p;
    const long total = (long)npro * mas;
    for (long idx = (long)blockIdx.x * blockDim.x + threadIdx.x; idx < total;
         idx += (long)gridDim.x * blockDim.x) {
        long m = idx / mas;
        long rem = idx - m * mas;
        long i = rem / (long)nlines_p;
        long j = rem - i * (long)nlines_p;
        stage[m * mas + i * (long)nlines_p + j] = src[m * (long)nlines_p + i * (long)npage + j];
    }
}

//   subroutine hip_gather_memcpy2d_push(dst, dst_off, dpitch, src, stage, npro, nmax_p, nlines_p, npage, width, height)
// dst = peer window base (0-based dst_off applied C-side); dpitch = per-peer stride (doubles); src = strided
// source; stage = contiguous staging (= the recv slot pattern); width = mas (doubles/peer), height = npro.
extern "C" void hip_gather_memcpy2d_push(double* dst, int64_t dst_off, int64_t dpitch,
                                         const double* src, double* stage,
                                         int npro, int nmax_p, int nlines_p, int npage,
                                         int64_t width, int64_t height)
{
    if (npro <= 0 || nmax_p <= 0 || nlines_p <= 0 || width <= 0 || height <= 0) return;
    if (g_push_stream == nullptr) {
        hipError_t cerr = hipStreamCreate(&g_push_stream);
        if (cerr != hipSuccess) {
            fprintf(stderr, "[hip_gather_memcpy2d_push] hipStreamCreate failed: %s\n",
                    hipGetErrorString(cerr));
            return;
        }
    }
    long total = (long)npro * (long)nmax_p * (long)nlines_p;
    int threads = 256;
    int blocks = (int)((total + threads - 1) / threads);
    if (blocks < 1) blocks = 1;
    hipLaunchKernelGGL(hip_gather3d_kernel, dim3(blocks), dim3(threads), 0, g_push_stream,
                       stage, src, npro, nmax_p, nlines_p, npage);
    hipError_t err = hipMemcpy2DAsync(dst + dst_off, (size_t)dpitch * sizeof(double),
                                      stage, (size_t)width * sizeof(double),
                                      (size_t)width * sizeof(double),
                                      (size_t)height, hipMemcpyDefault, g_push_stream);
    if (err != hipSuccess) {
        fprintf(stderr, "[hip_gather_memcpy2d_push] hipMemcpy2DAsync failed: %s\n",
                hipGetErrorString(err));
    }
    err = hipStreamSynchronize(g_push_stream);   // drain gather+DMA; commits the coherent copy to MALL
    if (err != hipSuccess) {
        fprintf(stderr, "[hip_gather_memcpy2d_push] hipStreamSynchronize failed: %s\n",
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
