# SIMD Kernels in OpenBLAS

This document explains where SIMD (Single Instruction Multiple Data) kernels for algorithms are written in the OpenBLAS codebase.

## Quick Answer

**SIMD kernels are located in the `kernel/` directory**, organized by CPU architecture. Each architecture subdirectory contains optimized assembly and C code that uses SIMD instructions specific to that platform.

## Directory Structure

All SIMD-optimized kernels live under:

```
kernel/
├── simd/              # Universal intrinsics headers for cross-architecture SIMD
│   ├── intrin.h       # Main intrinsics header
│   ├── intrin_sse.h   # x86 SSE intrinsics
│   ├── intrin_avx.h   # x86 AVX intrinsics
│   ├── intrin_avx512.h # x86 AVX-512 intrinsics
│   └── intrin_neon.h  # ARM NEON intrinsics
│
├── x86/               # 32-bit x86 (SSE, SSE2, SSE3, SSSE3, AVX)
├── x86_64/            # 64-bit x86 (SSE → AVX-512, AMX)
├── arm/               # 32-bit ARM (NEON, VFPv2, VFPv3)
├── arm64/             # 64-bit ARM (NEON/ASIMD, SVE, SME)
├── power/             # IBM POWER (VSX, VSX2, VSX3, ISA 3.0)
├── riscv64/           # RISC-V 64-bit (RVV - RISC-V Vector Extension)
├── loongarch64/       # LoongArch (LSX, LASX)
├── mips/              # MIPS 32-bit (MSA)
├── mips64/            # MIPS 64-bit (MSA)
├── zarch/             # IBM z/Architecture (Vector Facility)
├── sparc/             # SPARC (VIS)
├── ia64/              # Intel Itanium (legacy)
├── alpha/             # DEC Alpha (legacy, from GotoBLAS)
├── e2k/               # Elbrus
├── csky/              # T-Head C-SKY
└── generic/           # Plain C fallback kernels (no SIMD)
```

## SIMD Technologies by Architecture

| Architecture | SIMD Extensions | Example Files |
|--------------|----------------|---------------|
| **x86-64** | SSE, SSE2, SSE3, SSSE3, SSE4.1, AVX, AVX2, FMA3, AVX-512 | `dgemm_kernel_4x8_haswell.S` (AVX2)<br>`sgemm_kernel_16x4_skylakex.S` (AVX-512) |
| **ARM64** | NEON (Advanced SIMD), SVE, SME | `sgemm_kernel_sve_v2x8.S` (SVE)<br>`dgemm_kernel_8x4_thunderx2t99.S` (NEON) |
| **ARM** | NEON, VFPv3, VFPv2 | `sgemm_kernel_4x4_vfpv3.S` |
| **POWER** | VSX, VSX2, VSX3 | `dgemm_kernel_power10.S` |
| **RISC-V** | RVV (RISC-V Vector 0.7.1, 1.0) | `sgemm_kernel_rvv_v1.c` |
| **LoongArch64** | LSX, LASX | Files with `_lsx` or `_lasx` suffixes |
| **MIPS** | MSA (MIPS SIMD Architecture) | Files with `_msa` suffix |

## File Naming Conventions

SIMD kernel files follow these patterns:

### Assembly Files (`.S`)
```
{precision}{function}_kernel_{params}_{simd_variant}.S
```

**Examples:**
- `dgemm_kernel_4x8_haswell.S` - Double-precision GEMM for Haswell (AVX2)
- `sgemm_kernel_16x4_skylakex.S` - Single-precision GEMM for Skylake-X (AVX-512)
- `sgemm_kernel_sve_v2x8.S` - Single-precision GEMM using ARM SVE
- `daxpy_sse.S` - Double-precision AXPY using SSE

### C Files with Intrinsics (`.c`)
```
{precision}{function}_{simd_variant}.c
```

**Examples:**
- `dgemm_kernel_power10.c` - Uses POWER10 VSX intrinsics
- `sgemm_kernel_rvv_v1.c` - Uses RISC-V Vector intrinsics
- `axpy_sve.c` - ARM SVE intrinsics

### Precision Prefixes
- `s` = Single precision (float)
- `d` = Double precision (double)
- `c` = Complex single precision
- `z` = Complex double precision
- `sb` = BFloat16

## How Kernels Are Selected

Each CPU variant has a **`KERNEL.<CPU>`** configuration file that maps function roles to specific implementations.

### Example: `kernel/x86_64/KERNEL.HASWELL`
```makefile
DGEMMKERNEL    = dgemm_kernel_4x8_haswell.S
SGEMMKERNEL    = sgemm_kernel_8x4_haswell.S
DTRMMKERNEL    = dtrmm_kernel_4x8_haswell.c
DAXPYKERNEL    = daxpy_sse.S
DDOTKERNEL     = ddot_sse2.S
```

This file tells the build system:
- For double-precision GEMM on Haswell, use `dgemm_kernel_4x8_haswell.S` (AVX2 assembly)
- For single-precision GEMM on Haswell, use `sgemm_kernel_8x4_haswell.S`
- For AXPY, use the SSE version
- And so on...

### Finding Active Kernels

To find which kernel is used for a specific CPU:

1. Look at `kernel/{arch}/KERNEL.{CPU}` 
2. Find the variable for your operation (e.g., `DGEMMKERNEL`)
3. The value points to the source file

**Example workflow:**
```bash
# For Haswell dgemm kernel:
$ cat kernel/x86_64/KERNEL.HASWELL | grep DGEMMKERNEL
DGEMMKERNEL    =  dgemm_kernel_4x8_haswell.S

# The kernel source is at:
$ ls kernel/x86_64/dgemm_kernel_4x8_haswell.S
```

## Universal Intrinsics (`kernel/simd/`)

The `kernel/simd/` directory provides architecture-agnostic SIMD intrinsics wrappers:

```
kernel/simd/
├── intrin.h          # Master header, includes architecture-specific headers
├── intrin_sse.h      # x86 SSE wrapper
├── intrin_avx.h      # x86 AVX wrapper
├── intrin_avx512.h   # x86 AVX-512 wrapper
└── intrin_neon.h     # ARM NEON wrapper
```

These headers provide a **unified interface** so that some kernels can be written once and compiled for multiple architectures. They define common operations like:
- `v_load` / `v_store` - Load/store vectors
- `v_add` / `v_mul` - Arithmetic operations
- `v_fma` - Fused multiply-add

**Example usage:**
```c
#include "kernel/simd/intrin.h"

// This code can compile for SSE, AVX, NEON, etc.
void example_kernel(float *a, float *b, float *c, int n) {
    for (int i = 0; i < n; i += V_WIDTH) {
        v_type va = v_load(&a[i]);
        v_type vb = v_load(&b[i]);
        v_type vc = v_mul(va, vb);
        v_store(&c[i], vc);
    }
}
```

## Dynamic Architecture Selection

When built with `DYNAMIC_ARCH=1`, OpenBLAS includes **multiple CPU-specific kernel sets** in a single library and selects the appropriate one at runtime based on CPUID detection:

```
Library startup
    └─ driver/others/init.c
           └─ driver/others/dynamic.c (or dynamic_arm64.c, etc.)
                  └─ Reads CPUID → selects gotoblas_<CPU> function table
                         └─ All kernels use this dispatch table
```

**Example:** On an x86-64 system with AVX2, the library will:
1. Detect "Haswell" family CPU via CPUID
2. Select the `gotoblas_HASWELL` function table
3. Use Haswell-optimized kernels (AVX2) for all operations

## Writing New SIMD Kernels

To add a new SIMD kernel:

1. **Create the kernel file** in `kernel/{arch}/`
   - Use assembly (`.S`) for maximum performance
   - Or use intrinsics in C (`.c`) for portability
   
2. **Add to `KERNEL.<CPU>` file**
   ```makefile
   DGEMMKERNEL = my_new_dgemm_kernel.S
   ```

3. **Optional: Use universal intrinsics** from `kernel/simd/intrin.h` for cross-platform code

4. **Tune blocking parameters** in `param.h`:
   - `GEMM_P`, `GEMM_Q`, `GEMM_R` - Cache blocking sizes
   - `GEMM_UNROLL_M`, `GEMM_UNROLL_N` - Register blocking

5. **Test thoroughly** with `make test` and `make lapack-test`

## Common BLAS Operations and Their Kernels

| Operation | Description | Typical Kernel Variables |
|-----------|-------------|-------------------------|
| **GEMM** | Matrix-matrix multiply | `SGEMMKERNEL`, `DGEMMKERNEL`, `CGEMMKERNEL`, `ZGEMMKERNEL` |
| **GEMV** | Matrix-vector multiply | `SGEMVNKERNEL`, `DGEMVTKERNEL` |
| **AXPY** | y = αx + y | `SAXPYKERNEL`, `DAXPYKERNEL` |
| **DOT** | Dot product | `SDOTKERNEL`, `DDOTKERNEL` |
| **TRSM** | Triangular solve | `STRSMKERNEL_LN`, `DTRSMKERNEL_RU` |
| **COPY** | Copy vectors | `SCOPYKERNEL`, `DCOPYKERNEL` |
| **IAMAX** | Index of max element | `ISAMAXKERNEL`, `IDAMAXKERNEL` |

Each has precision variants (S/D/C/Z) and sometimes transpose/operation variants.

## References

- **Architecture overview:** See `docs/architecture.md` for the full three-layer architecture
- **Developer guide:** See `docs/developers.md` for kernel development details
- **Goto paper:** [Anatomy of High-Performance Matrix Multiplication](http://www.cs.utexas.edu/~flame/web/FLAMEPublications.html) - describes the algorithm that GEMM kernels implement

## Summary

!!! info "Key Takeaways"
    - **Location:** SIMD kernels are in `kernel/{architecture}/`
    - **Languages:** Assembly (`.S`) or C with intrinsics (`.c`)
    - **Selection:** CPU-specific `KERNEL.<CPU>` files map operations to implementations
    - **Cross-platform:** `kernel/simd/intrin*.h` provides universal intrinsics
    - **Examples:** `dgemm_kernel_4x8_haswell.S`, `sgemm_kernel_sve_v2x8.S`, `axpy_sse.S`
