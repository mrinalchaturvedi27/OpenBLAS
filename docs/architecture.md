# OpenBLAS Architectural Analysis

This document provides a comprehensive architectural overview of OpenBLAS for developers
who are new to the codebase and want to understand how it is organized, how the components
relate to one another, and how a typical computation flows from API call to hardware.

## High-Level Overview

OpenBLAS is a fork of GotoBLAS2 1.13 and implements:

- **BLAS** (Basic Linear Algebra Subprograms) — Levels 1, 2, and 3
- **CBLAS** — The standardized C interface to BLAS
- **LAPACK** — Higher-level linear algebra routines
- **LAPACKE** — The C interface to LAPACK

The library is designed around a three-layer architecture:

```
┌─────────────────────────────────────────────────────────────┐
│  User code                                                  │
│   (Fortran BLAS API / CBLAS C API / LAPACK / LAPACKE)       │
└──────────────────────┬──────────────────────────────────────┘
                       │  public API
┌──────────────────────▼──────────────────────────────────────┐
│  Interface Layer  (interface/)                              │
│   Argument validation, type dispatch, SMP threshold logic   │
└──────────────────────┬──────────────────────────────────────┘
                       │
┌──────────────────────▼──────────────────────────────────────┐
│  Driver Layer  (driver/level2, driver/level3, driver/others)│
│   Algorithmic skeleton (Goto algorithm), threading,         │
│   memory management, LAPACK routines                        │
└──────────────────────┬──────────────────────────────────────┘
                       │
┌──────────────────────▼──────────────────────────────────────┐
│  Kernel Layer  (kernel/<arch>/)                             │
│   Architecture-optimized inner loops (ASM, intrinsics, C)  │
└─────────────────────────────────────────────────────────────┘
```

## Directory Structure

```
OpenBLAS/
├── interface/          BLAS & CBLAS public entry points
│   ├── lapack/         LAPACK interface wrappers
│   └── netlib/         Netlib-compatible wrappers
│
├── driver/             Portable algorithmic logic
│   ├── level2/         Level-2 BLAS drivers (GEMV, GER, SYMV, …)
│   ├── level3/         Level-3 BLAS drivers (GEMM, TRSM, SYRK, …)
│   ├── mapper/         (reserved/experimental)
│   └── others/         Memory, threading, dynamic dispatch, utilities
│
├── kernel/             CPU-specific inner loops
│   ├── x86_64/         x86-64 (SSE2 → AVX-512 / AMX)
│   ├── arm64/          AArch64 / Neon / SVE
│   ├── arm/            ARMv5–ARMv7
│   ├── power/          IBM POWER (POWER6–POWER10)
│   ├── riscv64/        RISC-V 64-bit
│   ├── zarch/          IBM z/Architecture
│   ├── mips/ mips64/   MIPS 32 & 64
│   ├── ia64/           Intel Itanium (legacy)
│   ├── alpha/          DEC Alpha (legacy, from GotoBLAS)
│   ├── sparc/          SPARC
│   ├── loongarch64/    Loongson LoongArch
│   ├── csky/           T-Head C-SKY
│   ├── e2k/            Elbrus
│   ├── simd/           Architecture-agnostic Universal Intrinsics helpers
│   └── generic/        Plain C fallback kernels (used by all architectures)
│
├── lapack/             Optimized replacements for selected LAPACK routines
│   ├── getrf/ getf2/   LU factorization
│   ├── potrf/ potf2/   Cholesky factorization
│   ├── trtri/ trti2/   Triangular inverse
│   ├── trtrs/          Triangular solve
│   ├── laswp/          Row permutation
│   ├── lauu2/ lauum/   Upper triangular product
│   └── getrs/          LU back-substitution
│
├── lapack-netlib/      Full reference LAPACK (from Netlib), used when no
│                       optimized version exists
├── relapack/           Recursive LAPACK (optional, BUILD_RELAPACK=1)
├── reference/          Fortran BLAS reference implementation (sanity-check only)
│
├── exports/            Scripts to generate the library's public symbol table
├── cmake/              CMake build system modules
├── benchmark/          Simple standalone benchmarks
├── test/               BLAS test programs
├── ctest/              CBLAS test programs
├── utest/              OpenBLAS-specific regression tests
├── cpp_thread_test/    C++ threading integration tests
│
├── common.h            Master include: pulls in all common_*.h headers
├── common_*.h          Per-topic headers (param, thread, level1–3, arch, …)
├── cblas.h             Public CBLAS header
├── Makefile            Top-level build entry point
├── Makefile.system     Architecture detection, variable defaults
├── Makefile.rule       User-configurable build options
├── Makefile.prebuild   Pre-build step (runs getarch, c_check, f_check)
├── Makefile.$ARCH      Architecture-specific compiler flags & cache sizes
├── getarch.c           CPU identification utility (compiled for host)
├── getarch_2nd.c       Second-pass arch info (compiled for target)
├── cpuid_*.c           Per-architecture CPUID helpers
└── param.h             GEMM tuning parameters (block sizes, unroll factors)
```

## Layer Details

### Interface Layer (`interface/`)

Each BLAS function has one source file in `interface/`.  The same file is compiled
multiple times with different preprocessor macros (`SINGLE`, `DOUBLE`, `COMPLEX`,
`COMPLEX16`, `BFLOAT16`) to produce the `s`, `d`, `c`, `z`, `sb` prefix variants.

Responsibilities:

- Parse and validate arguments (calls `xerbla` on error).
- Decode the `transa`/`transb`/`uplo`/`side`/`diag` character arguments.
- Fill a `blas_arg_t` structure and decide whether to call the single-threaded
  or multi-threaded execution path based on problem size vs.
  `SMP_THRESHOLD_MIN / GEMM_MULTITHREAD_THRESHOLD`.
- For single-threaded path: call the appropriate driver function directly.
- For multi-threaded path: call `GEMM_THREAD` / the BLAS server.

Key files:

| File | Operation |
|------|-----------|
| `gemm.c` | Matrix–matrix multiply (GEMM) — most important file |
| `gemv.c` | Matrix–vector multiply |
| `trsm.c`, `trsv.c` | Triangular solve |
| `syrk.c`, `syr2k.c` | Symmetric rank-k update |
| `dot.c`, `axpy.c` | Level-1 vector operations |

### Driver Layer (`driver/`)

#### `driver/level3/`

Contains the **Goto algorithm** implementation (`level3.c`) used by all
Level-3 operations.  The algorithm partitions the matrix computation into
blocks sized `GEMM_P × GEMM_Q × GEMM_R`, packing panels into contiguous
buffers to maximise cache reuse:

1. Pack a column panel of `B` into a contiguous buffer (`OCOPY`).
2. For each row panel of `A`, pack it (`ICOPY`), then call the inner kernel
   (`GEMM_KERNEL`) to accumulate the result into `C`.

`gemm.c` in this directory includes either `level3.c` (single-threaded) or
`level3_thread.c` (multi-threaded) depending on the `THREADED_LEVEL3` macro.

#### `driver/level2/`

Template drivers for Level-2 operations (GEMV, GER, SYMV, SPMV, TRMV, TRSV,
…).  Thread-dispatch variants have a `_thread` suffix.

#### `driver/others/`

Infrastructure components:

| File | Purpose |
|------|---------|
| `blas_server.c` | pthreads-based SMP work queue |
| `blas_server_omp.c` | OpenMP SMP backend |
| `blas_server_win32.c` | Windows threading backend |
| `dynamic.c` | x86/x86-64 runtime CPU dispatch (`DYNAMIC_ARCH`) |
| `dynamic_arm64.c` | ARM64 runtime dispatch |
| `dynamic_power.c` | POWER runtime dispatch |
| `dynamic_riscv64.c` | RISC-V runtime dispatch |
| `dynamic_zarch.c` | z/Architecture runtime dispatch |
| `memory.c` | Aligned memory allocator for packing buffers |
| `init.c` | Library initialisation (`__attribute__((constructor))`) |
| `parameter.c` | Runtime tuning knobs |
| `openblas_set_num_threads.c` | `openblas_set_num_threads()` |
| `openblas_get_num_threads.c` | `openblas_get_num_threads()` |

### Kernel Layer (`kernel/<arch>/`)

Each architecture directory contains:

- **`KERNEL.<CPU>`** — a list of `make` variable assignments that map abstract
  kernel roles (e.g. `DGEMMKERNEL`, `DGEMM_INCOPY`, `DTRMM_KERNEL`) to
  specific source files.
- **`KERNEL`** — default (oldest / most-conservative) kernel selection file.
- **`KERNEL.generic`** — C-only fallback used when no optimised kernel exists.
- Assembly (`.S`), C with intrinsics (`.c`), or pure C (`.c`) kernels for each
  precision × transpose variant.

For `DYNAMIC_ARCH=1` builds, multiple CPU-specific kernel sets are compiled and
linked into a single library; the correct set is selected at runtime by
`dynamic.c` using CPUID.

## Preprocessor-Based Polymorphism

OpenBLAS uses heavy compile-time polymorphism through macros rather than
runtime polymorphism.  The same `interface/gemm.c` file is compiled ~10 times
with different flag combinations:

| Macro | Meaning |
|-------|---------|
| `SINGLE` / `DOUBLE` | Real single / double precision |
| `COMPLEX` | Complex numbers (paired with `SINGLE` or `DOUBLE`) |
| `BFLOAT16` | 16-bit brain float |
| `XDOUBLE` | Extended (80-bit) double |
| `GEMM3M` | Winograd 3-multiply complex GEMM algorithm |
| `INTERFACE64` | 64-bit (ILP64) integer ABI |

The macro `BLASFUNC(name)` expands to the correct symbol name
(e.g. `dgemm_`, `DGEMM`, `cblas_dgemm`) depending on whether the Fortran,
CBLAS, or C interface is being compiled.

## Data-Type Mapping

| Precision | Real type | Complex type | Fortran prefix | CBLAS prefix |
|-----------|-----------|--------------|----------------|--------------|
| Half (BF16) | `bfloat16` | — | `sb` | `cblas_sb` |
| Single | `float` | `float[2]` | `s` / `c` | `cblas_s` / `cblas_c` |
| Double | `double` | `double[2]` | `d` / `z` | `cblas_d` / `cblas_z` |
| Extended | `long double` | `long double[2]` | `q` / `x` | — |

## Threading Model

OpenBLAS supports three threading backends, selected at build time:

| Backend | Macro | Description |
|---------|-------|-------------|
| pthreads (default) | `USE_THREAD` | POSIX threads, portable |
| OpenMP | `USE_OPENMP` | Uses `#pragma omp parallel` via `blas_server_omp.c` |
| Windows threads | `OS_WINNT` | `CreateThread` via `blas_server_win32.c` |

The threading logic:

1. The interface layer decides whether to parallelise based on problem size.
2. It calls the **BLAS server** (`blas_server.c`), which dispatches work items
   to a pool of pre-created worker threads.
3. Each worker calls the same single-threaded kernel function with a different
   sub-range of the output matrix.
4. Workers synchronise on a spin-lock barrier before returning to the caller.

Key structures are defined in `common_thread.h`:

- `blas_arg_t` — packaged arguments for a kernel invocation.
- `blas_pool_t` — thread pool state.
- `gotoblas_t` — per-CPU function-pointer table (used for `DYNAMIC_ARCH`).

## CPU Detection and Dynamic Dispatch

The CPU detection pipeline runs **at build time**:

```
Makefile.prebuild
    └─ compiles & runs getarch (using cpuid_<arch>.c)
           └─ writes Makefile.conf  (sets TARGET, CORE, ARCH, ...)
                   └─ Makefile.conf included by Makefile.system
                          └─ selects kernel/$(ARCH)/KERNEL.$(CPU)
```

For `DYNAMIC_ARCH=1`, a second path operates **at runtime**:

```
Library initialisation (driver/others/init.c)
    └─ calls gotoblas_init()  (driver/others/dynamic.c)
           └─ reads CPUID → selects gotoblas_<CPU> function table
                  └─ fills global gotoblas pointer used by all kernels
```

## LAPACK Integration

OpenBLAS ships two sources of LAPACK:

1. **`lapack-netlib/`** — full Fortran reference implementation, compiled as-is.
2. **`lapack/`** — hand-optimised C replacements for performance-critical routines
   (LU/Cholesky factorisation, triangular solve/inverse).

The optimised routines in `lapack/` override the Netlib equivalents by being
placed earlier in the link order.  They call back into BLAS (specifically GEMM
and TRSM) to leverage the highly tuned kernels.

Optionally, **ReLAPACK** (`relapack/`) provides recursive LAPACK implementations
that achieve better cache efficiency for large problems by splitting work into
smaller tiles that are solved recursively.

## Key Entry Points for Reverse Engineering

| Starting point | What to look at next |
|----------------|---------------------|
| `cblas.h` | Public C API — all function signatures |
| `interface/gemm.c` | GEMM interface — representative of all Level-3 entry points |
| `driver/level3/level3.c` | Goto algorithm — the core computation loop |
| `kernel/x86_64/KERNEL.HASWELL` | Kernel selection for a concrete CPU |
| `kernel/x86_64/dgemm_kernel_4x8_haswell.S` | Inner assembly loop |
| `driver/others/dynamic.c` | Runtime CPU dispatch |
| `driver/others/blas_server.c` | Threading internals |
| `common_param.h` | `gotoblas_t` struct — the per-CPU dispatch table |
| `param.h` | GEMM blocking parameters per CPU |
