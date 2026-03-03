# Linear Algebra Algorithms

This page lists all linear algebra algorithms implemented in OpenBLAS along with their source directories.

---

## BLAS Level 1 — Vector Operations

These routines operate on individual vectors.

| Algorithm | Description | Directory |
|-----------|-------------|-----------|
| `asum` | Sum of absolute values of a vector | `interface/` |
| `axpy` | Vector addition: y = α·x + y | `interface/` |
| `axpby` | Vector addition: y = α·x + β·y | `interface/` |
| `copy` | Copy a vector | `interface/` |
| `dot` | Dot product of two vectors | `interface/` |
| `dsdot` | Double-precision dot product of single-precision vectors | `interface/` |
| `sdsdot` | Single-precision dot product with scalar accumulation | `interface/` |
| `imax` | Index of the element with the maximum absolute value | `interface/` |
| `nrm2` | Euclidean norm of a vector | `interface/` |
| `rot` | Apply a Givens rotation to two vectors | `interface/` |
| `rotg` | Generate a Givens rotation | `interface/` |
| `rotm` | Apply a modified Givens rotation | `interface/` |
| `rotmg` | Generate a modified Givens rotation | `interface/` |
| `scal` | Scale a vector by a scalar | `interface/` |
| `sum` | Sum of all elements of a vector | `interface/` |
| `swap` | Swap two vectors | `interface/` |

Kernel implementations for each routine exist per architecture under `kernel/<arch>/`
(e.g. `kernel/x86_64/`, `kernel/arm64/`, `kernel/generic/`, etc.).

---

## BLAS Level 2 — Matrix-Vector Operations

These routines perform matrix-vector products and rank-1 or rank-2 updates.

| Algorithm | Description | Directory |
|-----------|-------------|-----------|
| `gbmv` | General banded matrix-vector multiply | `interface/`, `driver/level2/` |
| `gemv` | General matrix-vector multiply | `interface/`, `driver/level2/` |
| `ger` | Rank-1 update of a general matrix | `interface/`, `driver/level2/` |
| `sbmv` / `sbgemv` | Symmetric banded matrix-vector multiply | `interface/`, `driver/level2/` |
| `spmv` | Packed symmetric matrix-vector multiply | `interface/`, `driver/level2/` |
| `spr` | Rank-1 update of a packed symmetric matrix | `interface/`, `driver/level2/` |
| `spr2` | Rank-2 update of a packed symmetric matrix | `interface/`, `driver/level2/` |
| `symv` | Symmetric matrix-vector multiply | `interface/`, `driver/level2/` |
| `syr` | Rank-1 update of a symmetric matrix | `interface/`, `driver/level2/` |
| `syr2` | Rank-2 update of a symmetric matrix | `interface/`, `driver/level2/` |
| `tbmv` | Triangular banded matrix-vector multiply | `interface/`, `driver/level2/` |
| `tbsv` | Triangular banded system solve | `interface/`, `driver/level2/` |
| `tpmv` | Packed triangular matrix-vector multiply | `interface/`, `driver/level2/` |
| `tpsv` | Packed triangular system solve | `interface/`, `driver/level2/` |
| `trmv` | Triangular matrix-vector multiply | `interface/`, `driver/level2/` |
| `trsv` | Triangular system solve | `interface/`, `driver/level2/` |

Complex-valued variants (prefixed with `z`) for most of the routines above, plus Hermitian
variants (`zhbmv`, `zhemv`, `zher`, `zher2`, `zhpmv`, `zhpr`, `zhpr2`), are implemented
in `interface/`.

---

## BLAS Level 3 — Matrix-Matrix Operations

These routines perform matrix-matrix products and related operations.

| Algorithm | Description | Directory |
|-----------|-------------|-----------|
| `gemm` | General matrix-matrix multiply: C = α·A·B + β·C | `interface/`, `driver/level3/` |
| `gemm_batch` | Batched general matrix-matrix multiply | `interface/`, `driver/level3/` |
| `gemm_batch_strided` | Strided batched general matrix-matrix multiply | `interface/`, `driver/level3/` |
| `gemmt` | General matrix-matrix multiply with triangular result | `interface/`, `driver/level3/` |
| `symm` | Symmetric matrix-matrix multiply | `interface/`, `driver/level3/` |
| `syrk` | Symmetric rank-k update | `interface/`, `driver/level3/` |
| `syr2k` | Symmetric rank-2k update | `interface/`, `driver/level3/` |
| `trsm` | Triangular system solve with multiple right-hand sides | `interface/`, `driver/level3/` |

Complex-valued variants (prefixed with `z`) and Hermitian variants (`hemm`, `herk`, `her2k`)
are included via the complex interface in `interface/`.

---

## Extension Functions

Non-standard routines provided by OpenBLAS for additional functionality.

| Algorithm | Description | Directory |
|-----------|-------------|-----------|
| `geadd` | Matrix addition: B = α·A + β·B | `interface/` |
| `imatcopy` | In-place matrix transposition / scaling | `interface/` |
| `omatcopy` | Out-of-place matrix transposition / scaling | `interface/` |
| `bf16dot` | Dot product of BFloat16 vectors | `interface/` |
| `tobf16` | Convert single/double-precision to BFloat16 | `interface/` |
| `bf16to` | Convert BFloat16 to single/double-precision | `interface/` |
| `sbgemmt` | BFloat16 general matrix-matrix multiply with triangular result | `interface/` |

---

## Optimised LAPACK Routines

OpenBLAS provides its own optimised implementations of a selected set of LAPACK
routines. The kernel code lives in `lapack/<routine>/` while the high-level
interface (callable from C or Fortran) lives in `interface/lapack/`.

| Algorithm | Description | Kernel directory | Interface directory |
|-----------|-------------|-----------------|---------------------|
| `getf2` | LU factorisation (unblocked) | `lapack/getf2/` | `interface/lapack/` |
| `getrf` | LU factorisation (blocked, parallel) | `lapack/getrf/` | `interface/lapack/` |
| `getrs` | Solve a system using an LU factorisation | `lapack/getrs/` | `interface/lapack/` |
| `gesv` | Complete linear system solve (LU) | — | `interface/lapack/` |
| `laswp` | Apply row permutations to a matrix | `lapack/laswp/` | `interface/lapack/` |
| `laed3` | Compute the secular equation (divide-and-conquer EVD helper) | `lapack/laed3/` | `interface/lapack/` |
| `lauu2` | Compute U·U^T or L^T·L — unblocked | `lapack/lauu2/` | `interface/lapack/` |
| `lauum` | Compute U·U^T or L^T·L — blocked, parallel | `lapack/lauum/` | `interface/lapack/` |
| `potf2` | Cholesky factorisation (unblocked) | `lapack/potf2/` | `interface/lapack/` |
| `potrf` | Cholesky factorisation (blocked, parallel) | `lapack/potrf/` | `interface/lapack/` |
| `potri` | Matrix inverse via Cholesky factorisation | — | `interface/lapack/` |
| `trti2` | Triangular matrix inverse (unblocked) | `lapack/trti2/` | `interface/lapack/` |
| `trtri` | Triangular matrix inverse (blocked, parallel) | `lapack/trtri/` | `interface/lapack/` |
| `trtrs` | Triangular system solve with multiple right-hand sides | `lapack/trtrs/` | `interface/lapack/` |

All of the above have complex-valued (`z`-prefixed) counterparts implemented in the same directories.

`laswp` also contains architecture-specific kernel optimisations under
`lapack/laswp/<arch>/` (e.g. `lapack/laswp/x86_64/`, `lapack/laswp/arm64/`).

---

## ReLAPACK — Recursive LAPACK Algorithms

OpenBLAS bundles [ReLAPACK](https://github.com/HPAC/ReLAPACK), a collection of
recursive variants of LAPACK compute kernels that often outperform the standard
blocked implementations.  All sources live under `relapack/src/`.

| Algorithm | Description |
|-----------|-------------|
| `[sdcz]getrf` | Recursive LU factorisation |
| `[sdcz]potrf` | Recursive Cholesky factorisation |
| `[sdcz]trtri` | Recursive triangular matrix inverse |
| `[sdcz]trsyl` / `[sdcz]trsyl_rec2` | Recursive triangular Sylvester equation solver |
| `[sdcz]lauum` | Recursive computation of U·U^T or L^T·L |
| `[sd]sytrf` / `[sdcz]sytrf_rec2` | Recursive symmetric/Hermitian indefinite factorisation |
| `[sd]sytrf_rook` / `[sdcz]sytrf_rook_rec2` | Recursive rook-pivoting symmetric/Hermitian factorisation |
| `[cz]hetrf` / `[cz]hetrf_rec2` | Recursive Hermitian indefinite factorisation |
| `[cz]hetrf_rook` / `[cz]hetrf_rook_rec2` | Recursive rook-pivoting Hermitian factorisation |
| `[sdcz]gbtrf` | Recursive LU factorisation of a banded matrix |
| `[sdcz]pbtrf` | Recursive Cholesky factorisation of a banded matrix |
| `[sd]sygst` / `[cz]hegst` | Recursive reduction to standard symmetric eigenproblem |
| `[sdcz]tgsyl` | Recursive triangular generalised Sylvester equation solver |
| `[sdcz]gemmt` | Recursive triangular matrix-matrix multiply accumulate |

Prefixes: `s` = single, `d` = double, `c` = complex single, `z` = complex double.

---

## Kernel Implementations

The low-level compute kernels are implemented per CPU architecture and live under
`kernel/<arch>/`.  Available architecture directories are:

| Directory | Architecture |
|-----------|-------------|
| `kernel/generic/` | Portable C reference kernels |
| `kernel/x86_64/` | 64-bit x86 (SSE, AVX, AVX-512, …) |
| `kernel/x86/` | 32-bit x86 |
| `kernel/arm64/` | 64-bit ARM (NEON, SVE, …) |
| `kernel/arm/` | 32-bit ARM |
| `kernel/power/` | IBM POWER / PowerPC (VSX, …) |
| `kernel/zarch/` | IBM Z / s390x |
| `kernel/riscv64/` | RISC-V 64-bit (RVV, …) |
| `kernel/mips64/` | 64-bit MIPS |
| `kernel/mips/` | 32-bit MIPS |
| `kernel/ia64/` | Intel Itanium |
| `kernel/alpha/` | Alpha |
| `kernel/loongarch64/` | LoongArch 64-bit |
| `kernel/csky/` | C-SKY |
| `kernel/e2k/` | Elbrus |
| `kernel/sparc/` | SPARC |

Each architecture directory contains optimised kernels for GEMM, GEMV, dot, axpy,
copy, scal, and other Level-1/2/3 operations.
