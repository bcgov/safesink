<p align="center">
  <img src="man/figures/logo.png" width="480">
</p>

# SafeSink

**SafeSink** provides a numerically stable implementation of the entropy‑regularized
optimal transport (Sinkhorn) algorithm with safeguards against a common but subtle
source of error: **misaligned marginals and cost matrices**.

Many optimal transport implementations assume that the order of the marginal
vectors exactly matches the row and column ordering of the cost matrix. If these
are misaligned, the solver will typically return a result **without any warning**,
even though the resulting transport plan is incorrect.

SafeSink provides a wrapper (`sinkhorn_aligned()`) that verifies and enforces
alignment using the names of the marginals and cost matrix.

---

## Installation

Install from GitHub using `pak`:

```r
pak::pak("bcgov/safesink")
```

---

## The alignment problem

Optimal transport solvers implicitly assume:

```
a[i] corresponds to row i of C
b[j] corresponds to column j of C
```

If the marginals are named but appear in a different order, this assumption can
silently break.

Example:

```r
a <- c(A = 0.4, B = 0.6)
b <- c(B = 0.5, A = 0.5)

C <- matrix(
  c(0,1,
    1,0),
  nrow = 2,
  dimnames = list(c("A","B"), c("A","B"))
)
```

Many OT solvers will treat the vectors as if their order matches the matrix and
produce an incorrect transport plan with no warning.

---

## Safe alignment

SafeSink prevents this problem by verifying and aligning inputs before solving:

```r
library(safesink)

sol <- sinkhorn_aligned(
  a,
  b,
  C,
  epsilon = 0.5,
  solver = sinkhorn_log
)
```

`sinkhorn_aligned()`:

- verifies that marginals match the cost matrix
- reorders the cost matrix to match marginal ordering
- delegates the computation to the chosen solver

---

## Core functions

| Function | Purpose |
|--------|--------|
| `sinkhorn_log()` | Log‑domain Sinkhorn solver |
| `sinkhorn_aligned()` | Safe wrapper ensuring correct marginal alignment |
| `check_transport()` | Validate marginal constraints and numerical sanity |
| `compare_transport()` | Compare two transport plans |
| `kl_score()` | KL divergence between observed and predicted transition matrices |
| `make_safe()` | Construct numerically safe marginal distributions |

---

## Example

```r
sol <- sinkhorn_log(a, b, C, epsilon = 0.5)

check_transport(sol$plan, a, b)
```

---

## Design goals

SafeSink focuses on:

* **Numerical stability** via a log‑domain Sinkhorn implementation
* **Safety** by preventing silent marginal misalignment
* **Diagnostics** for validating and comparing transport solutions

The package is intentionally lightweight and designed primarily as research
infrastructure for experiments involving optimal transport models.
