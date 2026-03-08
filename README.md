# SafeSink

SafeSink provides a numerically stable implementation of the
entropy-regularized optimal transport (Sinkhorn) algorithm.

The package focuses on robustness and diagnostics for research
applications comparing alternative transport solvers.

## Main functions

- `sinkhorn_log()` – log-domain Sinkhorn solver
- `sinkhorn_aligned()` – safe wrapper aligning marginals with cost matrix
- `check_transport()` – verify marginal constraints
- `compare_transport()` – compare transport plans
- `kl_score()` – KL divergence between transition matrices
- `make_safe()` – construct numerically safe marginals

## Example

```r
library(safesink)

sol <- sinkhorn_log(a, b, C, epsilon = 0.5)

check_transport(sol$plan, a, b)
