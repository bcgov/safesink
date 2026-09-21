library(testthat)

# ---- fixtures ---------------------------------------------------------------
# A small 3x3 setup with a known structure.
set.seed(1)
P_obs <- matrix(c(0.5, 0.2, 0.1,
                  0.1, 0.4, 0.2,
                  0.0, 0.1, 0.4), nrow = 3, byrow = TRUE)
P_hat <- matrix(c(0.4, 0.3, 0.1,
                  0.1, 0.3, 0.3,
                  0.05, 0.05, 0.4), nrow = 3, byrow = TRUE)
# Independence baseline: outer product of margins of P_obs (a genuine independence fit)
po_n  <- P_obs / sum(P_obs)
P_ind <- outer(rowSums(po_n), colSums(po_n))
# A symmetric cost matrix (zero diagonal, as the metrics require)
C <- matrix(c(0, 1, 2,
              1, 0, 1,
              2, 1, 0), nrow = 3, byrow = TRUE)

# =============================================================================
test_that("cwtvd_decompose row/col sum to total and match cwtvd()", {
  d <- cwtvd_decompose(P_obs, P_hat, C)
  # decomposition is exact: cells sum to total, rows and cols sum to total
  expect_equal(sum(d$cell), d$total)
  expect_equal(sum(d$row),  d$total)
  expect_equal(sum(d$col),  d$total)
  # total matches the scalar cwtvd() on the same inputs
  expect_equal(d$total, cwtvd(P_obs, P_hat, C))
})

test_that("perfect fit gives zero divergence, every row", {
  # model = observed -> zero everywhere
  expect_equal(cwtvd_decompose(P_obs, P_obs, C)$total, 0)
  expect_equal(unname(rowwise_score(P_obs, P_obs, "cwtvd", C = C, normalize = "conditional")),
               rep(0, nrow(P_obs)))
  expect_equal(unname(rowwise_score(P_obs, P_obs, "kl", normalize = "conditional")),
               rep(0, nrow(P_obs)))
})

test_that("relative improvement is 0 when model == baseline, 1 when model perfect", {
  # model == baseline -> zero improvement over baseline
  ri_eq <- rowwise_rel_improvement(P_obs, P_ind, P_ind, metric = "kl",
                                   normalize = "conditional")
  expect_equal(unname(ri_eq[!is.na(ri_eq)]),
               rep(0, sum(!is.na(ri_eq))))
  # model == observed (perfect) -> improvement 1 (score_model 0, so (base-0)/base = 1)
  ri_perf <- rowwise_rel_improvement(P_obs, P_obs, P_ind, metric = "kl",
                                     normalize = "conditional")
  expect_equal(unname(ri_perf[!is.na(ri_perf)]),
               rep(1, sum(!is.na(ri_perf))))
})

test_that("conditional normalization is invariant to row mass (the key placebo property)", {
  # Scaling one origin's total outflow must NOT change its conditional score,
  # because conditional scoring normalizes each row. This is the property that
  # keeps the exposure gradient free of a mobility-volume confound.
  P_obs2 <- P_obs; P_hat2 <- P_hat
  P_obs2[1, ] <- P_obs2[1, ] * 10   # inflate origin 1's mobility volume 10x
  P_hat2[1, ] <- P_hat2[1, ] * 10

  s1 <- rowwise_score(P_obs,  P_hat,  "kl", normalize = "conditional")
  s2 <- rowwise_score(P_obs2, P_hat2, "kl", normalize = "conditional")
  expect_equal(s1, s2)                      # conditional: unchanged

  # joint, by contrast, SHOULD change (mass-weighted): sanity that the two differ
  j1 <- rowwise_score(P_obs,  P_hat,  "kl", normalize = "joint")
  j2 <- rowwise_score(P_obs2, P_hat2, "kl", normalize = "joint")
  expect_false(isTRUE(all.equal(j1, j2)))
})

test_that("joint rowwise_score matches the decompose row output exactly", {
  expect_equal(rowwise_score(P_obs, P_hat, "kl", normalize = "joint"),
               kl_decompose(P_obs, P_hat)$row)
  expect_equal(rowwise_score(P_obs, P_hat, "cwtvd", C = C, normalize = "joint"),
               cwtvd_decompose(P_obs, P_hat, C)$row)
})

test_that("no-outflow origins return NA under conditional scoring", {
  P0 <- P_obs; P0[2, ] <- 0            # origin 2 has no observed moves
  s <- rowwise_score(P0, P_hat, "kl", normalize = "conditional")
  expect_true(is.na(s[2]))
  expect_false(any(is.na(s[-2])))
})

test_that("rel_improvement NAs origins the baseline already fits (eps guard)", {
  # Construct a baseline that perfectly fits origin 3's (conditional) row:
  P_base <- P_ind
  P_base[3, ] <- P_obs[3, ]            # baseline exactly matches obs row 3
  ri <- rowwise_rel_improvement(P_obs, P_hat, P_base, metric = "kl",
                                normalize = "conditional")
  expect_true(is.na(ri[3]))           # baseline score ~0 -> NA, not Inf
})

test_that("cwtvd requires a cost matrix", {
  expect_error(rowwise_score(P_obs, P_hat, "cwtvd", normalize = "conditional"),
               "C required")
  expect_error(rowwise_score(P_obs, P_hat, "cwtvd", normalize = "joint"),
               "C required")
})

test_that("kl and cwtvd give same ORDERING is NOT assumed (they can differ)", {
  # Not an equality test: just document that the two metrics are distinct objects.
  # (Guards against someone accidentally making cwtvd delegate to kl.)
  sk <- rowwise_score(P_obs, P_hat, "kl", normalize = "conditional")
  sc <- rowwise_score(P_obs, P_hat, "cwtvd", C = C, normalize = "conditional")
  expect_false(isTRUE(all.equal(sk, sc)))
})
