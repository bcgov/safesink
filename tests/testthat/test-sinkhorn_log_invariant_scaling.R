test_that("Sinkhorn is invariant to scaling of a and b", {

  set.seed(1)

  n <- 4
  m <- 6
  eps <- .5
  a <- runif(n)
  b <- runif(m)
  C <- matrix(runif(n * m), n, m)

  sol1 <- sinkhorn_log(a, b, C, eps)$plan
  sol2 <- sinkhorn_log(2*a, 3*b, C, eps)$plan

  expect_lt(max(abs(sol1 - sol2)), 1e-8)
})
