test_that("Sinkhorn preserves marginals", {

  set.seed(1)

  n <- 6
  m <- 5

  a <- runif(n)
  b <- runif(m)

  C <- matrix(runif(n*m), n, m)

  sol <- sinkhorn_log(a, b, C, epsilon = 0.5)

  P <- sol$plan

  a_norm <- a / sum(a)
  b_norm <- b / sum(b)

  expect_true(abs(sum(P) - 1) < 1e-10)
  expect_true(max(abs(rowSums(P) - a_norm)) < 1e-6)
  expect_true(max(abs(colSums(P) - b_norm)) < 1e-6)

})
