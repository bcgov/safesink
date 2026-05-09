test_that("sinkhorn_log agrees with POT", {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    skip("reticulate not installed")
  }

  reticulate::py_require("POT")

  if (!reticulate::py_module_available("ot")) {
    skip("Python POT not available")
  }

  set.seed(1)

  n <- 5
  m <- 5

  a <- runif(n)
  b <- runif(m)
  C <- matrix(runif(n * m), n, m)

  sol_r <- sinkhorn_log(a, b, C, epsilon = 0.5)
  P_r <- sol_r$plan

  ot <- reticulate::import("ot")
  P_py <- reticulate::py_to_r(
    ot$sinkhorn(a / sum(a), b / sum(b), C, 0.5)
  )

  expect_lt(max(abs(P_r - P_py)), 1e-6)
})

