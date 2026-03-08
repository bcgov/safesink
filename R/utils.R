#' Numerically stable log-sum-exp
#'
#' Computes \eqn{\log\left(\sum_i e^{x_i}\right)} in a numerically stable
#' way by subtracting the maximum element before exponentiation.
#'
#' This helper avoids overflow when working with large values in
#' log-space calculations, as required by the log-domain Sinkhorn
#' implementation.
#'
#' @param x Numeric vector.
#'
#' @return Scalar numeric value representing `log(sum(exp(x)))`.
#'
#' @keywords internal
log_sum_exp <- function(x) {
  stopifnot(is.numeric(x))
  xmax <- max(x)
  if (is.infinite(xmax)) return(-Inf)
  xmax + log(sum(exp(x - xmax)))
}

#' Row-wise log-sum-exp
#'
#' Computes the log-sum-exp of each row of a matrix. Uses
#' `matrixStats::rowLogSumExps()` when the matrixStats package is
#' available, otherwise falls back to a base R implementation.
#'
#' This function is used internally by the log-domain Sinkhorn solver
#' for stable row updates.
#'
#' @param M Numeric matrix.
#'
#' @return Numeric vector containing the log-sum-exp of each row.
#'
#' @keywords internal
row_lse <- function(M) {
  stopifnot(is.matrix(M))
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowLogSumExps(M)
  } else {
    apply(M, 1, log_sum_exp)
  }
}

#' Column-wise log-sum-exp
#'
#' Computes the log-sum-exp of each column of a matrix. Uses
#' `matrixStats::colLogSumExps()` when the matrixStats package is
#' available, otherwise falls back to a base R implementation.
#'
#' This function is used internally by the log-domain Sinkhorn solver
#' for stable column updates.
#'
#' @param M Numeric matrix.
#'
#' @return Numeric vector containing the log-sum-exp of each column.
#'
#' @keywords internal
col_lse <- function(M) {
  stopifnot(is.matrix(M))
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colLogSumExps(M)
  } else {
    apply(M, 2, log_sum_exp)
  }
}
