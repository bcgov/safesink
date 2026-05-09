#' Validate and enforce alignment of transport inputs
#'
#' Ensures that origin and destination marginals are correctly aligned
#' with the rows and columns of a transport matrix.
#'
#' This function performs *index validation and reordering only*.
#'
#' It is intended for use immediately prior to Sinkhorn-type solvers
#' to guarantee consistent ordering between marginals and matrix inputs.
#'
#' @param a Named numeric vector representing the origin marginal distribution.
#' @param b Named numeric vector representing the destination marginal distribution.
#' @param M Numeric transport cost matrix where rows correspond to the
#'   support of `a` and columns correspond to the support of `b`.
#'   Rectangular matrices are allowed.
#' @return A list containing:
#' \describe{
#'   \item{a}{Reordered origin marginal.}
#'   \item{b}{Reordered destination marginal.}
#'   \item{M}{Reordered matrix aligned to `a` and `b`.}
#' }
#'
#' @export
align_transport_inputs <- function(a, b, M) {
  
  stopifnot(is.numeric(a), is.numeric(b))
  stopifnot(!is.null(names(a)), !is.null(names(b)))
  stopifnot(!is.null(rownames(M)), !is.null(colnames(M)))
  stopifnot(is.matrix(M))
  
  stopifnot(setequal(names(a), rownames(M)))
  stopifnot(setequal(names(b), colnames(M)))
  
  M <- M[names(a), names(b), drop = FALSE]
  
  stopifnot(identical(names(a), rownames(M)))
  stopifnot(identical(names(b), colnames(M)))
  
  list(
    a = a,
    b = b,
    M = M
  )
}
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
