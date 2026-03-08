#' Construct a numerically safe probability vector
#'
#' Floors small values in a numeric vector and renormalizes the result so that
#' the entries sum to one. This is useful for constructing numerically stable
#' marginal distributions for optimal transport solvers.
#'
#' Each element of `x` is replaced by `max(x_i, eps)` and the resulting vector
#' is then normalized to sum to one.
#'
#' @param x Numeric vector.
#' @param eps Small positive floor applied to each element of `x`.
#'
#' @return Numeric vector with all entries at least `eps` and total sum equal to
#' one.
#'
#' @export
make_safe <- function(x, eps = 1e-12) {
  stopifnot(is.numeric(x))
  stopifnot(all(is.finite(x)))
  x <- pmax(x, eps)
  x / sum(x)
}
