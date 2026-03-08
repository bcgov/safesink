#' Align marginals with cost matrix before solving OT
#'
#' Ensures that the origin and destination marginals are correctly aligned
#' with the rows and columns of the cost matrix before calling a Sinkhorn
#' solver. This prevents subtle errors caused by mismatched ordering of
#' marginals and the cost matrix.
#'
#' The function checks that the names of `a` match the row names of `C`
#' and that the names of `b` match the column names of `C`. The cost
#' matrix is then reordered to match the ordering of the marginals before
#' passing the inputs to the specified solver.
#'
#' This wrapper is useful when working with named marginals, where the
#' ordering of elements may not match the ordering of the cost matrix.
#'
#' @param a Named numeric vector representing the origin marginal distribution.
#' @param b Named numeric vector representing the destination marginal distribution.
#' @param C Cost matrix whose row and column names correspond to the elements of
#' `a` and `b`.
#' @param epsilon Entropy regularization parameter passed to the solver.
#' @param solver Function implementing a Sinkhorn solver. The function must
#' accept arguments `(a, b, C, epsilon, ...)` and return a list containing at
#' least a transport matrix `plan`.
#' @param ... Additional arguments passed to the solver.
#'
#' @return The result returned by the specified solver.
#'
#' @export
sinkhorn_aligned <- function(a, b, C, epsilon, solver, ...) {

  stopifnot(!is.null(names(a)))
  stopifnot(!is.null(names(b)))
  stopifnot(!is.null(rownames(C)))
  stopifnot(!is.null(colnames(C)))
  stopifnot(setequal(names(a), rownames(C)))
  stopifnot(setequal(names(b), colnames(C)))

  C <- C[names(a), names(b), drop = FALSE]

  stopifnot(identical(names(a), rownames(C)))
  stopifnot(identical(names(b), colnames(C)))

  solver(a, b, C, epsilon, ...)
}
