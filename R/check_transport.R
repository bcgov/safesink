#' Check transport plan diagnostics
#'
#' Performs basic diagnostic checks on a transport matrix produced by a
#' Sinkhorn solver. The function verifies that the transport plan satisfies
#' the marginal constraints and reports several numerical diagnostics.
#'
#' The following quantities are reported:
#'
#' * total transported mass
#' * maximum row marginal error
#' * maximum column marginal error
#' * deviation of total mass from one
#' * minimum entry in the transport matrix
#' * presence of negative entries
#' * presence of non-finite values
#'
#' @param P Transport matrix.
#' @param a Origin marginal vector.
#' @param b Destination marginal vector.
#' @param tol Numerical tolerance used when checking for negative entries.
#'
#' @return A data frame containing diagnostic statistics. The result is
#' printed and returned invisibly.
#'
#' @export
check_transport <- function(P, a, b, tol = 1e-8) {

  a <- a / sum(a)
  b <- b / sum(b)

  row_err <- max(abs(rowSums(P) - a))
  col_err <- max(abs(colSums(P) - b))
  mass_err <- abs(sum(P) - 1)

  out <- data.frame(
    mass      = sum(P),
    row_err   = row_err,
    col_err   = col_err,
    mass_err  = mass_err,
    min_value = min(P),
    any_neg   = any(P < -tol),
    any_nan   = any(!is.finite(P))
  )

  print(out)

  invisible(out)
}
