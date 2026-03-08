#' Compare two transport plans
#'
#' Computes simple diagnostic statistics comparing two transport
#' matrices. This is useful for validating alternative solvers
#' (e.g., comparing SafeSink results against other implementations
#' such as Python POT).
#'
#' The following statistics are reported:
#'
#' * maximum absolute difference
#' * root mean squared error (RMSE)
#' * correlation between entries
#' * difference in total transported mass
#'
#' @param P1 First transport matrix.
#' @param P2 Second transport matrix.
#'
#' @return A data frame containing comparison diagnostics.
#' The result is printed and also returned invisibly.
#'
#' @importFrom stats cor
#'
#' @export
compare_transport <- function(P1, P2) {
  stopifnot(is.matrix(P1), is.matrix(P2))
  stopifnot(all(dim(P1) == dim(P2)))

  cor_val <- suppressWarnings(cor(as.vector(P1), as.vector(P2)))

  out <- data.frame(
    max_abs_diff = max(abs(P1 - P2)),
    rmse = sqrt(mean((P1 - P2)^2)),
    correlation = cor_val,
    mass_diff = abs(sum(P1) - sum(P2))
  )

  print(out)

  invisible(out)
}
