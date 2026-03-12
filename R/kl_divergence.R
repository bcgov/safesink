#' Kullback-Leibler divergence between transition matrices
#'
#' Computes KL divergence between an observed transition matrix and a
#' predicted matrix. Both matrices are normalized to sum to one.
#'
#' @param P_obs Observed transition matrix
#' @param P_hat Predicted transition matrix
#' @param tiny Small floor value used to avoid log(0)
#'
#' @return Scalar KL divergence
#' @export
kl_score <- function(P_obs, P_hat, tiny = 1e-15) {

  # Normalize distributions
  P_obs <- P_obs / sum(P_obs)
  P_hat <- P_hat / sum(P_hat)

  # Bound away from zero
  P_hat <- pmax(P_hat, tiny)

  # Renormalize after flooring
  P_hat <- P_hat / sum(P_hat)

  idx <- P_obs > 0

  sum(P_obs[idx] * log(P_obs[idx] / P_hat[idx]))
}

#' KL divergence decomposition
#'
#' Decomposes KL divergence between an observed transition matrix and
#' a predicted matrix into cell, row, and column contributions.
#'
#' @param P_obs Observed transition matrix
#' @param P_hat Predicted transition matrix
#' @param tiny Small floor value used to avoid log(0)
#'
#' @return A list containing:
#' \describe{
#'   \item{total}{Total KL divergence}
#'   \item{cell}{Matrix of cell contributions}
#'   \item{row}{Row contributions}
#'   \item{col}{Column contributions}
#' }
#' #' @examples
#' P <- matrix(c(.2,.3,.1,.4),2)
#' Phat <- matrix(c(.25,.25,.25,.25),2)
#' kl_score(P, Phat)
#' @export
kl_decompose <- function(P_obs, P_hat, tiny = 1e-15) {

  P_obs <- P_obs / sum(P_obs)
  P_hat <- P_hat / sum(P_hat)

  P_hat <- pmax(P_hat, tiny)
  P_hat <- P_hat / sum(P_hat)

  cell <- matrix(0, nrow(P_obs), ncol(P_obs))

  idx <- P_obs > 0
  cell[idx] <- P_obs[idx] * log(P_obs[idx] / P_hat[idx])

  list(
    total = sum(cell),
    cell = cell,
    row = rowSums(cell),
    col = colSums(cell)
  )
}


