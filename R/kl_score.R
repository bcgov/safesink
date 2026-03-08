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
