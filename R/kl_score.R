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
