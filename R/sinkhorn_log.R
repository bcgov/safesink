#' Entropy-regularized optimal transport via log-domain Sinkhorn
#'
#' Computes the entropy-regularized optimal transport plan between two
#' marginal distributions using a numerically stable log-domain
#' implementation of the Sinkhorn algorithm.
#'
#' The solver finds a transport matrix \eqn{P} that minimizes
#'
#' \deqn{
#' \sum_{ij} C_{ij} P_{ij} +
#' \varepsilon \sum_{ij} P_{ij}(\log P_{ij} - 1)
#' }
#'
#' subject to the marginal constraints
#'
#' \deqn{
#' \sum_j P_{ij} = a_i, \qquad
#' \sum_i P_{ij} = b_j
#' }
#'
#' The algorithm is implemented in log space to avoid numerical
#' underflow when costs are large or the entropy parameter is small.
#'
#' @param a Numeric vector representing the origin marginal distribution.
#' @param b Numeric vector representing the destination marginal distribution.
#' @param C Numeric cost matrix.
#' @param epsilon Positive scalar entropy regularization parameter.
#' @param max_iter Maximum number of Sinkhorn iterations.
#' @param tol Convergence tolerance for marginal error.
#' @param verbose Logical; if `TRUE`, prints iteration diagnostics.
#'
#' @return A list containing:
#' \describe{
#'   \item{plan}{Optimal transport matrix.}
#'   \item{log_u}{Log dual potentials for origin marginals.}
#'   \item{log_v}{Log dual potentials for destination marginals.}
#'   \item{iterations}{Number of Sinkhorn iterations performed.}
#'   \item{converged}{Logical indicating whether convergence was achieved.}
#' }
#'
#' @export
sinkhorn_log <- function(a, b, C, epsilon,
                         max_iter = 5000,
                         tol = 1e-9,
                         verbose = FALSE) {

  stopifnot(is.matrix(C))
  stopifnot(is.numeric(a), is.numeric(b))
  stopifnot(length(a) == nrow(C))
  stopifnot(length(b) == ncol(C))
  stopifnot(epsilon > 0)
  stopifnot(sum(a) > 0, sum(b) > 0)

  # normalize marginals
  a <- a / sum(a)
  b <- b / sum(b)

  n <- length(a)
  m <- length(b)

  log_a <- ifelse(a > 0, log(a), -Inf)
  log_b <- ifelse(b > 0, log(b), -Inf)

  logK <- -C / epsilon

  log_u <- rep(0, n)
  log_v <- rep(0, m)

  converged <- FALSE
  iter <- max_iter

  for (i in seq_len(max_iter)) {

    row_terms <- logK + matrix(log_v, n, m, byrow = TRUE)
    lse_rows <- row_lse(row_terms)

    log_u_new <- log_a - lse_rows
    log_u_new[a == 0] <- -Inf
    log_u <- log_u_new

    col_terms <- logK + matrix(log_u, n, m, byrow = FALSE)
    lse_cols <- col_lse(col_terms)

    log_v_new <- log_b - lse_cols
    log_v_new[b == 0] <- -Inf
    log_v <- log_v_new

    if (i %% 50 == 0) {

      logP <- outer(log_u, log_v, "+") + logK
      P <- exp(logP)

      err <- max(
        max(abs(rowSums(P) - a)),
        max(abs(colSums(P) - b))
      )

      if (verbose) {
        cat("Iter", i, "err:", err, "\n")
      }

      if (err < tol) {
        converged <- TRUE
        iter <- i
        break
      }
    }
  }

  logP <- outer(log_u, log_v, "+") + logK
  P <- exp(logP)

  list(
    plan = P,
    log_u = log_u,
    log_v = log_v,
    iterations = iter,
    converged = converged
  )
}
