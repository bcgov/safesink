#' Row-wise (per-origin) fit scoring for the graded arm
#'
#' These helpers sit on top of the metric primitives in kl_divergence.R
#' (\code{kl_score}, \code{kl_decompose}) and CWTVD.R (\code{cwtvd},
#' \code{cwtvd_decompose}). They contain no metric math of their own: the
#' conditional branch normalizes each origin's destination row and delegates the
#' actual scoring to the metric primitive, so any change to a metric's
#' convention (flooring, normalization) propagates here automatically.

#' Row-wise fit score, conditional or joint
#'
#' Returns a per-origin (row) fit score for either metric.
#'
#' The \strong{conditional} version normalizes each origin's destination row to
#' sum to one BEFORE scoring, so each occupation's score reflects only how well
#' the model predicts \emph{its own} destination distribution, independent of how
#' much total mobility mass it carries. This is the object to grade on exposure
#' for the placebo: it isolates fit quality from mobility volume, so the exposure
#' gradient is not contaminated by high-exposure occupations happening to be
#' high-mobility.
#'
#' The \strong{joint} version reproduces the mass-weighted row contribution from
#' \code{kl_decompose()} / \code{cwtvd_decompose()} (each row's share of the total
#' divergence).
#'
#' @param P_obs,P_hat Transition matrices (counts or probabilities), same dims.
#' @param metric "kl" or "cwtvd".
#' @param C Cost matrix (required for cwtvd, ignored for kl). For the conditional
#'   CWTVD of origin \code{i}, row \code{C[i, ]} is used as the destination costs.
#' @param normalize "conditional" (per-origin, volume-independent; default) or
#'   "joint" (mass-weighted row contribution).
#' @param tiny Floor passed through to the KL primitive.
#'
#' @return Numeric vector of length \code{nrow(P_obs)}: per-origin fit score
#'   (0 = perfect fit for that origin, higher = worse). NA for origins with no
#'   observed outflow under conditional normalization.
#' @seealso \code{\link{kl_score}}, \code{\link{cwtvd}},
#'   \code{\link{kl_decompose}}, \code{\link{cwtvd_decompose}}
#' @export
rowwise_score <- function(P_obs, P_hat, metric = c("kl", "cwtvd"),
                          C = NULL, normalize = c("conditional", "joint"),
                          tiny = 1e-15) {
  metric    <- match.arg(metric)
  normalize <- match.arg(normalize)
  P_obs <- as.matrix(P_obs)
  P_hat <- as.matrix(P_hat)

  if (metric == "cwtvd" && is.null(C)) stop("C required for cwtvd")

  if (normalize == "joint") {
    if (metric == "kl") {
      return(kl_decompose(P_obs, P_hat, tiny = tiny)$row)
    } else {
      return(cwtvd_decompose(P_obs, P_hat, C)$row)
    }
  }

  # conditional: normalize each origin's destination row, then delegate scoring
  # to the metric primitive on that single row.
  n  <- nrow(P_obs)
  m  <- ncol(P_obs)
  ro <- rowSums(P_obs)
  rh <- rowSums(P_hat)
  out <- numeric(n)

  for (i in seq_len(n)) {
    if (ro[i] <= 0) { out[i] <- NA_real_; next }        # no destination dist to score
    po <- P_obs[i, ] / ro[i]
    ph <- if (rh[i] > 0) P_hat[i, ] / rh[i] else rep(1 / m, m)  # empty model row -> uniform

    if (metric == "kl") {
      # kl_score renormalizes and floors internally; pass the single row.
      out[i] <- kl_score(po, ph, tiny = tiny)
    } else {
      # cwtvd expects matrices; pass 1-row matrices with this origin's cost row.
      out[i] <- cwtvd(matrix(po, nrow = 1),
                      matrix(ph, nrow = 1),
                      matrix(C[i, ], nrow = 1))
    }
  }
  out
}

#' Row-wise relative improvement over a baseline (e.g. independence)
#'
#' The placebo quantity: for each origin, how much better the fitted metric model
#' predicts its destination distribution than a baseline (independence) model, as
#' a proportion of the baseline's misfit.
#'
#' \deqn{rel_i = (score^{base}_i - score^{model}_i) / score^{base}_i}
#'
#' Positive = the model improves on the baseline for that origin. Grade THIS
#' vector on Eloundou for the placebo. Pre-window the graded relationship should
#' be flat (no purchase); post-window the graded-arm signal is this improvement
#' DROPPING at high exposure (fit degrades where LLM exposure is high). Use the
#' same \code{metric} and \code{normalize} pre- and post-window so the two are the
#' identical estimand.
#'
#' @param P_obs Observed matrix.
#' @param P_hat_model Fitted metric model's predicted matrix.
#' @param P_hat_base Baseline (independence) predicted matrix.
#' @param metric,C,normalize,tiny As in \code{rowwise_score()}.
#' @param eps Floor on the baseline score; origins whose baseline score is at or
#'   below \code{eps} (baseline already fits them near-perfectly) are returned NA
#'   rather than dividing by ~0.
#'
#' @return Numeric vector length \code{nrow(P_obs)}: per-origin relative
#'   improvement. NA where the origin has no outflow, either score is NA, or the
#'   baseline score is non-positive / below \code{eps}.
#' @seealso \code{\link{rowwise_score}}
#' @export
rowwise_rel_improvement <- function(P_obs, P_hat_model, P_hat_base,
                                    metric = c("kl", "cwtvd"), C = NULL,
                                    normalize = c("conditional", "joint"),
                                    tiny = 1e-15, eps = 1e-12) {
  metric    <- match.arg(metric)
  normalize <- match.arg(normalize)
  s_model <- rowwise_score(P_obs, P_hat_model, metric, C, normalize, tiny)
  s_base  <- rowwise_score(P_obs, P_hat_base,  metric, C, normalize, tiny)
  out <- (s_base - s_model) / s_base
  out[is.na(s_base) | is.na(s_model) | s_base <= eps] <- NA_real_
  out
}
