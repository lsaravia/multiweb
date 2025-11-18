#' Estimate the dynamic stability threshold of an ecological network
#'
#' This function estimates a stability threshold for an ecological network
#' using maximum eigenvalue simulations derived from the [calc_QSS()]
#' function. The goal is to approximate the dynamical stability of a system
#' *without requiring knowledge of the actual interaction strengths or the real
#' Jacobian matrix*.
#'
#' ## Conceptual Background
#'
#' Real ecosystems persist over ecological timescales, which implies that their
#' underlying community matrices (Jacobians) must be locally stable—that is,
#' they produce a negative real part of the maximum eigenvalue. Because
#' interaction strengths and self-regulation parameters are rarely known,
#' stability must be inferred indirectly.
#'
#' The [calc_QSS()] function provides such an indirect method by
#' generating random Jacobian matrices whose *sign structure* is determined by
#' the interaction network. Interaction magnitudes are drawn from a uniform
#' distribution whose limits are defined by the parameters \code{negative},
#' \code{positive}, and \code{selfDamping} (Assumed with default values in this function).
#'
#' The QSS value is the proportion of simulated matrices that are locally
#' stable (i.e., have negative maximum eigenvalues). A higher QSS implies that
#' the network topology is more compatible with stable dynamics.
#'
#' ## Purpose of this Function
#'
#' Because the true amount of self-regulation in a real system is unknown,
#' this function searches for the value of \code{selfDamping} for which a
#' target proportion of negative eigenvalues (typically 5%) is produced. This
#' "stability threshold" represents the minimal level of self-regulation
#' required for the network to behave as a persistently stable system.
#'
#' This approach follows the idea that:
#'
#' * Real food webs must be dynamically stable.
#' * Stability requires some level of intraspecific self-regulation.
#' * The amount of self-regulation needed *depends on network structure*.
#'
#' By interpolating the crossing point where the QSS proportion reaches a
#' target value (e.g., 0.05), we obtain an estimate of stability that can be
#' compared across networks—in the same spirit as structural stability,
#' modularity–stability relationships, and feasibility analysis.
#'
#'
#'
#' For details, see References below.
#'
#' ## Parameters
#'
#' @param g An \code{igraph} ecological network.
#' @param nsim Integer. Number of QSS simulations per self-damping value.
#' @param target_prob Numeric. Target proportion of negative eigenvalues
#'   (default = 0.05).
#' @param sd_start Numeric. Starting value of \code{selfDamping} (default = -1).
#' @param sd_step Numeric. Increment used to decrease \code{selfDamping}
#'   (default = -0.1).
#' @param sd_min Numeric. Minimum (most negative) value allowed during search.
#' @param parallel Logical. Whether to run QSS simulations in parallel.
#' @param verbose Logical. Print diagnostic messages.
#'
#' ## Returns
#'
#' A list with:
#'   \item{threshold}{Interpolated \code{selfDamping} value where the target
#'         proportion is reached.}
#'   \item{sd_lower}{Lower bracketing self-damping value.}
#'   \item{sd_upper}{Upper bracketing self-damping value.}
#'   \item{prob_lower}{Proportion of negative eigenvalues at \code{sd_lower}.}
#'   \item{prob_upper}{Proportion of negative eigenvalues at \code{sd_upper}.}
#'   \item{sd_values}{All tested self-damping values.}
#'   \item{probs}{Corresponding proportions of negative eigenvalues.}
#'   \item{qss_interp}{Distribution of maximum eigenvalues at the interpolated
#'         threshold value.}
#'
#' ## References
#'
#' Allesina, S. & Pascual, M. (2008). Network structure, predator–prey modules,
#' and stability in large food webs. *Theoretical Ecology*, **1**, 55–64.
#'
#' Grilli, J., Rogers, T. & Allesina, S. (2016). Modularity and stability in
#' ecological communities. *Nature Communications*, **7**, 12031.
#'
#' Monteiro, A. B. & Del Bianco Faria, L. (2017). Causal relationships between
#' population stability and food-web topology. *Functional Ecology*, **31**,
#' 1294–1300.
#'
#' Borrelli, J. J. (2015). Selection against instability: stable subgraphs are
#' most frequent in empirical food webs. *Oikos*, **124**, 1583–1588.
#'
#' @examples
#' \dontrun{
#' g <- generate_niche(40, 0.1)
#' result <- calc_stability_threshold(g, nsim = 500)
#' plot_stability_curve(result)
#'
#' # calculate the mean of the negative max eigenvalues at the threshold
#' mean(result$qss_raw_interp[result$qss_raw_interp < 0])
#'
#' }
#'
#' @export
calc_stability_threshold <- function(
    g,
    nsim = 1000,
    target_prob = 0.05,
    sd_start = -1,
    sd_step = -0.1,
    sd_min = -20,
    parallel = FALSE,
    verbose = TRUE
) {

  proportions <- c()
  sd_values <- c()

  sd_current <- sd_start
  last_below <- NULL

  if (parallel) {
    library(future.apply)
    plan(multisession)
  }

  repeat {
    if (verbose) cat("Testing selfDamping =", sd_current, "...\n")

    # run QSS with or without parallel backend
    if (parallel) {
      qss_raw <- unlist(
        future_lapply(1:nsim, function(i)
          calc_QSS(g, nsim = 1, returnRaw = TRUE, selfDamping = sd_current)
        )
      )
    } else {
      qss_raw <- calc_QSS(
        g,
        nsim = nsim,
        returnRaw = TRUE,
        selfDamping = sd_current
      )
    }

    prop_neg <- mean(qss_raw < 0)
    if (verbose) cat("  Proportion of negative values:", round(prop_neg, 3), "\n")

    proportions <- c(proportions, prop_neg)
    sd_values <- c(sd_values, sd_current)

    # Check threshold crossing
    if (prop_neg < target_prob) {
      last_below <- list(sd = sd_current, p = prop_neg)
    } else {

      if (!is.null(last_below)) {
        sd1 <- last_below$sd
        p1  <- last_below$p
        sd2 <- sd_current
        p2  <- prop_neg

        # Linear interpolation for threshold
        slope <- (p2 - p1) / (sd2 - sd1)
        sd_interp <- sd1 + (target_prob - p1) / slope

        if (verbose) {
          cat("Threshold crossed between", sd1, "and", sd2, "\n")
          cat("Interpolated selfDamping =", round(sd_interp, 4), "\n")
        }

        # NEW: compute QSS distribution at interpolated threshold
        qss_raw_interp <- calc_QSS(
          g,
          nsim = nsim,
          returnRaw = TRUE,
          selfDamping = sd_interp
        )

        return(list(
          threshold      = sd_interp,
          sd_lower       = sd1,
          sd_upper       = sd2,
          prob_lower     = p1,
          prob_upper     = p2,
          sd_values      = sd_values,
          probs          = proportions,
          qss_last       = qss_raw,
          qss_raw_interp = qss_raw_interp
        ))
      } else {
        stop("Target proportion exceeded before any self-damping value below it.")
      }
    }

    sd_current <- sd_current + sd_step
    if (sd_current < sd_min)
      stop("Reached minimum selfDamping without finding target proportion.")
  }
}

#' Plot stability threshold search curve
#'
#' @param result Output list from [calc_stability_threshold()]}.
#'
#' @export
plot_stability_curve <- function(result) {
  df <- data.frame(
    selfDamping = result$sd_values,
    prob_neg = result$probs
  )

  ggplot(df, aes(selfDamping, prob_neg)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    geom_vline(xintercept = result$threshold, linetype = "dotted", color = "blue") +
    labs(
      x = "Self-damping parameter",
      y = "Proportion of negative eigenvalues",
      title = "Stability Threshold Search",
      subtitle = paste("Interpolated threshold =", round(result$threshold, 4))
    ) +
    theme_minimal(base_size = 14)
}
