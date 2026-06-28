#' Lu's Agreement Coefficient for Subjective Judgments
#'
#' @description
#' Computes Lu's (1971) coefficient of agreement for subjective judgments
#' provided by multiple judges using ordered categorical ratings.
#'
#' @details
#' Consider a rating design with \eqn{n} subjects (e.g., items, objectives)
#' and \eqn{m} judges, where each judge assigns an ordered categorical rating
#' (e.g., 1 to 4, 1 to 5) to each subject. Let \eqn{X_{ij}} denote the rating
#' given by judge \eqn{j} to subject \eqn{i}. Lu's (1971) approach proceeds
#' in three main steps:
#'
#' \enumerate{
#'   \item \strong{Empirical probability scale.}
#'   All ratings are pooled across subjects and judges to obtain the empirical
#'   proportions \eqn{p_k} for each response category \eqn{k = 1, \dots, t}.
#'   These proportions define a cumulative distribution \eqn{F(k)} and each
#'   category is mapped to a point \eqn{y_k} on the unit interval corresponding
#'   to the midpoint of its probability interval:
#'   \deqn{
#'     y_k = F(k) - \frac{p_k}{2},
#'   }
#'   yielding a transformed score \eqn{Y_{ij}} for each rating.
#'
#'   \item \strong{Within-subject variance of ratings.}
#'   For each subject \eqn{i}, the mean of the transformed ratings is
#'   \eqn{\bar{Y}_i} and the within-subject sum of squares is
#'   \deqn{
#'     SS_{w,i} = \sum_{j=1}^m (Y_{ij} - \bar{Y}_i)^2.
#'   }
#'   The pooled within-subject variance is then
#'   \deqn{
#'     S_w^2 = \frac{\sum_{i=1}^n SS_{w,i}}{n (m - 1)}.
#'   }
#'
#'   \item \strong{Maximum-disagreement reference (maximum entropy).}
#'   As a reference for complete disagreement under a maximum-entropy model,
#'   Lu assumes all \eqn{t} categories are equally probable,
#'   \eqn{p_k = 1/t}. On the same probability scale, categories are located at
#'   equally spaced midpoints:
#'   \deqn{
#'     y_k^{*} = \frac{k - 0.5}{t}, \quad k = 1, \dots, t.
#'   }
#'   The variance of this discrete uniform distribution on the probability scale
#'   is
#'   \deqn{
#'     S_{w0}^2 = \sum_{k=1}^t (y_k^{*} - \mu_0)^2 \, p_k,
#'   }
#'   where \eqn{\mu_0} is the mean of \eqn{y_k^{*}} and \eqn{p_k = 1/t}.
#' }
#'
#' Lu's coefficient of agreement is defined as
#' \deqn{
#'   A = 1 - \frac{S_w^2}{S_{w0}^2},
#' }
#' so that \eqn{A = 1} indicates perfect agreement, \eqn{A = 0} indicates
#' agreement no better than expected under the maximum-entropy reference
#' (i.e., essentially random use of categories), and \eqn{A < 0} indicates
#' less agreement than expected under that reference model.
#'
#' In addition to the magnitude of agreement, Lu's formulation leads to a
#' chi-square test of the null hypothesis of agreement equal to that expected
#' under the maximum-entropy model, that is \eqn{H_0: A = 0}. Let
#' \eqn{S_w^2} be the observed within-subject variance and \eqn{S_{w0}^2} the
#' reference variance under maximum entropy. Lu shows that, under \eqn{H_0},
#' the statistic
#' \deqn{
#'   T = \frac{S_w^2}{S_{w0}^2} \, n (m - 1)
#' }
#' is approximately distributed as a chi-square with
#' \eqn{df = n (m - 1)} degrees of freedom. A one-sided lower-tail p-value
#' \eqn{p = P(\chi^2_{df} \le T)} evaluates whether the observed agreement
#' is significantly greater than that expected by chance:
#' small values of \eqn{T} (and hence small p-values) indicate
#' \eqn{S_w^2 < S_{w0}^2} and thus \eqn{A > 0}.
#'
#' When \code{group} is provided, the function computes Lu's coefficient
#' separately for each group of judges (e.g., internal vs. external experts)
#' and, if \code{ci = TRUE}, uses a simple MOVER-type (Method Of Variance
#' Estimates Recovery) approach to form confidence intervals for the differences
#' between group-specific coefficients based solely on the endpoints of the
#' individual confidence intervals.
#'
#' Bootstrap confidence intervals for Lu's coefficient are obtained by
#' resampling subjects with replacement and recomputing the coefficient in
#' each bootstrap sample. The \code{ci.method} argument controls how the
#' bootstrap confidence interval is constructed.
#'
#' @param data A numeric matrix or data frame of size \eqn{n \times m},
#'   where rows correspond to subjects (e.g., items, objectives) and columns
#'   correspond to judges. Entries should be ordered categorical ratings
#'   (numeric, character, or factors). All columns are coerced to a common
#'   set of ordered categories.
#' @param group Optional vector indicating groups of judges. Must have length
#'   equal to the number of columns (judges) in \code{data}. It can be a
#'   factor, character, or numeric vector and is internally converted to a
#'   factor. If \code{NULL} (default), all judges are treated as a single
#'   group.
#' @param ci Logical. If \code{TRUE}, bootstrap confidence intervals are
#'   computed for Lu's coefficient (and for group differences when
#'   \code{group} is not \code{NULL}). Default is \code{FALSE}.
#' @param ci.method Character string specifying the bootstrap confidence
#'   interval method. Options are:
#'   \itemize{
#'     \item \code{"perc"} Percentile interval (default).
#'     \item \code{"norm"} Normal-theory bootstrap interval.
#'     \item \code{"bca"} Bias-corrected and accelerated interval (BCa),
#'           using the \pkg{boot} package.
#'   }
#'   Ignored when \code{ci = FALSE}.
#' @param conf.level Confidence level for the intervals, typically between
#'   0 and 1. Default is \code{0.95}.
#' @param nd Integer. Number of decimal places used to round the reported
#'   estimates and interval endpoints. Default is \code{3}.
#' @param na.rm Logical. If \code{TRUE} (default), only subjects with complete
#'   data across all judges (no missing ratings in the row) are retained.
#'   If \code{FALSE} and any missing values are present, the function stops
#'   with an error.
#' @param B Integer. Number of bootstrap resamples used to compute confidence
#'   intervals when \code{ci = TRUE}. Default is \code{1000}.
#'
#' @return
#' If \code{group} is \code{NULL}, returns an object of class
#' \code{"LuAgree"} with components:
#' \itemize{
#'   \item \code{results} A data frame with one row and columns:
#'     \describe{
#'       \item{\code{A}}{Lu's agreement coefficient (rounded to \code{nd} decimals).}
#'       \item{\code{lwr.ci}}{Lower bound of the bootstrap confidence interval
#'         for \code{A} (or \code{NA} if \code{ci = FALSE}).}
#'       \item{\code{upr.ci}}{Upper bound of the bootstrap confidence interval
#'         for \code{A} (or \code{NA} if \code{ci = FALSE}).}
#'       \item{\code{chi.square}}{Value of Lu's chi-square statistic
#'         \eqn{T = (S_w^2/S_{w0}^2)\,n(m-1)}.}
#'       \item{\code{df}}{Degrees of freedom for the chi-square test,
#'         \eqn{df = n (m - 1)}.}
#'       \item{\code{p.value}}{One-sided lower-tail p-value
#'         \eqn{P(\chi^2_{df} \le T)} for testing \eqn{H_0: A = 0}
#'         against \eqn{H_1: A > 0}.}
#'     }
#'   \item \code{n_subjects} Number of subjects used in the analysis.
#'   \item \code{n_judges} Number of judges.
#'   \item \code{categories} The ordered rating categories used.
#'   \item \code{p_empirical} Empirical category probabilities (información
#'         auxiliar; no se imprime por defecto).
#'   \item \code{call} The matched function call.
#' }
#'
#' If \code{group} is not \code{NULL}, returns a list of class
#' \code{"LuAgree"} with components:
#' \itemize{
#'   \item \code{group_results} A data frame summarising Lu's coefficient
#'         by group of judges, with columns:
#'         \describe{
#'           \item{\code{group}}{Group label.}
#'           \item{\code{A}}{Lu's agreement coefficient in that group
#'             (rounded).}
#'           \item{\code{lwr.ci}}{Lower bound of the bootstrap confidence
#'             interval for \code{A} (or \code{NA} if \code{ci = FALSE}).}
#'           \item{\code{upr.ci}}{Upper bound of the bootstrap confidence
#'             interval for \code{A} (or \code{NA} if \code{ci = FALSE}).}
#'           \item{\code{chi.square}}{Lu's chi-square statistic for that group.}
#'           \item{\code{df}}{Degrees of freedom for the group-specific test.}
#'           \item{\code{p.value}}{One-sided lower-tail p-value for testing
#'             \eqn{H_0: A = 0} in that group.}
#'           \item{\code{n_subjects}}{Number of subjects.}
#'           \item{\code{n_judges}}{Number of judges in the group.}
#'         }
#'   \item \code{comparisons} A data frame of pairwise comparisons between
#'         groups (when \code{ci = TRUE}), including the estimated difference
#'         in coefficients and a MOVER-type confidence interval:
#'         \code{group1}, \code{group2}, \code{diff}, \code{lwr.ci},
#'         \code{upr.ci}. If \code{ci = FALSE} or there is only one group,
#'         this component is \code{NULL}.
#'   \item \code{n_subjects} Number of subjects used in the analysis.
#'   \item \code{categories} Ordered rating categories for the full set of
#'         judges.
#'   \item \code{p_empirical} Empirical category probabilities for the full
#'         panel of judges.
#'   \item \code{call} The matched function call.
#' }
#'
#' @references
#' Lu, K. H. (1971). A measure of agreement among subjective judgments.
#' \emph{Educational and Psychological Measurement}, 31(1), 57-68.
#' https://doi.org/10.1177/001316447103100105
#'
#' Rovinelli, Richard J.; Hambleton, Ronald K. (1976). On the Use of Content
#' Specialists in the Assessment of Criterion-Referenced Test Item Validity.Paper
#' presented at the Annual Meeting of the American Educational Research
#' Association (60th, San Francisco, California, April 19-23, 1976).
#' https://files.eric.ed.gov/fulltext/ED121845.pdf
#'
#' @examples
#' set.seed(123)
#' # 10 subjects, 5 judges, ratings 1-4
#' dat <- matrix(sample(1:4, size = 10 * 5, replace = TRUE),
#'               nrow = 10, ncol = 5)
#'
#' # Single group of judges
#' res1 <- LuAgree(dat)
#' res1
#'
#' # With bootstrap percentile confidence interval
#' \donttest{
#' res2 <- LuAgree(dat, ci = TRUE, B = 500, ci.method = "perc")
#' res2
#' }
#'
#' # Two groups of judges (first 3 vs last 2)
#' group_vec <- c("internal", "internal", "internal", "external", "external")
#' \donttest{
#' res3 <- LuAgree(dat, group = group_vec, ci = TRUE, B = 500)
#' res3
#' }
#'
#' @export
LuAgree <- function(data,
                    group      = NULL,
                    ci         = FALSE,
                    ci.method  = c("perc", "norm", "bca"),
                    conf.level = 0.95,
                    nd         = 3,
                    na.rm      = TRUE,
                    B          = 1000) {

  # --- 0. Argument checks -------------------------------------------------

  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (ncol(data) < 2L) {
    stop("`data` must have at least 2 judges (columns).")
  }

  if (!is.logical(ci) || length(ci) != 1L) {
    stop("`ci` must be a single logical value (TRUE/FALSE).")
  }

  ci.method <- match.arg(ci.method)

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      conf.level <= 0 || conf.level >= 1) {
    stop("`conf.level` must be a single number in (0, 1).")
  }

  nd <- as.integer(nd)
  if (is.na(nd) || nd < 0L) {
    stop("`nd` must be a non-negative integer.")
  }

  if (!is.logical(na.rm) || length(na.rm) != 1L) {
    stop("`na.rm` must be a single logical value (TRUE/FALSE).")
  }

  B <- as.integer(B)
  if (ci && (is.na(B) || B <= 0L)) {
    stop("`B` must be a positive integer when `ci = TRUE`.")
  }

  # --- 1. Handle missing data: complete cases only if na.rm = TRUE --------

  if (na.rm) {
    cc <- stats::complete.cases(data)
    data <- data[cc, , drop = FALSE]
  } else {
    if (anyNA(data)) {
      stop("Missing values found in `data`. Set `na.rm = TRUE` to ",
           "use only complete cases.")
    }
  }

  n <- nrow(data)
  m <- ncol(data)

  if (n < 1L) {
    stop("No valid subjects after applying `na.rm` / complete case filtering.")
  }
  if (m < 2L) {
    stop("`data` must have at least 2 judges (columns).")
  }

  # --- 2. Internal helper: compute Lu's A for a single panel of judges -----

  .Lu_core <- function(df) {
    if (!is.data.frame(df)) df <- as.data.frame(df)

    n_i <- nrow(df)
    m_i <- ncol(df)

    # Collect all non-missing values
    vals <- unlist(df, use.names = FALSE)
    vals <- vals[!is.na(vals)]

    if (length(vals) == 0L) {
      stop("All entries are NA; cannot compute Lu's coefficient.")
    }

    # Determine ordered categories
    if (is.factor(vals)) {
      levs <- levels(vals)
    } else {
      levs <- sort(unique(vals))
    }

    t_cat <- length(levs)
    if (t_cat < 2L) {
      stop("There must be at least 2 rating categories to compute agreement.")
    }

    # Coerce each column to an ordered factor with the same levels
    for (j in seq_len(m_i)) {
      colj <- df[[j]]
      df[[j]] <- factor(colj, levels = levs, ordered = TRUE)
    }

    # Empirical frequencies and probabilities
    all_fac <- unlist(df, use.names = FALSE)
    all_fac <- all_fac[!is.na(all_fac)]

    tab <- table(all_fac)
    p_emp <- rep(0, t_cat)
    names(p_emp) <- as.character(levs)
    p_emp[names(tab)] <- as.numeric(tab)
    p_emp <- p_emp / sum(p_emp)

    # Empirical probability midpoints y_k
    cemp <- cumsum(p_emp)
    y_emp <- cemp - p_emp / 2

    # Map ratings to Y matrix
    Y <- matrix(NA_real_, nrow = n_i, ncol = m_i)
    for (k in seq_along(levs)) {
      Y[df == levs[k]] <- y_emp[k]
    }

    # Check for any remaining NA in Y
    if (anyNA(Y)) {
      stop("Internal error: NA values found in transformed ratings (Y).")
    }

    # Within-subject variance S_w^2 (pooled)
    Sw2_num <- 0
    for (i in seq_len(n_i)) {
      row_i <- Y[i, ]
      mu_i  <- mean(row_i)
      ss_i  <- sum((row_i - mu_i)^2)
      Sw2_num <- Sw2_num + ss_i
    }
    Sw2 <- Sw2_num / (n_i * (m_i - 1L))

    # Maximum-entropy reference variance S_{w0}^2
    p_unif <- rep(1 / t_cat, t_cat)
    y_unif <- (seq_len(t_cat) - 0.5) / t_cat
    mu0    <- sum(y_unif * p_unif)
    Sw0_2  <- sum((y_unif - mu0)^2 * p_unif)

    A_val <- 1 - Sw2 / Sw0_2

    list(
      A           = A_val,
      Sw2         = Sw2,
      Sw0_2       = Sw0_2,
      categories  = levs,
      p_empirical = p_emp,
      n_subjects  = n_i,
      n_judges    = m_i
    )
  }

  # --- 3. Internal helper: bootstrap CI for a single panel -----------------

  .Lu_boot <- function(df, conf.level, B, ci.method) {
    n_i <- nrow(df)
    theta_hat <- .Lu_core(df)$A

    # Bootstrap replicates
    thetas <- numeric(B)
    for (b in seq_len(B)) {
      idx <- sample.int(n_i, n_i, replace = TRUE)
      thetas[b] <- .Lu_core(df[idx, , drop = FALSE])$A
    }

    alpha <- 1 - conf.level

    ci.method <- match.arg(ci.method, c("perc", "norm", "bca"))

    if (ci.method == "perc") {
      # Percentile CI
      lwr <- stats::quantile(thetas, probs = alpha / 2, na.rm = TRUE)
      upr <- stats::quantile(thetas, probs = 1 - alpha / 2, na.rm = TRUE)

    } else if (ci.method == "norm") {
      # Normal-theory CI from bootstrap mean/sd
      se_boot <- stats::sd(thetas, na.rm = TRUE)
      z <- stats::qnorm(1 - alpha / 2)
      lwr <- theta_hat - z * se_boot
      upr <- theta_hat + z * se_boot

    } else {  # "bca"
      if (!requireNamespace("boot", quietly = TRUE)) {
        warning("Package 'boot' is required for ci.method = 'bca'. ",
                "Falling back to percentile interval.")
        lwr <- stats::quantile(thetas, probs = alpha / 2, na.rm = TRUE)
        upr <- stats::quantile(thetas, probs = 1 - alpha / 2, na.rm = TRUE)
      } else {
        boot_fun <- function(data, indices) {
          d <- data[indices, , drop = FALSE]
          .Lu_core(d)$A
        }
        boot_obj <- boot::boot(df,
                               statistic = function(dat, ind) boot_fun(dat, ind),
                               R = B)
        ci_obj <- boot::boot.ci(boot_obj, conf = conf.level, type = "bca")
        if (is.null(ci_obj$bca)) {
          warning("BCa interval not available; falling back to percentile CI.")
          lwr <- stats::quantile(thetas, probs = alpha / 2, na.rm = TRUE)
          upr <- stats::quantile(thetas, probs = 1 - alpha / 2, na.rm = TRUE)
        } else {
          lwr <- ci_obj$bca[4]
          upr <- ci_obj$bca[5]
        }
      }
    }

    list(
      estimate    = theta_hat,
      lwr         = as.numeric(lwr),
      upr         = as.numeric(upr),
      boot.values = thetas
    )
  }

  # --- 4. No group: one panel of judges -----------------------------------

  if (is.null(group)) {
    core <- .Lu_core(data)

    # Chi-square test (Lu, 1971): T = (Sw2/Sw0_2) * n(m-1)
    df_chi <- core$n_subjects * (core$n_judges - 1L)
    Tobs   <- (core$Sw2 / core$Sw0_2) * df_chi
    pval   <- stats::pchisq(Tobs, df = df_chi, lower.tail = TRUE)

    if (ci) {
      boot_res <- .Lu_boot(data,
                           conf.level = conf.level,
                           B          = B,
                           ci.method  = ci.method)
      results_df <- data.frame(
        A          = round(core$A, nd),
        lwr.ci     = round(boot_res$lwr, nd),
        upr.ci     = round(boot_res$upr, nd),
        chi.square = round(Tobs, nd),
        df         = df_chi,
        p.value    = pval
      )
    } else {
      results_df <- data.frame(
        A          = round(core$A, nd),
        lwr.ci     = NA_real_,
        upr.ci     = NA_real_,
        chi.square = round(Tobs, nd),
        df         = df_chi,
        p.value    = pval
      )
    }

    res <- list(
      results     = results_df,
      n_subjects  = core$n_subjects,
      n_judges    = core$n_judges,
      categories  = core$categories,
      p_empirical = core$p_empirical,
      call        = match.call()
    )

    class(res) <- "LuAgree"
    return(res)
  }

  # --- 5. With groups: multiple panels of judges ---------------------------

  if (length(group) != m) {
    stop("Length of `group` must match the number of judges (columns) in `data`.")
  }
  group <- as.factor(group)

  g_levels <- levels(group)
  G <- length(g_levels)

  if (G < 1L) {
    stop("`group` must define at least one group of judges.")
  }

  # Global categories & probabilities (for info)
  global_core <- .Lu_core(data)

  # Per-group results
  A_vec    <- numeric(G)
  lwr_vec  <- rep(NA_real_, G)
  upr_vec  <- rep(NA_real_, G)
  m_vec    <- integer(G)
  chi_vec  <- numeric(G)
  df_vec   <- integer(G)
  pval_vec <- numeric(G)

  for (g_idx in seq_along(g_levels)) {
    g <- g_levels[g_idx]
    cols_g <- which(group == g)
    df_g   <- data[, cols_g, drop = FALSE]
    m_vec[g_idx] <- ncol(df_g)

    core_g <- .Lu_core(df_g)
    A_vec[g_idx] <- core_g$A

    # Chi-square test for this group
    df_chi_g      <- core_g$n_subjects * (core_g$n_judges - 1L)
    Tobs_g        <- (core_g$Sw2 / core_g$Sw0_2) * df_chi_g
    chi_vec[g_idx]  <- Tobs_g
    df_vec[g_idx]   <- df_chi_g
    pval_vec[g_idx] <- stats::pchisq(Tobs_g, df = df_chi_g, lower.tail = TRUE)

    if (ci) {
      boot_g <- .Lu_boot(df_g,
                         conf.level = conf.level,
                         B          = B,
                         ci.method  = ci.method)
      lwr_vec[g_idx] <- boot_g$lwr
      upr_vec[g_idx] <- boot_g$upr
    }
  }

  group_results <- data.frame(
    group      = g_levels,
    A          = round(A_vec, nd),
    lwr.ci     = if (ci) round(lwr_vec, nd) else NA_real_,
    upr.ci     = if (ci) round(upr_vec, nd) else NA_real_,
    chi.square = round(chi_vec, nd),
    df         = df_vec,
    p.value    = pval_vec,
    n_subjects = rep(n, G),
    n_judges   = m_vec,
    stringsAsFactors = FALSE
  )

  # Pairwise comparisons via simple MOVER (only if ci = TRUE and G >= 2)
  comparisons <- NULL
  if (ci && G >= 2L) {
    comp_list <- list()
    k <- 1L

    for (i in 1:(G - 1L)) {
      for (j in (i + 1L):G) {
        theta1 <- A_vec[i]
        theta2 <- A_vec[j]
        L1 <- lwr_vec[i]; U1 <- upr_vec[i]
        L2 <- lwr_vec[j]; U2 <- upr_vec[j]

        diff_ij <- theta1 - theta2

        # Simple MOVER for difference
        Ldiff <- diff_ij - sqrt((theta1 - L1)^2 + (U2 - theta2)^2)
        Udiff <- diff_ij + sqrt((U1 - theta1)^2 + (theta2 - L2)^2)

        comp_list[[k]] <- data.frame(
          group1 = g_levels[i],
          group2 = g_levels[j],
          diff   = round(diff_ij, nd),
          lwr.ci = round(Ldiff, nd),
          upr.ci = round(Udiff, nd),
          stringsAsFactors = FALSE
        )
        k <- k + 1L
      }
    }
    comparisons <- do.call(rbind, comp_list)
  }

  res <- list(
    group_results = group_results,
    comparisons   = comparisons,
    n_subjects    = n,
    categories    = global_core$categories,
    p_empirical   = global_core$p_empirical,
    call          = match.call()
  )

  class(res) <- "LuAgree"
  return(res)
}

#' @export
print.LuAgree <- function(x, ...) {
  cat("\nLu's agreement coefficient (LuAgree)\n\n")

  if (!is.null(x$group_results)) {
    # Grouped case
    cat("By group of judges:\n")
    print(x$group_results, row.names = FALSE, ...)
    if (!is.null(x$comparisons)) {
      cat("\nPairwise comparisons (MOVER-type CIs):\n")
      print(x$comparisons, row.names = FALSE, ...)
    }
  } else if (!is.null(x$results)) {
    # Single-panel case
    cat("Overall agreement:\n")
    print(x$results, row.names = FALSE, ...)
    cat("\nNumber of subjects:", x$n_subjects,
        "\nNumber of judges:  ", x$n_judges, "\n")
  } else {
    NextMethod()
  }

  invisible(x)
}
