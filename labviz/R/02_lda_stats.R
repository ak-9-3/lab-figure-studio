# ── LDA statistics ────────────────────────────────────────────────────────
# Both common assay layouts reduce to ONE model: a binomial GLM with a
# complementary log-log link.
#
#   "cells"  every seeded cell is a trial, organoids are successes.
#            cloglog(p) = b            ->  frequency f = exp(b)   (per cell)
#
#   "wells"  every well is a trial, scored +/-, seeded with `Dose` cells.
#            p(+) = 1 - exp(-f * dose)
#            cloglog(p) = log(dose) + b  ->  same f = exp(b)
#            i.e. classic ELDA (Hu & Smyth 2009) written as an offset.
#
# So the fitting, the CIs and the group comparisons are shared code, and
# only the figure differs between the two layouts.

lda_fit <- function(df, groups, mode = c("cells", "wells"), conf = 0.95) {
  mode <- match.arg(mode)
  df <- df[df$Group %in% groups, , drop = FALSE]
  df$Group <- factor(df$Group, levels = groups)

  if (mode == "cells") {
    trials <- df$CellCount; events <- df$Organoids; off <- rep(0, nrow(df))
  } else {
    trials <- df$Tested;    events <- df$Positive;  off <- log(df$Dose)
  }
  bad <- which(is.finite(events) & is.finite(trials) & events > trials)
  if (length(bad)) {
    i <- bad[1]
    lab <- if (mode == "cells") c("Organoids formed", "Cells seeded")
           else                 c("Wells positive", "Wells tested")
    stop(sprintf(
      paste0("'%s' is bigger than '%s' in %d of %d rows.\n",
             "Row %d says %s grew from %s.\n",
             "Those two columns are almost certainly mapped the wrong way round."),
      lab[1], lab[2], length(bad), length(events),
      i, format(events[i]), format(trials[i])), call. = FALSE)
  }

  md  <- data.frame(events = events, fails = trials - events,
                    Group = droplevels(df$Group), off = off)
  # A single group has no contrasts to estimate — fit the intercept alone.
  form <- if (nlevels(md$Group) < 2) cbind(events, fails) ~ 1
          else                       cbind(events, fails) ~ Group
  fit <- stats::glm(stats::update(form, . ~ . + offset(off)),
                    family = stats::binomial(link = "cloglog"), data = md)

  b <- stats::coef(fit); V <- stats::vcov(fit); z <- stats::qnorm(1 - (1 - conf) / 2)

  est <- do.call(rbind, lapply(seq_along(groups), function(i) {
    L <- numeric(length(b)); L[1] <- 1; if (i > 1) L[i] <- 1
    e  <- sum(L * b)
    se <- sqrt(drop(t(L) %*% V %*% L))
    data.frame(Group = groups[i], n_rows = sum(df$Group == groups[i]),
               Frequency = exp(e), CI_low = exp(e - z * se), CI_high = exp(e + z * se),
               SE_log = se, stringsAsFactors = FALSE)
  }))
  est$`1 in` <- 1 / est$Frequency
  est$Group  <- factor(est$Group, levels = groups)

  # Goodness of fit: does a single frequency per group explain all doses?
  gof <- tryCatch(
    stats::pchisq(fit$deviance, fit$df.residual, lower.tail = FALSE),
    error = function(e) NA_real_)

  list(fit = fit, est = est, groups = groups, mode = mode,
       gof_p = gof, deviance = fit$deviance, df_resid = fit$df.residual)
}

# Wald test on any pair of groups, on the log-frequency scale.
lda_pairwise <- function(lda, pairs, adjust = "none") {
  if (!length(pairs)) return(NULL)
  b <- stats::coef(lda$fit); V <- stats::vcov(lda$fit); g <- lda$groups
  out <- do.call(rbind, lapply(pairs, function(pr) {
    i <- match(pr[1], g); j <- match(pr[2], g)
    # contrast is B - A, so Ratio reads as "B relative to A"
    L <- numeric(length(b)); if (j > 1) L[j] <- L[j] + 1; if (i > 1) L[i] <- L[i] - 1
    d  <- sum(L * b)                       # log fold-change in frequency
    se <- sqrt(drop(t(L) %*% V %*% L))
    zv <- d / se
    data.frame(A = pr[1], B = pr[2],
               Ratio = exp(d), CI_low = exp(d - 1.96 * se), CI_high = exp(d + 1.96 * se),
               z = zv, p = 2 * stats::pnorm(-abs(zv)), stringsAsFactors = FALSE)
  }))
  out$p_adj <- stats::p.adjust(out$p, method = adjust)
  out$stars <- p_stars(out$p_adj)
  out
}

# Observed per-row rate, used for the dots on the bar chart.
lda_observed <- function(df, mode = c("cells", "wells")) {
  mode <- match.arg(mode)
  if (mode == "cells") {
    data.frame(Group = df$Group, rate = df$Organoids / df$CellCount,
               stringsAsFactors = FALSE)
  } else {
    # per-dose MLE from the negative fraction: f = -log(1 - p+) / dose
    p <- pmin(df$Positive / df$Tested, 1 - 1e-9)
    data.frame(Group = df$Group, rate = -log(1 - p) / df$Dose,
               stringsAsFactors = FALSE)
  }
}

# Which group pairs to compare, shared by the app and the exported script.
# When the groups form panels (see group_panels), "vs ref" means "vs the
# reference treatment inside each panel" — N16-ATRA against N16-DMSO, never
# against Ad16-DMSO. "all" adds every pair inside a panel plus the same
# treatment across panels.
lda_pairs <- function(groups, ref, mode = c("none", "ref", "all"),
                      panels = group_panels(groups)) {
  mode <- match.arg(mode)
  if (length(groups) < 2 || is.null(ref) || !ref %in% groups) return(list())
  if (mode == "none") return(list())
  if (is.null(panels))
    return(switch(mode,
      ref = lapply(setdiff(groups, ref), function(x) c(ref, x)),
      all = utils::combn(groups, 2, simplify = FALSE)))

  ref_lab <- as.character(panels$Label[panels$Group == ref])
  out <- list()
  for (pn in levels(panels$Panel)) {
    in_p <- panels$Group[panels$Panel == pn]
    if (mode == "ref") {
      r <- in_p[as.character(panels$Label[panels$Panel == pn]) == ref_lab]
      r <- if (length(r)) r[1] else in_p[1]
      out <- c(out, lapply(setdiff(in_p, r), function(x) c(r, x)))
    } else if (length(in_p) >= 2) {
      out <- c(out, utils::combn(in_p, 2, simplify = FALSE))
    }
  }
  if (mode == "all")
    for (lb in levels(panels$Label)) {
      same <- panels$Group[panels$Label == lb]
      if (length(same) >= 2) out <- c(out, utils::combn(same, 2, simplify = FALSE))
    }
  out
}
