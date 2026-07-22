
    # ══════════════════════════════════════════════════════════════════════════
    # snpPGS plot builders (ggplot2)
    #
    # Every plot is a ggplot object printed at the end of the render function
    # (return TRUE on success, FALSE to leave the image blank). ggplot2 is
    # supplied by the jamovi runtime and is only ever referenced namespaced
    # (ggplot2::…), exactly like snpStats' .render_ld_plot — it is intentionally
    # NOT a declared dependency (see CLAUDE.md).
    #
    # Panels are laid out with facets (mode = Weighted/Unweighted; ROC/calibration
    # polytomous also facet by comparison), replacing the old par(mfrow=…) juggling.
    # ROC/calibration use coord_equal() instead of base asp = 1, which letterboxes
    # on a small device instead of throwing `invalid value for graphical parameter
    # "pin"`. Image size (.pgs_plot_size) is set from the plotWidth/plotHeight
    # options in .init (before the first render, so it is not upscaled and blurry)
    # and re-asserted here for the data-dependent facet count.
    # ══════════════════════════════════════════════════════════════════════════

    # Shared palette (groups / curves / comparisons) and theme.
    .pgs_pal <- c("#2980B9", "#C0392B", "#27AE60", "#8E44AD", "#E67E22", "#16A085")

    # Total image size for a facet grid of nc columns × nr rows. plotWidth/Height
    # are the total for up to two side-by-side panels; a wider grid (polytomous
    # ROC/calibration) scales by half-widths so each panel keeps its size. The
    # SAME formula is used in .init (predicted nc/nr) and in each render function,
    # so the image is sized before the first render — no blurry upscale-then-resize.
    .pgs_plot_size <- function(opts, nc = 1L, nr = 1L) {
      pw <- opts$plotWidth  %||% 600
      ph <- opts$plotHeight %||% 400
      w  <- if (nc <= 2L) pw else pw * nc / 2
      c(w = round(w), h = round(ph * max(1L, nr)))
    }

    .pgs_theme <- function(base_size = 12) {
      ggplot2::theme_minimal(base_size = base_size) +
        ggplot2::theme(
          plot.title       = ggplot2::element_text(face = "bold",
                                                    size = base_size + 1),
          plot.subtitle    = ggplot2::element_text(colour = "#555555",
                                                    size = base_size - 2),
          strip.text       = ggplot2::element_text(face = "bold"),
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_line(colour = "#ECECEC"),
          legend.position  = "bottom",
          legend.title     = ggplot2::element_blank(),
          plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
          panel.background = ggplot2::element_rect(fill = "white", colour = NA))
    }

    # ════════════════════════════════════════════════════════════════════════
    # plotDist — score distribution, one facet per scoring mode.
    #   No / continuous response: overall distribution.
    #   Binary response: one density/histogram per group, overlaid, with group
    #   mean lines and a Mann-Whitney annotation (two groups).
    # ════════════════════════════════════════════════════════════════════════
    plotDist <- function(image, cache, opts) {

      all_scores <- cache$all_scores
      resp       <- cache$resp
      respCol    <- cache$respCol
      if (is.null(all_scores) || length(all_scores) == 0) return(FALSE)

      plot_type <- opts$distPlotType                 # density | histogram | both
      nb        <- opts$histBreaks
      bins      <- if (!is.null(nb) && !is.na(nb) && nb >= 2) as.integer(nb) else 30L

      is_binary <- !is.null(resp) && (is.factor(resp) ||
                    length(unique(resp[!is.na(resp)])) == 2)
      # Group order follows the response's own level order (data order of
      # appearance / user-defined for a factor); without this the character
      # `group` column makes ggplot order the legend/colours alphabetically.
      grp_levels <- if (is_binary)
        levels(droplevels(factor(resp[!is.na(resp)]))) else "Overall"

      seg <- list(); mean_rows <- list(); ann <- list()
      for (ml in names(all_scores)) {
        scores <- all_scores[[ml]]
        valid  <- !is.na(scores)
        if (!is.null(resp)) valid <- valid & !is.na(resp)
        sc <- scores[valid]
        rs <- if (!is.null(resp)) resp[valid] else NULL
        if (length(sc) < 2) next

        grp <- if (is_binary) as.character(rs) else rep("Overall", length(sc))
        seg[[length(seg) + 1L]] <-
          data.frame(mode = ml, group = grp, score = sc, stringsAsFactors = FALSE)
        for (g in unique(grp))
          mean_rows[[length(mean_rows) + 1L]] <-
            data.frame(mode = ml, group = g, mean = mean(sc[grp == g]),
                       stringsAsFactors = FALSE)

        if (is_binary) {
          lvls <- levels(factor(rs))
          if (length(lvls) == 2) {
            mw <- tryCatch(wilcox.test(sc[grp == lvls[2]], sc[grp == lvls[1]],
                                       exact = FALSE), error = function(e) NULL)
            if (!is.null(mw)) {
              p_fmt <- if (mw$p.value < 0.001) "p < 0.001"
                       else paste0("p = ", round(mw$p.value, 3))
              ann[[length(ann) + 1L]] <- data.frame(
                mode = ml, label = paste0("Mann-Whitney  ", p_fmt),
                stringsAsFactors = FALSE)
            }
          }
        }
      }

      df <- do.call(rbind, seg)
      if (is.null(df) || nrow(df) == 0) return(FALSE)
      mean_df <- do.call(rbind, mean_rows)
      ann_df  <- if (length(ann)) do.call(rbind, ann) else NULL
      # Honor the response level order in the legend/colour mapping.
      df$group      <- factor(df$group,      levels = grp_levels)
      mean_df$group <- factor(mean_df$group, levels = grp_levels)
      n_modes <- length(unique(df$mode))
      { sz <- .pgs_plot_size(opts, n_modes, 1L); image$setSize(sz[["w"]], sz[["h"]]) }

      lab_fun <- if (is_binary && !is.null(respCol) && nchar(respCol) > 0)
                   function(g) paste0(respCol, "=", g) else function(g) g

      p <- ggplot2::ggplot(df, ggplot2::aes(x = score, colour = group, fill = group))
      if (plot_type %in% c("histogram", "both"))
        p <- p + ggplot2::geom_histogram(
          ggplot2::aes(y = ggplot2::after_stat(density)),
          bins = bins, position = "identity", alpha = 0.30, colour = NA)
      if (plot_type %in% c("density", "both"))
        p <- p + ggplot2::geom_density(alpha = 0.18, linewidth = 0.9)
      p <- p +
        ggplot2::geom_vline(data = mean_df,
          ggplot2::aes(xintercept = mean, colour = group),
          linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = .pgs_pal, labels = lab_fun) +
        ggplot2::scale_fill_manual(values = .pgs_pal, labels = lab_fun) +
        ggplot2::labs(x = "PGS Score", y = "Density") +
        .pgs_theme()
      if (n_modes > 1) p <- p + ggplot2::facet_wrap(~ mode, scales = "free")
      if (!is.null(ann_df))
        p <- p + ggplot2::geom_text(data = ann_df, inherit.aes = FALSE,
          ggplot2::aes(x = Inf, y = Inf, label = label),
          hjust = 1.05, vjust = 1.5, size = 3, colour = "#555555")
      if (!is_binary)
        p <- p + ggplot2::theme(legend.position = "none")

      print(p)
      TRUE
    }

    # ════════════════════════════════════════════════════════════════════════
    # plotStrat — PGS vs continuous response scatter with a linear fit.
    # (Binary response is shown by plotDist.)
    # ════════════════════════════════════════════════════════════════════════
    plotStrat <- function(image, cache, opts) {

      all_scores <- cache$all_scores
      resp       <- cache$resp
      respCol    <- cache$respCol
      if (is.null(all_scores) || is.null(resp)) return(FALSE)
      is_binary <- is.factor(resp) || length(unique(resp[!is.na(resp)])) == 2
      if (is_binary) return(FALSE)

      pts <- list(); line_rows <- list(); ann <- list()
      for (ml in names(all_scores)) {
        d <- data.frame(mode = ml, pgs = all_scores[[ml]],
                        resp = as.numeric(resp), stringsAsFactors = FALSE)
        d <- d[complete.cases(d), ]
        if (nrow(d) < 3) next
        pts[[length(pts) + 1L]] <- d

        fit <- tryCatch(lm(resp ~ pgs, data = d), error = function(e) NULL)
        if (!is.null(fit)) {
          xs  <- seq(min(d$pgs), max(d$pgs), length.out = 100)
          line_rows[[length(line_rows) + 1L]] <- data.frame(
            mode = ml, pgs = xs,
            resp = as.numeric(predict(fit, data.frame(pgs = xs))),
            stringsAsFactors = FALSE)
          r2 <- summary(fit)$r.squared
          pv <- coef(summary(fit))[2, 4]
          p_fmt <- if (pv < 0.001) "p < 0.001" else paste0("p = ", round(pv, 3))
          ann[[length(ann) + 1L]] <- data.frame(mode = ml,
            label = paste0("R² = ", round(r2, 3), "   ", p_fmt),
            stringsAsFactors = FALSE)
        }
      }

      df <- do.call(rbind, pts)
      if (is.null(df) || nrow(df) == 0) return(FALSE)
      line_df <- if (length(line_rows)) do.call(rbind, line_rows) else NULL
      ann_df  <- if (length(ann)) do.call(rbind, ann) else NULL
      n_modes <- length(unique(df$mode))
      { sz <- .pgs_plot_size(opts, n_modes, 1L); image$setSize(sz[["w"]], sz[["h"]]) }

      y_lab <- if (!is.null(respCol) && nchar(respCol) > 0) respCol else "Response"
      p <- ggplot2::ggplot(df, ggplot2::aes(pgs, resp)) +
        ggplot2::geom_point(colour = .pgs_pal[1], alpha = 0.40, size = 1.1)
      if (!is.null(line_df))
        p <- p + ggplot2::geom_line(data = line_df, colour = .pgs_pal[2],
                                    linewidth = 1)
      p <- p +
        ggplot2::labs(x = "PGS Score", y = y_lab) +
        .pgs_theme() + ggplot2::theme(legend.position = "none")
      if (n_modes > 1) p <- p + ggplot2::facet_wrap(~ mode, scales = "free")
      if (!is.null(ann_df))
        p <- p + ggplot2::geom_text(data = ann_df, inherit.aes = FALSE,
          ggplot2::aes(x = -Inf, y = Inf, label = label),
          hjust = -0.05, vjust = 1.5, size = 3, colour = "#555555")

      print(p)
      TRUE
    }

    # ════════════════════════════════════════════════════════════════════════
    # plotForest — OR (binary) / β (continuous) per percentile category, using
    # the same logistic / linear fits as .fillPercentileTable. Facet per mode;
    # reference category drawn as a diamond, others as point + CI whiskers.
    # (Polytomous response is skipped.)
    # ════════════════════════════════════════════════════════════════════════
    plotForest <- function(image, cache, opts) {

      all_scores <- cache$all_scores
      resp       <- cache$resp
      covs       <- cache$covs
      if (is.null(all_scores) || is.null(resp)) return(FALSE)

      n_lvls    <- length(unique(resp[!is.na(resp)]))
      resp_type <- if (!is.factor(resp) && n_lvls > 5) "continuous"
                   else if (n_lvls == 2)               "binary"
                   else if (n_lvls > 2)                "polytomous"
                   else                                return(FALSE)
      resp_levels <- if (resp_type == "polytomous")
        levels(droplevels(factor(resp[!is.na(resp)]))) else character(0)

      breaks_num <- suppressWarnings(as.numeric(
        strsplit(trimws(opts$percentileBreaks), ",")[[1]]))
      breaks_num <- sort(unique(breaks_num[!is.na(breaks_num) &
                                           breaks_num > 0 & breaks_num < 100]))
      if (length(breaks_num) == 0) breaks_num <- c(20, 40, 60, 80, 90, 95)

      ref_opt  <- opts$pgsRefCategory
      has_covs <- !is.null(covs) && ncol(covs) > 0

      make_labels <- function(brks) {
        n <- length(brks); lbl <- character(n + 1)
        lbl[1] <- paste0("<P", brks[1])
        for (i in seq_len(n - 1))
          lbl[i + 1] <- paste0("P", brks[i], "–P", brks[i + 1])
        lbl[n + 1] <- paste0(">P", brks[n]); lbl
      }
      cat_labels <- make_labels(breaks_num)
      n_cats     <- length(cat_labels)

      null_v <- if (resp_type == "continuous") 0 else 1
      x_lab  <- if (resp_type == "continuous") "β vs reference (95% CI)"
                else "Odds Ratio (95% CI)"

      rows <- list()
      for (ml in names(all_scores)) {
        scores  <- all_scores[[ml]]
        cuts    <- quantile(scores, breaks_num / 100, na.rm = TRUE)
        cat_idx <- pmin(pmax(findInterval(scores, cuts) + 1L, 1L), n_cats)
        ref_idx <- switch(ref_opt, lowest = 1L, highest = n_cats,
          middle = as.integer(pmin(pmax(
            findInterval(median(scores, na.rm = TRUE), cuts) + 1L, 1L), n_cats)))

        df <- if (has_covs)
          data.frame(cat = cat_idx, resp = resp, covs, check.names = FALSE)
        else data.frame(cat = cat_idx, resp = resp)
        df <- df[complete.cases(df), ]
        if (nrow(df) < n_cats) next
        df$cat <- relevel(factor(df$cat, levels = seq_len(n_cats),
                                 labels = cat_labels), ref = cat_labels[ref_idx])
        cov_terms <- if (has_covs) safe_rhs(names(covs)) else ""

        # One category-indexed est/lo/hi block per contrast (reference category
        # holds the null value). Emitted with a `comparison` label so polytomous
        # responses facet by contrast (mode × comparison) like the ROC plot.
        add_block <- function(comparison, est, lo, hi) {
          rows[[length(rows) + 1L]] <<- data.frame(
            mode = ml, comparison = comparison,
            category = factor(cat_labels, levels = cat_labels),
            est = est, lo = lo, hi = hi, is_ref = seq_len(n_cats) == ref_idx,
            stringsAsFactors = FALSE)
        }
        # Pull a "cat<label>"-named coefficient vector into category order.
        pull <- function(cf, cis, transform, null) {
          est <- rep(NA_real_, n_cats); lo <- est; hi <- est
          est[ref_idx] <- null; lo[ref_idx] <- null; hi[ref_idx] <- null
          for (ci in seq_len(n_cats)) {
            if (ci == ref_idx) next
            nm <- paste0("cat", cat_labels[ci])
            if (is.null(cf) || !nm %in% rownames(cf)) next
            b <- cf[nm, 1]
            v <- if (!is.null(cis) && nm %in% rownames(cis)) cis[nm, ]
                 else c(b - 1.96 * cf[nm, 2], b + 1.96 * cf[nm, 2])
            est[ci] <- transform(b); lo[ci] <- transform(v[1]); hi[ci] <- transform(v[2])
          }
          list(est = est, lo = lo, hi = hi)
        }

        if (resp_type == "binary") {
          df$resp <- factor(df$resp)
          frm <- if (has_covs) as.formula(paste("resp ~ cat +", cov_terms)) else resp ~ cat
          fit <- tryCatch(glm(frm, data = df, family = binomial()),
                          error = function(e) NULL)
          if (is.null(fit)) next
          r <- pull(coef(summary(fit)),
                    tryCatch(confint.default(fit), error = function(e) NULL), exp, 1)
          add_block("", r$est, r$lo, r$hi)

        } else if (resp_type == "continuous") {
          df$resp <- as.numeric(df$resp)
          frm <- if (has_covs) as.formula(paste("resp ~ cat +", cov_terms)) else resp ~ cat
          fit <- tryCatch(lm(frm, data = df), error = function(e) NULL)
          if (is.null(fit)) next
          r <- pull(coef(summary(fit)),
                    tryCatch(confint(fit), error = function(e) NULL), identity, 0)
          add_block("", r$est, r$lo, r$hi)

        } else {  # polytomous — one OR block per non-reference response level
          df$resp <- factor(df$resp)
          frm <- if (has_covs) as.formula(paste("resp ~ cat +", cov_terms)) else resp ~ cat
          fit <- tryCatch(nnet::multinom(frm, data = df, trace = FALSE),
                          error = function(e) NULL)
          if (is.null(fit)) next
          cf_mat <- coef(fit)                                  # (K-1) × p
          vc     <- tryCatch(vcov(fit), error = function(e) NULL)
          for (lv in resp_levels[-1]) {
            est <- rep(NA_real_, n_cats); lo <- est; hi <- est
            est[ref_idx] <- 1; lo[ref_idx] <- 1; hi[ref_idx] <- 1
            for (ci in seq_len(n_cats)) {
              if (ci == ref_idx) next
              nm <- paste0("cat", cat_labels[ci])
              if (!lv %in% rownames(cf_mat) || !nm %in% colnames(cf_mat)) next
              b  <- cf_mat[lv, nm]
              vn <- paste0(lv, ":", nm)
              vd <- if (!is.null(vc) && vn %in% rownames(vc)) vc[vn, vn] else NA_real_
              se <- if (!is.na(vd) && vd > 0) sqrt(vd) else NA_real_
              est[ci] <- exp(b)
              if (!is.na(se)) { lo[ci] <- exp(b - 1.96 * se); hi[ci] <- exp(b + 1.96 * se) }
            }
            add_block(paste0(lv, " vs ", resp_levels[1]), est, lo, hi)
          }
        }
      }

      df <- do.call(rbind, rows)
      if (is.null(df) || nrow(df) == 0) return(FALSE)
      n_modes <- length(unique(df$mode))
      n_comps <- length(unique(df$comparison))
      if (n_comps == 1L) { nc <- n_modes; nr <- 1L }
      else               { nc <- n_comps; nr <- n_modes }
      { sz <- .pgs_plot_size(opts, nc, nr); image$setSize(sz[["w"]], sz[["h"]]) }

      pt   <- df[!df$is_ref, ]
      refp <- df[df$is_ref, ]
      sub  <- if (resp_type == "polytomous")
        paste0("Each level vs ", resp_levels[1],
               if (has_covs) paste0("; adjusted for ", paste(names(covs), collapse = ", ")))
      else if (has_covs)
        paste0("Adjusted for: ", paste(names(covs), collapse = ", ")) else NULL

      p <- ggplot2::ggplot(df, ggplot2::aes(x = est, y = category)) +
        ggplot2::geom_vline(xintercept = null_v, linetype = "dashed",
                            colour = "#888888")
      if (nrow(pt) > 0)
        p <- p +
          ggplot2::geom_errorbarh(data = pt,
            ggplot2::aes(xmin = lo, xmax = hi), height = 0.22,
            colour = .pgs_pal[1], linewidth = 0.8, na.rm = TRUE) +
          ggplot2::geom_point(data = pt, size = 2.4, colour = "#2C3E50",
                              na.rm = TRUE)
      if (nrow(refp) > 0)
        p <- p + ggplot2::geom_point(data = refp, shape = 23, size = 3.4,
          fill = .pgs_pal[2], colour = .pgs_pal[2], na.rm = TRUE)
      p <- p +
        # Pin the y order to the percentile sequence (lowest at top, matching the
        # percentile table). Every geom uses a subset (pt/refp) as its own data, so
        # ggplot would otherwise train the discrete scale from observed values in
        # layer order — appending the reference category (absent from pt) last and
        # breaking the order.
        ggplot2::scale_y_discrete(limits = rev(cat_labels)) +
        ggplot2::labs(subtitle = sub, x = x_lab, y = NULL) +
        .pgs_theme() + ggplot2::theme(legend.position = "none")
      if (n_comps == 1L && n_modes > 1) p <- p + ggplot2::facet_wrap(~ mode)
      else if (n_comps > 1)             p <- p + ggplot2::facet_grid(mode ~ comparison)

      print(p)
      TRUE
    }

    # ════════════════════════════════════════════════════════════════════════
    # plotROC — ROC curve(s). Binary: one facet per mode. Polytomous: each
    # non-reference level vs the reference, faceted mode × comparison. When
    # covariates are present, PGS+covariates and covariates-only curves are
    # overlaid. AUC via the trapezoidal rule; a note flags AUC < 0.5.
    # ════════════════════════════════════════════════════════════════════════
    plotROC <- function(image, cache, opts) {

      all_scores <- cache$all_scores
      resp       <- cache$resp
      covs       <- cache$covs
      if (is.null(all_scores) || is.null(resp)) return(FALSE)

      lvls_all <- levels(droplevels(factor(resp[!is.na(resp)])))
      if (length(lvls_all) < 2L) return(FALSE)
      has_covs  <- !is.null(covs) && ncol(covs) > 0
      ref_lv    <- lvls_all[1]
      comp_lvls <- lvls_all[-1]
      n_comps   <- length(comp_lvls)

      roc_curve <- function(score, label) {
        thr <- sort(unique(score), decreasing = TRUE)
        tpr <- vapply(thr, function(t) mean(score[label == 1] >= t), numeric(1))
        fpr <- vapply(thr, function(t) mean(score[label == 0] >= t), numeric(1))
        tpr <- c(0, tpr, 1); fpr <- c(0, fpr, 1)
        list(fpr = fpr, tpr = tpr,
             auc = sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2))
      }

      curves <- list(); ann <- list()
      for (ml in names(all_scores)) {
        scores <- all_scores[[ml]]
        for (comp_lv in comp_lvls) {
          comp_lab <- if (n_comps == 1L) "" else paste0(comp_lv, " vs ", ref_lv)
          mask <- !is.na(resp) & (as.character(resp) == ref_lv |
                                  as.character(resp) == comp_lv)
          if (sum(mask) < 10L) next
          df <- data.frame(pgs = scores[mask], resp = resp[mask])
          if (has_covs) df <- cbind(df, covs[mask, , drop = FALSE])
          df <- df[complete.cases(df), ]
          if (nrow(df) < 10L) next
          label01 <- as.integer(as.character(df$resp) == comp_lv)
          df$resp_fac <- factor(as.character(df$resp), levels = c(ref_lv, comp_lv))

          add_curve <- function(name, fpr, tpr, auc) {
            curves[[length(curves) + 1L]] <<- data.frame(
              mode = ml, comparison = comp_lab, curve = name,
              fpr = fpr, tpr = tpr, stringsAsFactors = FALSE)
            ann[[length(ann) + 1L]] <<- data.frame(
              mode = ml, comparison = comp_lab, curve = name, auc = auc,
              stringsAsFactors = FALSE)
          }

          # PGS curve from the univariate logistic P(event | PGS), not the raw
          # score. The fitted probability orients discrimination to the event, so
          # the curve is invariant to which level is the reference — matching the
          # covariate/adjusted curves (model-based, already reference-invariant).
          # For a positively-associated PGS this is the identical empirical ROC as
          # the raw score; only the sign of the reference↔event choice differs.
          fit_pgs <- tryCatch(glm(resp_fac ~ pgs, data = df, family = binomial()),
                              error = function(e) NULL)
          if (!is.null(fit_pgs)) {
            r1 <- tryCatch(roc_curve(predict(fit_pgs, type = "response"), label01),
                           error = function(e) NULL)
            if (!is.null(r1)) add_curve("PGS", r1$fpr, r1$tpr, r1$auc)
          }

          if (has_covs) {
            cov_terms   <- safe_rhs(names(covs))
            fit_adj <- tryCatch(glm(
              as.formula(paste("resp_fac ~ pgs +", cov_terms)),
              data = df, family = binomial()), error = function(e) NULL)
            if (!is.null(fit_adj)) {
              r2 <- tryCatch(roc_curve(predict(fit_adj, type = "response"), label01),
                             error = function(e) NULL)
              if (!is.null(r2)) add_curve("PGS + covariates", r2$fpr, r2$tpr, r2$auc)
            }
            fit_cov <- tryCatch(glm(
              as.formula(paste("resp_fac ~", cov_terms)),
              data = df, family = binomial()), error = function(e) NULL)
            if (!is.null(fit_cov)) {
              r3 <- tryCatch(roc_curve(predict(fit_cov, type = "response"), label01),
                             error = function(e) NULL)
              if (!is.null(r3)) add_curve("Covariates only", r3$fpr, r3$tpr, r3$auc)
            }
          }
        }
      }

      df <- do.call(rbind, curves)
      if (is.null(df) || nrow(df) == 0) return(FALSE)
      ann_df <- do.call(rbind, ann)
      # One AUC label block per facet, curves listed on separate lines.
      ann_df <- ann_df[order(ann_df$mode, ann_df$comparison, ann_df$curve), ]
      lab_df <- do.call(rbind, lapply(
        split(ann_df, list(ann_df$mode, ann_df$comparison), drop = TRUE),
        function(g) data.frame(mode = g$mode[1], comparison = g$comparison[1],
          label = paste(sprintf("%s  AUC = %.3f", g$curve, g$auc), collapse = "\n"),
          stringsAsFactors = FALSE)))

      n_modes <- length(unique(df$mode))
      if (n_comps == 1L) { nc <- n_modes; nr <- 1L }
      else               { nc <- n_comps; nr <- n_modes }
      { sz <- .pgs_plot_size(opts, nc, nr); image$setSize(sz[["w"]], sz[["h"]]) }

      curve_lvls <- c("PGS", "PGS + covariates", "Covariates only")
      df$curve <- factor(df$curve, levels = curve_lvls)
      curve_cols <- stats::setNames(.pgs_pal[1:3], curve_lvls)
      curve_ltys <- c("PGS" = "solid", "PGS + covariates" = "22",
                      "Covariates only" = "42")
      p <- ggplot2::ggplot(df, ggplot2::aes(fpr, tpr, colour = curve,
                                            linetype = curve)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = "#AAAAAA") +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::scale_colour_manual(values = curve_cols) +
        ggplot2::scale_linetype_manual(values = curve_ltys) +
        ggplot2::labs(
          subtitle = if (n_comps > 1) paste0("Each category vs ", ref_lv) else NULL,
          x = "1 − Specificity (FPR)", y = "Sensitivity (TPR)") +
        ggplot2::geom_text(data = lab_df, inherit.aes = FALSE,
          ggplot2::aes(x = 1, y = 0, label = label),
          hjust = 1, vjust = 0, size = 3, colour = "#333333") +
        .pgs_theme()
      if (n_comps == 1L && n_modes > 1) p <- p + ggplot2::facet_wrap(~ mode)
      else if (n_comps > 1) p <- p + ggplot2::facet_grid(mode ~ comparison)

      print(p)
      TRUE
    }

    # ════════════════════════════════════════════════════════════════════════
    # plotCalibration — observed vs predicted event rate by decile of predicted
    # probability, with a loess smooth and a Hosmer-Lemeshow annotation. Same
    # faceting as plotROC.
    # ════════════════════════════════════════════════════════════════════════
    plotCalibration <- function(image, cache, opts) {

      all_scores <- cache$all_scores
      resp       <- cache$resp
      covs       <- cache$covs
      if (is.null(all_scores) || is.null(resp)) return(FALSE)

      lvls_all <- levels(droplevels(factor(resp[!is.na(resp)])))
      if (length(lvls_all) < 2L) return(FALSE)
      has_covs  <- !is.null(covs) && ncol(covs) > 0
      ref_lv    <- lvls_all[1]
      comp_lvls <- lvls_all[-1]
      n_comps   <- length(comp_lvls)

      pts <- list(); smooth <- list(); ann <- list()
      for (ml in names(all_scores)) {
        scores <- all_scores[[ml]]
        for (comp_lv in comp_lvls) {
          comp_lab <- if (n_comps == 1L) "" else paste0(comp_lv, " vs ", ref_lv)
          mask <- !is.na(resp) & (as.character(resp) == ref_lv |
                                  as.character(resp) == comp_lv)
          if (sum(mask) < 20L) next
          df <- data.frame(pgs = scores[mask], resp = resp[mask])
          if (has_covs) df <- cbind(df, covs[mask, , drop = FALSE])
          df <- df[complete.cases(df), ]
          if (nrow(df) < 20L) next
          label01   <- as.integer(as.character(df$resp) == comp_lv)
          df$resp_fac <- factor(as.character(df$resp), levels = c(ref_lv, comp_lv))
          cov_terms <- if (has_covs) safe_rhs(names(covs)) else ""
          frm <- if (has_covs) as.formula(paste("resp_fac ~ pgs +", cov_terms))
                 else resp_fac ~ pgs
          fit <- tryCatch(glm(frm, data = df, family = binomial()),
                          error = function(e) NULL)
          if (is.null(fit)) next

          pred    <- predict(fit, type = "response")
          n_bins  <- max(min(10L, floor(nrow(df) / 15L)), 3L)
          qb <- unique(quantile(pred, seq(0, 1, 1 / n_bins), na.rm = TRUE))
          if (length(qb) < 2L) next
          bin_idx   <- cut(pred, breaks = qb, include.lowest = TRUE, labels = FALSE)
          obs_rate  <- tapply(label01, bin_idx, mean,   na.rm = TRUE)
          mean_pred <- tapply(pred,    bin_idx, mean,   na.rm = TRUE)
          bin_n     <- tapply(label01, bin_idx, length)

          pts[[length(pts) + 1L]] <- data.frame(
            mode = ml, comparison = comp_lab,
            mean_pred = as.numeric(mean_pred), obs_rate = as.numeric(obs_rate),
            n = as.numeric(bin_n), stringsAsFactors = FALSE)

          if (length(mean_pred) >= 5L) {
            lo_fit <- tryCatch(loess(obs_rate ~ mean_pred,
              weights = as.numeric(bin_n), span = 0.9), error = function(e) NULL)
            if (!is.null(lo_fit)) {
              xs <- seq(min(mean_pred, na.rm = TRUE),
                        max(mean_pred, na.rm = TRUE), length.out = 100)
              ys <- tryCatch(as.numeric(predict(lo_fit, xs)),
                             error = function(e) NULL)
              if (!is.null(ys))
                smooth[[length(smooth) + 1L]] <- data.frame(
                  mode = ml, comparison = comp_lab, mean_pred = xs, obs_rate = ys,
                  stringsAsFactors = FALSE)
            }
          }

          exp_ev  <- tapply(pred,    bin_idx, sum, na.rm = TRUE)
          obs_ev  <- tapply(label01, bin_idx, sum, na.rm = TRUE)
          exp_ne  <- as.numeric(bin_n) - exp_ev
          obs_ne  <- as.numeric(bin_n) - obs_ev
          ok      <- exp_ev > 0 & exp_ne > 0
          if (sum(ok) >= 2L) {
            hl <- sum((obs_ev[ok] - exp_ev[ok])^2 / exp_ev[ok] +
                      (obs_ne[ok] - exp_ne[ok])^2 / exp_ne[ok])
            hl_df <- sum(ok) - 2L
            hl_p  <- pchisq(hl, df = hl_df, lower.tail = FALSE)
            p_fmt <- if (hl_p < 0.001) "p < 0.001" else paste0("p = ", round(hl_p, 3))
            ann[[length(ann) + 1L]] <- data.frame(mode = ml, comparison = comp_lab,
              label = sprintf("Hosmer-Lemeshow χ²(%d) = %.2f  %s",
                              hl_df, hl, p_fmt), stringsAsFactors = FALSE)
          }
        }
      }

      df <- do.call(rbind, pts)
      if (is.null(df) || nrow(df) == 0) return(FALSE)
      smooth_df <- if (length(smooth)) do.call(rbind, smooth) else NULL
      ann_df    <- if (length(ann)) do.call(rbind, ann) else NULL

      n_modes <- length(unique(df$mode))
      if (n_comps == 1L) { nc <- n_modes; nr <- 1L }
      else               { nc <- n_comps; nr <- n_modes }
      { sz <- .pgs_plot_size(opts, nc, nr); image$setSize(sz[["w"]], sz[["h"]]) }

      p <- ggplot2::ggplot(df, ggplot2::aes(mean_pred, obs_rate)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = "#AAAAAA")
      if (!is.null(smooth_df))
        p <- p + ggplot2::geom_line(data = smooth_df, colour = .pgs_pal[1],
                                    linewidth = 0.9)
      p <- p +
        ggplot2::geom_point(ggplot2::aes(size = n), shape = 21,
          fill = "#2C3E50", colour = "#2C3E50", alpha = 0.65, na.rm = TRUE) +
        ggplot2::scale_size_area(max_size = 5, guide = "none") +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(
          subtitle = if (n_comps > 1) paste0("Each category vs ", ref_lv) else NULL,
          x = "Mean predicted probability", y = "Observed event rate") +
        .pgs_theme() + ggplot2::theme(legend.position = "none")
      if (!is.null(ann_df))
        p <- p + ggplot2::geom_text(data = ann_df, inherit.aes = FALSE,
          ggplot2::aes(x = 0, y = 1, label = label),
          hjust = 0, vjust = 1, size = 3, colour = "#555555")
      if (n_comps == 1L && n_modes > 1) p <- p + ggplot2::facet_wrap(~ mode)
      else if (n_comps > 1) p <- p + ggplot2::facet_grid(mode ~ comparison)

      print(p)
      TRUE
    }
