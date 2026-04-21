#' @importFrom R6 R6Class
#' @import jmvcore
#' @importFrom genetics genotype allele LD
#' @importFrom haplo.stats setupGeno hapl.em haplo.glm haplo.glm.control
#' @import ggplot2
source("R/snp_helpers.R")

snpLDHaploClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
  "snpLDHaploClass",
  inherit = snpLDHaploBase,
  private = list(

    # Private storage for the LD heatmap render function
    .ld_store  = NULL,
    .ld_nms    = NULL,
    .ld_metric = NULL,

    .init = function() {
      self$results$ldGroup$setVisible(FALSE)
      self$results$haploGroup$setVisible(FALSE)
    },

    .run = function() {
      data           <- self$data
      opts           <- self$options
      response_var   <- opts$response
      snp_vars       <- opts$snps
      covariate_vars <- opts$covariates

      run_ldAnalysis      <- isTRUE(opts$ldAnalysis)
      run_ldMatrix        <- isTRUE(opts$ldMatrix)
      run_ldPlot          <- isTRUE(opts$ldPlot)
      run_haploFreq       <- isTRUE(opts$haploFreq)
      run_haploAssoc      <- isTRUE(opts$haploAssoc)
      run_haploInteraction <- isTRUE(opts$haploInteraction)
      run_subpop          <- isTRUE(opts$subpop)

      # ── Validate: need SNPs ──────────────────────────────────────
      if (length(snp_vars) == 0) {
        self$results$validationMsg$setContent(
          "<p style='color:red;'>Please add at least one SNP variable.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return()
      }
      val      <- validate_snp_vars(snp_vars, data)
      snp_vars <- val$valid_snps
      if (nchar(val$bad_html) > 0) {
        self$results$validationMsg$setContent(val$bad_html)
        self$results$validationMsg$setVisible(TRUE)
      } else {
        self$results$validationMsg$setVisible(FALSE)
      }

      # ── Validate: need ≥2 SNPs ───────────────────────────────────
      if (length(snp_vars) < 2) {
        self$results$validationMsg$setContent(
          "<p style='color:red;'>LD and haplotype analyses require at least 2 SNPs.</p>")
        self$results$validationMsg$setVisible(TRUE)
        return()
      }

      # ── Prepare response / covariates ────────────────────────────
      response_raw  <- if (!is.null(response_var) && response_var != "")
                         data[[response_var]] else NULL
      response_type <- detect_response_type(response_raw, opts$responseType)
      response      <- prepare_response(response_raw, response_type)
      cov_df        <- prepare_covariates(data, covariate_vars)

      # if (run_subpop && (is.null(response_raw) || response_type == "quantitative"))
      #   run_subpop <- FALSE
      if (run_haploInteraction && (is.null(cov_df) || ncol(cov_df) == 0)) {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>Haplotype \u00D7 covariate interaction requires at least one covariate.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_haploInteraction <- FALSE
      }
      if (run_haploAssoc && is.null(response_raw)) {
        self$results$validationMsg$setContent(
          "<p style='color:orange;'>Haplotype association requires a response variable.</p>")
        self$results$validationMsg$setVisible(TRUE)
        run_haploAssoc <- FALSE
      }

      # ── Visibility ───────────────────────────────────────────────
      self$results$ldGroup$setVisible(run_ldAnalysis || run_ldMatrix || run_ldPlot)
      self$results$haploGroup$setVisible(run_haploFreq || run_haploAssoc || run_haploInteraction)

      # ── Complete-case mask ───────────────────────────────────────
      n_rows        <- nrow(data)
      complete_mask <- rep(TRUE, n_rows)
      if (!is.null(response))              complete_mask <- complete_mask & !is.na(response)
      if (!is.null(cov_df) && ncol(cov_df) > 0)
        complete_mask <- complete_mask & complete.cases(cov_df)

      # ── Parse genotypes ──────────────────────────────────────────
      geno_list <- list()
      for (snp_nm in snp_vars) {
        snp_raw     <- data[[snp_nm]]
        user_levels <- get_snp_level_order(snp_raw)
        geno_obj    <- parse_genotype(snp_raw, user_levels)
        if (!is.null(geno_obj)) geno_list[[snp_nm]] <- geno_obj
      }
      if (length(geno_list) < 2) return()

      # ── LD analysis ──────────────────────────────────────────────
      if (run_ldAnalysis || run_ldMatrix || run_ldPlot)
        private$.run_ld(geno_list, opts, run_ldAnalysis, run_ldMatrix, run_ldPlot)

      # ── Haplotype analysis ───────────────────────────────────────
      if (run_haploFreq || run_haploAssoc || run_haploInteraction)
        private$.run_haplo(geno_list, data, response, response_raw, response_type,
                           cov_df, opts, run_haploFreq, run_haploAssoc,
                           run_haploInteraction, run_subpop, complete_mask)
    },

    # ── LD ───────────────────────────────────────────────────────────────────
    .run_ld = function(geno_list, opts, run_ldAnalysis, run_ldMatrix, run_ldPlot) {
      nms   <- names(geno_list)
      n     <- length(nms)
      pairs <- combn(nms, 2, simplify = FALSE)

      ld_store <- list()
      for (pair in pairs) {
        key    <- paste(pair, collapse = "___")
        ld_res <- tryCatch(genetics::LD(geno_list[[pair[1]]], geno_list[[pair[2]]]),
                           error = function(e) NULL)
        if (!is.null(ld_res)) ld_store[[key]] <- ld_res
      }

      if (run_ldAnalysis) {
        tbl <- self$results$ldGroup$ldTable
        for (pair in pairs) {
          key    <- paste(pair, collapse = "___")
          ld_res <- ld_store[[key]]
          if (is.null(ld_res)) next
          tbl$addRow(rowKey = paste(pair, collapse="_"), values = list(
            snp1   = pair[1], snp2 = pair[2],
            r2     = round(ld_res$`r`^2,  3),
            Dprime = round(ld_res$`D'`,   3),
            D      = round(ld_res$`D`,    3),
            pval   = ld_res$`P-value`))
        }
      }

      if (run_ldMatrix) {
        mtbl   <- self$results$ldGroup$ldMatrixTable
        metric <- opts$ldMetric
        for (nm in nms) {
          safe_nm <- gsub("[^A-Za-z0-9_]","_",nm)
          mtbl$addColumn(name = safe_nm, title = nm, type = "text")
        }
        upper_mat <- matrix("", n, n, dimnames = list(nms, nms))
        lower_mat <- matrix("", n, n, dimnames = list(nms, nms))
        diag(upper_mat) <- nms; diag(lower_mat) <- nms

        for (pair in pairs) {
          key    <- paste(pair, collapse = "___")
          ld_res <- ld_store[[key]]
          if (is.null(ld_res)) next
          p_val  <- ld_res$`P-value`
          p_str  <- if (!is.na(p_val)) { if (p_val < 0.001) "< .001" else sprintf("%.3f",p_val) } else ""
          up_val <- switch(metric,
            Dprime = sprintf("%.3f", round(ld_res$`D'`,  3)),
            r2     = sprintf("%.3f", round(ld_res$`r`^2, 3)),
            D      = sprintf("%.3f", round(ld_res$`D`,   3)))
          upper_mat[pair[1], pair[2]] <- up_val
          lower_mat[pair[2], pair[1]] <- p_str
        }
        for (i in seq_len(n)) {
          row_vals <- list(snp = nms[i])
          for (j in seq_len(n)) {
            safe_nm <- gsub("[^A-Za-z0-9_]","_",nms[j])
            row_vals[[safe_nm]] <- if(i==j) nms[i] else if(j>i) upper_mat[i,j] else lower_mat[i,j]
          }
          mtbl$addRow(rowKey = paste0("row_",i), values = row_vals)
        }
        metric_label <- switch(metric, Dprime="D'", r2="r²", D="D")
        mtbl$setNote(key="layout",
                     note=paste0("Upper triangle: ", metric_label,
                                 ". Lower triangle: P-value. Diagonal: SNP name."))
      }

      if (run_ldPlot) {
        private$.ld_store  <- ld_store
        private$.ld_nms    <- nms
        private$.ld_metric <- opts$ldMetric
        self$results$ldGroup$ldPlotImage$setState(
          list(ld_store = ld_store, nms = nms, metric = opts$ldMetric))
      }
    },

    .render_ld_plot = function(image, ggtheme, theme, ...) {
      state <- image$state
      if (is.null(state)) return(FALSE)
      ld_store <- state$ld_store; nms <- state$nms; metric <- state$metric; n <- length(nms)

      metric_label <- switch(metric, Dprime="D'", r2="r²", D="D")
      df_rows <- list()
      for (i in seq_len(n)) for (j in seq_len(n)) {
        val <- if (i==j) 1.0 else {
          key <- paste(c(nms[min(i,j)], nms[max(i,j)]), collapse="___")
          ld_res <- ld_store[[key]]
          if (!is.null(ld_res)) switch(metric,
            Dprime = abs(as.numeric(ld_res$`D'`)),
            r2     = as.numeric(ld_res$`r`)^2,
            D      = abs(as.numeric(ld_res$`D`))) else NA_real_
        }
        df_rows[[length(df_rows)+1L]] <- data.frame(
          SNP1  = factor(nms[i], levels=rev(nms)),
          SNP2  = factor(nms[j], levels=nms),
          value = val, stringsAsFactors=FALSE)
      }
      df <- do.call(rbind, df_rows)

      p_mat <- matrix(NA_real_, n, n, dimnames=list(nms,nms))
      for (pk in names(ld_store)) {
        parts <- strsplit(pk,"___")[[1]]
        pv    <- ld_store[[pk]]$`P-value`
        p_mat[parts[1],parts[2]] <- pv; p_mat[parts[2],parts[1]] <- pv
      }

      df$label <- ""
      for (k in seq_len(nrow(df))) {
        i_nm <- as.character(df$SNP1[k]); j_nm <- as.character(df$SNP2[k])
        i_idx <- which(nms==i_nm); j_idx <- which(nms==j_nm)
        if (i_idx > j_idx) {
          pv <- p_mat[i_nm, j_nm]
          df$label[k] <- if (!is.na(pv)) { if(pv<0.001) "<.001" else sprintf("%.3f",pv) } else ""
        } else if (i_idx < j_idx) {
          key <- paste(c(nms[min(i_idx,j_idx)], nms[max(i_idx,j_idx)]), collapse="___")
          ld_res <- ld_store[[key]]
          if (!is.null(ld_res)) {
            raw <- switch(metric, r2=ld_res$`r`^2, Dprime=ld_res$`D'`, D=ld_res$`D`)
            df$label[k] <- sprintf("%.3f", round(as.numeric(raw),3))
          }
        } else { df$label[k] <- i_nm }
      }

      colour_label <- switch(metric, Dprime="|D'|", r2="r²", D="|D|")
      p <- ggplot2::ggplot(df, ggplot2::aes(x=SNP2, y=SNP1, fill=value)) +
        ggplot2::geom_tile(colour="white", linewidth=0.5) +
        ggplot2::geom_text(ggplot2::aes(label=label), size=3, colour="grey10") +
        ggplot2::scale_fill_gradientn(
          colours  = c("#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b"),
          limits   = c(0,1), na.value="grey85", name=colour_label) +
        ggplot2::scale_x_discrete(position="bottom") +
        ggplot2::labs(title=paste0("LD Heatmap  •  upper: ",metric_label," | lower: p-value"),
                      x=NULL, y=NULL) +
        ggplot2::theme_minimal(base_size=11) +
        ggplot2::theme(
          axis.text.x=ggplot2::element_text(angle=45,hjust=1,vjust=1),
          axis.text.y=ggplot2::element_text(hjust=1),
          panel.grid=ggplot2::element_blank(),
          legend.position="right",
          plot.title=ggplot2::element_text(size=11,face="bold",
                                           margin=ggplot2::margin(b=8)))
      print(p); TRUE
    },

    # ── Haplotypes ───────────────────────────────────────────────────────────
    .run_haplo = function(geno_list, data, response, response_raw, response_type,
                      cov_df, opts, run_haploFreq, run_haploAssoc,
                      run_haploInteraction, run_subpop, complete_mask) {
    
      # ── Common Data Prep (Logic remains the same) ──
      snp_names   <- names(geno_list)
      allele_mat  <- do.call(cbind, lapply(snp_names, function(nm) genetics::allele(geno_list[[nm]])))
      geno_setup  <- tryCatch(haplo.stats::setupGeno(allele_mat, locus.label = snp_names), 
                              error = function(e) NULL)
      
      if (is.null(geno_setup)) return()

      # Metadata needed by sub-functions
      u_alleles <- attr(geno_setup, "unique.alleles")

      # Missing Management
      # Count only rows newly excluded at this stage: rows passing complete_mask
      # that have at least all SNP missing.  Rows already excluded by
      # complete_mask (response/covariate missing) are NOT double-counted.
      snp_miss_mask <- apply(is.na(allele_mat), 1, all)
      keep   <- complete_mask & !snp_miss_mask
      n_miss <- sum(snp_miss_mask & complete_mask)

      # ── Dispatch to specialized methods ──
      if (run_haploFreq) {
          private$.compute_haplo_freqs(geno_setup, response_raw, response_type, keep, 
                                      n_miss, opts, run_subpop, snp_names, u_alleles)
      }

      if (run_haploAssoc && !is.null(response)) {
          private$.compute_haplo_assoc(geno_setup, response, response_type, cov_df, keep, 
            n_miss, opts, snp_names, u_alleles)
      }

      if (run_haploInteraction && !is.null(cov_df) && !is.null(response)) {
          private$.compute_haplo_interaction(geno_setup, response, response_type, cov_df, keep, 
            n_miss, opts, snp_names, u_alleles)
      }
    },

    .compute_haplo_freqs = function(geno_setup, response_raw, response_type, keep, 
                                    n_miss, opts, run_subpop, snp_names, u_alleles) {
      
      tbl <- self$results$haploGroup$haploFreqTable
      tbl$setTitle("<b>Haplotype Frequencies</b>")
      do_strat_haplo <- isTRUE(run_subpop) && !is.null(response_raw) &&
                        identical(response_type, "binary")
      grp_levels_haplo <- levels(response_raw[keep])

      if (do_strat_haplo) {
        tbl$addColumn(name="freq_g0", title=as.character(grp_levels_haplo[1]),
                      type="number", format="zto")
        tbl$addColumn(name="freq_g1", title=as.character(grp_levels_haplo[2]),
                      type="number", format="zto")
      }

      em_all <- tryCatch(
        haplo.stats::haplo.em(subset_geno(geno_setup, keep), locus.label=snp_names),
        error=function(e) NULL)
      if (!is.null(em_all)) {
        freqs    <- em_all$hap.prob
        rare_sum <- 0
        em_grp   <- list()
        grp_freq <- list()
        if (do_strat_haplo) {
          for (lvl in grp_levels_haplo) {
            keep_lvl <- keep & as.character(response_raw)==lvl
            if (sum(keep_lvl) < 5) next
            em_grp[[lvl]] <- tryCatch(
              haplo.stats::haplo.em(subset_geno(geno_setup, keep_lvl), locus.label=snp_names),
              error=function(e) NULL)
          }
          grp_freq <- lapply(em_grp, function(em_g) {
            if (is.null(em_g)) return(list())
            setNames(as.list(round(em_g$hap.prob,3)),
                      sapply(seq_len(nrow(em_g$haplotype)), function(j)
                        decode_haplo_row(as.numeric(em_g$haplotype[j,]), u_alleles)))
          })
          grp_levels <- levels(as.factor(response_raw[keep]))
          
          # Toggle visibility and update titles based on actual group names
          tbl$getColumn('freq_g0')$setVisible(TRUE)
          tbl$getColumn('freq_g0')$setTitle(as.character(grp_levels[1]))
          
          tbl$getColumn('freq_g1')$setVisible(TRUE)
          tbl$getColumn('freq_g1')$setTitle(as.character(grp_levels[2]))
        } else {
            # Hide them if the user unchecks the option
            tbl$getColumn('freq_g0')$setVisible(FALSE)
            tbl$getColumn('freq_g1')$setVisible(FALSE)
        }
        # Sort indices by frequency descending (rare rows are excluded in the loop)
        sorted_idx <- order(freqs, decreasing = TRUE)
        for (i in sorted_idx) {
          if (freqs[i] < opts$haploFreqMin) { rare_sum <- rare_sum + freqs[i]; next }
          label    <- decode_haplo_row(as.numeric(em_all$haplotype[i,]), u_alleles)
          row_vals <- list(haplotype=label, freq=round(freqs[i],3))
          if (do_strat_haplo) {
            row_vals$freq_g0 <- grp_freq[[grp_levels_haplo[1]]][[label]] %||% NA_real_
            row_vals$freq_g1 <- grp_freq[[grp_levels_haplo[2]]][[label]] %||% NA_real_
          }
          tbl$addRow(rowKey=paste0("f",i), values=row_vals)
        }
        if (rare_sum > 0) {
          row_vals <- list(haplotype=paste0("Rare (<",opts$haploFreqMin,")"),
                            freq=round(rare_sum,3))
          if (do_strat_haplo) {
            em0 <- em_grp[[grp_levels_haplo[1]]]; em1 <- em_grp[[grp_levels_haplo[2]]]
            rare_g0 <- if (!is.null(em0))
              round(sum(em0$hap.prob[em0$hap.prob < opts$haploFreqMin]),3) else NA_real_
            rare_g1 <- if (!is.null(em1))
              round(sum(em1$hap.prob[em1$hap.prob < opts$haploFreqMin]),3) else NA_real_
            row_vals$freq_g0 <- if(!is.na(rare_g0) && rare_g0>0) rare_g0 else NA_real_
            row_vals$freq_g1 <- if(!is.na(rare_g1) && rare_g1>0) rare_g1 else NA_real_
          }
          tbl$addRow(rowKey="rare_freq", values=row_vals)
        }
      }
      if (n_miss > 0)
        tbl$setNote(note=paste0(n_miss," observation(s) with missing data excluded."),
                    key="missing_snp")
      else
        tbl$setNote(note=NULL, key="missing_snp")
    },

    .compute_haplo_assoc = function(geno_setup, response, response_type, cov_df, keep, 
                                    n_miss, opts, snp_names, u_alleles) {

      # ── Haplotype association ─────────────────────────────────────

        family     <- if (response_type == "binary") "binomial" else "gaussian"
        y_sub      <- if (response_type == "binary") {
          as.numeric(as.factor(response[keep])) - 1L
        } else {
          response[keep]
        }

        m_model      <- data.frame(y = y_sub)
        m_model$geno <- subset_geno(geno_setup, keep)
        if (!is.null(cov_df)) {
          m_model    <- cbind(m_model, cov_df[keep, , drop = FALSE])
          formula_str <- paste("y ~ geno +", paste(names(cov_df), collapse = " + "))
        } else {
          formula_str <- "y ~ geno"
        }

        haplo_fit <- tryCatch(
          haplo.stats::haplo.glm(
            as.formula(formula_str),
            family    = family,
            data      = m_model,
            na.action = na.geno.keep,
            control   = haplo.stats::haplo.glm.control(
                          haplo.freq.min = opts$haploFreqMin)
          ),
          error = function(e) {
            self$results$validationMsg$setContent(
              paste0("<b>Haplotype GLM error:</b> ", e$message))
            NULL
          }
        )

        if (!is.null(haplo_fit)) {
          tbl <- self$results$haploGroup$haploAssocTable
          tbl$setTitle("<b>Haplotype Association</b>")


          # Set effect column title to match response type
          tbl$getColumn("effect")$setTitle(
            if (response_type == "binary") "OR" else "β")

          # ── LRT: overall haplotype association ─────────────────────
          # Fit a null model (covariates only, no haplotype term) using plain
          # glm/lm — haplo.glm requires a geno term so cannot be used here.
          # Compare deviances against haplo_fit directly.
          null_formula_str <- if (!is.null(cov_df) && ncol(cov_df) > 0)
            paste("y ~", paste(names(cov_df), collapse = " + "))
          else
            "y ~ 1"

          haplo_null_fit <- tryCatch(
            if (family == "binomial")
              glm(as.formula(null_formula_str), family = binomial(), data = m_model)
            else
              lm(as.formula(null_formula_str), data = m_model),
            error = function(e) NULL
          )
          # debug

          p_lrt_assoc <- NA_real_
          if (!is.null(haplo_null_fit)) {
            dev_diff <- deviance(haplo_null_fit) - haplo_fit$deviance
            df_diff  <- (haplo_fit$df.null - haplo_fit$df.residual) - (haplo_null_fit$df.null - haplo_null_fit$df.residual)
            p_lrt_assoc <- if (!is.na(dev_diff) && !is.na(df_diff) && df_diff > 0)
              pchisq(dev_diff, df = df_diff, lower.tail = FALSE)
            else NA_real_
          }
          tbl$setNote(
            note = paste0("Likelihood ratio test for overall haplotype association: P = ",
                          format.pval(p_lrt_assoc, digits = 3)),
            key  = "lrt_assoc")

          # ── Decode haplotype label from haplo.unique row ────────────
          # haplo.unique stores allele characters directly (e.g. "C", "T", "A"),
          # one per locus — confirmed from diagnostic output.
          label_from_unique_row <- function(row_vec) {
            paste(as.character(row_vec), collapse = "-")
          }

          # ── Pull coefficients and CIs ───────────────────────────────
          # haplo.glm model matrix column names for haplotype terms are stored
          # in haplo_fit$haplo.names.  The actual rownames of coef() follow the
          # convention "geno" + haplo.names (the model frame column is "geno").
          # We match positionally rather than by name to be robust to separator
          # differences across haplo.stats versions.
          coef_sum <- tryCatch(summary(haplo_fit)$coefficients, error = function(e) NULL)
          ci_mat   <- tryCatch(confint(haplo_fit, level = opts$ciWidth / 100),
                              error = function(e) NULL)

          # summary(haplo_fit)$coefficients columns are: coef | SE | t.stat | pval
          # (haplo.glm uses its own summary method, not summary.glm)
          # Modern haplo.stats uses "SE" (capital); guard against both spellings.
          haplo_rows <- if (!is.null(coef_sum)) {
            grep("^geno", rownames(coef_sum))
          } else integer(0)
          se_col <- if (!is.null(coef_sum) && "SE" %in% colnames(coef_sum)) "SE" else "se"

          # Helper: get beta, se, pval, ci for haplo-term at position pos
          # (1-based within haplo_rows, matching order of haplo.common).
          get_stats <- function(pos) {
            row_idx <- if (!is.na(pos) && pos >= 1L && pos <= length(haplo_rows))
                        haplo_rows[pos] else NA_integer_
            if (is.na(row_idx) || is.null(coef_sum) ||
                row_idx < 1L || row_idx > nrow(coef_sum)) {
              return(list(beta = NA_real_, se = NA_real_, pval = NA_real_,
                          ci_lo = NA_real_, ci_hi = NA_real_))
            }
            rn   <- rownames(coef_sum)[row_idx]
            beta <- coef_sum[row_idx, "coef"]
            se   <- coef_sum[row_idx, se_col]
            pval <- coef_sum[row_idx, "pval"]
            # CI from ci_mat if available, otherwise Wald ± z * se
            if (!is.null(ci_mat) && rn %in% rownames(ci_mat)) {
              ci_lo <- ci_mat[rn, 1]
              ci_hi <- ci_mat[rn, 2]
            } else {
              z     <- qnorm(1 - (1 - opts$ciWidth / 100) / 2)
              ci_lo <- beta - z * se
              ci_hi <- beta + z * se
            }
            list(beta = beta, se = se, pval = pval, ci_lo = ci_lo, ci_hi = ci_hi)
          }

          make_row <- function(label, freq, stats) {
            b  <- stats$beta
            lo <- stats$ci_lo
            hi <- stats$ci_hi
            list(
              haplotype = label,
              freq      = round(freq, 4),
              effect    = if (response_type == "binary") exp(b)  else b,
              ciLow     = if (response_type == "binary") exp(lo) else lo,
              ciHigh    = if (response_type == "binary") exp(hi) else hi,
              pval      = stats$pval
            )
          }

          # ── Reference haplotype (OR = 1 by definition) ────────────
          base_idx   <- haplo_fit$haplo.base
          base_label <- label_from_unique_row(haplo_fit$haplo.unique[base_idx, ])
          base_freq  <- haplo_fit$haplo.freq[base_idx]
          tbl$addRow(rowKey = "base", values = list(
            haplotype   = paste0(base_label, " (Ref)"),
            freq      = round(base_freq, 4),
            effect    = if (response_type == "binary") 1.0 else 0.0,
            ciLow     = '',
            ciHigh    = '',
            pval      = ''
          ))

          # ── Common haplotypes ──────────────────────────────────────
          # haplo.common: integer index vector into haplo.unique rows,
          # in the same order as haplo.names / the GLM coefficients.
          # Sort by frequency descending so most-frequent haplotypes appear first.
          common_idx   <- haplo_fit$haplo.common
          common_freqs <- haplo_fit$haplo.freq[common_idx]
          sorted_j     <- order(common_freqs, decreasing = TRUE)
          for (j in sorted_j) {
            h_idx   <- common_idx[j]
            h_label <- label_from_unique_row(haplo_fit$haplo.unique[h_idx, ])
            h_freq  <- haplo_fit$haplo.freq[h_idx]
            stats   <- get_stats(j)
            tbl$addRow(rowKey = paste0("h", j),
                      values = make_row(h_label, h_freq, stats))
          }

          # ── Rare combined term ─────────────────────────────────────
          has_rare <- isTRUE(haplo_fit$haplo.rare.term) ||
                      (length(haplo_fit$haplo.rare) > 0)
          if (has_rare) {
            rare_freq <- sum(haplo_fit$haplo.freq[haplo_fit$haplo.rare])
            # rare term is the last haplotype coefficient
            stats     <- get_stats(length(common_idx) + 1L)
            tbl$addRow(rowKey = "rare",
                      values = make_row(
                        paste0("Rare (<", opts$haploFreqMin, ")"),
                        rare_freq, stats))
          }

                # ── Add covariate note ─────────────────────────────
          note_key <- "covariates"
          if (!is.null(cov_df) && ncol(cov_df) > 0) {
            cov_names <- names(cov_df)
            cov_names <- sapply(cov_names, function(x) {
              if (!is.null(self$data[[x]])) {
                attr(self$data[[x]], "label") %||% x
              } else x
            })
            note_txt <- paste0("Model adjusted for: ", paste(cov_names, collapse = ", "))
            tbl$setNote(note = note_txt, key = note_key)
          } else {
            tbl$setNote(note = NULL, key = note_key)
          }

          # ── Missing note
          if (n_miss > 0)
            tbl$setNote(
              note = paste0(n_miss,
                            " observation(s) with missing data excluded."),
              key  = "missing_snp")
          else
            tbl$setNote(note = NULL, key = "missing_snp")

        }
    },

    .compute_haplo_interaction_ = function(geno_setup, response, response_type, cov_df, keep,
                                          n_miss, opts, snp_names, u_alleles) {
      # ── Haplotype × covariate interaction ─────────────────────────

      int_var  <- names(cov_df)[1]   # first covariate is always the interaction term
      adj_vars <- setdiff(names(cov_df), int_var)

      tbl_int <- self$results$haploGroup$haploInteractionTable
      tbl_int$getColumn("effect")$setTitle(
        if (response_type == "binary") "OR" else "\u03B2")

      # ── Interaction type ──────────────────────────────────────────
      # "multiplicative"        y ~ geno * int_var   (default)
      # "conditional_on_haplo"  y ~ geno / int_var   (covar effect within each haplotype)
      # "conditional_on_covar"  y ~ int_var / geno   (haplotype effect within each covar level)
      int_type <- if (is.null(opts$haploInteractionType)) "multiplicative"
                  else opts$haploInteractionType

      # Human-readable formula token for table title / note
      formula_token <- switch(int_type,
        multiplicative        = paste0("Haplotype \u00D7 ", int_var),
        conditional_on_covar  = paste0(int_var, " | Haplotype"),
        conditional_on_haplo  = paste0("Haplotype | ", int_var))
      tbl_int$setTitle(paste0("<b>Interaction: ", formula_token, "</b>"))

      if (length(adj_vars) > 0){
        note_parts <- paste0(". Adjusted for: ", paste(adj_vars, collapse = ", "))
        tbl_int$setNote(note = note_parts, key = "intcov")
      }

      # ── Missing note ──────────────────────────────────────────────
      if (n_miss > 0)
        tbl_int$setNote(
          note = paste0(n_miss, " observation(s) with missing data excluded."),
          key  = "missing_snp")
      else
        tbl_int$setNote(note = NULL, key = "missing_snp")

      family_int <- if (response_type == "binary") "binomial" else "gaussian"
      y_int <- if (response_type == "binary")
        as.numeric(as.factor(response[keep])) - 1L
      else
        response[keep]

      m_int      <- data.frame(y = y_int)
      m_int$geno <- subset_geno(geno_setup, keep)
      if (!is.null(cov_df))
        m_int <- cbind(m_int, cov_df[keep, , drop = FALSE])

      adj_part <- if (length(adj_vars) > 0)
        paste("+", paste(adj_vars, collapse = "+")) else ""

      # ── Build formulae depending on parameterisation ──────────────
      # The "additive" formula (no interaction) is always used for the LRT.
      formula_add_str <- paste("y ~ geno +", int_var, adj_part)

      formula_fit_str <- switch(int_type,
        multiplicative       = paste("y ~ geno *",   int_var, adj_part),
        conditional_on_haplo = paste("y ~ geno /",   int_var, adj_part),
        conditional_on_covar = paste("y ~", int_var, "/ geno", adj_part)
      )

      haplo_fit_int <- tryCatch(
        haplo.stats::haplo.glm(
          as.formula(formula_fit_str),
          family    = family_int,
          data      = m_int,
          na.action = na.geno.keep,
          control   = haplo.stats::haplo.glm.control(haplo.freq.min = opts$haploFreqMin)
        ),
        error = function(e) {
          self$results$validationMsg$setContent(
            paste0("<b>Haplotype interaction GLM error:</b> ", e$message))
          NULL
        }
      )
      haplo_fit_add <- tryCatch(
        haplo.stats::haplo.glm(
          as.formula(formula_add_str),
          family    = family_int,
          data      = m_int,
          na.action = na.geno.keep,
          control   = haplo.stats::haplo.glm.control(haplo.freq.min = opts$haploFreqMin)
        ),
        error = function(e) NULL
      )

      p_inter_haplo <- NA_real_   # initialised here so note is always settable

      if (!is.null(haplo_fit_int) && !is.null(haplo_fit_add)) {

        # ── LRT: fitted model vs additive baseline ────────────────────
        # anova.haplo.glm can fail on mismatched EM sets; use deviance diff directly.
        dev_diff <- haplo_fit_add$deviance - haplo_fit_int$deviance
        df_diff  <- haplo_fit_add$df.residual - haplo_fit_int$df.residual
        p_inter_haplo <- if (!is.na(dev_diff) && !is.na(df_diff) && df_diff > 0)
          pchisq(dev_diff, df = df_diff, lower.tail = FALSE)
        else NA_real_

        coef_sum_int <- tryCatch(summary(haplo_fit_int)$coefficients,
                                 error = function(e) NULL)
        ci_int       <- tryCatch(confint(haplo_fit_int, level = opts$ciWidth / 100),
                                 error = function(e) NULL)
        se_col_int   <- if (!is.null(coef_sum_int) &&
                            "SE" %in% colnames(coef_sum_int)) "SE" else "se"

        if (!is.null(coef_sum_int)) {
          all_rows_int <- rownames(coef_sum_int)

          # ── Shared helper: decode haplo.unique row → allele label ──
          decode_haplo_label <- function(row_vec)
            paste(as.character(row_vec), collapse = "-")

          rare_label <- paste0("Rare (<", opts$haploFreqMin, ")")

          # ── Build coef-name → display-label map ───────────────────
          # haplo.glm names geno main-effect terms as "geno.N" or "geno.rare";
          # nested terms take the form "geno.N:int_varLevel" (multiplicative /
          # conditional_on_haplo) or "int_varLevel:geno.N" (conditional_on_covar).
          base_idx   <- haplo_fit_int$haplo.base
          base_label <- decode_haplo_label(haplo_fit_int$haplo.unique[base_idx, ])

          geno_main_rows <- grep("^geno[^:]+$", all_rows_int, value = TRUE)
          raw_to_label   <- character(0)

          for (rn in geno_main_rows) {
            suffix <- sub("^geno\\.", "", rn)
            display_label <- if (grepl("^[0-9]+$", suffix)) {
              idx <- as.integer(suffix)
              if (!is.na(idx) && idx >= 1L && idx <= nrow(haplo_fit_int$haplo.unique))
                decode_haplo_label(haplo_fit_int$haplo.unique[idx, ])
              else paste0("Haplotype ", suffix)
            } else if (grepl("rare", suffix, ignore.case = TRUE)) {
              rare_label
            } else {
              suffix
            }
            raw_to_label[rn] <- display_label

            # Nested terms: "geno.N:int_varLevel" (geno / covar direction)
            inter_rns <- grep(paste0("^", rn, ":"), all_rows_int, value = TRUE)
            for (irn in inter_rns)
              raw_to_label[irn] <- paste0(
                display_label, " | ",
                sub(paste0("^", rn, ":"), "", irn))
          }

          # For conditional_on_covar (y ~ int_var / geno), the model parameterises
          # as: intercept = base haplotype (reference covariate level),
          #     int_varLevel rows = int_varLevel WITHIN base haplotype,
          #     int_varLevel:geno.N rows = haplotype N effect within int_varLevel.
          if (int_type == "conditional_on_covar") {

            # Covariate main-effect rows (e.g. "SEXMale") — these are the effect
            # of that covariate level WITHIN the base haplotype, so label them
            # "SEXMale:<base_haplotype>"
            covar_main_rows <- grep(paste0("^", int_var, "[^:]*$"), all_rows_int, value = TRUE)
            for (crn in covar_main_rows) {
              covar_level <- sub(paste0("^", int_var), "", crn)
              raw_to_label[crn] <- paste0(int_var, covar_level, ":", base_label)

              # Nested "int_varLevel:geno.N" rows — haplotype N within that level
              nested_rns <- grep(paste0("^", crn, ":geno"), all_rows_int, value = TRUE)
              for (nrn in nested_rns) {
                geno_part <- sub(paste0("^", crn, ":"), "", nrn)
                g_suffix  <- sub("^geno\\.", "", geno_part)
                g_label   <- if (grepl("^[0-9]+$", g_suffix)) {
                  idx <- as.integer(g_suffix)
                  if (!is.na(idx) && idx >= 1L && idx <= nrow(haplo_fit_int$haplo.unique))
                    decode_haplo_label(haplo_fit_int$haplo.unique[idx, ])
                  else paste0("Haplotype ", g_suffix)
                } else if (grepl("rare", g_suffix, ignore.case = TRUE)) {
                  rare_label
                } else g_suffix
                raw_to_label[nrn] <- paste0(int_var, covar_level, ":", g_label)
              }
            }
          }

          # ── Classify coefficient rows ──────────────────────────────
          # "main" geno rows: always present, represent haplotype main effects
          # "nested" rows: contain ":" — these are the conditional/interaction terms
          # Rows for int_var itself (covariate main effect) are present in
          # multiplicative and conditional_on_covar models but not shown.
          main_rows   <- grep("^geno[^:]+$", all_rows_int)
          nested_rows <- grep(":", all_rows_int)  # all terms with ":"

          # ── Sort main rows by frequency descending ─────────────────
          get_freq <- function(r) {
            suffix <- sub("^geno\\.", "", all_rows_int[r])
            if (grepl("^[0-9]+$", suffix)) {
              idx <- as.integer(suffix)
              if (!is.na(idx) && idx >= 1L && idx <= length(haplo_fit_int$haplo.freq))
                return(haplo_fit_int$haplo.freq[idx])
            }
            0  # rare / unknown → sort last
          }
          main_freqs      <- sapply(main_rows, get_freq)
          sorted_main_pos <- main_rows[order(main_freqs, decreasing = TRUE)]

          # ── Build show_rows: each main row followed by its nested rows ──
          show_rows <- integer(0)
          for (mr in sorted_main_pos) {
            show_rows <- c(show_rows, mr)
            rn_mr     <- all_rows_int[mr]
            child_rows <- nested_rows[grepl(paste0("^", rn_mr, ":"),
                                            all_rows_int[nested_rows])]
            if (length(child_rows) > 0)
              show_rows <- c(show_rows, child_rows)
          }
          # For conditional_on_covar the nested rows are "covarLevel:geno.N";
          # they won't be caught by the prefix match above — append them now.
          show_rows <- c(show_rows, setdiff(nested_rows, show_rows))

          # ── Helper: extract beta / CI for one row ──────────────────
          get_row_stats <- function(r) {
            raw_nm <- all_rows_int[r]
            beta   <- coef_sum_int[r, "coef"]
            pval   <- coef_sum_int[r, "pval"]
            if (!is.null(ci_int) && raw_nm %in% rownames(ci_int)) {
              ci_lo <- ci_int[raw_nm, 1]; ci_hi <- ci_int[raw_nm, 2]
            } else {
              z     <- qnorm(1 - (1 - opts$ciWidth / 100) / 2)
              se    <- coef_sum_int[r, se_col_int]
              ci_lo <- beta - z * se; ci_hi <- beta + z * se
            }
            list(beta = beta, pval = pval, ci_lo = ci_lo, ci_hi = ci_hi)
          }

          # ── Reference row ──────────────────────────────────────────
          # For multiplicative and conditional_on_covar: base haplotype is reference.
          # For conditional_on_haplo (geno / covar): the base haplotype row itself
          # is the "intercept within base haplotype" and is not in coef_sum_int;
          # we still emit a Ref row for orientation.
          # (base_idx / base_label already computed above for label map)

          tbl_int$addRow(
            rowKey = "base",
            values = list(
              term   = paste0(base_label, " (Ref)"),
              effect = if (response_type == "binary") 1.0 else 0.0,
              ciLow  = '', ciHigh = '', pval = ''
            )
          )

          # ── Emit rows ─────────────────────────────────────────────
          for (r in show_rows) {
            raw_nm <- all_rows_int[r]
            label  <- raw_to_label[raw_nm]

            # Fallback label if map lookup missed this row
            if (is.na(label) || length(label) == 0) {
              label <- raw_nm   # raw coefficient name is always better than nothing
            }

            st <- get_row_stats(r)
            tbl_int$addRow(
              rowKey = paste0("hi", r),
              values = list(
                term   = label,
                effect = if (response_type == "binary") exp(st$beta)  else st$beta,
                ciLow  = if (response_type == "binary") exp(st$ci_lo) else st$ci_lo,
                ciHigh = if (response_type == "binary") exp(st$ci_hi) else st$ci_hi,
                pval   = st$pval
              )
            )
          }
        }
      }

      # ── LRT note (set regardless of fit outcome) ──────────────────
      if (!is.na(p_inter_haplo))
        tbl_int$setNote(
          note = paste0("Likelihood ratio test for interaction (vs additive): P = ",
                        format.pval(p_inter_haplo, digits = 3)),
          key  = "lrt_inter")
      else
        tbl_int$setNote(note = NULL, key = "lrt_inter")
    },


    .compute_haplo_interaction = function(geno_setup, response, response_type, cov_df, keep,
                                          n_miss, opts, snp_names, u_alleles) {
      
      int_var  <- names(cov_df)[1]   # first covariate is the interaction term
      adj_vars <- setdiff(names(cov_df), int_var)

      tbl_int <- self$results$haploGroup$haploInteractionTable
      tbl_int$getColumn("effect")$setTitle(
        if (response_type == "binary") "OR" else "\u03B2")

      int_type <- if (is.null(opts$haploInteractionType)) "multiplicative"
                  else opts$haploInteractionType

      # Set Title
      formula_token <- switch(int_type,
        multiplicative        = paste0("Haplotype \u00D7 ", int_var),
        conditional_on_covar  = paste0(int_var, " | Haplotype"),
        conditional_on_haplo  = paste0("Haplotype | ", int_var))
      tbl_int$setTitle(paste0("<b>Interaction: ", formula_token, "</b>"))

      if (length(adj_vars) > 0){
        note_parts <- paste0(". Adjusted for: ", paste(adj_vars, collapse = ", "))
        tbl_int$setNote(note = note_parts, key = "intcov")
      }

      # ── Missing note ──────────────────────────────────────────────
      if (n_miss > 0)
        tbl_int$setNote(
          note = paste0(n_miss, " observation(s) with missing data excluded."),
          key  = "missing_snp")
      else
        tbl_int$setNote(note = NULL, key = "missing_snp")

      family_int <- if (response_type == "binary") "binomial" else "gaussian"
      y_int <- if (response_type == "binary")
        as.numeric(as.factor(response[keep])) - 1L
      else
        response[keep]

      m_int      <- data.frame(y = y_int)
      m_int$geno <- subset_geno(geno_setup, keep)
      if (!is.null(cov_df))
        m_int <- cbind(m_int, cov_df[keep, , drop = FALSE])

      adj_part <- if (length(adj_vars) > 0)
        paste("+", paste(adj_vars, collapse = "+")) else ""

      # ── Build formulae depending on parameterisation ──────────────
      # The "additive" formula (no interaction) is always used for the LRT.
      formula_add_str <- paste("y ~ geno +", int_var, adj_part)

      formula_fit_str <- switch(int_type,
        multiplicative       = paste("y ~ geno *",   int_var, adj_part),
        conditional_on_haplo = paste("y ~ geno /",   int_var, adj_part),
        conditional_on_covar = paste("y ~", int_var, "/ geno", adj_part)
      )

      haplo_fit_int <- tryCatch(
        haplo.stats::haplo.glm(
          as.formula(formula_fit_str),
          family    = family_int,
          data      = m_int,
          na.action = na.geno.keep,
          control   = haplo.stats::haplo.glm.control(haplo.freq.min = opts$haploFreqMin)
        ),
        error = function(e) {
          self$results$validationMsg$setContent(
            paste0("<b>Haplotype interaction GLM error:</b> ", e$message))
          NULL
        }
      )
      haplo_fit_add <- tryCatch(
        haplo.stats::haplo.glm(
          as.formula(formula_add_str),
          family    = family_int,
          data      = m_int,
          na.action = na.geno.keep,
          control   = haplo.stats::haplo.glm.control(haplo.freq.min = opts$haploFreqMin)
        ),
        error = function(e) NULL
      )


      if (!is.null(haplo_fit_int) && !is.null(haplo_fit_add)) {

        p_inter_haplo <- NA_real_   # initialised here so note is always settable

        # ── LRT: fitted model vs additive baseline ────────────────────
        # anova.haplo.glm can fail on mismatched EM sets; use deviance diff directly.
        dev_diff <- haplo_fit_add$deviance - haplo_fit_int$deviance
        df_diff  <- haplo_fit_add$df.residual - haplo_fit_int$df.residual
        p_inter_haplo <- if (!is.na(dev_diff) && !is.na(df_diff) && df_diff > 0)
          pchisq(dev_diff, df = df_diff, lower.tail = FALSE)
        else NA_real_

        if (!is.na(p_inter_haplo))
          tbl_int$setNote(
            note = paste0("Likelihood ratio test for interaction (vs additive): P = ",
                          format.pval(p_inter_haplo, digits = 3)),
            key  = "lrt_inter")
        else
          tbl_int$setNote(note = NULL, key = "lrt_inter")

        coef_sum_int <- tryCatch(summary(haplo_fit_int)$coefficients, error = function(e) NULL)
        ci_int       <- tryCatch(confint(haplo_fit_int, level = opts$ciWidth / 100), error = function(e) NULL)
        se_col_int   <- if (!is.null(coef_sum_int) && "SE" %in% colnames(coef_sum_int)) "SE" else "se"

        if (!is.null(coef_sum_int)) {
          all_rows_int <- rownames(coef_sum_int)
          decode_haplo_label <- function(row_vec) paste(as.character(row_vec), collapse = "-")
          rare_label <- paste0("Rare (<", opts$haploFreqMin, ")")
          
          base_idx   <- haplo_fit_int$haplo.base
          base_label <- decode_haplo_label(haplo_fit_int$haplo.unique[base_idx, ])
          
          # We need to find the reference level of the interaction variable
          # haplo.glm usually leaves the first level as the baseline in the intercept
          int_var_factor <- as.factor(cov_df[[int_var]])
          ref_covar_lvl  <- levels(int_var_factor)[1]

          raw_to_label <- character(0)

          # 1. Map geno main suffixes (geno.1, geno.rare) to Allele Strings
          for (rn in all_rows_int) {
            # Extract geno part from strings like "geno.12" or "SEXFemale:geno.12"
            if (grepl("geno\\.", rn)) {
              g_match  <- regmatches(rn, regexpr("geno\\.([0-9]+|rare)", rn))
              g_suffix <- sub("geno\\.", "", g_match)
              
              g_label <- if (grepl("^[0-9]+$", g_suffix)) {
                idx <- as.integer(g_suffix)
                decode_haplo_label(haplo_fit_int$haplo.unique[idx, ])
              } else { rare_label }
              
              # Replace the geno component in the raw name with the label
              if (grepl(":", rn)) {
                # Handle nested terms: e.g., "SEXFemale:geno.12" -> "SEXFemale:G-C-A"
                raw_to_label[rn] <- sub("geno\\.([0-9]+|rare)", g_label, rn)
              } else {
                # Handle main geno terms: e.g., "geno.12" -> "G-C-A"
                raw_to_label[rn] <- g_label
              }
            }
          }

          # 2. Specifically handle Covariate Main Effects for "conditional models"
          covar_main_rows <- grep(paste0("^", int_var, "[^:]*$"), all_rows_int, value = TRUE)
          if (int_type == "conditional_on_covar") {
            # These appear as "SEXFemale" but represent "SEXFemale:BaseHaplotype"
            for (crn in covar_main_rows) {
               if (crn != covar_main_rows[1])
                 raw_to_label[crn] <- paste0(crn, ":", base_label)
              else
                raw_to_label[crn] <- crn
            }
          } else if (int_type == "conditional_on_haplo") {
            # These appear as "SEXFemale" but represent "BaseHaplotype:SEXFemale"
            for (crn in covar_main_rows) {
               raw_to_label[crn] <- paste0(base_label, ":", crn)
            }
          } else {
            for (crn in covar_main_rows) {
               raw_to_label[crn] <- crn
            }
          }
          
          
          # 3. Identify rows to show
          # here add option to show all terms
          show_mask <- grepl("geno", all_rows_int) | grepl(int_var, all_rows_int)           
          show_rows <- which(show_mask)

          # --- Emit Rows ---
          
          # Reference Row: Now specifically labeled with Ref Covariate
          tbl_int$addRow(
            rowKey = "base",
            values = list(
              term   = paste0(base_label, " (Ref)"),
              effect = if (response_type == "binary") 1.0 else 0.0,
              ciLow  = '', ciHigh = '', pval = ''
            )
          )

          # ── Helper: extract beta / CI for one row ──────────────────
          get_row_stats <- function(r) {
            raw_nm <- all_rows_int[r]
            beta   <- coef_sum_int[r, "coef"]
            pval   <- coef_sum_int[r, "pval"]
            if (!is.null(ci_int) && raw_nm %in% rownames(ci_int)) {
              ci_lo <- ci_int[raw_nm, 1]; ci_hi <- ci_int[raw_nm, 2]
            } else {
              z     <- qnorm(1 - (1 - opts$ciWidth / 100) / 2)
              se    <- coef_sum_int[r, se_col_int]
              ci_lo <- beta - z * se; ci_hi <- beta + z * se
            }
            list(beta = beta, pval = pval, ci_lo = ci_lo, ci_hi = ci_hi)
          }

          for (r in show_rows) {
            raw_nm <- all_rows_int[r]
            label  <- raw_to_label[raw_nm] %||% raw_nm
            st     <- get_row_stats(r) # Uses logic from your original code
            
            tbl_int$addRow(
              rowKey = paste0("hi", r),
              values = list(
                term   = label,
                effect = if (response_type == "binary") exp(st$beta) else st$beta,
                ciLow  = if (response_type == "binary") exp(st$ci_lo) else st$ci_lo,
                ciHigh = if (response_type == "binary") exp(st$ci_hi) else st$ci_hi,
                pval   = st$pval
              )
            )
          }
        }
      }
      # ── LRT note (set regardless of fit outcome) ──────────────────
      if (int_type == "conditional_on_haplo") {
        tbl_int$setNote(
          note = paste0("This model misses terms reference haplotype:covariate; interpret with caution."),
          key  = "cond_hap")
      }
    }
  )
)
