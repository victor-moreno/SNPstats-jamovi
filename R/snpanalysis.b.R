#' @importFrom R6 R6Class
#' @import jmvcore
#' @import genetics
#' @import haplo.stats

# ── Helpers ────────────────────────────────────────────────────────────────────

#' Detect if a character vector looks like diploid genotypes (A/B format)
is_snp_column <- function(x) {
  vals <- unique(na.omit(as.character(x)))
  if (length(vals) == 0 || length(vals) > 3) return(FALSE)
  all(grepl("^[A-Za-z0-9]+/[A-Za-z0-9]+$", vals))
}

#' Parse a SNP column into a genetics::genotype object
#' Returns NULL if parsing fails.
parse_genotype <- function(x) {
  tryCatch(
    genetics::genotype(as.character(x), sep = "/"),
    error = function(e) NULL
  )
}

#' Determine reference genotype (most frequent homozygote)
get_ref_genotype <- function(geno) {
  if (is.null(geno)) return(NULL)
  sm <- summary(geno)
  gf <- sm$genotype.freq
  # Homozygotes: allele1 == allele2
  alleles <- rownames(gf)
  is_homoz <- sapply(alleles, function(g) {
    parts <- strsplit(g, "/")[[1]]
    length(parts) == 2 && parts[1] == parts[2]
  })
  homoz_gf <- gf[is_homoz, , drop = FALSE]
  if (nrow(homoz_gf) == 0) return(alleles[1])
  rownames(homoz_gf)[which.max(homoz_gf[, "Count"])]
}

#' Reorder genotype frequency table: ref homozygote first, then het, then alt
reorder_geno <- function(gf, ref) {
  alleles <- rownames(gf)
  na_row  <- alleles == "NA"
  other   <- alleles[!na_row]
  # ref first, then hets, then other homozygotes
  is_homoz <- sapply(other, function(g) {
    parts <- strsplit(g, "/")[[1]]
    length(parts) == 2 && parts[1] == parts[2]
  })
  ordered <- c(
    ref,
    other[!is_homoz & other != ref],
    other[is_homoz  & other != ref]
  )
  ordered <- unique(ordered[ordered %in% other])
  final <- c(ordered, alleles[na_row])
  gf[final[final %in% alleles], , drop = FALSE]
}

#' Encode SNP under a given genetic model as a numeric/factor vector
encode_model <- function(geno_char, ref, model) {
  # alleles of reference homozygote
  ref_allele <- strsplit(ref, "/")[[1]][1]
  # count ref alleles per genotype
  dosage <- sapply(geno_char, function(g) {
    if (is.na(g) || g == "NA") return(NA_integer_)
    parts <- strsplit(g, "/")[[1]]
    sum(parts == ref_allele)
  })

  switch(model,
    codominant = {
      # factor with ref as first level
      lvls <- c(ref,
                unique(geno_char[geno_char != ref & !is.na(geno_char)]))
      factor(geno_char, levels = lvls)
    },
    dominant = {
      # 0 = ref homozygote, 1 = het or alt homozygote
      ifelse(is.na(dosage), NA_integer_, as.integer(dosage < 2))
    },
    recessive = {
      # 0 = ref + het, 1 = alt homozygote
      ifelse(is.na(dosage), NA_integer_, as.integer(dosage == 0))
    },
    overdominant = {
      # 0 = homozygote (either), 1 = heterozygote
      is_het <- sapply(geno_char, function(g) {
        if (is.na(g)) return(NA)
        parts <- strsplit(g, "/")[[1]]
        length(unique(parts)) > 1
      })
      as.integer(is_het)
    },
    logadditive = {
      # 0/1/2 copies of alt allele
      2L - as.integer(dosage)
    }
  )
}

#' Fit association model for one SNP under one genetic model
#' Returns list(effect, ci_low, ci_high, pval, global_p, comparison, aic)
fit_model <- function(snp_enc, response, covariates_df, model_name,
                      response_type, ci_width) {
  alpha <- 1 - ci_width / 100

  df <- data.frame(resp = response, snp = snp_enc)
  if (!is.null(covariates_df) && ncol(covariates_df) > 0) {
    df <- cbind(df, covariates_df)
    cov_formula <- paste("+", paste(names(covariates_df), collapse = "+"))
  } else {
    cov_formula <- ""
  }

  formula_full <- as.formula(paste("resp ~ snp", cov_formula))
  formula_null <- as.formula(paste("resp ~ 1", cov_formula))

  tryCatch({
    if (response_type == "binary") {
      fit_full <- glm(formula_full, data = df, family = binomial())
      fit_null <- glm(formula_null, data = df, family = binomial())
    } else {
      fit_full <- lm(formula_full, data = df)
      fit_null <- lm(formula_null, data = df)
    }

    # Global LRT p-value
    lrt <- tryCatch(
      anova(fit_null, fit_full, test = "LRT"),
      error = function(e) NULL
    )
    global_p <- if (!is.null(lrt)) lrt[2, "Pr(>Chi)"] else NA_real_

    coefs <- summary(fit_full)$coefficients
    snp_rows <- grep("^snp", rownames(coefs))

    if (length(snp_rows) == 0) {
      return(NULL)
    }

    ci <- tryCatch(
      confint(fit_full, level = ci_width / 100)[snp_rows, , drop = FALSE],
      error = function(e) matrix(NA, nrow = length(snp_rows), ncol = 2)
    )

    results <- lapply(seq_along(snp_rows), function(i) {
      row <- snp_rows[i]
      beta <- coefs[row, "Estimate"]
      pval <- coefs[row, ifelse(response_type == "binary",
                                 "Pr(>|z|)", "Pr(>|t|)")]
      ci_lo <- ci[i, 1]
      ci_hi <- ci[i, 2]

      if (response_type == "binary") {
        list(
          effect   = exp(beta),
          ci_low   = exp(ci_lo),
          ci_high  = exp(ci_hi),
          pval     = pval,
          global_p = global_p,
          comparison = sub("^snp", "", rownames(coefs)[row])
        )
      } else {
        list(
          effect   = beta,
          ci_low   = ci_lo,
          ci_high  = ci_hi,
          pval     = pval,
          global_p = global_p,
          comparison = sub("^snp", "", rownames(coefs)[row])
        )
      }
    })
    results
  }, error = function(e) NULL)
}


# ── Main analysis class ────────────────────────────────────────────────────────

snpAnalysisClass <- if (requireNamespace("jmvcore", quietly=TRUE)) R6::R6Class(
  "snpAnalysisClass",
  inherit = snpAnalysisBase,
  private = list(

    .init = function() {
      # Hide all optional groups immediately — shown only when options are on
      self$results$covDescGroup$setVisible(FALSE)
      self$results$ldGroup$setVisible(FALSE)
      self$results$haploGroup$setVisible(FALSE)

      # Initialise the per-SNP array when SNPs are assigned
      snp_names <- self$options$snps
      if (length(snp_names) == 0) return()

      arr <- self$results$snpResults
      for (nm in snp_names) {
        arr$addItem(key = nm)
      }
    },

    .run = function() {
      data      <- self$data
      opts      <- self$options
      response_var   <- opts$response
      snp_vars       <- opts$snps
      covariate_vars <- opts$covariates

      # ── Validation ──────────────────────────────────────────────
      self$results$validationMsg$setContent("")

      # Nothing assigned yet — keep results panel empty and silent
      if (length(snp_vars) == 0) return()

      needs_response <- opts$snpAssoc || opts$haploAssoc || opts$subpop || opts$covDesc
      if (needs_response && (is.null(response_var) || response_var == "")) {
        self$results$validationMsg$setContent(paste0(
          "<b>Please assign a response variable</b> (required for association, ",
          "stratification, and covariate descriptives)."))
        return()
      }

      response_raw <- if (!is.null(response_var) && response_var != "") data[[response_var]] else NULL

      # Check SNP columns
      bad_snps <- character(0)
      for (v in snp_vars) {
        if (!is_snp_column(data[[v]])) bad_snps <- c(bad_snps, v)
      }
      if (length(bad_snps) > 0) {
        self$results$validationMsg$setContent(paste0(
          "<b>Warning:</b> The following columns do not appear to contain ",
          "diploid genotypes (X/Y format): ",
          paste(bad_snps, collapse = ", "),
          ". They will be skipped."))
        snp_vars <- setdiff(snp_vars, bad_snps)
      }

      # ── Determine response type ──────────────────────────────────
      response_type <- opts$responseType
      if (!is.null(response_raw) && response_type == "auto") {
        n_unique <- length(unique(na.omit(response_raw)))
        if (n_unique <= 2) {
          response_type <- "binary"
        } else if (is.numeric(response_raw)) {
          response_type <- "quantitative"
        } else {
          response_type <- "binary"
        }
      } else if (is.null(response_raw)) {
        response_type <- "binary"  # default, won't be used
      }

      # Prepare response
      if (!is.null(response_raw)) {
        if (response_type == "binary") {
          response <- as.integer(as.factor(response_raw)) - 1L
          response[is.na(response_raw)] <- NA_integer_
        } else {
          response <- as.numeric(response_raw)
        }
      } else {
        response <- NULL
      }

      # ── Covariates ───────────────────────────────────────────────
      if (length(covariate_vars) > 0) {
        cov_df <- data[, covariate_vars, drop = FALSE]
        # Convert factors appropriately
        for (v in covariate_vars) {
          if (!is.numeric(cov_df[[v]])) {
            cov_df[[v]] <- as.factor(cov_df[[v]])
          }
        }
      } else {
        cov_df <- NULL
      }

      # ── Show/hide optional result groups ────────────────────────
      show_cov_desc <- isTRUE(opts$covDesc) && length(covariate_vars) > 0
      self$results$covDescGroup$setVisible(show_cov_desc)
      self$results$ldGroup$setVisible(isTRUE(opts$ldAnalysis))
      self$results$haploGroup$setVisible(
        isTRUE(opts$haploFreq) || isTRUE(opts$haploAssoc))

      # ── Covariate descriptives ───────────────────────────────────
      if (show_cov_desc && !is.null(cov_df)) {
        private$.run_cov_desc(cov_df, response_raw, response_type)
      }

      # ── Per-SNP analyses ─────────────────────────────────────────
      arr <- self$results$snpResults
      geno_list <- list()  # for LD/haplotype later

      for (snp_nm in snp_vars) {
        snp_raw  <- data[[snp_nm]]
        geno_obj <- parse_genotype(snp_raw)
        if (is.null(geno_obj)) next

        geno_list[[snp_nm]] <- geno_obj
        item <- arr$get(key = snp_nm)
        snp_summary <- summary(geno_obj)
        n_typed  <- snp_summary$n.typed
        n_total  <- snp_summary$n.total
        pct      <- round(n_typed / n_total * 100, 1)

        item$typingRate$setContent(sprintf(
          "<b>Typed samples:</b> %d / %d (%.1f%%)", n_typed, n_total, pct))

        ref <- get_ref_genotype(geno_obj)

        # Allele frequencies
        if (opts$allFreq) {
          private$.fill_allele_freq(item$allFreqTable, snp_summary,
                                     snp_nm, response_raw, opts$subpop,
                                     response_type, snp_raw)
        }

        # Genotype frequencies
        if (opts$genoFreq) {
          private$.fill_geno_freq(item$genoFreqTable, snp_summary, ref,
                                   snp_raw, response, response_type, opts$subpop,
                                   response_raw)
        }

        # HWE
        if (opts$hweTest) {
          private$.fill_hwe(item$hweTable, geno_obj, snp_nm,
                             response_raw, opts$subpop)
        }

        # Association
        if (opts$snpAssoc) {
          private$.fill_assoc(item$assocTable, snp_raw, ref, response,
                               cov_df, response_type, opts)
        }
      }

      # ── LD analysis ──────────────────────────────────────────────
      if (opts$ldAnalysis && length(geno_list) >= 2) {
        private$.run_ld(geno_list, opts)
      }

      # ── Haplotype analysis ───────────────────────────────────────
      if ((opts$haploFreq || opts$haploAssoc) && length(geno_list) >= 2) {
        private$.run_haplo(geno_list, data, response, response_type,
                            cov_df, opts)
      }
    },

    # ── Covariate descriptives ────────────────────────────────────
    .run_cov_desc = function(cov_df, response_raw, response_type) {
      tbl <- self$results$covDescGroup$covDescTable
      for (v in names(cov_df)) {
        col <- cov_df[[v]]
        if (is.factor(col) || is.character(col)) {
          col <- as.factor(col)
          lvl_counts <- table(col, useNA = "no")
          for (lvl in names(lvl_counts)) {
            tbl$addRow(rowKey = paste0(v, "_", lvl), values = list(
              variable = v,
              level    = lvl,
              n        = as.integer(lvl_counts[lvl]),
              stat     = paste0(round(lvl_counts[lvl] / sum(lvl_counts) * 100, 1), "%")
            ))
          }
        } else {
          mn  <- mean(col, na.rm = TRUE)
          sdv <- sd(col,   na.rm = TRUE)
          tbl$addRow(rowKey = v, values = list(
            variable = v,
            level    = "",
            n        = sum(!is.na(col)),
            stat     = sprintf("%.2f ± %.2f", mn, sdv)
          ))
        }
      }
    },

    # ── Allele frequencies ────────────────────────────────────────
    .fill_allele_freq = function(tbl, sm, snp_nm, response_raw, subpop,
                                  response_type, snp_raw) {
      af <- sm$allele.freq

      # Stratify by binary response only
      do_strat <- isTRUE(subpop) && !is.null(response_raw) &&
                  identical(response_type, "binary")
      grp_levels <- character(0)
      if (do_strat) {
        grp_levels <- sort(unique(na.omit(as.character(response_raw))))
        if (length(grp_levels) < 2L || length(grp_levels) > 5L)
          do_strat <- FALSE
      }

      if (do_strat) {
        for (g in grp_levels) {
          safe <- gsub("[^A-Za-z0-9]", "_", g)
          tbl$addColumn(name   = paste0("count_", safe),
                        title  = paste0("n (", g, ")"),
                        type   = "integer")
          tbl$addColumn(name   = paste0("prop_", safe),
                        title  = paste0("Prop (", g, ")"),
                        type   = "number",
                        format = "zto")
        }
      }

      allele_names <- rownames(af)
      for (i in seq_len(sm$nallele)) {
        al <- allele_names[i]
        row_vals <- list(
          allele = al,
          count  = as.integer(af[i, "Count"]),
          prop   = round(af[i, "Proportion"], 4)
        )
        if (do_strat) {
          for (g in grp_levels) {
            safe   <- gsub("[^A-Za-z0-9]", "_", g)
            mask_g <- !is.na(response_raw) & as.character(response_raw) == g
            # split every diploid genotype into its two alleles
            all_alleles <- unlist(strsplit(
              as.character(snp_raw)[mask_g & !is.na(snp_raw)], "/"))
            n_al  <- sum(all_alleles == al, na.rm = TRUE)
            n_tot <- length(all_alleles)
            row_vals[[paste0("count_", safe)]] <- as.integer(n_al)
            row_vals[[paste0("prop_",  safe)]] <-
              if (n_tot > 0L) round(n_al / n_tot, 4) else NA_real_
          }
        }
        tbl$addRow(rowKey = al, values = row_vals)
      }
    },

    # ── Genotype frequencies ──────────────────────────────────────
    .fill_geno_freq = function(tbl, sm, ref, snp_raw, response,
                                response_type, subpop, response_raw) {
      gf <- sm$genotype.freq
      gf <- tryCatch(reorder_geno(gf, ref), error = function(e) gf)

      # Stratify by binary response only
      do_strat <- isTRUE(subpop) && !is.null(response_raw) &&
                  identical(response_type, "binary")
      grp_levels <- character(0)
      if (do_strat) {
        grp_levels <- sort(unique(na.omit(as.character(response_raw))))
        if (length(grp_levels) < 2L || length(grp_levels) > 5L)
          do_strat <- FALSE
      }

      # Add dynamic stratification columns before any rows
      if (do_strat) {
        for (g in grp_levels) {
          safe <- gsub("[^A-Za-z0-9]", "_", g)
          tbl$addColumn(name   = paste0("count_", safe),
                        title  = paste0("n (", g, ")"),
                        type   = "integer")
          tbl$addColumn(name   = paste0("prop_", safe),
                        title  = paste0("Prop (", g, ")"),
                        type   = "number",
                        format = "zto")
        }
      }

      # Quantitative: make mean±SE column visible
      if (response_type == "quantitative") {
        tbl$getColumn("responseStat")$setVisible(TRUE)
      }

      snp_char <- as.character(snp_raw)
      for (i in seq_len(nrow(gf))) {
        geno_name <- rownames(gf)[i]
        if (geno_name == "NA") next          # skip missing-genotype row

        cnt  <- as.integer(gf[i, "Count"])
        prop <- if (is.na(gf[i, "Proportion"])) NA_real_
                else round(gf[i, "Proportion"], 4)

        resp_stat <- ""
        if (response_type == "quantitative") {
          mask <- snp_char == geno_name & !is.na(response)
          if (sum(mask) > 0) {
            mn <- mean(response[mask], na.rm = TRUE)
            se <- sd(response[mask],   na.rm = TRUE) / sqrt(sum(mask))
            resp_stat <- sprintf("%.2f (%.2f)", mn, se)
          }
        }

        row_vals <- list(
          genotype     = geno_name,
          count        = cnt,
          prop         = prop,
          responseStat = resp_stat
        )

        if (do_strat) {
          for (g in grp_levels) {
            safe  <- gsub("[^A-Za-z0-9]", "_", g)
            mask_g <- !is.na(response_raw) & as.character(response_raw) == g
            n_g    <- sum(mask_g & snp_char == geno_name, na.rm = TRUE)
            n_tot  <- sum(mask_g & !is.na(snp_raw) & snp_char != "NA",
                          na.rm = TRUE)
            row_vals[[paste0("count_", safe)]] <- as.integer(n_g)
            row_vals[[paste0("prop_",  safe)]] <-
              if (n_tot > 0L) round(n_g / n_tot, 4) else NA_real_
          }
        }

        tbl$addRow(rowKey = geno_name, values = row_vals)
      }
    },

    # ── Hardy-Weinberg test ───────────────────────────────────────
    .fill_hwe = function(tbl, geno_obj, snp_nm, response_raw, subpop) {
      # Overall
      hw <- tryCatch(genetics::HWE.exact(geno_obj), error = function(e) NULL)
      if (is.null(hw)) return()

      st <- hw$statistic
      pr <- hw$parameter
      tbl$addRow(rowKey = "All", values = list(
        group = "All subjects",
        n11   = as.integer(st["N11"]),
        n12   = as.integer(st["N12"]),
        n22   = as.integer(st["N22"]),
        n1    = as.integer(pr["N1"]),
        n2    = as.integer(pr["N2"]),
        pval  = hw$p.value
      ))

      # Stratified by response
      if (subpop && !is.null(response_raw)) {
        lvls <- unique(na.omit(as.character(response_raw)))
        if (length(lvls) <= 5) {
          for (lvl in lvls) {
            mask <- !is.na(response_raw) & as.character(response_raw) == lvl
            hw_sub <- tryCatch(
              genetics::HWE.exact(geno_obj[mask]),
              error = function(e) NULL)
            if (is.null(hw_sub)) next
            st2 <- hw_sub$statistic
            pr2 <- hw_sub$parameter
            tbl$addRow(rowKey = lvl, values = list(
              group = paste0("Response = ", lvl),
              n11   = as.integer(st2["N11"]),
              n12   = as.integer(st2["N12"]),
              n22   = as.integer(st2["N22"]),
              n1    = as.integer(pr2["N1"]),
              n2    = as.integer(pr2["N2"]),
              pval  = hw_sub$p.value
            ))
          }
        }
      }
    },

    # ── SNP association ───────────────────────────────────────────
    .fill_assoc = function(tbl, snp_raw, ref, response, cov_df,
                            response_type, opts) {
      models <- c()
      if (opts$modelCodominant)   models <- c(models, "codominant")
      if (opts$modelDominant)     models <- c(models, "dominant")
      if (opts$modelRecessive)    models <- c(models, "recessive")
      if (opts$modelOverdominant) models <- c(models, "overdominant")
      if (opts$modelLogAdditive)  models <- c(models, "logadditive")

      model_labels <- c(
        codominant   = "Codominant",
        dominant     = "Dominant",
        recessive    = "Recessive",
        overdominant = "Overdominant",
        logadditive  = "Log-additive"
      )

      row_key <- 0L
      for (mdl in models) {
        snp_enc <- encode_model(as.character(snp_raw), ref, mdl)
        res_list <- fit_model(snp_enc, response, cov_df, mdl,
                               response_type, opts$ciWidth)
        if (is.null(res_list)) next

        first_row <- TRUE
        for (res in res_list) {
          row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            model      = if (first_row) model_labels[mdl] else "",
            comparison = res$comparison,
            effect     = res$effect,
            ciLow      = res$ci_low,
            ciHigh     = res$ci_high,
            pval       = res$pval,
            globalP    = if (first_row) res$global_p else NA_real_
          ))
          first_row <- FALSE
        }
      }
    },

    # ── Linkage disequilibrium ────────────────────────────────────
    .run_ld = function(geno_list, opts) {
      tbl   <- self$results$ldGroup$ldTable
      nms   <- names(geno_list)
      pairs <- combn(nms, 2, simplify = FALSE)

      for (pair in pairs) {
        g1 <- geno_list[[pair[1]]]
        g2 <- geno_list[[pair[2]]]
        ld_res <- tryCatch(genetics::LD(g1, g2), error = function(e) NULL)
        if (is.null(ld_res)) next

        row_vals <- list(
          snp1   = pair[1],
          snp2   = pair[2],
          D      = round(ld_res$`D`,    4),
          Dprime = round(ld_res$`D'`,   4),
          r      = round(ld_res$`r`,    4),
          pval   = ld_res$`P-value`
        )
        tbl$addRow(rowKey = paste(pair, collapse = "_"), values = row_vals)
      }
    },

    # ── Haplotype analysis ────────────────────────────────────────
    .run_haplo = function(geno_list, data, response, response_type,
                           cov_df, opts) {
      # Build allele matrix for haplo.stats (2 cols per SNP)
      snp_names <- names(geno_list)
      allele_list <- lapply(snp_names, function(nm) {
        g <- geno_list[[nm]]
        genetics::allele(g)  # matrix of allele1, allele2
      })
      allele_mat <- do.call(cbind, allele_list)

      # Build a lookup: for each SNP, what are the actual allele labels per code?
      # haplo.stats encodes alleles as integers; we recover labels from allele matrix
      get_allele_labels <- function(allele_cols) {
        vals <- unique(na.omit(c(allele_cols)))
        sort(vals)
      }
      snp_allele_labels <- lapply(allele_list, function(al) get_allele_labels(al))

      geno_setup <- tryCatch(
        haplo.stats::setupGeno(allele_mat),
        error = function(e) NULL
      )
      if (is.null(geno_setup)) return()

      # Get the allele label lookup from setupGeno attributes
      # attr(geno_setup, "unique.alleles") is a list per locus
      unique_alleles <- attr(geno_setup, "unique.alleles")

      # Helper: decode numeric haplotype row to allele string
      decode_haplotype <- function(hap_row, unique_alleles) {
        parts <- character(length(hap_row))
        for (i in seq_along(hap_row)) {
          code <- hap_row[i]
          locus_alleles <- unique_alleles[[i]]
          if (!is.na(code) && code >= 1 && code <= length(locus_alleles)) {
            parts[i] <- locus_alleles[code]
          } else {
            parts[i] <- as.character(code)
          }
        }
        paste(parts, collapse = "-")
      }

      # Complete cases mask
      keep <- if (!is.null(response)) !is.na(response) else rep(TRUE, nrow(allele_mat))
      if (!is.null(cov_df)) {
        keep <- keep & complete.cases(cov_df)
      }

      # Haplotype frequencies
      em_res <- NULL
      hap_labels_full <- NULL  # labels for all haplotypes above threshold

      if (opts$haploFreq || opts$haploAssoc) {
        em_res <- tryCatch(
          haplo.stats::haplo.em(geno_setup[keep, ]),
          error = function(e) NULL
        )
        if (!is.null(em_res)) {
          hap_df <- as.data.frame(em_res$haplotype)
          hap_df$freq <- em_res$hap.prob

          # Decode haplotype labels using actual allele names
          if (!is.null(unique_alleles)) {
            hap_labels_full <- apply(hap_df[, seq_len(ncol(hap_df) - 1), drop = FALSE],
                                     1, decode_haplotype, unique_alleles = unique_alleles)
          } else {
            # Fallback: use allele labels from genetics object
            hap_labels_full <- apply(hap_df[, seq_len(ncol(hap_df) - 1), drop = FALSE],
                                     1, paste, collapse = "-")
          }

          if (opts$haploFreq) {
            tbl <- self$results$haploGroup$haploFreqTable
            is_rare <- hap_df$freq < opts$haploFreqMin
            rare_freq_sum <- sum(hap_df$freq[is_rare], na.rm = TRUE)

            for (i in seq_len(nrow(hap_df))) {
              if (is_rare[i]) next
              tbl$addRow(rowKey = paste0("hap_", i), values = list(
                haplotype = hap_labels_full[i],
                freq      = round(hap_df$freq[i], 4)
              ))
            }
            if (rare_freq_sum > 0) {
              tbl$addRow(rowKey = "rare_combined", values = list(
                haplotype = "Rare (combined)",
                freq      = round(rare_freq_sum, 4)
              ))
            }
          }
        }
      }

      # Haplotype association
      if (opts$haploAssoc && !is.null(response)) {
        family <- if (response_type == "binary") "binomial" else "gaussian"
        y_sub  <- response[keep]
        x_sub  <- if (!is.null(cov_df)) as.data.frame(cov_df[keep, , drop = FALSE]) else NULL

        # haplo.glm needs the geno matrix subset and y as numeric
        geno_sub <- geno_setup[keep, , drop = FALSE]

        haplo_fit <- tryCatch({
          if (is.null(x_sub) || ncol(x_sub) == 0) {
            haplo.stats::haplo.glm(
              y_sub ~ geno_sub,
              family         = family,
              haplo.effect   = "additive",
              haplo.freq.min = opts$haploFreqMin,
              data           = data.frame(y_sub = y_sub)
            )
          } else {
            cov_names <- names(x_sub)
            fit_data  <- cbind(data.frame(y_sub = y_sub), x_sub)
            cov_formula <- paste(cov_names, collapse = " + ")
            full_formula <- as.formula(paste("y_sub ~ geno_sub +", cov_formula))
            haplo.stats::haplo.glm(
              full_formula,
              family         = family,
              haplo.effect   = "additive",
              haplo.freq.min = opts$haploFreqMin,
              data           = fit_data
            )
          }
        }, error = function(e) {
          # Capture error for user feedback
          self$results$validationMsg$setContent(paste0(
            "<b>Haplotype association error:</b> ", conditionMessage(e)))
          NULL
        })

        if (!is.null(haplo_fit)) {
          tbl  <- self$results$haploGroup$haploAssocTable
          coef_mat <- summary(haplo_fit)$coefficients
          ci   <- tryCatch(
            confint(haplo_fit, level = opts$ciWidth / 100),
            error = function(e) matrix(NA_real_, nrow = nrow(coef_mat), ncol = 2,
                                       dimnames = list(rownames(coef_mat), c("2.5 %", "97.5 %")))
          )

          haplo_rows <- grep("^haplo\\.", rownames(coef_mat))

          # Build haplotype code -> label map from em_res if available
          haplo_label_map <- list()
          if (!is.null(em_res) && !is.null(hap_labels_full)) {
            hap_df <- as.data.frame(em_res$haplotype)
            hap_df$freq <- em_res$hap.prob
            for (i in seq_len(nrow(hap_df))) {
              # haplo.glm uses "haplo.X.X.X" notation matching haplotype codes
              code_str <- paste(as.integer(hap_df[i, seq_len(ncol(hap_df) - 1)]), collapse = ".")
              key <- paste0("haplo.", code_str)
              haplo_label_map[[key]] <- hap_labels_full[i]
            }
          }

          for (i in haplo_rows) {
            row_nm <- rownames(coef_mat)[i]
            beta   <- coef_mat[i, "Estimate"]
            pv     <- coef_mat[i, ncol(coef_mat)]

            ci_lo <- if (row_nm %in% rownames(ci)) ci[row_nm, 1] else NA_real_
            ci_hi <- if (row_nm %in% rownames(ci)) ci[row_nm, 2] else NA_real_

            if (response_type == "binary") {
              eff   <- exp(beta)
              ci_lo <- if (!is.na(ci_lo)) exp(ci_lo) else NA_real_
              ci_hi <- if (!is.na(ci_hi)) exp(ci_hi) else NA_real_
            } else {
              eff <- beta
            }

            # Decode label
            label <- if (!is.null(haplo_label_map[[row_nm]])) {
              haplo_label_map[[row_nm]]
            } else {
              sub("^haplo\\.", "", row_nm)
            }

            # Get frequency for this haplotype
            hap_freq <- NA_real_
            if (!is.null(em_res) && !is.null(haplo_label_map[[row_nm]])) {
              idx <- which(hap_labels_full == haplo_label_map[[row_nm]])
              if (length(idx) > 0) hap_freq <- em_res$hap.prob[idx[1]]
            }

            tbl$addRow(rowKey = row_nm, values = list(
              haplotype = label,
              freq      = hap_freq,
              effect    = eff,
              ciLow     = ci_lo,
              ciHigh    = ci_hi,
              pval      = pv
            ))
          }
        }
      }
    }
  )
)  # end R6Class
