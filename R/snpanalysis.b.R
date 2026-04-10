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

snpAnalysisClass <- R6::R6Class(
  "snpAnalysisClass",
  inherit = jmvcore::Analysis,
  private = list(

    .init = function() {
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
      if (is.null(response_var) || response_var == "") {
        self$results$validationMsg$setContent(
          "<b>Please assign a response variable.</b>")
        return()
      }
      if (length(snp_vars) == 0) {
        self$results$validationMsg$setContent(
          "<b>Please assign at least one SNP variable.</b>")
        return()
      }

      response_raw <- data[[response_var]]

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
      if (response_type == "auto") {
        n_unique <- length(unique(na.omit(response_raw)))
        if (n_unique <= 2) {
          response_type <- "binary"
        } else if (is.numeric(response_raw)) {
          response_type <- "quantitative"
        } else {
          response_type <- "binary"  # fallback for categorical with >2 levels
        }
      }

      # Prepare response
      if (response_type == "binary") {
        response <- as.integer(as.factor(response_raw)) - 1L
        response[is.na(response_raw)] <- NA_integer_
      } else {
        response <- as.numeric(response_raw)
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

      # ── Covariate descriptives ───────────────────────────────────
      if (opts$covDesc && !is.null(cov_df)) {
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
                                     snp_nm, response_raw, opts$subpop)
        }

        # Genotype frequencies
        if (opts$genoFreq) {
          private$.fill_geno_freq(item$genoFreqTable, snp_summary, ref,
                                   snp_raw, response, response_type, opts$subpop)
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
    .fill_allele_freq = function(tbl, sm, snp_nm, response_raw, subpop) {
      af <- sm$allele.freq
      for (i in seq_len(sm$nallele)) {
        tbl$addRow(rowKey = rownames(af)[i], values = list(
          allele = rownames(af)[i],
          count  = as.integer(af[i, "Count"]),
          prop   = round(af[i, "Proportion"], 4)
        ))
      }
    },

    # ── Genotype frequencies ──────────────────────────────────────
    .fill_geno_freq = function(tbl, sm, ref, snp_raw, response,
                                response_type, subpop) {
      gf <- sm$genotype.freq
      gf <- tryCatch(reorder_geno(gf, ref), error = function(e) gf)

      for (i in seq_len(nrow(gf))) {
        geno_name <- rownames(gf)[i]
        cnt  <- as.integer(gf[i, "Count"])
        prop <- if (is.na(gf[i, "Proportion"])) NA_real_
                else round(gf[i, "Proportion"], 4)

        resp_stat <- ""
        if (response_type == "quantitative" && geno_name != "NA") {
          mask <- as.character(snp_raw) == geno_name & !is.na(response)
          if (sum(mask) > 0) {
            mn  <- mean(response[mask], na.rm = TRUE)
            se  <- sd(response[mask],  na.rm = TRUE) / sqrt(sum(mask))
            resp_stat <- sprintf("%.2f (%.2f)", mn, se)
          }
        }

        row_vals <- list(
          genotype     = geno_name,
          count        = cnt,
          prop         = prop,
          responseStat = resp_stat
        )
        tbl$addRow(rowKey = geno_name, values = row_vals)
      }

      # Show responseStat column only for quantitative response
      if (response_type == "quantitative") {
        tbl$getColumn("responseStat")$setVisible(TRUE)
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
          D      = if (opts$ldD)      round(ld_res$`D`,       4) else NA,
          Dprime = if (opts$ldDprime) round(ld_res$`D'`,      4) else NA,
          r      = if (opts$ldR)      round(ld_res$`r`,       4) else NA,
          pval   = if (opts$ldPval)   ld_res$`P-value`           else NA
        )
        tbl$addRow(rowKey = paste(pair, collapse = "_"), values = row_vals)
      }
    },

    # ── Haplotype analysis ────────────────────────────────────────
    .run_haplo = function(geno_list, data, response, response_type,
                           cov_df, opts) {
      # Build allele matrix for haplo.stats (2 cols per SNP)
      snp_names <- names(geno_list)
      allele_mat <- do.call(cbind, lapply(snp_names, function(nm) {
        g <- geno_list[[nm]]
        al <- genetics::allele(g)  # returns matrix of allele1, allele2
        al
      }))

      geno_setup <- tryCatch(
        haplo.stats::setupGeno(allele_mat),
        error = function(e) NULL
      )
      if (is.null(geno_setup)) return()

      # Complete cases mask
      keep <- !is.na(response)
      if (!is.null(cov_df)) {
        keep <- keep & complete.cases(cov_df)
      }

      # Haplotype frequencies
      if (opts$haploFreq) {
        em_res <- tryCatch(
          haplo.stats::haplo.em(geno_setup[keep, ]),
          error = function(e) NULL
        )
        if (!is.null(em_res)) {
          tbl <- self$results$haploGroup$haploFreqTable
          hap_df <- as.data.frame(em_res$haplotype)
          hap_df$freq <- em_res$hap.prob

          # Label haplotypes
          hap_labels <- apply(hap_df[, snp_names, drop = FALSE], 1,
                               paste, collapse = "-")

          # Filter rare
          is_rare <- hap_df$freq < opts$haploFreqMin
          for (i in seq_len(nrow(hap_df))) {
            label <- if (is_rare[i]) "rare" else hap_labels[i]
            tbl$addRow(rowKey = as.character(i), values = list(
              haplotype = label,
              freq      = round(hap_df$freq[i], 4)
            ))
          }
        }
      }

      # Haplotype association
      if (opts$haploAssoc) {
        family <- if (response_type == "binary") "binomial" else "gaussian"
        y_sub  <- response[keep]
        x_sub  <- if (!is.null(cov_df)) cov_df[keep, , drop = FALSE] else NULL

        haplo_fit <- tryCatch({
          if (is.null(x_sub)) {
            haplo.stats::haplo.glm(
              y_sub ~ geno_setup[keep, ],
              family  = family,
              haplo.effect = "additive",
              haplo.freq.min = opts$haploFreqMin
            )
          } else {
            haplo.stats::haplo.glm(
              y_sub ~ geno_setup[keep, ] + as.matrix(x_sub),
              family  = family,
              haplo.effect = "additive",
              haplo.freq.min = opts$haploFreqMin
            )
          }
        }, error = function(e) NULL)

        if (!is.null(haplo_fit)) {
          tbl  <- self$results$haploGroup$haploAssocTable
          coef <- summary(haplo_fit)$coefficients
          ci   <- tryCatch(confint(haplo_fit), error = function(e)
                   matrix(NA, nrow = nrow(coef), ncol = 2))

          haplo_rows <- grep("^haplo", rownames(coef))
          for (i in haplo_rows) {
            nm   <- rownames(coef)[i]
            beta <- coef[i, "Estimate"]
            pv   <- coef[i, ncol(coef)]
            if (response_type == "binary") {
              eff <- exp(beta)
              ci_lo <- exp(ci[nm, 1])
              ci_hi <- exp(ci[nm, 2])
            } else {
              eff   <- beta
              ci_lo <- ci[nm, 1]
              ci_hi <- ci[nm, 2]
            }
            tbl$addRow(rowKey = nm, values = list(
              haplotype = sub("^haplo\\.", "", nm),
              freq      = NA_real_,
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
)
