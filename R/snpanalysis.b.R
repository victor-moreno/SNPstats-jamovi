#' @importFrom R6 R6Class
#' @import jmvcore
#' @import genetics
#' @import haplo.stats
#' @import ggplot2

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

  # remove missings
  df <- df[complete.cases(df), , drop = FALSE]

  formula_full <- as.formula(paste("resp ~ snp", cov_formula))
  formula_null <- as.formula(paste("resp ~ 1", cov_formula))

  tryCatch({
    if (response_type == "binary") {
      fit_full <- glm(formula_full, data = df, family = binomial())
      fit_null <- glm(formula_null, data = df, family = binomial())
      lrtest <-'Chisq'
      lrtest_label <- 'Pr(>Chi)'
    } else {
      fit_full <- lm(formula_full, data = df)
      fit_null <- lm(formula_null, data = df)
      lrtest <-'F'
      lrtest_label <- 'Pr(>F)'
    }

    # Global LRT p-value
    lrt <- tryCatch(
      anova(fit_null, fit_full, test = lrtest),
      error = function(e) NULL
    )
    global_p <- if (!is.null(lrt)) lrt[2, lrtest_label] else NA_real_
    
    aic_val <- AIC(fit_full)

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
          aic      = aic_val,
          comparison = sub("^snp", "", rownames(coefs)[row])
        )
      } else {
        list(
          effect   = beta,
          ci_low   = ci_lo,
          ci_high  = ci_hi,
          pval     = pval,
          global_p = global_p,
          aic      = aic_val,
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
      self$results$snpSummaryTablesGroup$setVisible(isTRUE(self$options$snpSummary))
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

      needs_response <- opts$snpAssoc || opts$haploAssoc || opts$subpop
      if (needs_response && (is.null(response_var) || response_var == "")) {
        self$results$validationMsg$setContent(paste0(
          "<b>Please assign a response variable</b> (required for association or ",
          "stratification)."))
        return()
      } else{
          self$results$validationMsg$setContent(NULL)
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
      self$results$ldGroup$setVisible(isTRUE(opts$ldAnalysis) || isTRUE(opts$ldMatrix) || isTRUE(opts$ldPlot))
      self$results$haploGroup$setVisible(
        isTRUE(opts$haploFreq) || isTRUE(opts$haploAssoc))

      # ── Covariate descriptives ───────────────────────────────────
      if (show_cov_desc && !is.null(cov_df)) {
        private$.run_cov_desc(cov_df, response_raw, response_type)
      }

      # ── SNP summary table ─────────────────────────────────────────
      if (isTRUE(opts$snpSummary)) {
        private$.fill_snp_summary(data, snp_vars, response_raw, response_type, opts$subpop)
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
      needs_ld <- (isTRUE(opts$ldAnalysis) || isTRUE(opts$ldMatrix) || isTRUE(opts$ldPlot))
      if (needs_ld && length(geno_list) >= 2) {
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
          first_row = TRUE
          for (lvl in names(lvl_counts)) {
            tbl$addRow(rowKey = paste0(v, "_", lvl), values = list(
              variable = ifelse(first_row,v,''),
              level    = lvl,
              n        = as.integer(lvl_counts[lvl]),
              stat     = paste0(round(lvl_counts[lvl] / sum(lvl_counts) * 100, 1), "%")
            ))
            first_row = FALSE
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

    # ── SNP summary table ─────────────────────────────────────────
    .fill_snp_summary = function(data, snp_vars, response_raw, response_type, subpop) {
      tbl <- self$results$snpSummaryTablesGroup$snpSummaryTable

      do_strat <- isTRUE(subpop) &&
                  !is.null(response_raw) &&
                  identical(response_type, "binary")

      grp_levels <- character(0)
      if (do_strat) {
        grp_levels <- sort(unique(na.omit(as.character(response_raw))))
        if (length(grp_levels) != 2L) do_strat <- FALSE
      }

      row_key <- 0L

      for (snp_nm in snp_vars) {
        snp_raw  <- data[[snp_nm]]
        geno_obj <- parse_genotype(snp_raw)
        if (is.null(geno_obj)) next

        sm  <- summary(geno_obj)
        ref <- get_ref_genotype(geno_obj)

        # Derive allele labels: A = major (ref allele), B = minor
        af_all       <- sm$allele.freq
        allele_nms   <- rownames(af_all)
        ref_allele   <- strsplit(ref, "/")[[1]][1]
        alt_allele   <- allele_nms[allele_nms != ref_allele]
        alt_allele   <- if (length(alt_allele) > 0) alt_allele[1] else "?"
        alleles_label <- paste0(ref_allele, "/", alt_allele)  # e.g. "C/T"

        # Helper: compute summary stats for a subset mask (NULL = all)
        compute_row <- function(mask) {
          snp_sub  <- if (is.null(mask)) snp_raw  else snp_raw[mask]
          geno_sub <- if (is.null(mask)) geno_obj else parse_genotype(snp_sub)
          if (is.null(geno_sub)) return(NULL)

          sm_sub <- summary(geno_sub)

          # N typed (non-missing)
          n_typed <- sm_sub$n.typed

          # Allele frequencies → MAF defined as freq of B (alt/minor) allele
          af <- sm_sub$allele.freq
          props <- af[, "Proportion"]
          # B allele is alt_allele (defined in outer scope via <<- / parent env)
          maf <- if (alt_allele %in% rownames(af)) {
            af[alt_allele, "Proportion"]
          } else if (length(props) >= 2) {
            min(props, na.rm = TRUE)
          } else NA_real_

          # Genotype counts in ref/het/alt order
          gf <- sm_sub$genotype.freq
          gf <- tryCatch(reorder_geno(gf, ref), error = function(e) gf)
          # Remove NA row if present
          gf <- gf[rownames(gf) != "NA", , drop = FALSE]
          counts <- as.integer(gf[, "Count"])
          geno_str <- if (length(counts) == 3) {
            paste(counts, collapse = " / ")
          } else if (length(counts) == 2) {
            paste(c(counts, 0L), collapse = " / ")
          } else {
            paste(counts, collapse = " / ")
          }

          # HWE exact test
          hwe_p <- tryCatch({
            hw <- genetics::HWE.exact(geno_sub)
            hw$p.value
          }, error = function(e) NA_real_)

          list(n = n_typed, maf = maf, genoCounts = geno_str, hwePval = hwe_p)
        }

        # Overall row — show alleles label here only
        res_all <- compute_row(NULL)
        if (!is.null(res_all)) {
          row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            snp        = snp_nm,
            alleles    = alleles_label,
            group      = if (do_strat) "All" else "",
            n          = res_all$n,
            maf        = round(res_all$maf, 4),
            genoCounts = res_all$genoCounts,
            hwePval    = res_all$hwePval
          ))
        }

        # Stratified rows — alleles blank (already shown on All row)
        if (do_strat) {
          for (lvl in grp_levels) {
            mask_lvl <- !is.na(response_raw) & as.character(response_raw) == lvl
            res_lvl  <- compute_row(mask_lvl)
            if (is.null(res_lvl)) next
            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              snp        = "",
              alleles    = "",
              group      = as.character(lvl),
              n          = res_lvl$n,
              maf        = round(res_lvl$maf, 4),
              genoCounts = res_lvl$genoCounts,
              hwePval    = res_lvl$hwePval
            ))
          }
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
                        title  = paste0("N (", g, ")"),
                        type   = "integer")
          tbl$addColumn(name   = paste0("prop_", safe),
                        title  = paste0("% (", g, ")"),
                        type   = "number",
                        format = "dp=1")
        }
      }

      allele_names <- rownames(af)
      for (i in seq_len(sm$nallele)) {
        al <- allele_names[i]
        row_vals <- list(
          allele = al,
          count  = as.integer(af[i, "Count"]),
          prop   = round(af[i, "Proportion"]*100, 1)
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
              if (n_tot > 0L) round(n_al / n_tot * 100, 1) else NA_real_
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

      if (do_strat) {
        for (g in grp_levels) {
          safe <- gsub("[^A-Za-z0-9]", "_", g)
          tbl$addColumn(name   = paste0("count_", safe),
                        title  = paste0("N (", g, ")"),
                        type   = "integer")
          tbl$addColumn(name   = paste0("prop_", safe),
                        title  = paste0("% (", g, ")"),
                        type   = "number",
                        format = "dp=1")
        }
      }

      if (response_type == "quantitative") {
        tbl$getColumn("responseStat")$setVisible(TRUE)
      }

      snp_char <- as.character(snp_raw)
      for (i in seq_len(nrow(gf))) {
        geno_name <- rownames(gf)[i]
        if (geno_name == "NA") next

        cnt  <- as.integer(gf[i, "Count"])
        prop <- if (is.na(gf[i, "Proportion"])) NA_real_
                else round(gf[i, "Proportion"]*100, 1)

        resp_stat <- ""
        if (response_type == "quantitative") {
          mask <- (snp_char == geno_name) & !is.na(response)
          # FIX: Added na.rm = TRUE to prevent crash on missing SNP values
          if (sum(mask, na.rm = TRUE) > 0) {
            mn <- mean(response[mask], na.rm = TRUE)
            se <- sd(response[mask],   na.rm = TRUE) / sqrt(sum(mask, na.rm = TRUE))
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
            safe   <- gsub("[^A-Za-z0-9]", "_", g)
            mask_g <- !is.na(response_raw) & as.character(response_raw) == g
            n_g    <- sum(mask_g & snp_char == geno_name, na.rm = TRUE)
            n_tot  <- sum(mask_g & !is.na(snp_raw) & snp_char != "NA", na.rm = TRUE)
            row_vals[[paste0("count_", safe)]] <- as.integer(n_g)
            row_vals[[paste0("prop_",  safe)]] <-
              if (n_tot > 0L) round(n_g / n_tot * 100, 1) else NA_real_
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
              pval  = hw_sub$p.value
            ))
          }
        }
      }
    },

    # ── Linkage disequilibrium ────────────────────────────────────
    .run_ld = function(geno_list, opts) {
      nms   <- names(geno_list)
      n     <- length(nms)
      pairs <- combn(nms, 2, simplify = FALSE)

      # Compute all pairwise LD results once
      ld_store <- list()   # keyed by "snp1___snp2"
      for (pair in pairs) {
        key    <- paste(pair, collapse = "___")
        ld_res <- tryCatch(genetics::LD(geno_list[[pair[1]]], geno_list[[pair[2]]]),
                           error = function(e) NULL)
        if (!is.null(ld_res)) ld_store[[key]] <- ld_res
      }

      # ── Pairwise table ─────────────────────────────────────────
      if (isTRUE(opts$ldAnalysis)) {
        tbl <- self$results$ldGroup$ldTable
        for (pair in pairs) {
          key    <- paste(pair, collapse = "___")
          ld_res <- ld_store[[key]]
          if (is.null(ld_res)) next
          tbl$addRow(rowKey = paste(pair, collapse = "_"), values = list(
            snp1   = pair[1],
            snp2   = pair[2],
            r2     = round(ld_res$`r`^2,  3),
            Dprime = round(ld_res$`D'`,   3),
            D      = round(ld_res$`D`,    3),
            pval   = ld_res$`P-value`
          ))
        }
      }

      # ── LD matrix ─────────────────────────────────────────────
      if (isTRUE(opts$ldMatrix)) {
        mtbl   <- self$results$ldGroup$ldMatrixTable
        metric <- opts$ldMetric   # "r2", "Dprime", or "D"

        # Add one column per SNP (beyond the row-label column already defined)
        for (nm in nms) {
          safe_nm <- gsub("[^A-Za-z0-9_]", "_", nm)
          mtbl$addColumn(name = safe_nm, title = nm, type = "text")
        }

        # Build n×n value matrices
        upper_mat <- matrix("", n, n, dimnames = list(nms, nms))
        lower_mat <- matrix("", n, n, dimnames = list(nms, nms))
        diag(upper_mat) <- nms
        diag(lower_mat) <- nms

        for (pair in pairs) {
          key    <- paste(pair, collapse = "___")
          ld_res <- ld_store[[key]]
          if (is.null(ld_res)) next

          p_val <- ld_res$`P-value`
          p_str <- if (!is.na(p_val)) {
            if (p_val < 0.001) "< .001" else sprintf("%.3f", p_val)
          } else ""

          up_val <- switch(metric,
            Dprime = sprintf("%.3f", round(ld_res$`D'`, 3)),
            r2      = sprintf("%.3f", round(ld_res$`r`^2,   3)),
            D      = sprintf("%.3f", round(ld_res$`D`,   3))
          )

          upper_mat[pair[1], pair[2]] <- up_val
          upper_mat[pair[2], pair[1]] <- ""          # will be lower
          lower_mat[pair[1], pair[2]] <- ""
          lower_mat[pair[2], pair[1]] <- p_str
        }

        for (i in seq_len(n)) {
          row_vals <- list(snp = nms[i])
          for (j in seq_len(n)) {
            safe_nm <- gsub("[^A-Za-z0-9_]", "_", nms[j])
            if (i == j) {
              row_vals[[safe_nm]] <- nms[i]
            } else if (j > i) {
              row_vals[[safe_nm]] <- upper_mat[i, j]
            } else {
              row_vals[[safe_nm]] <- lower_mat[i, j]
            }
          }
          mtbl$addRow(rowKey = paste0("row_", i), values = row_vals)
        }
        # Add footnote explaining the layout
        metric_label <- switch(metric,
          Dprime = "D\'", r2 = "r²", D = "D")
        mtbl$setNote(
          key  = "layout",
          note = paste0("Upper triangle: ", metric_label,
                        ". Lower triangle: P-value. Diagonal: SNP name."))
      }

      # ── Store LD data for heatmap plot ──────────────────────────
      if (isTRUE(opts$ldPlot)) {
        private$.ld_store  <- ld_store
        private$.ld_nms    <- nms
        private$.ld_metric <- opts$ldMetric
        self$results$ldGroup$ldPlotImage$setState(list(
          ld_store = ld_store, nms = nms, metric = opts$ldMetric))
      }
    },

    # ── LD heatmap render ─────────────────────────────────────────
    .render_ld_plot = function(image, ggtheme, theme, ...) {
      state <- image$state
      if (is.null(state)) return(FALSE)

      ld_store <- state$ld_store
      nms      <- state$nms
      metric   <- state$metric
      n        <- length(nms)

      # Build a data frame for ggplot (full symmetric matrix for heatmap)
      metric_label <- switch(metric, Dprime = "D'", r2 = "r²", D = "D")
      df_rows <- list()
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          if (i == j) {
            val <- 1.0
          } else {
            key <- if (i < j) paste(c(nms[i], nms[j]), collapse = "___")
                   else        paste(c(nms[j], nms[i]), collapse = "___")
            ld_res <- ld_store[[key]]
            val <- if (!is.null(ld_res)) {
              switch(metric,
                Dprime = abs(as.numeric(ld_res$`D'`)),
                r2     = as.numeric(ld_res$`r`)^2,
                D      = abs(as.numeric(ld_res$`D`))
              )
            } else NA_real_
          }
          df_rows[[length(df_rows) + 1L]] <- data.frame(
            SNP1  = factor(nms[i], levels = rev(nms)),
            SNP2  = factor(nms[j], levels = nms),
            value = val,
            stringsAsFactors = FALSE
          )
        }
      }
      df <- do.call(rbind, df_rows)

      # Value string for lower-triangle annotation (p-values)
      p_mat <- matrix(NA_real_, n, n, dimnames = list(nms, nms))
      for (pair_key in names(ld_store)) {
        parts <- strsplit(pair_key, "___")[[1]]
        pv    <- ld_store[[pair_key]]$`P-value`
        p_mat[parts[1], parts[2]] <- pv
        p_mat[parts[2], parts[1]] <- pv
      }

      df$label <- ""
      for (k in seq_len(nrow(df))) {
        i_nm <- as.character(df$SNP1[k])
        j_nm <- as.character(df$SNP2[k])
        i_idx <- which(nms == i_nm)
        j_idx <- which(nms == j_nm)
        if (i_idx > j_idx) {   # lower triangle → p-value
          pv <- p_mat[i_nm, j_nm]
          df$label[k] <- if (!is.na(pv)) {
            if (pv < 0.001) "<.001" else sprintf("%.3f", pv)
          } else ""
        } else if (i_idx < j_idx) {  # upper triangle → metric value
          key <- paste(c(nms[min(i_idx, j_idx)], nms[max(i_idx, j_idx)]),
                       collapse = "___")
          ld_res <- ld_store[[key]]
          if (!is.null(ld_res)) {
            raw <- switch(metric,
              r2      = ld_res$`r`^2,
              Dprime = ld_res$`D'`,
              D      = ld_res$`D`
            )
            df$label[k] <- sprintf("%.3f", round(as.numeric(raw), 3))
          }
        } else {
          df$label[k] <- i_nm   # diagonal
        }
      }

      colour_label <- switch(metric,
        Dprime = "|D'|", r2 = "r²", D = "|D|")

      # Use scale_fill_gradientn with explicit continuous colours to avoid
      # jamovi theme interfering with the fill scale via ggPalette()
      p <- ggplot2::ggplot(df, ggplot2::aes(x = SNP2, y = SNP1, fill = value)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
        ggplot2::geom_text(ggplot2::aes(label = label),
                           size = 3, colour = "grey10") +
        ggplot2::scale_fill_gradientn(
          colours  = c("#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"),
          limits   = c(0, 1),
          na.value = "grey85",
          name     = colour_label
        ) +
        ggplot2::scale_x_discrete(position = "bottom") +
        ggplot2::scale_y_discrete() +
        ggplot2::labs(
          title = paste0("LD Heatmap  •  upper: ", metric_label,
                         "  |  lower: p-value"),
          x = NULL, y = NULL
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
          axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y     = ggplot2::element_text(hjust = 1),
          panel.grid      = ggplot2::element_blank(),
          legend.position = "right",
          plot.title      = ggplot2::element_text(size = 11, face = "bold",
                                                  margin = ggplot2::margin(b = 8))
        )

      print(p)
      TRUE
    },

    # ── SNP association ───────────────────────────────────────────
    .fill_assoc = function(tbl, snp_raw, ref, response, cov_df,
                            response_type, opts) {

#      tbl <- self$results$assocGroup$assocTable

      # ── Dynamic column header ─────────────────────────
      effect_col <- tbl$getColumn("effect")
      if (response_type == "binary") {
        effect_col$setTitle("OR")
      } else {
        effect_col$setTitle("\u03B2")
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

      note_key <- "lrt"
      lrt_note <- NULL
      
      row_key <- 0L
      for (mdl in models) {
        snp_enc <- encode_model(as.character(snp_raw), ref, mdl)
        res_list <- fit_model(snp_enc, response, cov_df, mdl,
                               response_type, opts$ciWidth)
        if (is.null(res_list)) next

        if (mdl == "codominant" && !is.null(res_list) && length(res_list) > 0) {
          gp <- res_list[[1]]$global_p
          if (!is.na(gp)) {
            lrt_note <- paste0(
              "Codominant model: likelihood ratio test P = ",
              format.pval(gp, digits = 3)
            )
          }
        }

        # ── Set / remove LRT note ──────────────────────────
        if (!is.null(lrt_note)) {
          tbl$setNote(note = lrt_note, key = note_key)
        } else {
          tbl$setNote(note = NULL, key = note_key)
        }

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
            AIC = if (first_row && !is.nan(res$aic)) res$aic else ""
          ))
          first_row <- FALSE
        }
      }
    },

    .run_haplo = function(geno_list, data, response, response_type,
                           cov_df, opts) {
      
      snp_names <- names(geno_list)
      allele_list <- lapply(snp_names, function(nm) genetics::allele(geno_list[[nm]]))
      allele_mat  <- do.call(cbind, allele_list)
      
      geno_setup <- tryCatch(
        haplo.stats::setupGeno(allele_mat, locus.label = snp_names),
        error = function(e) NULL
      )
      if (is.null(geno_setup)) return()
      
      u_alleles <- attr(geno_setup, "unique.alleles")

      decode_haplo_row <- function(codes, label_list) {
        # Check if this is the "rare" indicator from haplo.df
        if (all(codes == "*") || any(codes == "*")) return("Rare (combined)")
        
        parts <- character(length(codes))
        for (i in seq_along(codes)) {
          idx <- as.numeric(codes[i])
          parts[i] <- if (!is.na(idx)) label_list[[i]][idx] else "?"
        }
        paste(parts, collapse = "-")
      }

      keep <- if (!is.null(response)) !is.na(response) else rep(TRUE, nrow(allele_mat))
      if (!is.null(cov_df)) keep <- keep & complete.cases(cov_df)

      # 1. Haplotype Frequencies (EM)
      if (opts$haploFreq) {
        em_res <- tryCatch(
          haplo.stats::haplo.em(geno_setup[keep, ], locus.label = snp_names),
          error = function(e) NULL
        )
        
        if (!is.null(em_res)) {
          tbl <- self$results$haploGroup$haploFreqTable
          freqs <- em_res$hap.prob
          rare_sum <- 0
          
          for (i in seq_along(freqs)) {
            if (freqs[i] < opts$haploFreqMin) {
              rare_sum <- rare_sum + freqs[i]
              next
            }
            label <- decode_haplo_row(as.numeric(em_res$haplotype[i, ]), u_alleles)
            tbl$addRow(rowKey = paste0("f", i), values = list(
              haplotype = label,
              freq      = round(freqs[i], 4)
            ))
          }
          if (rare_sum > 0) {
            tbl$addRow(rowKey = "rare_freq", values = list(
              haplotype = "Rare (combined)",
              freq      = round(rare_sum, 4)
            ))
          }
        }
      }

 # 2. Association Table (Fixes Empty Table Issue)
    if (opts$haploAssoc && !is.null(response)) {
        family <- if (response_type == "binary") binomial else gaussian
        y_sub  <- if (response_type == "binary") as.numeric(as.factor(response[keep])) - 1 else response[keep]
        
        # Merge geno_setup and covariates into one data frame
        # We name the geno_setup column explicitly to match the formula
        m_model <- data.frame(y = y_sub)
        m_model$geno <- geno_setup
        if (!is.null(cov_df)) {
            m_model <- cbind(m_model, cov_df[keep, , drop = FALSE])
            cov_names <- names(cov_df)
            formula_str <- paste("y ~ geno +", paste(cov_names, collapse = " + "))
        } else {
            formula_str <- "y ~ geno"
        }
        haplo_fit <- tryCatch({
          haplo.stats::haplo.glm(
            as.formula(formula_str),
            family   = family,
            data     = m_model,
            na.action = "na.geno.keep", 
            control  = haplo.stats::haplo.glm.control(haplo.freq.min = opts$haploFreqMin)
          )
        }, error = function(e) {
          self$results$validationMsg$setContent(paste0("Haplotype GLM Error: ", e$message))
          NULL
        })
        if (!is.null(haplo_fit)) {
          tbl <- self$results$haploGroup$haploAssocTable

          # DEBUG DUMP
          dbg <- c()
          dbg <- c(dbg, paste0("<b>names(haplo_fit):</b> ", paste(names(haplo_fit), collapse=", ")))
          dbg <- c(dbg, paste0("<b>coef names:</b> ", paste(names(coefficients(haplo_fit)), collapse=", ")))
          dbg <- c(dbg, paste0("<b>haplo.common:</b> ", paste(haplo_fit$haplo.common, collapse=", ")))
          dbg <- c(dbg, paste0("<b>haplo.base:</b> ",   paste(haplo_fit$haplo.base,   collapse=", ")))
          dbg <- c(dbg, paste0("<b>haplo.rare:</b> ",   paste(haplo_fit$haplo.rare,   collapse=", ")))
          dbg <- c(dbg, paste0("<b>haplo.freq (len):</b> ", length(haplo_fit$haplo.freq)))
          dbg <- c(dbg, paste0("<b>names(haplo.freq):</b> ", paste(names(haplo_fit$haplo.freq), collapse=", ")))
          dbg <- c(dbg, paste0("<b>names(haplo_fit$haplo):</b> ", paste(names(haplo_fit$haplo), collapse=", ")))
          if (!is.null(haplo_fit$var.mat))
            dbg <- c(dbg, paste0("<b>var.mat rownames:</b> ", paste(rownames(haplo_fit$var.mat), collapse=", ")))
          hd_str <- tryCatch({
            hd2 <- haplo.stats::haplo.df(haplo_fit)
            paste0("dim=", nrow(hd2), "x", ncol(hd2), " cols:", paste(colnames(hd2), collapse="|"), "<br>",
                   paste(capture.output(print(head(hd2))), collapse="<br>"))
          }, error = function(e) paste("haplo.df ERROR:", e$message))
          dbg <- c(dbg, paste0("<b>haplo.df():</b><br>", hd_str))
          self$results$validationMsg$setContent(paste(dbg, collapse="<br>"))
          # END DEBUG

          # Using haplo.df correctly extracts coefficients regardless of naming convention
          h_df <- haplo.stats::haplo.df(haplo_fit)
          ci   <- tryCatch(confint(haplo_fit, level = opts$ciWidth/100), error = function(e) NULL)
          
          n_snps <- length(snp_names)
          for (i in 1:nrow(h_df)) {
            codes <- h_df[i, 1:n_snps]
            label <- decode_haplo_row(codes, u_alleles)
            beta  <- h_df[i, "Estimate"]
            
            row_vals <- list(
              haplotype = label,
              freq      = h_df[i, "Hap-Freq"],
              pval      = h_df[i, "p-value"]
            )

            if (!is.na(beta)) {
              row_vals$effect <- if (response_type == "binary") exp(beta) else beta
              # Map CI by identifying the correct coefficient row name
              target <- if (any(codes == "*")) "geno_setup.rare" else paste0("geno_setup.", paste(as.numeric(codes), collapse="."))
              if (!is.null(ci) && target %in% rownames(ci)) {
                row_vals$ciLow  <- if (response_type == "binary") exp(ci[target, 1]) else ci[target, 1]
                row_vals$ciHigh <- if (response_type == "binary") exp(ci[target, 2]) else ci[target, 2]
              }
            } else {
              # This is the reference group
              row_vals$effect <- if (response_type == "binary") 1.0 else 0.0
              row_vals$haplotype <- paste0(label, " (Ref)")
            }
            tbl$addRow(rowKey = paste0("assoc", i), values = row_vals)
          }
        }
      }
    },

    # ── Private LD storage for plot render ────────────────────────
    .ld_store  = NULL,
    .ld_nms    = NULL,
    .ld_metric = NULL
  )
)  # end R6Class
