# ── snp_compute.R ─────────────────────────────────────────────────────────────
#
# Pure-R computation layer for SNPstats.
# No jamovi dependencies — every function takes plain R objects and returns
# plain R lists / data frames.  Both the jamovi class (snpStats_b.R) and the
# exported snpStats() function (snpStats.R) call these functions so that
# results are guaranteed to be identical in both contexts.
#
# Source order:
#   snp_helpers.R   (utilities: parse_genotype, encode_model, fit_model, …)
#   snp_compute.R   (this file)
#   snpStats_b.R    (jamovi class — calls snp_prepare() + compute_*())
#   snpStats.R      (exported function — calls snp_prepare() + compute_*())
# ──────────────────────────────────────────────────────────────────────────────

# ── Formatting helpers ────────────────────────────────────────────────────────

#' Round a p-value to a display string.
#' Values < 0.001 are shown as "< 0.001"; otherwise 3 significant figures.
fmt_pval <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return('')
  vapply(x, function(p) {
    if (!is.numeric(p)) return(as.character(p))
    if (p < 0.001)     return("< 0.001")
    format.pval(p, digits = 3, eps = 0.001, nsmall = 3, scientific = FALSE)
  }, '')
}

#' Round an effect / CI value (OR or beta) to 3 decimal places.
fmt3 <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return('')
  formatC(round(as.numeric(x), 3), format = "f", flag = "#", digits = 3)
}



# ══════════════════════════════════════════════════════════════════════════════
# snp_prepare — single entry point for all preprocessing
# ══════════════════════════════════════════════════════════════════════════════

#' Validate, clean, and pre-process all inputs.
#'
#' Returns a named list ("prep object") consumed by every compute_* function.
#' Call this once, then pass the result to whichever compute_* functions you
#' need.  The jamovi .run() method and the exported snpStats() function both
#' delegate to this function so their preprocessing is identical.
#'
#' @param data           data.frame
#' @param snps           character vector of SNP column names
#' @param response       character scalar — column name, or NULL
#' @param covariates     character vector of covariate column names
#' @param response_type  "auto" | "binary" | "quantitative" | "categorical"
#' @param rm_snp_missing logical — exclude rows missing any SNP before fitting
#'
#' @return Named list with elements:
#'   $data, $snp_vars, $snp_data (per-SNP parsed objects),
#'   $response_var, $response_raw, $response_type, $response_enc,
#'   $cov_df, $complete_mask, $n_rows, $warnings
snp_prepare <- function(data, snps, response = NULL, covariates = NULL,
                        response_type = "auto", rm_snp_missing = FALSE) {

  # ── Response ─────────────────────────────────────────────────────────
  response_var <- if (!is.null(response) && nchar(response) > 0) response else NULL
  response_raw <- if (!is.null(response_var)) data[[response_var]] else NULL
  rtype        <- detect_response_type(response_raw, response_type)
  response_enc <- prepare_response(response_raw, rtype)

  # ── Covariates ────────────────────────────────────────────────────────
  cov_df <- prepare_covariates(data, covariates %||% character(0))
  if (is.null(cov_df) && !is.null(response_raw))
    cov_df <- data.frame(row.names = seq_len(nrow(data)))

  # ── Complete-case mask ────────────────────────────────────────────────
  n_rows        <- nrow(data)
  complete_mask <- rep(TRUE, n_rows)
  if (!is.null(response_enc))                 complete_mask <- complete_mask & !is.na(response_enc)
  if (!is.null(cov_df) && ncol(cov_df) > 0)  complete_mask <- complete_mask & complete.cases(cov_df)

  # ── SNP columns ─────────────────────────────────────────────
  if (length(snps) == 0) {
    val<- list(bad_html = "No SNPs specified — skipping SNP processing")
    snp_data <- character(0)
    snp_vars <- character(0)
  } else {
    # ── Validate SNP columns ─────────────────────────────────────────────
    val      <- validate_snp_vars(snps, data)
    snp_vars <- val$valid_snps

    # ── Optional: additionally exclude rows missing any SNP ──────────────
    if (isTRUE(rm_snp_missing) && length(snp_vars) > 0) {
      snp_mat <- as.data.frame(
        lapply(data[, snp_vars, drop = FALSE],
              function(col) clean_null_alleles(as.character(col))),
        stringsAsFactors = FALSE)
      complete_mask <- complete_mask & complete.cases(snp_mat)
    }

    # ── Per-SNP: clean → parse → subset to complete cases ────────────────
    snp_data <- lapply(setNames(snp_vars, snp_vars), function(nm) {
      raw         <- data[[nm]]
      user_levels <- get_snp_level_order(raw)           # from original factor
      clean       <- clean_null_alleles(as.character(raw))
      snp_mask    <- complete_mask & !is.na(clean)      # per-SNP complete-case
      clean_cc    <- clean[snp_mask]
      geno_cc     <- parse_genotype(clean_cc, user_levels)
      if (is.null(geno_cc)) return(NULL)
      ref         <- get_ref_genotype(geno_cc, user_levels)
      summary_cc  <- summary(geno_cc)
      list(
        raw         = raw,
        clean       = clean,
        user_levels = user_levels,
        snp_mask    = snp_mask,
        clean_cc    = clean_cc,
        geno_cc     = geno_cc,
        ref         = ref,
        summary_cc  = summary_cc,
        n_typed     = sum(snp_mask),
        n_missing   = sum(!is.na(clean) & complete_mask & !snp_mask) +
                      sum(is.na(clean)  & complete_mask)
      )
    })
    snp_data <- Filter(Negate(is.null), snp_data)
  } 

  list(
    data          = data,
    snp_vars      = names(snp_data),
    snp_data      = snp_data,
    response_var  = response_var,
    response_raw  = response_raw,
    response_type = rtype,
    response_enc  = response_enc,
    cov_df        = cov_df,
    complete_mask = complete_mask,
    n_rows        = n_rows,
    warnings      = val$bad_html
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Descriptive compute functions
# ══════════════════════════════════════════════════════════════════════════════

#' Compute allele frequency table for one SNP.
#'
#' @param snp_nm  SNP name (key into prep$snp_data)
#' @param prep    prep object from snp_prepare()
#' @param subpop  logical — stratify by response group
#' @param show_missing  logical — include a missing-count row
#'
#' @return data.frame with columns: allele, overall, stat_g0..stat_gN (uniform)
compute_allele_freq <- function(snp_nm, prep, subpop = FALSE,
                                show_missing = FALSE) {
  sd          <- prep$snp_data[[snp_nm]]
  sm          <- sd$summary_cc
  af          <- sm$allele.freq
  user_levels <- sd$user_levels
  rtype       <- prep$response_type

  # Use full-length response so masks align with sd$clean / prep$complete_mask
  response_raw_full <- prep$response_raw          # length = n_rows
  response_raw_cc   <- if (!is.null(response_raw_full))
                         response_raw_full[sd$snp_mask] else NULL  # cc-subset

  # Allele display order: ref allele first
  allele_nms <- rownames(af)[rownames(af) != "NA"]
  ref_allele <- strsplit(sd$ref, "/", fixed = TRUE)[[1]][1]
  if (!is.null(user_levels)) {
    for (g in user_levels) {
      p <- strsplit(g, "/", fixed = TRUE)[[1]]
      if (length(p) == 2 && p[1] == p[2]) { ref_allele <- p[1]; break }
    }
  }
  if (ref_allele %in% allele_nms)
    allele_nms <- c(ref_allele, setdiff(allele_nms, ref_allele))

  do_strat   <- isTRUE(subpop) && !is.null(response_raw_cc) &&
                rtype %in% c("binary", "categorical")
  grp_levels <- if (do_strat) levels(as.factor(response_raw_full)) else character(0)
  n_grp      <- length(grp_levels)

  # Per-group allele totals (denominator = 2 × typed in group)
  # Uses cc-subset: response_raw_cc and sd$clean_cc are both snp_mask-aligned
  alleles_split_cc <- if (do_strat) strsplit(sd$clean_cc, "/", fixed = TRUE) else NULL
  resp_chr_cc      <- if (do_strat) as.character(response_raw_cc) else NULL

  grp_allele_totals <- if (do_strat)
    sapply(grp_levels, function(lvl) {
      idx <- resp_chr_cc == lvl & !is.na(resp_chr_cc)
      length(unlist(alleles_split_cc[idx]))
    }) else NULL

  make_row <- function(allele_lbl, overall_str, grp_strs = rep("", n_grp)) {
    row <- list(allele = allele_lbl, overall = overall_str)
    for (j in seq_len(n_grp)) row[[paste0("stat_g", j-1L)]] <- grp_strs[j]
    as.data.frame(row, stringsAsFactors = FALSE)
  }

  rows <- lapply(allele_nms, function(al) {
    if (!al %in% rownames(af)) return(NULL)
    count <- as.integer(af[al, "Count"])
    prop  <- round(af[al, "Proportion"] * 100, 1)
    grp_strs <- if (do_strat) {
      vapply(grp_levels, function(lvl) {
        idx   <- resp_chr_cc == lvl & !is.na(resp_chr_cc)
        n_al  <- sum(unlist(alleles_split_cc[idx]) == al, na.rm = TRUE)
        fmt_cat(n_al, grp_allele_totals[[lvl]])
      }, character(1))
    } else rep("", 0L)
    make_row(al, fmt_catpct(count, prop), grp_strs)
  })
  rows <- Filter(Negate(is.null), rows)

  if (isTRUE(show_missing) && sd$n_missing > 0) {
    n_elig <- sum(prep$complete_mask)
    # Per-group missing: rows in complete_mask that have NA SNP for this group
    grp_miss <- if (do_strat) {
      clean_full   <- sd$clean
      resp_chr_full <- as.character(response_raw_full)
      # per-group total eligible = complete_mask rows in that group
      grp_eligible <- sapply(grp_levels, function(lvl)
        sum(prep$complete_mask & !is.na(resp_chr_full) &
            resp_chr_full == lvl, na.rm = TRUE))
      vapply(grp_levels, function(lvl) {
        n_lv <- sum(is.na(clean_full) & prep$complete_mask &
                    !is.na(resp_chr_full) & resp_chr_full == lvl, na.rm = TRUE)
        fmt_cat(n_lv, grp_eligible[[lvl]])
      }, character(1))
    } else rep("", 0L)
    rows <- c(rows, list(make_row("Missing", fmt_cat(sd$n_missing, n_elig), grp_miss)))
  }

  do.call(rbind, rows)
}


#' Compute genotype frequency table for one SNP.
#'
#' @return data.frame with columns: genotype, overall, stat_g0..stat_gN (uniform),
#'         response_stat (mean±SE for quantitative)
compute_geno_freq <- function(snp_nm, prep, subpop = FALSE,
                              show_missing = FALSE) {
  sd           <- prep$snp_data[[snp_nm]]
  sm           <- sd$summary_cc
  rtype        <- prep$response_type
  response_cc  <- if (!is.null(prep$response_enc)) prep$response_enc[sd$snp_mask] else NULL

  # Full-length response for missing-row alignment; cc-subset for typed-sample counts
  response_raw_full <- prep$response_raw
  response_raw_cc   <- if (!is.null(response_raw_full))
                         response_raw_full[sd$snp_mask] else NULL

  geno_obj <- sd$geno_cc
  snp_chr  <- as.character(geno_obj)   # normalised, matches gf rownames

  gf <- tryCatch(
    reorder_geno(sm$genotype.freq, sd$ref, sd$user_levels),
    error = function(e) sm$genotype.freq)
  gf <- gf[rownames(gf) != "NA", , drop = FALSE]

  do_strat   <- isTRUE(subpop) && !is.null(response_raw_cc) &&
                rtype %in% c("binary", "categorical")
  grp_levels <- if (do_strat) levels(as.factor(response_raw_full)) else character(0)
  n_grp      <- length(grp_levels)

  # Denominators: typed observations per group (snp_chr not NA within group)
  resp_chr_cc <- if (do_strat) as.character(response_raw_cc) else NULL
  strat_totals <- if (do_strat)
    sapply(grp_levels, function(lvl)
      sum(resp_chr_cc == lvl & !is.na(resp_chr_cc) & !is.na(snp_chr), na.rm = TRUE))
  else NULL

  make_row <- function(geno_lbl, overall_str, resp_stat = "",
                       grp_strs = rep("", n_grp)) {
    row <- list(genotype = geno_lbl, overall = overall_str,
                response_stat = resp_stat)
    for (j in seq_len(n_grp)) row[[paste0("stat_g", j-1L)]] <- grp_strs[j]
    as.data.frame(row, stringsAsFactors = FALSE)
  }

  rows <- lapply(seq_len(nrow(gf)), function(i) {
    geno  <- rownames(gf)[i]
    if (geno == "NA") return(NULL)
    count <- as.integer(gf[i, "Count"])
    prop  <- gf[i, "Proportion"] * 100

    # Quantitative response stat
    resp_stat <- ""
    if (rtype == "quantitative" && !is.null(response_cc)) {
      mask   <- snp_chr == geno & !is.na(snp_chr) & !is.na(response_cc)
      n_mask <- sum(mask)
      if (n_mask > 0) {
        resp_num  <- as.numeric(response_cc)
        resp_stat <- sprintf("%.2f (%.2f)",
          mean(resp_num[mask], na.rm = TRUE),
          sd(resp_num[mask], na.rm = TRUE) / sqrt(n_mask))
      }
    }

    grp_strs <- if (do_strat) {
      vapply(grp_levels, function(lvl) {
        idx   <- resp_chr_cc == lvl & !is.na(resp_chr_cc)
        n_g   <- sum(idx & snp_chr == geno, na.rm = TRUE)
        fmt_cat(n_g, strat_totals[[lvl]])
      }, character(1))
    } else rep("", 0L)

    make_row(geno, fmt_catpct(count, prop), resp_stat, grp_strs)
  })
  rows <- Filter(Negate(is.null), rows)

  if (isTRUE(show_missing) && sd$n_missing > 0) {
    n_elig        <- sum(prep$complete_mask)
    clean_full    <- sd$clean
    resp_chr_full <- if (!is.null(response_raw_full))
                       as.character(response_raw_full) else NULL

    grp_miss <- if (do_strat) {
      # Denominator: all complete-case rows in that group (typed + missing SNP)
      grp_eligible <- sapply(grp_levels, function(lvl)
        sum(prep$complete_mask & !is.na(resp_chr_full) &
            resp_chr_full == lvl, na.rm = TRUE))
      vapply(grp_levels, function(lvl) {
        n_lv <- sum(is.na(clean_full) & prep$complete_mask &
                    !is.na(resp_chr_full) & resp_chr_full == lvl, na.rm = TRUE)
        fmt_cat(n_lv, grp_eligible[[lvl]])
      }, character(1))
    } else rep("", 0L)

    rows <- c(rows, list(make_row("Missing",
                                  fmt_cat(sd$n_missing, n_elig),
                                  "", grp_miss)))
  }

  do.call(rbind, rows)
}


#' Compute Hardy-Weinberg test for one SNP, optionally stratified.
#'
#' @return list with elements:
#'   $col_labels  character(3) — actual genotype labels for n11/n12/n22
#'   $rows        data.frame: group, n11, n12, n22, missing, pval
compute_hwe <- function(snp_nm, prep, subpop = FALSE,
                        show_missing = FALSE) {
  sd           <- prep$snp_data[[snp_nm]]
  geno_obj     <- sd$geno_cc
  response_raw <- if (!is.null(prep$response_raw)) prep$response_raw[sd$snp_mask] else NULL

  hw <- tryCatch(genetics::HWE.exact(geno_obj), error = function(e) NULL)
  if (is.null(hw)) return(NULL)

  # Column labels from genotype.freq (user order via reorder_geno)
  get_ordered_counts <- function(go) {
    gf <- tryCatch(
      reorder_geno(summary(go)$genotype.freq, sd$ref, sd$user_levels),
      error = function(e) summary(go)$genotype.freq)
    gf <- gf[rownames(gf) != "NA", , drop = FALSE]
    list(labels = rownames(gf), counts = as.integer(gf[, "Count"]))
  }

  info       <- get_ordered_counts(geno_obj)
  col_labels <- if (length(info$labels) == 3) info$labels else c("AA", "AB", "BB")

  rows <- list()
  miss_val <- if (isTRUE(show_missing)) sd$n_missing else NA_integer_
  rows[[1]] <- data.frame(
    group   = "All subjects",
    n11     = info$counts[1L], n12 = info$counts[2L], n22 = info$counts[3L],
    missing = miss_val, pval = hw$p.value,
    stringsAsFactors = FALSE)

  if (isTRUE(subpop) && !is.null(response_raw)) {
    lvls <- levels(as.factor(response_raw))
    if (length(lvls) <= 5) {
      for (lvl in lvls) {
        mask <- as.character(response_raw) == lvl & !is.na(response_raw)
        if (sum(mask) == 0) next
        hw_sub <- tryCatch(genetics::HWE.exact(geno_obj[mask]), error = function(e) NULL)
        if (is.null(hw_sub)) next
        sub_info <- get_ordered_counts(geno_obj[mask])
        rows[[length(rows) + 1L]] <- data.frame(
          group   = lvl,
          n11     = sub_info$counts[1L], n12 = sub_info$counts[2L],
          n22     = sub_info$counts[3L],
          missing = NA_integer_, pval = hw_sub$p.value,
          stringsAsFactors = FALSE)
      }
    }
  }

  list(col_labels = col_labels, rows = do.call(rbind, rows))
}


#' Compute SNP summary table (one row per SNP, optionally per group).
#'
#' @return data.frame: snp, alleles, group, n, missing, maf, geno_counts, hwe_pval
compute_snp_summary <- function(prep, subpop = FALSE) {
  rtype        <- prep$response_type
  response_raw <- prep$response_raw
  do_strat     <- isTRUE(subpop) && !is.null(response_raw) &&
                  rtype %in% c("binary", "categorical")
  grp_levels   <- if (do_strat) levels(as.factor(response_raw)) else NULL

  rows <- lapply(prep$snp_vars, function(nm) {
    sd          <- prep$snp_data[[nm]]
    user_levels <- sd$user_levels
    geno_cc     <- sd$geno_cc
    snp_cc      <- sd$clean_cc
    resp_cc     <- if (do_strat) response_raw[sd$snp_mask] else NULL
    sm_cc       <- sd$summary_cc
    ref         <- sd$ref
    ref_allele  <- strsplit(ref, "/", fixed = TRUE)[[1]][1]
    af_all      <- sm_cc$allele.freq
    allele_nms  <- rownames(af_all)
    alt_allele  <- setdiff(allele_nms[allele_nms != "NA"], ref_allele)
    alt_allele  <- if (length(alt_allele)) alt_allele[1] else "?"
    alleles_lbl <- paste0(ref_allele, "/", alt_allele)
    n_total     <- prep$n_rows
    n_excl      <- n_total - sd$n_typed

    compute_row_stats <- function(g_obj, a_allele) {
      sm   <- summary(g_obj)
      af   <- sm$allele.freq
      props <- af[rownames(af) != "NA", "Proportion"]
      maf  <- if (a_allele %in% rownames(af)) af[a_allele, "Proportion"]
              else min(props, na.rm = TRUE)
      gf   <- tryCatch(reorder_geno(sm$genotype.freq, ref, user_levels),
                       error = function(e) sm$genotype.freq)
      gf   <- gf[rownames(gf) != "NA", , drop = FALSE]
      cnts <- as.integer(gf[, "Count"])
      geno_str <- switch(as.character(length(cnts)),
        "3" = paste(cnts, collapse = " / "),
        "2" = paste(c(cnts, 0L), collapse = " / "),
        paste(cnts, collapse = " / "))
      hwe <- tryCatch(genetics::HWE.exact(g_obj)$p.value, error = function(e) NA_real_)
      list(n = sm$n.typed, maf = round(maf, 4), geno_counts = geno_str, hwe_pval = hwe)
    }

    res_all <- compute_row_stats(geno_cc, alt_allele)
    out <- list(data.frame(
      snp = nm, alleles = alleles_lbl,
      group   = if (do_strat) "All" else "",
      n       = res_all$n,
      missing = if (n_excl > 0L) n_excl else NA_integer_,
      maf     = res_all$maf, geno_counts = res_all$geno_counts,
      hwe_pval = res_all$hwe_pval,
      stringsAsFactors = FALSE))

    if (do_strat) {
      resp_chr <- as.character(resp_cc)
      strat_totals <- table(factor(resp_cc, levels = grp_levels))
      for (lvl in grp_levels) {
        mask <- !is.na(resp_chr) & resp_chr == lvl
        if (sum(mask) == 0) next
        g_sub <- tryCatch(parse_genotype(snp_cc[mask], user_levels), error = function(e) NULL)
        if (is.null(g_sub)) next
        res_s  <- compute_row_stats(g_sub, alt_allele)
        n_excl_s <- max(0L, as.integer(strat_totals[lvl]) - res_s$n)
        out[[length(out) + 1L]] <- data.frame(
          snp = "", alleles = "", group = lvl,
          n = res_s$n,
          missing = if (n_excl_s > 0L) n_excl_s else NA_integer_,
          maf = res_s$maf, geno_counts = res_s$geno_counts,
          hwe_pval = res_s$hwe_pval,
          stringsAsFactors = FALSE)
      }
    }
    do.call(rbind, out)
  })

  do.call(rbind, Filter(Negate(is.null), rows))
}


#' Compute covariate descriptive table.
#'
#' Returns a uniformly structured data frame regardless of whether stratification
#' is active. Every row always has the same columns so rbind() never misaligns:
#'   variable, level, overall, stat_g0 … stat_gN (when do_strat), pval
#'
#' @return list: $table data.frame, $notes character
compute_cov_desc <- function(prep, subpop = FALSE) {
  cov_df       <- prep$cov_df
  response_raw <- prep$response_raw
  rtype        <- prep$response_type
  response_var <- prep$response_var

  if (is.null(cov_df)) return(NULL)

  do_strat   <- isTRUE(subpop) && !is.null(response_raw) &&
                rtype %in% c("binary", "categorical")
  valid_resp  <- if (!is.null(response_raw)) !is.na(response_raw)
                 else rep(TRUE, prep$n_rows)

  # Group metadata
  grp_levels <- NULL; mask_list <- NULL; totals <- NULL
  if (do_strat) {
    grp_fac    <- as.factor(response_raw)
    grp_levels <- levels(grp_fac)
    mask_list  <- lapply(grp_levels, function(l) valid_resp & as.character(grp_fac) == l)
    names(mask_list) <- grp_levels
    totals     <- sapply(mask_list, sum)
  }
  n_grp <- length(grp_levels)   # 0 when not stratified

  # ── Uniform row builder ────────────────────────────────────────────────────
  # Always produces a row with: variable, level, overall,
  # stat_g0..stat_g(n_grp-1) (empty string when not stratified), pval (NA when not)
  make_row <- function(variable, level, overall,
                       grp_stats = rep("", n_grp),
                       pval      = NA_real_) {
    row <- list(variable = as.character(variable),
                level    = as.character(level),
                overall  = as.character(overall))
    for (j in seq_len(n_grp))
      row[[paste0("stat_g", j - 1L)]] <- as.character(grp_stats[j])
    row$pval <- as.numeric(pval)
    as.data.frame(row, stringsAsFactors = FALSE)
  }

  # Compute per-group stat string for a logical mask (same length as prep$n_rows)
  grp_stats_for <- function(mask, fmt_fn = function(n, tot) fmt_cat(n, tot)) {
    if (!do_strat) return(rep("", 0L))
    vapply(grp_levels, function(lvl)
      fmt_fn(sum(mask & mask_list[[lvl]], na.rm = TRUE), totals[[lvl]]),
      character(1))
  }

  grp_stats_cont <- function(x_vec) {
    if (!do_strat) return(rep("", 0L))
    vapply(grp_levels, function(lvl)
      fmt_cont(x_vec[mask_list[[lvl]]]),
      character(1))
  }

  rows <- list()

  # ── Response variable rows ─────────────────────────────────────────────────
  if (!is.null(response_raw) && !is.null(rtype)) {
    if (rtype == "quantitative") {
      rows[[length(rows)+1L]] <- make_row(
        response_var, "Mean \u00B1 SD",
        fmt_cont(response_raw),
        grp_stats_cont(response_raw))
      mask <- !is.na(response_raw)
      rows[[length(rows)+1L]] <- make_row(
        "", "Valid", as.character(sum(mask)),
        grp_stats_for(mask, function(n, tot) as.character(n)))
      if (sum(!mask) > 0)
        rows[[length(rows)+1L]] <- make_row(
          "", "Missing", fmt_cat(sum(!mask), length(response_raw)),
          grp_stats_for(!mask))
    } else {
      mask <- valid_resp
      rows[[length(rows)+1L]] <- make_row(
        response_var, "Valid", as.character(sum(mask)),
        grp_stats_for(mask, function(n, tot) as.character(n)))
      if (sum(!valid_resp) > 0)
        rows[[length(rows)+1L]] <- make_row(
          "", "Missing", fmt_cat(sum(!valid_resp), length(response_raw)),
          grp_stats_for(!valid_resp))
    }
  }

  # ── Covariate rows ─────────────────────────────────────────────────────────
  if (ncol(cov_df) > 0) {
    for (v in names(cov_df)) {
      col    <- cov_df[[v]]
      n      <- length(col)
      n_miss <- sum(is.na(col))
      is_cat <- is.factor(col) || is.character(col)
      if (is_cat && !is.factor(col)) col <- factor(col)

      if (is_cat) {
        # p-value: chi-squared test across groups
        pval_cat <- if (do_strat) tryCatch({
          ct <- table(col[valid_resp], as.factor(response_raw)[valid_resp])
          suppressWarnings(chisq.test(ct)$p.value)
        }, error = function(e) NA_real_) else NA_real_

        first <- TRUE
        for (lvl in levels(col)) {
          mask <- !is.na(col) & col == lvl
          rows[[length(rows)+1L]] <- make_row(
            if (first) v else "", lvl,
            fmt_cat(sum(mask), n),
            grp_stats_for(mask),
            if (first) pval_cat else NA_real_)
          first <- FALSE
        }
        if (n_miss > 0)
          rows[[length(rows)+1L]] <- make_row(
            "", "Missing", fmt_cat(n_miss, n),
            grp_stats_for(is.na(col)))

      } else {
        # Continuous variable: t-test (2 groups) or ANOVA (>2)
        pval_cont <- if (do_strat) tryCatch({
          grps <- split(col[valid_resp], as.factor(response_raw)[valid_resp])
          if (length(grps) == 2)
            t.test(grps[[1]], grps[[2]])$p.value
          else
            summary(aov(col ~ as.factor(response_raw)))[[1]][["Pr(>F)"]][1]
        }, error = function(e) NA_real_) else NA_real_

        rows[[length(rows)+1L]] <- make_row(
          v, "Mean \u00B1 SD",
          fmt_cont(col),
          grp_stats_cont(col),
          pval_cont)

        if (n_miss > 0)
          rows[[length(rows)+1L]] <- make_row(
            "", "Missing", fmt_cat(n_miss, n),
            grp_stats_for(is.na(col)))
      }
    }
  }

  list(table = do.call(rbind, rows),
       grp_levels = grp_levels,
       notes = character(0))
}


# ══════════════════════════════════════════════════════════════════════════════
# Association compute functions
# ══════════════════════════════════════════════════════════════════════════════

# Internal helpers mirroring private methods in snpStatsClass

.geno_labels_for_model <- function(model, all_genos, ref) {
  het  <- all_genos[all_genos != ref & sapply(all_genos, function(g) {
    p <- strsplit(g, "/")[[1]]; length(p) == 2 && p[1] != p[2] })]
  hom2 <- all_genos[all_genos != ref & sapply(all_genos, function(g) {
    p <- strsplit(g, "/")[[1]]; length(p) == 2 && p[1] == p[2] })]
  switch(model,
    codominant   = all_genos,
    dominant     = c(ref, paste(c(het, hom2), collapse = "-")),
    recessive    = c(paste(c(ref, het), collapse = "-"), hom2[1]),
    overdominant = c(paste(c(ref, hom2), collapse = "-"), het[1]),
    logadditive  = all_genos)
}

.split_genos <- function(gl) {
  strsplit(gl, "-", fixed = TRUE)[[1]]
}

# NOTE: .bic_from_aic is kept as a legacy fallback only.
# fit_model() now returns BIC(fit_full) directly as res$bic, which is always
# correct.  The formula below was incorrect for quantitative (lm) responses
# because R's AIC/BIC for lm counts the residual variance (sigma^2) as an
# extra parameter (k_lm = p_betas + 1_intercept + 1_sigma), but the old
# code used k = 1 + n_cov + snp_df, missing the +1 for sigma.
.bic_from_aic <- function(aic_val, mdl, n_fit, n_cov, response_type = "binary") {
  snp_df <- c(codominant=2L, dominant=1L, recessive=1L, overdominant=1L, logadditive=1L)
  if (is.null(aic_val) || is.na(aic_val) || is.nan(aic_val)) return(NA_real_)
  # For lm (quantitative), R counts sigma^2 as an extra free parameter, so k is
  # one higher than for glm(binomial).  Adjust accordingly.
  sigma_extra <- if (!is.null(response_type) && response_type == "quantitative") 1L else 0L
  k <- 1L + n_cov + snp_df[[mdl]] + sigma_extra
  round(aic_val + k * (log(n_fit) - 2), 2)
}

.compute_stats <- function(geno_labels, snp_char, response, response_type,
                            response_raw = NULL) {
  split_genos <- .split_genos
  ref_al <- NULL
  for (lbl in geno_labels) {
    # geno_labels can contain compound dash-joined labels (e.g. "A/A-G/G" for
    # overdominant/recessive/dominant collapsed groups).  Split on "-" first to
    # obtain individual "X/Y" genotype strings before checking for homozygosity.
    individual_genos <- strsplit(lbl, "-", fixed = TRUE)[[1]]
    for (g in individual_genos) {
      parts <- strsplit(g, "/", fixed = TRUE)[[1]]
      if (length(parts) == 2 && parts[1] == parts[2]) { ref_al <- parts[1]; break }
    }
    if (!is.null(ref_al)) break
  }
  norm_snp_char <- function(sc) {
    if (is.null(ref_al)) return(sc)
    sapply(sc, function(g) {
      if (is.na(g)) return(NA_character_)
      p <- strsplit(g, "/", fixed = TRUE)[[1]]
      if (length(p) == 2 && p[1] != p[2] && p[2] == ref_al) paste0(p[2], "/", p[1])
      else g
    }, USE.NAMES = FALSE)
  }
  sc <- norm_snp_char(snp_char)

  if (response_type == "binary") {
    resp_grp <- if (!is.null(response_raw)) response_raw else response
    lv       <- levels(as.factor(resp_grp))
    if (length(lv) < 2) lv <- c(lv, "")
    n_col0   <- sum(resp_grp == lv[1] & !is.na(resp_grp))
    n_col1   <- sum(resp_grp == lv[2] & !is.na(resp_grp))
    stats0   <- character(length(geno_labels))
    stats1   <- character(length(geno_labels))
    for (i in seq_along(geno_labels)) {
      mask      <- sc %in% split_genos(geno_labels[i]) & !is.na(resp_grp)
      n0        <- sum(mask & resp_grp == lv[1])
      n1        <- sum(mask & resp_grp == lv[2])
      stats0[i] <- sprintf("%d (%.1f%%)", n0, if (n_col0 > 0) n0/n_col0*100 else 0)
      stats1[i] <- sprintf("%d (%.1f%%)", n1, if (n_col1 > 0) n1/n_col1*100 else 0)
    }
    list(s0 = stats0, s1 = stats1)
  } else if (response_type == "categorical") {
    resp_grp  <- if (!is.null(response_raw)) response_raw else response
    cats      <- levels(as.factor(resp_grp))
    n_cats    <- sapply(cats, function(c) sum(resp_grp == c & !is.na(resp_grp)))
    cat_stats <- lapply(cats, function(cat) {
      n_cat <- n_cats[[cat]]
      sapply(seq_along(geno_labels), function(i) {
        mask <- sc %in% split_genos(geno_labels[i]) & !is.na(resp_grp)
        n    <- sum(mask & resp_grp == cat)
        sprintf("%d (%.1f%%)", n, if (n_cat > 0) n/n_cat*100 else 0)
      })
    })
    names(cat_stats) <- cats
    stats_total <- sapply(seq_along(geno_labels), function(i) {
      mask <- sc %in% split_genos(geno_labels[i]) & !is.na(resp_grp)
      sprintf("%d", sum(mask))
    })
    list(s0 = stats_total, s1 = rep("", length(geno_labels)),
         by_cat = cat_stats, cats = cats)
  } else {
    stats0 <- character(length(geno_labels))
    for (i in seq_along(geno_labels)) {
      vals      <- as.numeric(response[sc %in% split_genos(geno_labels[i]) & !is.na(response)])
      stats0[i] <- if (length(vals) > 0) fmt_cont(vals) else "—"
    }
    list(s0 = stats0, s1 = rep("", length(geno_labels)))
  }
}


#' Compute association table results for one SNP.
#'
#' @param snp_nm  SNP name
#' @param prep    prep object from snp_prepare()
#' @param models  character vector of model names
#' @param ci_width  numeric CI width (e.g. 95)
#'
#' @return list:
#'   $snp          character
#'   $response_type character
#'   $col_titles   list(stat0, stat1, effect)
#'   $notes        character vector
#'   $rows         data.frame: model, genotype, stat0, stat1, effect, ci_low,
#'                             ci_high, pval, aic, bic
compute_assoc <- function(snp_nm, prep,
                          models   = c("codominant","dominant","recessive",
                                       "overdominant","logadditive"),
                          ci_width = 95) {
  sd           <- prep$snp_data[[snp_nm]]
  response_enc <- prep$response_enc[sd$snp_mask]
  response_raw <- if (!is.null(prep$response_raw)) prep$response_raw[sd$snp_mask] else NULL
  cov_df_cc    <- if (!is.null(prep$cov_df)) prep$cov_df[sd$snp_mask, , drop=FALSE] else NULL
  rtype        <- prep$response_type
  # Use the normalised genotype strings from geno_cc (het orientation is
  # ref-allele-first, consistent with user_levels and .geno_labels_for_model).
  # Using clean_cc here caused overdominant counts to be 0 for SNPs where the
  # raw data had het alleles in non-ref-first order (e.g. "G/A" vs "A/G").
  snp_char     <- as.character(sd$geno_cc)
  ref          <- sd$ref
  user_levels  <- sd$user_levels
  n_cov        <- if (!is.null(cov_df_cc)) ncol(cov_df_cc) else 0L
  n_miss       <- prep$n_rows - sd$n_typed

  all_genos <- if (!is.null(user_levels) && length(user_levels) > 0)
    user_levels
  else
    c(ref, setdiff(sort(unique(snp_char[!is.na(snp_char)])), ref))

  n_fit <- sum(!is.na(snp_char) & !is.na(response_enc) &
               (if (!is.null(cov_df_cc) && ncol(cov_df_cc)>0) complete.cases(cov_df_cc) else TRUE))

  model_labels <- c(codominant="Codominant", dominant="Dominant",
                    recessive="Recessive", overdominant="Overdominant",
                    logadditive="Log-additive")
  is_categorical <- rtype == "categorical"

  # column metadata
  col_titles <- list(
    effect = if (rtype %in% c("binary","categorical")) "OR" else "\u03B2",
    stat0  = switch(rtype,
      binary      = if (!is.null(response_raw)) levels(as.factor(response_raw))[1] else "Group 0",
      categorical = "N (%)",
      "Mean (SD)"),
    stat1  = switch(rtype,
      binary      = if (!is.null(response_raw)) levels(as.factor(response_raw))[2] else "Group 1",
      ""))

  notes <- character(0)
  if (!is.null(cov_df_cc) && ncol(cov_df_cc) > 0)
    notes <- c(notes, paste0("Adjusted for: ", paste(names(cov_df_cc), collapse=", ")))
  if (n_miss > 0)
    notes <- c(notes, paste0(n_miss, " observation(s) excluded."))
  if (is_categorical && !is.null(response_raw))
    notes <- c(notes, paste0("Reference category: \u2018",
      levels(as.factor(response_raw))[1], "\u2019"))

  all_rows  <- list()
  any_clamp <- FALSE

  for (mdl in models) {
    snp_enc  <- encode_model(snp_char, ref, mdl, user_levels)
    res_list <- fit_model(snp_enc, response_enc, cov_df_cc, mdl, rtype, ci_width)
    if (is.null(res_list)) next

    geno_labels <- .geno_labels_for_model(mdl, all_genos, ref)
    aic_val     <- { a <- res_list[[1]]$aic
                     if (!is.null(a) && !is.nan(a)) round(a, 2) else NA_real_ }
    bic_val     <- { b <- res_list[[1]]$bic
                     if (!is.null(b) && !is.nan(b)) round(b, 2) else NA_real_ }

    if (is_categorical) {
      cats <- unique(sapply(res_list, `[[`, "category"))
      st   <- .compute_stats(geno_labels, snp_char, response_enc, rtype, response_raw)
      for (cat in cats) {
        cat_res <- res_list[sapply(res_list, function(r) r$category == cat)]
        cat_sts <- if (!is.null(st$by_cat)) st$by_cat[[cat]] else rep("", length(geno_labels))
        n_cat   <- sum(as.character(response_raw) == cat, na.rm=TRUE)
        if (mdl == "logadditive") {
          res <- cat_res[[1]]
          all_rows[[length(all_rows)+1L]] <- data.frame(
            model    = paste0(model_labels[mdl], " \u2014 ", cat, " (n=", n_cat, ")"),
            genotype = "Per allele", stat0 = as.character(n_cat), stat1 = "",
            effect   = res$effect, ci_low = res$ci_low, ci_high = res$ci_high,
            pval     = res$pval,   aic    = aic_val,     bic    = bic_val,
            stringsAsFactors = FALSE)
          next
        }
        all_rows[[length(all_rows)+1L]] <- data.frame(
          model    = paste0(model_labels[mdl], " \u2014 ", cat, " (n=", n_cat, ")"),
          genotype = geno_labels[1],
          stat0    = if (length(cat_sts)>=1) cat_sts[1] else "",
          stat1    = "", effect=1., ci_low="", ci_high="",
          pval     = cat_res[[1]]$global_p, aic=aic_val, bic=bic_val,
          stringsAsFactors = FALSE)
        for (i in seq_along(cat_res)) {
          res <- cat_res[[i]]
          gl  <- if ((i+1)<=length(geno_labels)) geno_labels[i+1] else res$comparison
          if (is.na(res$effect)) any_clamp <- TRUE
          all_rows[[length(all_rows)+1L]] <- data.frame(
            model = "", genotype = gl,
            stat0 = if ((i+1)<=length(cat_sts)) cat_sts[i+1] else "",
            stat1 = "",
            effect=res$effect, ci_low=res$ci_low, ci_high=res$ci_high,
            pval=res$pval, aic="", bic="",
            stringsAsFactors = FALSE)
        }
      }
      next
    }

    # Binary / quantitative
    st <- .compute_stats(geno_labels, snp_char, response_enc, rtype, response_raw)
    if (mdl == "logadditive") {
      res <- res_list[[1]]
      if (rtype == "binary" && !is.null(response_raw)) {
        lv        <- levels(as.factor(response_raw))
        stat0_val <- sprintf("%d", sum(response_raw == lv[1], na.rm=TRUE))
        stat1_val <- sprintf("%d", sum(response_raw == lv[2], na.rm=TRUE))
      } else {
        stat0_val <- fmt_cont(as.numeric(response_enc))
        stat1_val <- ""
      }
      if (is.na(res$effect)) any_clamp <- TRUE
      all_rows[[length(all_rows)+1L]] <- data.frame(
        model="Log-additive", genotype="Per allele",
        stat0=stat0_val, stat1=stat1_val,
        effect=res$effect, ci_low=res$ci_low, ci_high=res$ci_high,
        pval=res$pval, aic=aic_val, bic=bic_val,
        stringsAsFactors = FALSE)
      next
    }
    pval_row1 <- if (mdl=="codominant") res_list[[1]]$global_p else NA_real_
    all_rows[[length(all_rows)+1L]] <- data.frame(
      model    = model_labels[mdl], genotype = geno_labels[1],
      stat0    = st$s0[1], stat1 = st$s1[1],
      effect   = if (rtype=="binary") 1. else 0.,
      ci_low   = NA_real_, ci_high = NA_real_,
      pval     = pval_row1, aic = aic_val, bic = bic_val,
      stringsAsFactors = FALSE)
    for (i in seq_along(res_list)) {
      res <- res_list[[i]]
      gl  <- if ((i+1)<=length(geno_labels)) geno_labels[i+1] else res$comparison
      if (is.na(res$effect)) any_clamp <- TRUE
      all_rows[[length(all_rows)+1L]] <- data.frame(
        model="", genotype=gl,
        stat0=if ((i+1)<=length(st$s0)) st$s0[i+1] else "—",
        stat1=if ((i+1)<=length(st$s1)) st$s1[i+1] else "",
        effect=res$effect, ci_low=res$ci_low, ci_high=res$ci_high,
        pval=res$pval, aic=NA_real_, bic=NA_real_,
        stringsAsFactors = FALSE)
    }
  }

  if (any_clamp)
    notes <- c(notes, "One or more OR/CI suppressed (shown as NA) due to complete or quasi-complete separation.")

  list(snp           = snp_nm,
       response_type = rtype,
       col_titles    = col_titles,
       notes         = notes,
       rows          = do.call(rbind, all_rows))
}


# ── Joint-mask geno_list helper ───────────────────────────────────────────────
#
# genetics::LD() and haplo.stats::setupGeno() both require every input vector
# to have the *same* length.  Each SNP's geno_cc is subset by its own per-SNP
# snp_mask (complete_mask & !is.na(clean)), so SNPs with different missingness
# patterns yield vectors of different lengths, causing silent NULL returns from
# LD() and "subscript too long" errors in subset_geno().
#
# This helper computes the row-wise intersection of all snp_masks (joint_mask),
# re-parses each SNP on that intersection, and returns equal-length genotype
# objects.  All three compute functions (compute_ld, compute_haplo_freq,
# compute_haplo_assoc) call it so the fix is in one place.
#
# @return list: $geno_list (named, equal-length), $joint_mask (logical, n_rows)
.make_joint_geno_list <- function(prep) {
  snp_vars   <- prep$snp_vars
  joint_mask <- Reduce(`&`, lapply(snp_vars,
                                   function(nm) prep$snp_data[[nm]]$snp_mask))
  geno_list  <- lapply(snp_vars, function(nm) {
    sd <- prep$snp_data[[nm]]
    parse_genotype(sd$clean[joint_mask], sd$user_levels)
  })
  names(geno_list) <- snp_vars
  geno_list <- Filter(Negate(is.null), geno_list)
  list(geno_list = geno_list, joint_mask = joint_mask)
}
# ─────────────────────────────────────────────────────────────────────────────


#' Compute pairwise LD for all SNP pairs.
#'
#' @return data.frame: snp1, snp2, r2, Dprime, D, pval
compute_ld <- function(prep, metric = "r2") {
  jg        <- .make_joint_geno_list(prep)
  geno_list <- jg$geno_list
  if (length(geno_list) < 2) return(NULL)

  nms   <- names(geno_list)
  pairs <- combn(nms, 2, simplify = FALSE)

  rows <- lapply(pairs, function(pair) {
    ld_res <- tryCatch(
      genetics::LD(geno_list[[pair[1]]], geno_list[[pair[2]]]),
      error = function(e) NULL)
    if (is.null(ld_res)) return(NULL)
    data.frame(
      snp1   = pair[1], snp2 = pair[2],
      r2     = round(ld_res$`r`^2, 3),
      Dprime = round(ld_res$`D'`,  3),
      D      = round(ld_res$`D`,   3),
      pval   = ld_res$`P-value`,
      stringsAsFactors = FALSE)
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}


#' Compute haplotype frequencies.
#'
#' @return list: $table data.frame (haplotype, freq, [group freqs]), $notes
compute_haplo_freq <- function(prep, subpop = FALSE, min_freq = 0.01) {
  jg        <- .make_joint_geno_list(prep)
  geno_list <- jg$geno_list
  if (length(geno_list) < 2) return(NULL)

  snp_names  <- names(geno_list)
  allele_mat <- do.call(cbind, lapply(snp_names, function(nm) genetics::allele(geno_list[[nm]])))
  geno_setup <- tryCatch(haplo.stats::setupGeno(allele_mat, locus.label=snp_names),
                         error = function(e) NULL)
  if (is.null(geno_setup)) return(NULL)

  u_alleles     <- attr(geno_setup, "unique.alleles")
  # geno_list is already on the joint mask so allele_mat has no cross-SNP
  # missingness; keep is all-TRUE.  response_raw is subset to joint_mask so
  # its length matches nrow(allele_mat).
  keep          <- rep(TRUE, nrow(allele_mat))
  response_raw  <- if (!is.null(prep$response_raw)) prep$response_raw[jg$joint_mask] else NULL

  do_strat   <- isTRUE(subpop) && !is.null(response_raw) && prep$response_type == "binary"
  grp_levels <- if (do_strat) levels(as.factor(response_raw)) else NULL

  em_all <- tryCatch(haplo.stats::haplo.em(subset_geno(geno_setup, keep),
                                           locus.label=snp_names),
                     error = function(e) NULL)
  if (is.null(em_all)) return(NULL)

  freqs      <- em_all$hap.prob
  grp_freq   <- list()
  if (do_strat) {
    for (lvl in grp_levels) {
      keep_lvl <- keep & !is.na(response_raw) & as.character(response_raw) == lvl
      if (sum(keep_lvl) < 5) next
      em_g <- tryCatch(haplo.stats::haplo.em(subset_geno(geno_setup, keep_lvl),
                                             locus.label=snp_names),
                       error=function(e) NULL)
      if (!is.null(em_g))
        grp_freq[[lvl]] <- setNames(
          as.list(round(em_g$hap.prob, 3)),
          sapply(seq_len(nrow(em_g$haplotype)), function(j)
            decode_haplo_row(as.numeric(em_g$haplotype[j,]), u_alleles)))
    }
  }

  sorted_idx <- order(freqs, decreasing=TRUE)
  rare_sum   <- 0
  rows <- lapply(sorted_idx, function(i) {
    if (freqs[i] < min_freq) { rare_sum <<- rare_sum + freqs[i]; return(NULL) }
    label    <- decode_haplo_row(as.numeric(em_all$haplotype[i,]), u_alleles)
    row      <- list(haplotype=label, freq=round(freqs[i],3))
    if (do_strat) for (lvl in grp_levels)
      row[[lvl]] <- grp_freq[[lvl]][[label]] %||% NA_real_
    as.data.frame(row, stringsAsFactors=FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (rare_sum > 0)
    rows[[length(rows)+1L]] <- as.data.frame(
      c(list(haplotype="Rare (combined)", freq=round(rare_sum,3)),
        if (do_strat) setNames(as.list(rep(NA_real_, length(grp_levels))), grp_levels)),
      stringsAsFactors=FALSE)

  list(table = do.call(rbind, rows),
       notes = paste0("Min frequency: ", min_freq,
                      ". N missing: ", sum(!jg$joint_mask & prep$complete_mask)))
}


#' Compute haplotype association.
#'
#' @return list: $table data.frame, $notes character
compute_haplo_assoc <- function(prep, ci_width = 95, min_freq = 0.01) {
  jg        <- .make_joint_geno_list(prep)
  geno_list <- jg$geno_list
  if (length(geno_list) < 2 || is.null(prep$response_enc)) return(NULL)

  snp_names  <- names(geno_list)
  allele_mat <- do.call(cbind, lapply(snp_names, function(nm) genetics::allele(geno_list[[nm]])))
  geno_setup <- tryCatch(haplo.stats::setupGeno(allele_mat, locus.label=snp_names),
                         error=function(e) NULL)
  if (is.null(geno_setup)) return(NULL)

  u_alleles <- attr(geno_setup, "unique.alleles")
  # geno_list is on the joint mask; allele_mat has no cross-SNP missingness.
  keep   <- rep(TRUE, nrow(allele_mat))
  n_miss <- sum(!jg$joint_mask & prep$complete_mask)

  response  <- prep$response_enc[jg$joint_mask]
  cov_df    <- if (!is.null(prep$cov_df) && nrow(prep$cov_df) > 0)
                 prep$cov_df[jg$joint_mask, , drop = FALSE] else prep$cov_df
  rtype     <- prep$response_type
  fam       <- if (rtype == "binary") binomial() else gaussian()
  cov_nms   <- if (!is.null(cov_df) && ncol(cov_df) > 0) names(cov_df) else character(0)

  em_all <- tryCatch(haplo.stats::haplo.em(subset_geno(geno_setup, keep),
                                           locus.label=snp_names),
                     error=function(e) NULL)
  if (is.null(em_all)) return(NULL)

  haplo_freq <- em_all$hap.prob
  haplo_names <- sapply(seq_len(nrow(em_all$haplotype)), function(j)
    decode_haplo_row(as.numeric(em_all$haplotype[j,]), u_alleles))
  rare_mask  <- haplo_freq < min_freq
  if (all(rare_mask)) return(NULL)

  fit <- tryCatch({
    haplo.stats::haplo.glm(
      response ~ 1,
      geno     = subset_geno(geno_setup, keep),
      family   = fam,
      haplo.freq.min = min_freq,
      control  = haplo.stats::haplo.glm.control(haplo.base = which(!rare_mask)[1]))
  }, error=function(e) NULL)
  if (is.null(fit)) return(NULL)

  coefs   <- summary(fit)$coefficients
  z_crit  <- qnorm(1 - (1 - ci_width/100)/2)
  hap_rows <- grep("^hap\\.", rownames(coefs))

  rows <- lapply(hap_rows, function(r) {
    nm   <- sub("^hap\\.", "", rownames(coefs)[r])
    beta <- coefs[r, "Estimate"]
    se   <- coefs[r, "Std. Error"]
    pval <- coefs[r, grep("Pr\\(", colnames(coefs))[1]]
    freq <- haplo_freq[match(nm, haplo_names)] %||% NA_real_
    data.frame(
      haplotype = nm, freq = round(freq, 3),
      effect    = if (rtype=="binary") .exp_or(beta) else beta,
      ci_low    = if (rtype=="binary") .exp_or(beta - z_crit*se) else beta - z_crit*se,
      ci_high   = if (rtype=="binary") .exp_or(beta + z_crit*se) else beta + z_crit*se,
      pval      = pval,
      stringsAsFactors = FALSE)
  })

  list(table = do.call(rbind, rows),
       notes = paste0(if (n_miss>0) paste0(n_miss," missing. "),
                      "Min frequency: ", min_freq))
}


