#' @importFrom R6 R6Class
#' @import jmvcore
#' @importFrom genetics genotype allele HWE.exact
source("R/snp_helpers.R")

snpDescClass <- if (requireNamespace("jmvcore", quietly = TRUE)) R6::R6Class(
  "snpDescClass",
  inherit = snpDescBase,
  private = list(

    .init = function() {
      self$results$covDescGroup$setVisible(FALSE)
      self$results$snpSummaryTablesGroup$setVisible(FALSE)
      snp_names <- self$options$snps
      if (length(snp_names) == 0) return()
      arr <- self$results$snpResults
      for (nm in snp_names) arr$addItem(key = nm)
    },

    .run = function() {
      data           <- self$data
      opts           <- self$options
      response_var   <- opts$response
      snp_vars       <- opts$snps
      covariate_vars <- opts$covariates

      run_snpSummary  <- isTRUE(opts$snpSummary)
      run_allFreq     <- isTRUE(opts$allFreq)
      run_genoFreq    <- isTRUE(opts$genoFreq)
      run_hweTest     <- isTRUE(opts$hweTest)
      run_subpop      <- isTRUE(opts$subpop)
      run_covDesc     <- isTRUE(opts$covDesc)
      run_showMissing <- isTRUE(opts$showMissing)
      run_rmSnpMissing <- isTRUE(opts$rmSnpMissing)

      # ── Validate SNPs ────────────────────────────────────────────
      if (length(snp_vars) == 0) {
        self$results$validationMsgSNP$setContent(
          "<p style='color:red;'>Please add at least one SNP variable.</p>")
        self$results$validationMsgSNP$setVisible(TRUE)
      } else {
        self$results$validationMsgSNP$setVisible(FALSE)
      }

      val <- validate_snp_vars(snp_vars, data)
      snp_vars <- val$valid_snps
      if (nchar(val$bad_html) > 0) {
        self$results$validationMsgGeno$setContent(val$bad_html)
        self$results$validationMsgGeno$setVisible(TRUE)
      } else {
        self$results$validationMsgGeno$setVisible(FALSE)
      }

      # ── Response / covariates ────────────────────────────────────
      response_raw <- if (!is.null(response_var) && response_var != "")
                        data[[response_var]] else NULL
      response_type <- detect_response_type(response_raw, opts$responseType)
      response      <- prepare_response(response_raw, response_type)

      cov_df        <- prepare_covariates(data, covariate_vars)
      # Create an empty but non-NULL cov_df if no covariates selected
      # This allows response to be shown even without covariates
      if (is.null(cov_df) && !is.null(response_raw)) {
        cov_df <- data.frame(row.names = seq_len(nrow(data)))
      }

      if (run_subpop && (is.null(response_raw) || response_type == "quantitative")) {
        run_subpop <- FALSE
      }

      # ── Visibility ───────────────────────────────────────────────
      self$results$covDescGroup$setVisible(run_covDesc && (!is.null(cov_df) || !is.null(response_raw)))
      self$results$snpSummaryTablesGroup$setVisible(run_snpSummary && length(snp_vars) > 0)

      # ── Covariate descriptives ───────────────────────────────────
      if (run_covDesc && !is.null(cov_df)) {
        private$.run_cov_desc(cov_df, response_raw, response_type, run_subpop, response_var, data, 
          snp_vars, run_rmSnpMissing)
      }

      # ── SNP summary table ────────────────────────────────────────
      if (run_snpSummary) {
        private$.fill_snp_summary(data, snp_vars, response_raw, response_type,
                                  run_subpop, cov_df)
      }

      # ── Complete-case mask (response + covariates) ───────────────
      n_rows        <- nrow(data)
      complete_mask <- rep(TRUE, n_rows)
      if (!is.null(response))              complete_mask <- complete_mask & !is.na(response)
      if (!is.null(cov_df) && ncol(cov_df) > 0)
        complete_mask <- complete_mask & complete.cases(cov_df)

      # ── Per-SNP descriptives ─────────────────────────────────────
      arr <- self$results$snpResults
      for (snp_nm in snp_vars) {
        snp_raw     <- data[[snp_nm]]

        user_levels <- get_snp_level_order(snp_raw)
        geno_obj    <- parse_genotype(snp_raw, user_levels)
        if (is.null(geno_obj)) next

        # For each SNP, we need complete cases: no missing in SNP, response, or covariates
        snp_complete_mask <- complete_mask & !is.na(snp_raw)
        
        # Count missing by category for reporting
        missing_snp_only <- is.na(snp_raw) & complete_mask
        missing_resp_cov <- !complete_mask & !is.na(snp_raw)
        missing_both <- !complete_mask & is.na(snp_raw)
        
        n_missing_snp_only <- sum(missing_snp_only)
        n_missing_resp_cov <- sum(missing_resp_cov)
        n_missing_both <- sum(missing_both)
        total_missing <- sum(!snp_complete_mask)
        
        # For stratified reporting, track missing by response level
        if (run_subpop && (response_type == "binary" || response_type == "categorical")) {
          resp_levels <- levels(as.factor(response_raw))
          n_miss_by_level <- sapply(resp_levels, function(lvl) {
            # Missing within this response level (including any missing reason)
            sum((is.na(snp_raw) | !complete_mask) & response_raw == lvl, na.rm=TRUE)
          })
          names(n_miss_by_level) <- resp_levels
        } else {
          n_miss_by_level <- NULL
        }

        snp_raw_cc  <- snp_raw[snp_complete_mask]
        geno_obj_cc <- parse_genotype(snp_raw_cc, user_levels)
        response_cc <- if (!is.null(response))     response[snp_complete_mask]     else NULL
        resp_raw_cc <- if (!is.null(response_raw)) response_raw[snp_complete_mask] else NULL
        if (is.null(geno_obj_cc)) next

        item    <- arr$get(key = snp_nm)
        n_typed <- sum(snp_complete_mask)
        n_total_eligible <- sum(complete_mask)  # cases with complete response+covariates
        
        item$typingRate$setContent(sprintf(
          "<b>Typed samples:</b> %d / %d (%.1f%%)", n_typed, n_total_eligible,
          if (n_total_eligible > 0) n_typed / n_total_eligible * 100 else 0))

        snp_summary_cc <- summary(geno_obj_cc)
        ref            <- get_ref_genotype(geno_obj_cc, user_levels)

        if (run_allFreq)
          private$.fill_allele_freq(item$allFreqTable, snp_summary_cc,
                                    snp_nm, resp_raw_cc, run_subpop,
                                    response_type, snp_raw_cc,
                                    run_showMissing, n_miss_by_level, 
                                    n_total_eligible, total_missing,
                                    user_levels = user_levels)
        if (run_genoFreq)
          private$.fill_geno_freq(item$genoFreqTable, snp_summary_cc, ref,
                                  snp_raw_cc, response_cc, response_type,
                                  run_subpop, resp_raw_cc,
                                  run_showMissing, n_miss_by_level, 
                                  n_total_eligible, total_missing,
                                  user_levels = user_levels)
        if (run_hweTest)
          private$.fill_hwe(item$hweTable, geno_obj_cc, snp_nm,
                            resp_raw_cc, run_subpop,
                            run_showMissing, n_miss_by_level, n_total_eligible, total_missing)
      }
    },

  .run_cov_desc = function(cov_df, response_raw, response_type,
    subpop = FALSE, response_var = NULL,
    data = NULL, snp_vars = NULL,
    rm_snp_missing = FALSE) {

      tbl <- self$results$covDescGroup$covDescTable

      has_covariates <- !is.null(cov_df) && ncol(cov_df) > 0

      n_removed_snp <- 0L
      if (isTRUE(rm_snp_missing) && !is.null(data) && !is.null(snp_vars) && length(snp_vars) > 0) {
        snp_cc <- complete.cases(data[, snp_vars, drop=FALSE])
        n_removed_snp <- sum(!snp_cc)
        if (n_removed_snp > 0) {
          cov_df <- cov_df[snp_cc, , drop=FALSE]
          if (!is.null(response_raw))
            response_raw <- response_raw[snp_cc]
        }
      }

      has_response <- !is.null(response_raw) && !is.null(response_type)
      is_binary <- has_response && response_type == "binary"
      is_cat_resp <- has_response && response_type == "categorical"
      is_cont_resp <- has_response && response_type == "quantitative"

      do_strat <- isTRUE(subpop) && (is_binary || is_cat_resp)

      valid_resp <- if (has_response) !is.na(response_raw) else rep(TRUE, nrow(cov_df))

      if (do_strat) {
        grp_fac  <- as.factor(response_raw)
        grp_lvls <- levels(grp_fac)

        mask_list <- lapply(grp_lvls, function(l) valid_resp & grp_fac == l)
        names(mask_list) <- grp_lvls
        totals <- sapply(mask_list, sum)

        for (i in seq_along(grp_lvls)) {
          nm <- paste0("stat_g", i-1)
          tbl$getColumn(nm)$setTitle(grp_lvls[i])
          tbl$getColumn(nm)$setVisible(TRUE)
        }
        tbl$getColumn("pval")$setVisible(TRUE)

        get_counts <- function(mask) sapply(mask_list, function(m) sum(mask & m))
      }

      has_cont <- FALSE

      # ---- Response block ----
      # ALWAYS show response description if response exists
      if (has_response) {

        if (is_cont_resp) {
          # For continuous response: show mean ± SD
          row_vals <- list(variable = response_var, level = "Mean \u00B1 SD",
                          stat_overall = fmt_cont(response_raw))
          tbl$addRow(rowKey = paste0(response_var, "_mean"), values = row_vals)

          # Show valid count
          mask <- !is.na(response_raw)
          row_vals <- list(variable = "", level = "Valid",
                          stat_overall = sum(mask))
          tbl$addRow(rowKey = paste0(response_var, "_valid"), values = row_vals)

          # Show missing count if any
          n_miss <- sum(!mask)
          if (n_miss > 0) {
            miss_vals <- list(variable = "", level = "Missing",
                              stat_overall = fmt_cat(n_miss, length(response_raw)))
            tbl$addRow(rowKey = paste0(response_var, "_missing"), values = miss_vals)
          }

        } else {
          # For categorical/binary response: show counts
          mask <- valid_resp
          row_vals <- list(variable = response_var, level = "Valid",
                          stat_overall = sum(mask))

          if (do_strat) {
            cnt <- get_counts(mask)
            for (i in seq_along(cnt))
              row_vals[[paste0("stat_g", i-1)]] <- fmt_cat(cnt[i], sum(mask))
          }

          tbl$addRow(rowKey = response_var, values = row_vals)

          n_miss <- sum(!valid_resp)
          if (n_miss > 0) {
            miss_vals <- list(variable = "", level = "Missing",
                            stat_overall = fmt_cat(n_miss, length(response_raw)))
            if (do_strat) {
              for (i in seq_along(grp_lvls)) {
                miss_vals[[paste0("stat_g", i-1)]] <- ''
              }
              miss_vals$pval <- ''
            }
            tbl$addRow(rowKey = paste0(response_var, "_missing"), values = miss_vals)
          }
        }

      }

      # ---- Covariates ----
      # Only show covariates if there are any
      if (!is.null(cov_df) && ncol(cov_df) > 0) {
        for (v in names(cov_df)) {

          col <- cov_df[[v]]
          is_cat <- is.factor(col) || is.character(col)
          if (is_cat && !is.factor(col)) col <- factor(col)

          n <- length(col)
          n_miss <- sum(is.na(col))

          if (is_cat) {

            lvls <- levels(col)

            pval <- if (do_strat) {
              tryCatch({
                ct <- table(col[valid_resp], grp_fac[valid_resp])
                suppressWarnings(chisq.test(ct)$p.value)
              }, error = function(e) NA_real_)
            } else NA_real_

            first <- TRUE
            for (lvl in lvls) {

              mask <- !is.na(col) & col == lvl

              row_vals <- list(
                variable = if (first) v else "",
                level = lvl,
                stat_overall = fmt_cat(sum(mask), n)
              )

              if (do_strat) {
                cnt <- get_counts(mask)
                for (i in seq_along(cnt))
                  row_vals[[paste0("stat_g", i-1)]] <- fmt_cat(cnt[i], totals[i])
                row_vals$pval <- if (first) pval else ''
              }

              tbl$addRow(rowKey = paste0(v, "_", lvl), values = row_vals)
              first <- FALSE
            }

            if (n_miss > 0) {
              mask <- is.na(col)
              row_vals <- list(variable = "", level = "Missing",
                              stat_overall = fmt_cat(n_miss, n))

              if (do_strat) {
                cnt <- get_counts(mask)
                for (i in seq_along(cnt))
                  row_vals[[paste0("stat_g", i-1)]] <- fmt_cat(cnt[i], totals[i])
                row_vals$pval <- ''
              }

              tbl$addRow(rowKey = paste0(v, "_missing"), values = row_vals)
            }

          } else {

            has_cont <- TRUE

            row_vals <- list(variable = v, level = "Mean \u00B1 SD",
                            stat_overall = fmt_cont(col))

            if (do_strat) {
              for (i in seq_along(mask_list)) {
                vals <- col[mask_list[[i]]]
                row_vals[[paste0("stat_g", i-1)]] <- fmt_cont(vals)
              }
              row_vals$pval <- tryCatch({
                groups <- split(col[valid_resp], grp_fac[valid_resp])
                if (length(groups) == 2)
                  t.test(groups[[1]], groups[[2]])$p.value
                else
                  summary(aov(col ~ grp_fac))[[1]][["Pr(>F)"]][1]
              }, error=function(e) NA_real_)
            }

            tbl$addRow(rowKey = v, values = row_vals)

            if (n_miss > 0) {
              mask <- is.na(col)
              row_vals <- list(variable = "", level = "Missing",
                              stat_overall = fmt_cat(n_miss, n))

              if (do_strat) {
                cnt <- get_counts(mask)
                for (i in seq_along(cnt))
                  row_vals[[paste0("stat_g", i-1)]] <- fmt_cat(cnt[i], totals[i])
                row_vals$pval <- ''
              }

              tbl$addRow(rowKey = paste0(v, "_missing"), values = row_vals)
            }
          }
        }
      }

      # ---- Notes ----
      notes <- c()

      if (n_removed_snp > 0)
        notes <- c(notes, paste0("Removed ", n_removed_snp, " cases with missing SNP values."))

      tbl$setNote(note = if (length(notes)) paste(notes, collapse=" ") else NULL,
                  key  = "cov_desc_note")
    },

    .fill_snp_summary = function(data, snp_vars, response_raw, response_type,
                             subpop, cov_df = NULL) {

      tbl <- self$results$snpSummaryTablesGroup$snpSummaryTable

      do_strat <- isTRUE(subpop) && !is.null(response_raw) && (
        response_type == "binary" || response_type == "categorical")
      has_response <- !is.null(response_raw)
      has_cov <- !is.null(cov_df) && ncol(cov_df) > 0
      n_total <- nrow(data)

      grp_levels <- if (do_strat) levels(response_raw) else NULL
      tbl$getColumn("group")$setVisible(do_strat)

      base_cc <- rep(TRUE, n_total)
      if (has_response) base_cc <- base_cc & !is.na(response_raw)
      if (has_cov) base_cc <- base_cc & complete.cases(cov_df)

      if (do_strat) {
        resp_base <- response_raw[base_cc]
        stratum_totals <- table(factor(resp_base, levels = grp_levels))
      }

      row_key <- 0L

      for (snp_nm in snp_vars) {

        snp_raw <- data[[snp_nm]]
        user_levels_sum <- get_snp_level_order(snp_raw)
        geno_obj <- parse_genotype(snp_raw, user_levels_sum)
        if (is.null(geno_obj)) next

        cc_mask <- base_cc & !is.na(snp_raw)
        n_cc <- sum(cc_mask)
        n_excluded <- n_total - n_cc

        snp_cc <- snp_raw[cc_mask]
        geno_cc <- parse_genotype(snp_cc, user_levels_sum)
        if (is.null(geno_cc)) next

        resp_cc <- if (has_response) response_raw[cc_mask] else NULL

        sm_cc <- summary(geno_cc)
        ref <- get_ref_genotype(geno_cc, user_levels_sum)

        af_all <- sm_cc$allele.freq
        allele_nms <- rownames(af_all)
        ref_allele <- strsplit(ref, "/")[[1]][1]
        alt_allele <- setdiff(allele_nms, ref_allele)
        alt_allele <- if (length(alt_allele)) alt_allele[1] else "?"
        alleles_label <- paste0(ref_allele, "/", alt_allele)

        compute_row <- function(mask = NULL) {
          snp_sub <- if (is.null(mask)) snp_cc else snp_cc[mask]
          geno_sub <- if (is.null(mask)) geno_cc else parse_genotype(snp_sub, user_levels_sum)
          if (is.null(geno_sub)) return(NULL)

          sm <- summary(geno_sub)
          af <- sm$allele.freq
          props <- af[, "Proportion"]

          maf <- if (alt_allele %in% rownames(af)) {
            af[alt_allele, "Proportion"]
          } else if (length(props) >= 2) {
            min(props, na.rm = TRUE)
          } else NA_real_

          gf <- sm$genotype.freq
          gf <- tryCatch(reorder_geno(gf, ref, user_levels_sum), error = function(e) gf)
          gf <- gf[rownames(gf) != "NA", , drop = FALSE]

          counts <- as.integer(gf[, "Count"])
          len <- length(counts)
          geno_str <- if (len == 3) paste(counts, collapse = " / ")
                      else if (len == 2) paste(c(counts, 0L), collapse = " / ")
                      else paste(counts, collapse = " / ")

          hwe <- tryCatch(genetics::HWE.exact(geno_sub)$p.value, error = function(e) NA_real_)

          list(n = sm$n.typed, maf = maf, genoCounts = geno_str, hwePval = hwe)
        }

        res_all <- compute_row()
        if (!is.null(res_all)) {
          row_key <- row_key + 1L
          tbl$addRow(rowKey = as.character(row_key), values = list(
            snp = snp_nm,
            alleles = alleles_label,
            group = if (do_strat) "All" else "",
            n = res_all$n,
            missing = if (n_excluded > 0L) n_excluded else '',
            maf = round(res_all$maf, 4),
            genoCounts = res_all$genoCounts,
            hwePval = res_all$hwePval
          ))
        }

        if (do_strat) {
          resp_cc_chr <- as.character(resp_cc)

          for (lvl in grp_levels) {
            mask <- !is.na(resp_cc_chr) & resp_cc_chr == lvl
            res <- compute_row(mask)
            if (is.null(res)) next

            n_excl <- max(0L, as.integer(stratum_totals[lvl] - res$n))

            row_key <- row_key + 1L
            tbl$addRow(rowKey = as.character(row_key), values = list(
              snp = "",
              alleles = "",
              group = lvl,
              n = res$n,
              missing = if (n_excl > 0L) n_excl else '',
              maf = round(res$maf, 4),
              genoCounts = res$genoCounts,
              hwePval = res$hwePval
            ))
          }
        }
      }

      if (has_response || has_cov) {
        parts <- c(if (has_cov) "covariates", if (has_response) "response")
        tbl$setNote(
          note = paste0("Complete cases used: rows missing any ",
                        paste(parts, collapse = " or "),
                        " or SNP value are excluded."),
          key = "missing_resp_cov"
        )
      } else {
        tbl$setNote(note = NULL, key = "missing_resp_cov")
      }
    },

    .fill_allele_freq = function(tbl, sm, snp_nm, response_raw, subpop,
                                  response_type, snp_raw, show_missing=FALSE,
                                  n_miss_by_level=NULL, n_total_eligible=0L, 
                                  total_missing=0L, user_levels=NULL) {

      af <- sm$allele.freq
      allele_names <- rownames(af)

      if (!is.null(user_levels) && length(user_levels)) {
        ref <- NULL
        for (g in user_levels) {
          p <- strsplit(g, "/", fixed=TRUE)[[1]]
          if (length(p)==2 && p[1]==p[2]) { ref <- p[1]; break }
        }
        if (!is.null(ref) && ref %in% allele_names)
          allele_names <- c(ref, setdiff(allele_names, ref))
      }

      do_strat <- isTRUE(subpop) && !is.null(response_raw) && (
        response_type == "binary" || response_type == "categorical")

      if (do_strat) {
        grp_levels <- levels(response_raw)
        resp_chr <- as.character(response_raw)
        
        # Pre-compute allele counts by stratum
        alleles_split <- strsplit(as.character(snp_raw), "/")
        
        for (i in seq_along(grp_levels)) {
          tbl$addColumn(name=paste0("stat_g", i-1),
                        title=grp_levels[i],
                        type="string")
        }
      }

      # Add allele frequency rows
      for (al in allele_names) {
        if (!al %in% rownames(af)) next

        count <- as.integer(af[al, "Count"])
        prop <- round(af[al, "Proportion"] * 100, 1)
        
        row_vals <- list(
          allele = al,
          stat = fmt_catpct(count, prop)
        )

        if (do_strat) {
          for (i in seq_along(grp_levels)) {
            lvl <- grp_levels[i]
            idx <- resp_chr == lvl
            all_al <- unlist(alleles_split[idx])
            n_al <- sum(all_al == al, na.rm=TRUE)
            n_tot <- length(all_al)
            row_vals[[paste0("stat_g", i-1)]] <- fmt_cat(n_al, n_tot)
          }
        }

        tbl$addRow(rowKey=al, values=row_vals)
      }

      # Add missing row if requested
      if (isTRUE(show_missing) && total_missing > 0L) {  
        miss_vals <- list(
          allele = "Missing", 
          stat = fmt_cat(total_missing, n_total_eligible + total_missing)
        )
        
        if (do_strat && !is.null(n_miss_by_level)) {
          for (j in seq_along(grp_levels)) {
            lvl <- grp_levels[j]
            miss_count <- if (lvl %in% names(n_miss_by_level)) n_miss_by_level[lvl] else 0
            miss_vals[[paste0("stat_g", j-1)]] <- fmt_catn(miss_count)
          }
        }
        tbl$addRow(rowKey="missing", values=miss_vals)
      }
    }, 

    .fill_geno_freq = function(tbl, sm, ref, snp_raw, response,
                                response_type, subpop, response_raw,
                                show_missing=FALSE, n_miss_by_level=NULL, 
                                n_total_eligible=0L, total_missing=0L,
                                user_levels=NULL) {

      if (response_type == "quantitative") {
        tbl$getColumn("responseStat")$setVisible(TRUE)
        if (!is.numeric(response)) response <- as.numeric(as.character(response))
      }

      gf <- tryCatch(reorder_geno(sm$genotype.freq, ref, user_levels), 
                     error=function(e) sm$genotype.freq)
      gf <- gf[rownames(gf) != "NA", , drop = FALSE]

      do_strat <- isTRUE(subpop) && !is.null(response_raw) && (
        response_type == "binary" || response_type == "categorical")

      if (do_strat) {
        grp_levels <- levels(response_raw)
        resp_chr <- as.character(response_raw)

        for (i in seq_along(grp_levels)) {
          tbl$addColumn(name=paste0("stat_g", i-1),
                        title=grp_levels[i],
                        type="string")
        }
      }

      # Add genotype frequency rows
      for (i in seq_len(nrow(gf))) {
        geno <- rownames(gf)[i]
        if (geno == "NA") next
     
        count <- as.integer(gf[i, "Count"])
        prop <- gf[i, "Proportion"] * 100
        
        row_vals <- list(
          genotype = geno, 
          stat = fmt_catpct(count, prop),
          responseStat = ""
        )

        # Add response statistics for quantitative traits
        if (response_type == "quantitative" && !is.null(response)) {
          mask <- snp_raw == geno & !is.na(response)
          n_mask <- sum(mask)
          if (n_mask > 0) {
            mn <- mean(response[mask], na.rm=TRUE)
            se <- sd(response[mask], na.rm=TRUE) / sqrt(n_mask)
            row_vals$responseStat <- sprintf("%.2f (%.2f)", mn, se)
          }
        }

        if (do_strat) {
          for (j in seq_along(grp_levels)) {
            lvl <- grp_levels[j]
            idx <- resp_chr == lvl
            n_g <- sum(idx & snp_raw == geno)
            n_tot <- sum(idx)
            row_vals[[paste0("stat_g", j-1)]] <- fmt_cat(n_g, n_tot)
          }
        }

        tbl$addRow(rowKey=geno, values=row_vals)
      }

      # Add missing row if requested
      if (isTRUE(show_missing) && total_missing > 0) {
        miss_vals <- list(
          genotype = "Missing", 
          stat = fmt_cat(total_missing, n_total_eligible + total_missing),
          responseStat = ""
        )
        
        if (do_strat && !is.null(n_miss_by_level)) {
          for (j in seq_along(grp_levels)) {
            lvl <- grp_levels[j]
            miss_count <- if (lvl %in% names(n_miss_by_level)) n_miss_by_level[lvl] else 0
            miss_vals[[paste0("stat_g", j-1)]] <- fmt_catn(miss_count)
          }
        }
        tbl$addRow(rowKey="missing", values=miss_vals)
      }
    },

    .fill_hwe = function(tbl, geno_obj, snp_nm, response_raw, subpop,
                          show_missing=FALSE, n_miss_by_level=NULL, 
                          n_total_eligible=0L, total_missing=0L) {

      tbl$getColumn("missing")$setVisible(isTRUE(show_missing))

      hw <- tryCatch(genetics::HWE.exact(geno_obj), error=function(e) NULL)
      if (is.null(hw)) return()

      add_row <- function(key, label, st, miss, p) {
        tbl$addRow(rowKey=key, values=list(
          group=label,
          n11=as.integer(st["N11"]),
          n12=as.integer(st["N12"]),
          n22=as.integer(st["N22"]),
          missing=miss,
          pval=p
        ))
      }

      st <- hw$statistic
      add_row("All", "All subjects", st,
              if (isTRUE(show_missing)) total_missing else '',
              hw$p.value)

      if (isTRUE(subpop) && !is.null(response_raw)) {
        lvls <- levels(response_raw)
        if (length(lvls) <= 5) {
          for (lvl in lvls) {
            mask <- response_raw == lvl & !is.na(response_raw)
            if (sum(mask) == 0) next
            hw_sub <- tryCatch(genetics::HWE.exact(geno_obj[mask]), error=function(e) NULL)
            if (is.null(hw_sub)) next
            st2 <- hw_sub$statistic
            miss_count <- if (isTRUE(show_missing) && !is.null(n_miss_by_level) && 
                            lvl %in% names(n_miss_by_level)) n_miss_by_level[lvl] else ''
            add_row(lvl, lvl, st2, miss_count, hw_sub$p.value)
          }
        }
      }
    }
  )
)