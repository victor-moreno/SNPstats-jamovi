# ── snp_compute.R ─────────────────────────────────────────────────────────────
#
# Pure-R computation layer for SNPstats.
# No jamovi dependencies — every function takes plain R objects and returns
# plain R lists / data frames.  Both the jamovi class (snpStats.b.R) and the
# exported snpStats() function call these functions so results are identical
# in both contexts.
#
# Sections (in dependency order):
#   1. Constants and operators
#   2. Formatting
#   3. Genotype parsing
#   4. SNP validation
#   5. Response / covariate preparation
#   6. Genetic model encoding
#   7. Formula and model fitting
#   8. Data preparation  (snp_prepare)
#   9. Descriptive statistics
#  10. Association analysis
#  11. Haplotype analysis
# ──────────────────────────────────────────────────────────────────────────────

#' @importFrom genetics genotype allele HWE.exact LD
#' @importFrom R6 R6Class
#' @import jmvcore


# ══════════════════════════════════════════════════════════════════════════════
# 1. Constants and operators
# ══════════════════════════════════════════════════════════════════════════════

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Regex that matches null-allele genotype codings (0/0, 0|0, 0>0, 00).
.NULL_ALLELE_PAT <- "^0[/|>]0$|^00$"

# Canonical display labels for the five genetic models.
.MODEL_LABELS <- c(codominant   = "Codominant",
                   dominant     = "Dominant",
                   recessive    = "Recessive",
                   overdominant = "Overdominant",
                   logadditive  = "Additive")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Formatting
# ══════════════════════════════════════════════════════════════════════════════

#' Round a p-value to a display string.
#' Values < 0.001 are shown as "< 0.001"; otherwise 3 significant figures.
fmt_pval <- function(x) {
  if (is.null(x) || length(x) == 0) return('')
  vapply(x, function(p) {
    if (is.na(p))       return('')
    if (!is.numeric(p)) return(as.character(p))
    if (p < 0.001)      return("< 0.001")
    format.pval(p, digits = 3, eps = 0.001, nsmall = 3, scientific = FALSE)
  }, '')
}

#' Round an effect / CI value to 3 decimal places.
fmt3 <- function(x) {
  # "" arises when rbind() coerces mixed numeric/character columns (e.g. ref-row ci_low="")
  if (is.null(x) || length(x) == 0 || is.na(x) || identical(x, "")) return('')
  formatC(round(as.numeric(x), 3), format = "f", flag = "#", digits = 3)
}

# Categorical: N (%)
fmt_cat   <- function(n, total) sprintf("%d (%.1f%%)", n, if (total > 0) 100*n/total else 0)
fmt_catpct <- function(n, pct)  sprintf("%d (%.1f%%)", n, pct)

# Continuous: mean ± SD. Empty / all-missing cells render as an em dash (—) so
# every table treats a zero-observation cell the same way (the stratified-by-
# genotype table used to call this directly and show a bare "NA"); a single
# observation keeps "mean ± NA" since the SD is undefined for n = 1.
fmt_cont <- function(x) {
  if (all(is.na(x))) return("—")
  sprintf("%.2f ± %.2f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

#' Sample skewness (bias-corrected).
skewness <- function(x) {
  x <- x[!is.na(x)]; n <- length(x)
  if (n < 3) return(NA_real_)
  x  <- x - mean(x)
  s2 <- sum(x ^ 2)
  if (s2 == 0) return(NA_real_)            # constant vector: skewness undefined
  y <- sqrt(n) * sum(x ^ 3) / (s2 ^ (3/2))
  y * ((1 - 1 / n)) ^ (3/2)
}


# ══════════════════════════════════════════════════════════════════════════════
# 3. Genotype parsing
# ══════════════════════════════════════════════════════════════════════════════

#' Split a normalised "A/B" genotype string into its two alleles.
split_alleles <- function(g) strsplit(g, "/", fixed = TRUE)[[1]]

#' Check that a vector of unique "A/B" genotype strings is biallelic.
#' Returns list($ok, $reason, $alleles).
check_biallelic <- function(vals) {
  pairs <- lapply(vals, split_alleles)
  bad   <- which(sapply(pairs, length) != 2)
  if (length(bad) > 0)
    return(list(ok = FALSE,
                reason = paste0("cannot split into two alleles: ",
                                paste(vals[bad], collapse = ", "))))
  alleles <- unique(unlist(pairs))
  if (length(alleles) > 2)
    return(list(ok = FALSE,
                reason = paste0("more than 2 alleles found (",
                                paste(sort(alleles), collapse = ", "),
                                "); only biallelic SNPs are supported")))
  a <- alleles[1]; b <- if (length(alleles) == 2) alleles[2] else alleles[1]
  valid    <- c(paste0(a,"/",a), paste0(a,"/",b),
                paste0(b,"/",a), paste0(b,"/",b))
  bad_geno <- vals[!vals %in% valid]
  if (length(bad_geno) > 0)
    return(list(ok = FALSE,
                reason = paste0("unexpected genotype(s): ",
                                paste(bad_geno, collapse = ", "))))
  list(ok = TRUE, reason = NULL, alleles = alleles)
}

#' Replace null-allele coded genotypes (0/0, 0|0, 0>0, 00) with NA.
clean_null_alleles <- function(x) {
  x_chr <- as.character(x)
  is_null <- !is.na(x_chr) & grepl(.NULL_ALLELE_PAT, x_chr, ignore.case = TRUE)
  x_chr[is_null] <- NA_character_
  x_chr
}

#' Detect genotype separator; returns NULL if not a valid biallelic SNP column.
detect_snp_sep <- function(x) {
  vals <- unique(na.omit(clean_null_alleles(as.character(x))))
  if (length(vals) == 0 || length(vals) > 10) return(NULL)
  for (sep in c("/", "|", ">")) {
    pat  <- paste0("^.+", if (sep == "|") "\\|" else sep, ".+$")
    if (all(grepl(pat, vals))) {
      norm <- if (sep == "/") vals else sub(sep, "/", vals, fixed = TRUE)
      if (!isTRUE(check_biallelic(norm)$ok)) return(NULL)
      return(sep)
    }
  }
  if (all(nchar(vals) == 2) && all(grepl("^[A-Za-z0-9]{2}$", vals))) {
    norm <- paste0(substr(vals, 1, 1), "/", substr(vals, 2, 2))
    if (!isTRUE(check_biallelic(norm)$ok)) return(NULL)
    return("")
  }
  NULL
}

#' Full biallelic check returning a reason string; used for validation messages.
snp_biallelic_check <- function(x) {
  vals <- unique(na.omit(clean_null_alleles(as.character(x))))
  sep  <- NULL
  for (s in c("/", "|", ">")) {
    pat <- paste0("^.+", if (s == "|") "\\|" else s, ".+$")
    if (all(grepl(pat, vals))) { sep <- s; break }
  }
  if (is.null(sep) && all(nchar(vals) == 2) &&
      all(grepl("^[A-Za-z0-9]{2}$", vals))) sep <- ""
  if (is.null(sep)) return(list(ok = FALSE, reason = "unrecognised format"))
  norm <- if (sep == "") paste0(substr(vals,1,1),"/",substr(vals,2,2)) else
          if (sep == "/") vals else sub(sep, "/", vals, fixed = TRUE)
  check_biallelic(norm)
}

#' Extract the user-defined genotype level order from a jamovi factor column.
get_snp_level_order <- function(x) {
  if (!is.factor(x)) return(NULL)
  lvls <- levels(x)
  if (length(lvls) == 0) return(NULL)
  lvls <- lvls[!grepl(.NULL_ALLELE_PAT, lvls, ignore.case = TRUE)]
  if (length(lvls) == 0) return(NULL)

  norm <- lvls
  for (sep in c("|", ">")) {
    pat <- paste0("^.+", if (sep == "|") "\\|" else sep, ".+$")
    if (all(grepl(pat, norm))) {
      norm <- sub(sep, "/", norm, fixed = TRUE); break
    }
  }
  if (all(nchar(norm) == 2) && all(grepl("^[A-Za-z0-9]{2}$", norm)))
    norm <- paste0(substr(norm, 1, 1), "/", substr(norm, 2, 2))
  if (!all(grepl("^.+/.+$", norm))) return(NULL)

  # Alphabetical order = auto-generated levels (Jamovi default), not user intent.
  if (identical(norm, sort(norm))) return(NULL)

  ref_allele <- NULL
  for (g in norm) {
    parts <- strsplit(g, "/", fixed = TRUE)[[1]]
    if (length(parts) == 2 && parts[1] == parts[2]) { ref_allele <- parts[1]; break }
  }
  if (!is.null(ref_allele)) {
    norm <- sapply(norm, function(g) {
      parts <- strsplit(g, "/", fixed = TRUE)[[1]]
      if (length(parts) == 2 && parts[1] != parts[2] && parts[2] == ref_allele)
        paste0(parts[2], "/", parts[1])
      else g
    }, USE.NAMES = FALSE)
  }

  is_het <- function(g) {
    parts <- strsplit(g, "/", fixed = TRUE)[[1]]
    length(parts) == 2 && parts[1] != parts[2]
  }
  if (length(norm) == 3) {
    is_het_vec <- sapply(norm, is_het)
    hom_levels <- norm[!is_het_vec]; het_levels <- norm[is_het_vec]
    if (length(het_levels) == 1 && length(hom_levels) == 2) {
      ref_hom <- if (!is.null(ref_allele))
        hom_levels[sapply(hom_levels, function(g) {
          parts <- strsplit(g, "/", fixed = TRUE)[[1]]
          length(parts) == 2 && parts[1] == ref_allele })]
      else hom_levels[1]
      alt_hom <- hom_levels[hom_levels != ref_hom[1]]
      if (length(ref_hom) == 1 && length(alt_hom) == 1)
        norm <- c(ref_hom, het_levels, alt_hom)
    }
  } else {
    if (length(norm) >= 2 && is_het(norm[1])) norm[c(1, 2)] <- norm[c(2, 1)]
  }
  norm
}

#' Normalise and parse a raw genotype vector via genetics::genotype().
parse_genotype <- function(x, user_levels = NULL) {
  x_chr <- clean_null_alleles(as.character(x))
  sep <- detect_snp_sep(x_chr)
  if (is.null(sep)) return(NULL)

  if (sep == "") {
    x_norm <- ifelse(is.na(x_chr), NA_character_,
                     paste0(substr(x_chr,1,1), "/", substr(x_chr,2,2)))
  } else if (sep == "/") {
    x_norm <- x_chr
  } else {
    x_norm <- ifelse(is.na(x_chr), NA_character_,
                     sub(sep, "/", x_chr, fixed = TRUE))
  }

  ref_allele <- NULL
  if (!is.null(user_levels)) {
    for (g in user_levels) {
      parts <- strsplit(g, "/", fixed = TRUE)[[1]]
      if (length(parts) == 2 && parts[1] == parts[2]) { ref_allele <- parts[1]; break }
    }
  }
  if (is.null(ref_allele)) {
    all_pairs   <- strsplit(as.character(na.omit(x_norm)), "/", fixed = TRUE)
    all_alleles <- unique(unlist(all_pairs))
    for (pr in all_pairs) {
      if (length(pr) == 2 && pr[1] == pr[2]) { ref_allele <- pr[1]; break }
    }
    if (is.null(ref_allele) && length(all_alleles) >= 1) ref_allele <- all_alleles[1]
  }
  if (!is.null(ref_allele)) {
    # Split every genotype once (was a per-row strsplit inside sapply). Same
    # rule: a heterozygote whose second allele is the reference is reordered
    # ref-first; everything else (incl. malformed / non-2-allele) is left as the
    # original string, matching the previous `else g` exactly.
    nn         <- which(!is.na(x_norm))
    gs         <- x_norm[nn]
    parts_list <- strsplit(gs, "/", fixed = TRUE)
    x_norm[nn] <- vapply(seq_along(gs), function(i) {
      parts <- parts_list[[i]]
      if (length(parts) == 2 && parts[1] != parts[2] && parts[2] == ref_allele)
        paste0(parts[2], "/", parts[1])
      else
        gs[i]
    }, character(1))
  }
  genetics::genotype(x_norm, sep = "/")
}

#' Determine reference genotype (user-specified first, then most-frequent homozygote).
get_ref_genotype <- function(geno, user_levels = NULL) {
  if (is.null(geno)) return(NULL)
  # A monomorphic SNP collapses genotype.freq to an unlabelled named vector
  # (c(Count, Proportion)) with no row names, so there is nothing to rank —
  # the single observed genotype is the reference.
  obs_genos <- unique(as.character(geno)[!is.na(as.character(geno))])
  if (!is.null(user_levels) && length(user_levels) > 0) {
    sm    <- summary(geno)
    gf0   <- sm$genotype.freq
    obs   <- if (is.matrix(gf0)) rownames(gf0)[rownames(gf0) != "NA"] else obs_genos
    for (lvl in user_levels) if (lvl %in% obs) return(lvl)
  }
  sm      <- summary(geno)
  gf      <- sm$genotype.freq
  if (!is.matrix(gf))
    return(if (length(obs_genos) >= 1L) obs_genos[1L] else NULL)
  alleles <- rownames(gf)
  is_hom  <- sapply(alleles, function(g) {
    p <- strsplit(g, "/")[[1]]; length(p) == 2 && p[1] == p[2] })
  homz_gf <- gf[is_hom, , drop = FALSE]
  if (nrow(homz_gf) == 0) return(alleles[1])
  rownames(homz_gf)[which.max(homz_gf[, "Count"])]
}

#' Reorder genotype frequency table: user order first, then ref/het/alt fallback.
reorder_geno <- function(gf, ref, user_levels = NULL) {
  alleles <- rownames(gf)
  na_row  <- alleles == "NA"
  other   <- alleles[!na_row]

  canon_key <- function(g) {
    parts <- strsplit(g, "/", fixed = TRUE)[[1]]
    if (length(parts) == 2) paste(sort(parts), collapse = "/") else g
  }

  if (!is.null(user_levels) && length(user_levels) > 0) {
    other_keys      <- setNames(sapply(other, canon_key), other)
    canon_to_actual <- setNames(names(other_keys), other_keys)
    ordered <- character(0)
    for (ul in user_levels) {
      ck <- canon_key(ul)
      if (ck %in% names(canon_to_actual)) {
        actual_nm <- canon_to_actual[[ck]]
        if (!actual_nm %in% ordered) ordered <- c(ordered, actual_nm)
      }
    }
    ordered <- c(ordered, other[!other %in% ordered])
  } else {
    is_hom  <- sapply(other, function(g) { p <- strsplit(g,"/")[[1]]; length(p)==2 && p[1]==p[2] })
    ordered <- c(ref,
                 other[!is_hom & other != ref],
                 other[ is_hom & other != ref])
    ordered <- unique(ordered[ordered %in% other])
  }

  final <- c(ordered, alleles[na_row])
  out   <- gf[final[final %in% alleles], , drop = FALSE]
  # genetics::summary()$genotype.freq can collapse to a named vector when only
  # one genotype is observed in a small stratum; ensure always a matrix.
  if (!is.matrix(out)) {
    out <- matrix(out, nrow = 1L,
                  dimnames = list(final[final %in% alleles], names(out)))
  }
  out
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. SNP validation
# ══════════════════════════════════════════════════════════════════════════════

#' Validate SNP variables; return $valid_snps and $bad_html.
validate_snp_vars <- function(snp_vars, data) {
  bad_snps <- character(0); bad_msgs <- character(0)
  for (v in snp_vars) {
    col <- data[[v]]
    if (all(is.na(clean_null_alleles(as.character(col))))) next
    chk <- snp_biallelic_check(col)
    if (!isTRUE(chk$ok)) {
      bad_snps <- c(bad_snps, v)
      bad_msgs <- c(bad_msgs, paste0("<b>", v, "</b>: ", chk$reason))
    }
  }
  html <- if (length(bad_snps) > 0)
    paste0("<p style='color:red;'>The following SNP columns were skipped ",
           "(accepted formats: A/B, A|B, A>B, or AB; exactly 2 alleles required):</p>",
           "<ul>", paste0("<li>", bad_msgs, "</li>", collapse = ""), "</ul>")
  else ""
  list(valid_snps = setdiff(snp_vars, bad_snps), bad_html = html)
}

#' Compute effect-allele frequency and exact HWE p-value for one SNP column.
snp_af_hwe <- function(col, effect_allele = NULL) {
  out <- list(effect_af = NA_real_, hwe_p = NA_real_)

  if (is.numeric(col)) {
    vals <- col[!is.na(col)]
    if (length(vals) == 0 || !all(vals %in% c(0, 1, 2))) return(out)
    out$effect_af <- mean(vals) / 2
    n0 <- sum(vals == 0); n1 <- sum(vals == 1); n2 <- sum(vals == 2)
    geno_obj <- tryCatch(
      genetics::genotype(rep(c("A/A", "A/B", "B/B"), c(n0, n1, n2)), sep = "/"),
      error = function(e) NULL)
    if (!is.null(geno_obj)) {
      hw <- tryCatch(genetics::HWE.exact(geno_obj), error = function(e) NULL)
      if (!is.null(hw)) out$hwe_p <- hw$p.value
    }
    return(out)
  }

  user_levels <- get_snp_level_order(col)
  geno_obj    <- tryCatch(parse_genotype(col, user_levels), error = function(e) NULL)
  if (is.null(geno_obj)) return(out)
  sm <- tryCatch(summary(geno_obj), error = function(e) NULL)
  if (is.null(sm)) return(out)

  af           <- sm$allele.freq
  allele_names <- rownames(af)
  ea <- toupper(trimws(as.character(effect_allele %||% "")))
  if (nchar(ea) == 0 || !ea %in% allele_names) {
    props <- af[allele_names != "NA", "Proportion", drop = TRUE]
    ea    <- allele_names[allele_names != "NA"][which.min(props)]
  }
  if (ea %in% allele_names)
    out$effect_af <- af[ea, "Proportion"]

  hw <- tryCatch(genetics::HWE.exact(geno_obj), error = function(e) NULL)
  if (!is.null(hw)) out$hwe_p <- hw$p.value
  out
}


# ══════════════════════════════════════════════════════════════════════════════
# 5. Response and covariate preparation
# ══════════════════════════════════════════════════════════════════════════════

#' Detect response type ("binary" / "quantitative" / "categorical" / "none").
detect_response_type <- function(response_raw, responseType_opt) {
  if (responseType_opt != "auto") return(responseType_opt)
  if (is.null(response_raw)) return("none")
  n_unique <- length(unique(na.omit(response_raw)))
  if (n_unique == 2) "binary"
  else if (is.numeric(response_raw)) "quantitative"
  else if (n_unique > 2 & n_unique <= 6) "categorical"
  else "none"
}

#' Prepare response as integer (binary), factor (categorical), or numeric (quantitative).
prepare_response <- function(response_raw, response_type) {
  if (is.null(response_raw) || is.null(response_type)) return(NULL)
  if (response_type == "binary") {
    r <- as.integer(as.factor(response_raw)) - 1L
    r[is.na(response_raw)] <- NA_integer_
    r
  } else if (response_type == "categorical") {
    as.factor(response_raw)
  } else {
    as.numeric(response_raw)
  }
}

#' Prepare covariate data frame (factor-encode character columns).
prepare_covariates <- function(data, covariate_vars) {
  if (length(covariate_vars) == 0) return(NULL)
  cov_df <- data[, covariate_vars, drop = FALSE]
  for (v in covariate_vars)
    if (!is.numeric(cov_df[[v]])) cov_df[[v]] <- as.factor(cov_df[[v]])
  cov_df
}


# ══════════════════════════════════════════════════════════════════════════════
# 6. Genetic model encoding
# ══════════════════════════════════════════════════════════════════════════════

#' Build ordered vector of all genotypes for a SNP, respecting user-defined
#' level order when available.
.all_genos_for_snp <- function(user_levels, snp_char, ref) {
  if (!is.null(user_levels) && length(user_levels) > 0)
    user_levels
  else {
    obs     <- sort(unique(snp_char[!is.na(snp_char)]))
    non_ref <- setdiff(obs, ref)
    is_hom  <- function(g) { p <- strsplit(g, "/", fixed = TRUE)[[1]]; length(p) == 2 && p[1] == p[2] }
    # A monomorphic SNP has no non-reference genotypes; vapply() (not sapply())
    # returns logical(0) rather than a list() here, so the subset stays valid.
    if (length(non_ref) == 0L) return(ref)
    hom <- vapply(non_ref, is_hom, logical(1))
    c(ref, non_ref[!hom], non_ref[hom])
  }
}

# Exponentiate a log-OR and clamp extreme values (|beta| > log(1e4)) to NA.
# Separation produces OR > 10000 or < 0.0001, which are uninformative.
.exp_or <- function(x) {
  threshold <- log(1e4)
  clamped   <- !is.na(x) & abs(x) > threshold
  v         <- exp(x)
  v[clamped] <- NA_real_
  v
}

#' Encode SNP under a given genetic model.
encode_model <- function(geno_char, ref, model, user_levels = NULL) {
  ref_allele <- strsplit(ref, "/")[[1]][1]
  dosage <- sapply(geno_char, function(g) {
    if (is.na(g) || g == "NA") return(NA_integer_)
    sum(strsplit(g, "/")[[1]] == ref_allele)
  })

  make_geno_factor <- function(values_01, ref_label, alt_label) {
    lbl <- ifelse(is.na(values_01), NA_character_,
                  ifelse(values_01 == 0L, ref_label, alt_label))
    factor(lbl, levels = c(ref_label, alt_label))
  }

  is_hom_fn <- function(g) { p <- strsplit(g, "/")[[1]]; length(p) == 2 && p[1] == p[2] }
  # vapply() returns logical(0) for an empty vector; sapply() returns list(),
  # which cannot index — so a monomorphic SNP (no non-reference genotypes) would
  # crash. Use this everywhere a homozygosity mask is built from observed genos.
  is_hom_vec <- function(v) if (length(v) == 0L) logical(0) else vapply(v, is_hom_fn, logical(1))
  obs  <- unique(na.omit(geno_char[geno_char != "NA"]))
  het  <- obs[obs != ref & !is_hom_vec(obs)]
  hom2 <- obs[obs != ref &  is_hom_vec(obs)]

  switch(model,
    codominant = {
      if (!is.null(user_levels) && length(user_levels) > 0) {
        obs_genos <- unique(geno_char[!is.na(geno_char) & geno_char != "NA"])
        lvls <- user_levels[user_levels %in% obs_genos]
        lvls <- c(lvls, obs_genos[!obs_genos %in% lvls])
      } else {
        non_ref <- unique(geno_char[geno_char != ref & !is.na(geno_char) & geno_char != "NA"])
        lvls <- c(ref, non_ref[!is_hom_vec(non_ref)], non_ref[is_hom_vec(non_ref)])
      }
      factor(geno_char, levels = lvls)
    },
    dominant = {
      alt_label <- paste(c(het, hom2), collapse = "-")
      make_geno_factor(as.integer(dosage < 2), ref, alt_label)
    },
    recessive = {
      ref_label <- paste(c(ref, het), collapse = "-")
      make_geno_factor(as.integer(dosage == 0), ref_label, hom2[1])
    },
    overdominant = {
      ref_label <- paste(c(ref, hom2), collapse = "-")
      is_het_01 <- as.integer(sapply(geno_char, function(g) {
        if (is.na(g) || g == "NA") NA_integer_
        else as.integer(length(unique(strsplit(g, "/")[[1]])) > 1)
      }))
      make_geno_factor(is_het_01, ref_label, het[1])
    },
    logadditive = 2L - as.integer(dosage)
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# 7. Formula and model fitting
# ══════════════════════════════════════════════════════════════════════════════

# Wrap a column name in backticks for safe formula interpolation.
safe_term <- function(x) paste0("`", gsub("`", "\\`", x, fixed = TRUE), "`")

# Escape a vector of names and collapse to a "+"-joined RHS string.
safe_rhs <- function(nms) paste(sapply(nms, safe_term), collapse = " + ")

# Evaluate `code` under a fixed RNG seed, restoring the previous RNG state on
# exit. haplo.glm's EM consumes R's RNG, so without this its haplotype effect
# estimates drift between identical calls; seeding makes them reproducible
# without disturbing the caller's RNG stream.
with_fixed_seed <- function(code, seed = 20240920L) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv))
           get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv))
            rm(".Random.seed", envir = .GlobalEnv))
  set.seed(seed)
  force(code)
}

# Evaluate `expr`, capturing any warnings it emits (so they neither print nor
# abort) alongside its value. Used to surface model-fit warnings (separation,
# non-convergence) as table notes instead of silently swallowing them.
# Returns list(value, warnings). If `expr` errors, value is NULL.
with_warnings <- function(expr) {
  warns <- character(0)
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) NULL),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = val, warnings = warns)
}

# Human-readable convergence / separation issues for a fitted model.
# `fit` is a glm/lm/multinom object (or NULL); `warnings` is the vector captured
# by with_warnings(). Returns a character vector (empty = clean fit).
fit_diagnostics <- function(fit, warnings = character(0)) {
  issues <- character(0)
  if (!is.null(fit)) {
    # glm/lm expose $converged (logical); nnet::multinom exposes $convergence
    # (0 = converged, 1 = maximum iterations reached).
    if (isFALSE(fit$converged))
      issues <- c(issues, "model did not converge")
    if (!is.null(fit$convergence) && !isTRUE(fit$convergence == 0))
      issues <- c(issues, "model did not converge (maximum iterations reached)")
  }
  if (any(grepl("fitted probabilities numerically 0 or 1", warnings, fixed = TRUE)))
    issues <- c(issues, "complete or quasi-complete separation (fitted probabilities of 0 or 1)")
  if (length(warnings) && any(grepl("did not converge", warnings, fixed = TRUE)))
    issues <- c(issues, "fitting algorithm reported non-convergence")
  unique(issues)
}

#' Fit association model for one SNP under one genetic model.
#' For categorical response, fits nnet::multinom and returns one result per category.
#' `null_fit` (binary/quantitative only) lets the caller supply the covariate-only
#' null model — identical across all genetic models — so it is fit once, not per model.
fit_model <- function(snp_enc, response, covariates_df, model_name,
                      response_type, ci_width, null_fit = NULL) {
  df <- data.frame(resp = response, snp = snp_enc)
  cov_formula <- ""
  if (!is.null(covariates_df) && ncol(covariates_df) > 0) {
    df          <- cbind(df, covariates_df)
    cov_formula <- paste("+", safe_rhs(names(covariates_df)))
  }
  df <- df[complete.cases(df), , drop = FALSE]
  formula_full <- as.formula(paste("resp ~ snp", cov_formula))
  formula_null <- as.formula(paste("resp ~ 1",   cov_formula))

  if (response_type == "categorical") {
    if (!requireNamespace("nnet", quietly = TRUE)) return(NULL)
    df$resp <- as.factor(df$resp)
    ff       <- with_warnings(nnet::multinom(formula_full, data = df, trace = FALSE))
    fit_full <- ff$value
    fit_null <- tryCatch(nnet::multinom(formula_null, data = df, trace = FALSE),
                         error = function(e) NULL)
    if (is.null(fit_full)) return(NULL)
    lrt      <- tryCatch(anova(fit_null, fit_full), error = function(e) NULL)
    global_p <- if (!is.null(lrt) && nrow(lrt) >= 2) lrt[2, "Pr(Chi)"] else NA_real_
    aic_val  <- AIC(fit_full)
    bic_val  <- BIC(fit_full)
    coefs    <- summary(fit_full)$coefficients
    ses      <- summary(fit_full)$standard.errors
    cats     <- rownames(coefs)
    snp_cols <- grep("^snp", colnames(coefs))
    if (length(snp_cols) == 0) return(NULL)
    z_crit   <- qnorm(1 - (1 - ci_width / 100) / 2)
    result <- list()
    for (cat in cats) {
      for (j in snp_cols) {
        beta  <- coefs[cat, j]; se <- ses[cat, j]
        ci_lo <- beta - z_crit * se; ci_hi <- beta + z_crit * se
        pval  <- 2 * (1 - pnorm(abs(beta / se)))
        result[[length(result) + 1L]] <- list(
          category   = cat,
          comparison = sub("^snp", "", colnames(coefs)[j]),
          effect     = .exp_or(beta), ci_low = .exp_or(ci_lo), ci_high = .exp_or(ci_hi),
          pval = pval, global_p = global_p, aic = aic_val, bic = bic_val,
          is_categorical = TRUE)
      }
    }
    attr(result, "diagnostics") <- fit_diagnostics(fit_full, ff$warnings)
    return(result)
  }

  tryCatch({
    if (response_type == "binary") {
      ff        <- with_warnings(glm(formula_full, data = df, family = binomial()))
      fit_full  <- ff$value
      fit_null  <- null_fit %||% glm(formula_null, data = df, family = binomial())
      lrtest <- "Chisq"; lrtest_label <- "Pr(>Chi)"; pval_col <- "Pr(>|z|)"
    } else {
      ff        <- with_warnings(lm(formula_full, data = df))
      fit_full  <- ff$value
      fit_null  <- null_fit %||% lm(formula_null, data = df)
      lrtest <- "F"; lrtest_label <- "Pr(>F)"; pval_col <- "Pr(>|t|)"
    }
    if (is.null(fit_full)) return(NULL)
    lrt      <- tryCatch(anova(fit_null, fit_full, test = lrtest), error = function(e) NULL)
    global_p <- if (!is.null(lrt)) lrt[2, lrtest_label] else NA_real_
    aic_val  <- AIC(fit_full)
    bic_val  <- BIC(fit_full)
    coefs    <- summary(fit_full)$coefficients
    snp_rows <- grep("^snp", rownames(coefs))
    if (length(snp_rows) == 0) return(NULL)
    ci <- tryCatch(
      confint(fit_full, level = ci_width / 100)[snp_rows, , drop = FALSE],
      error = function(e) matrix(NA, nrow = length(snp_rows), ncol = 2))
    res <- lapply(seq_along(snp_rows), function(i) {
      row  <- snp_rows[i]
      beta <- coefs[row, "Estimate"]; pval <- coefs[row, pval_col]
      ci_lo <- ci[i, 1]; ci_hi <- ci[i, 2]
      if (response_type == "binary")
        list(effect = .exp_or(beta), ci_low = .exp_or(ci_lo), ci_high = .exp_or(ci_hi),
             pval = pval, global_p = global_p, aic = aic_val, bic = bic_val,
             comparison = sub("^snp", "", rownames(coefs)[row]))
      else
        list(effect = beta, ci_low = ci_lo, ci_high = ci_hi,
             pval = pval, global_p = global_p, aic = aic_val, bic = bic_val,
             comparison = sub("^snp", "", rownames(coefs)[row]))
    })
    attr(res, "diagnostics") <- fit_diagnostics(fit_full, ff$warnings)
    res
  }, error = function(e) NULL)
}

#' Fit SNP × covariate interaction model.
fit_interaction_model <- function(snp_enc, response, covariates_df,
                                  interaction_var, model_name,
                                  response_type, ci_width,
                                  conditional = FALSE, cond_var = interaction_var) {
  df <- data.frame(resp = response, snp = snp_enc)
  adj_covs <- character(0)
  if (!is.null(covariates_df) && ncol(covariates_df) > 0) {
    df       <- cbind(df, covariates_df)
    adj_covs <- setdiff(names(covariates_df), interaction_var)
  }
  if (!(interaction_var %in% names(df))) return(NULL)
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 5) return(NULL)
  iv_safe  <- safe_term(interaction_var)
  adj_part <- if (length(adj_covs) > 0) paste("+", safe_rhs(adj_covs)) else ""

  if (conditional) {
    if (cond_var == "snp") {
      formula_fit <- as.formula(paste("resp ~ snp /", iv_safe, adj_part))
    } else {
      formula_fit <- as.formula(paste("resp ~", iv_safe, "/ snp", adj_part))
    }
    formula_add <- as.formula(paste("resp ~ snp +", iv_safe, adj_part))
    if (response_type == "binary") {
      ff      <- with_warnings(glm(formula_fit, data = df, family = binomial()))
      fit     <- ff$value
      fit_add <- glm(formula_add, data = df, family = binomial())
      pval_col <- "Pr(>|z|)"; lrtest <- "Chisq"; lrtest_label <- "Pr(>Chi)"
    } else {
      ff      <- with_warnings(lm(formula_fit, data = df))
      fit     <- ff$value
      fit_add <- lm(formula_add, data = df)
      pval_col <- "Pr(>|t|)"; lrtest <- "F"; lrtest_label <- "Pr(>F)"
    }
    if (is.null(fit)) return(NULL)
    lrt_cond <- tryCatch(anova(fit_add, fit, test = lrtest), error = function(e) NULL)
    p_inter  <- lrt_cond[2, lrtest_label]
    coefs    <- summary(fit)$coefficients
    ci_mat   <- tryCatch(confint(fit, level = ci_width / 100),
                         error = function(e) matrix(NA, nrow = nrow(coefs), ncol = 2))
    aic_val  <- AIC(fit); bic_val <- BIC(fit)
    all_rows <- rownames(coefs)
    inter_rows_idx <- grep(":", all_rows)
    snp_rows_idx   <- setdiff(grep("^snp", all_rows), inter_rows_idx)
    adj_rows_idx   <- setdiff(
      seq_along(all_rows),
      c(grep("^\\(Intercept\\)", all_rows),
        snp_rows_idx, inter_rows_idx,
        which(startsWith(all_rows, interaction_var))))
    if (length(inter_rows_idx) == 0) return(NULL)
    all_keep <- unique(c(snp_rows_idx, inter_rows_idx, adj_rows_idx))
    first_inter_done <- FALSE
    result <- lapply(all_keep, function(r) {
      term  <- all_rows[r]
      beta  <- coefs[r, "Estimate"]; pval <- coefs[r, pval_col]
      ci_lo <- ci_mat[r, 1]; ci_hi <- ci_mat[r, 2]
      is_inter_term <- r %in% inter_rows_idx
      row_type <- if (r %in% snp_rows_idx) "snp"
                  else if (is_inter_term)   "interaction"
                  else                      "adjustment"
      attach_p <- is_inter_term && !first_inter_done
      if (attach_p) first_inter_done <<- TRUE
      if (response_type == "binary")
        list(term = term, effect = .exp_or(beta), ci_low = .exp_or(ci_lo), ci_high = .exp_or(ci_hi),
             pval = pval, pval_interaction = if (attach_p) p_inter else NA_real_,
             aic = aic_val, bic = bic_val, row_type = row_type)
      else
        list(term = term, effect = beta, ci_low = ci_lo, ci_high = ci_hi,
             pval = pval, pval_interaction = if (attach_p) p_inter else NA_real_,
             aic = aic_val, bic = bic_val, row_type = row_type)
    })
    attr(result, "pval_interaction") <- p_inter
    attr(result, "diagnostics")      <- fit_diagnostics(fit, ff$warnings)
    result
  } else {
    formula_int  <- as.formula(paste("resp ~ snp *", iv_safe, adj_part))
    formula_main <- as.formula(paste("resp ~ snp +", iv_safe, adj_part))
    tryCatch({
      if (response_type == "binary") {
        ff       <- with_warnings(glm(formula_int,  data = df, family = binomial()))
        fit_int  <- ff$value
        fit_main <- glm(formula_main, data = df, family = binomial())
        pval_col <- "Pr(>|z|)"; lrtest <- "Chisq"; lrtest_label <- "Pr(>Chi)"
      } else {
        ff       <- with_warnings(lm(formula_int,  data = df))
        fit_int  <- ff$value
        fit_main <- lm(formula_main, data = df)
        pval_col <- "Pr(>|t|)"; lrtest <- "F"; lrtest_label <- "Pr(>F)"
      }
      if (is.null(fit_int)) return(NULL)
      lrt      <- tryCatch(anova(fit_main, fit_int, test = lrtest), error = function(e) NULL)
      p_inter  <- if (!is.null(lrt)) lrt[2, lrtest_label] else NA_real_
      aic_val  <- AIC(fit_int); bic_val <- BIC(fit_int)
      coefs    <- summary(fit_int)$coefficients
      ci       <- tryCatch(confint(fit_int, level = ci_width / 100),
                           error = function(e) matrix(NA, nrow = nrow(coefs), ncol = 2,
                                                      dimnames = list(rownames(coefs), c("lo","hi"))))
      all_rows   <- rownames(coefs)
      snp_rows   <- grep("^snp", all_rows)
      inter_rows <- which(
        (startsWith(all_rows, "snp") & grepl(":", all_rows, fixed = TRUE) &
           sapply(all_rows, function(x) any(grepl(interaction_var, strsplit(x, ":")[[1]], fixed = TRUE)))) |
        (startsWith(all_rows, interaction_var) & grepl(":snp", all_rows, fixed = TRUE)))
      covar_rows <- setdiff(which(startsWith(all_rows, interaction_var)), inter_rows)
      adj_rows   <- setdiff(seq_along(all_rows),
                            c(grep("^\\(Intercept\\)", all_rows),
                              snp_rows, inter_rows, covar_rows))
      keep_rows  <- unique(c(snp_rows, inter_rows, covar_rows, adj_rows))
      if (length(keep_rows) == 0) return(NULL)
      result <- lapply(keep_rows, function(r) {
        beta  <- coefs[r, "Estimate"]; pval <- coefs[r, pval_col]
        ci_lo <- ci[r, 1]; ci_hi <- ci[r, 2]
        is_inter <- r %in% inter_rows
        row_type <- if (r %in% snp_rows)       "snp"
                    else if (r %in% inter_rows) "interaction"
                    else if (r %in% covar_rows) "covariate"
                    else                        "adjustment"
        if (response_type == "binary")
          list(term = all_rows[r], effect = .exp_or(beta), ci_low = .exp_or(ci_lo), ci_high = .exp_or(ci_hi),
               pval = pval, pval_interaction = if (is_inter) p_inter else NA_real_,
               aic = aic_val, bic = bic_val, is_first = (r == keep_rows[1]), row_type = row_type)
        else
          list(term = all_rows[r], effect = beta, ci_low = ci_lo, ci_high = ci_hi,
               pval = pval, pval_interaction = if (is_inter) p_inter else NA_real_,
               aic = aic_val, bic = bic_val, is_first = (r == keep_rows[1]), row_type = row_type)
      })
      attr(result, "pval_interaction") <- p_inter
      attr(result, "diagnostics")      <- fit_diagnostics(fit_int, ff$warnings)
      result
    }, error = function(e) NULL)
  }
}


# ══════════════════════════════════════════════════════════════════════════════
# 8. Data preparation
# ══════════════════════════════════════════════════════════════════════════════

#' Validate, clean, and pre-process all inputs.
#'
#' Returns a named list ("prep object") consumed by every compute_* function.
#' Call this once, then pass the result to whichever compute_* functions you need.
#'
#' @return Named list: $data, $snp_vars, $snp_data, $response_var,
#'   $response_raw, $response_type, $response_enc, $cov_df,
#'   $complete_mask, $n_rows, $warnings
snp_prepare <- function(data, snps, response = NULL, covariates = NULL,
                        response_type = "auto", rm_snp_missing = FALSE) {

  response_var <- if (!is.null(response) && nchar(response) > 0) response else NULL
  response_raw <- if (!is.null(response_var)) data[[response_var]] else NULL
  # A character response carries no stored level order, so downstream as.factor()
  # calls (encoding, stratification, descriptives) would sort it alphabetically.
  # Normalize once here to a factor in order of appearance in the data, so the
  # reference is the first observed level; existing factors (ordinal or
  # user-ordered, as jamovi delivers nominal variables) keep their intended order.
  if (is.character(response_raw))
    response_raw <- factor(response_raw, levels = unique(response_raw[!is.na(response_raw)]))
  rtype        <- detect_response_type(response_raw, response_type)
  response_enc <- prepare_response(response_raw, rtype)

  cov_df <- prepare_covariates(data, covariates %||% character(0))
  if (is.null(cov_df) && !is.null(response_raw))
    cov_df <- data.frame(row.names = seq_len(nrow(data)))

  n_rows        <- nrow(data)
  complete_mask <- rep(TRUE, n_rows)
  if (!is.null(response_enc))                 complete_mask <- complete_mask & !is.na(response_enc)
  if (!is.null(cov_df) && ncol(cov_df) > 0)  complete_mask <- complete_mask & complete.cases(cov_df)

  if (length(snps) == 0) {
    val      <- list(bad_html = "No SNPs specified — skipping SNP processing")
    snp_data <- list()
    snp_vars <- character(0)
  } else {
    val      <- validate_snp_vars(snps, data)
    snp_vars <- val$valid_snps

    if (isTRUE(rm_snp_missing) && length(snp_vars) > 0) {
      snp_mat <- as.data.frame(
        lapply(data[, snp_vars, drop = FALSE],
               function(col) clean_null_alleles(as.character(col))),
        stringsAsFactors = FALSE)
      complete_mask <- complete_mask & complete.cases(snp_mat)
    }

    snp_data <- lapply(setNames(snp_vars, snp_vars), function(nm) {
      raw         <- data[[nm]]
      user_levels <- get_snp_level_order(raw)
      clean       <- clean_null_alleles(as.character(raw))
      snp_mask    <- complete_mask & !is.na(clean)
      clean_cc    <- clean[snp_mask]
      geno_cc     <- parse_genotype(clean_cc, user_levels)
      if (is.null(geno_cc)) return(NULL)
      ref        <- get_ref_genotype(geno_cc, user_levels)
      summary_cc <- summary(geno_cc)
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
        # eligible rows (complete response/covariates) whose SNP genotype is
        # missing; snp_mask == complete_mask & !is.na(clean), so this is exactly
        # the complement within the eligible set.
        n_missing   = sum(is.na(clean) & complete_mask)
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
# 9. Descriptive statistics
# ══════════════════════════════════════════════════════════════════════════════

#' Compute allele frequency table for one SNP.
#' @return data.frame: allele, overall, stat_g0..stat_gN
compute_allele_freq <- function(snp_nm, prep, subpop = FALSE, show_missing = FALSE) {
  sd          <- prep$snp_data[[snp_nm]]
  sm          <- sd$summary_cc
  af          <- sm$allele.freq
  user_levels <- sd$user_levels
  rtype       <- prep$response_type

  response_raw_full <- prep$response_raw
  response_raw_cc   <- if (!is.null(response_raw_full)) response_raw_full[sd$snp_mask] else NULL

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
        idx  <- resp_chr_cc == lvl & !is.na(resp_chr_cc)
        n_al <- sum(unlist(alleles_split_cc[idx]) == al, na.rm = TRUE)
        fmt_cat(n_al, grp_allele_totals[[lvl]])
      }, character(1))
    } else rep("", 0L)
    make_row(al, fmt_catpct(count, prop), grp_strs)
  })
  rows <- Filter(Negate(is.null), rows)

  if (isTRUE(show_missing) && sd$n_missing > 0) {
    n_elig <- sum(prep$complete_mask)
    grp_miss <- if (do_strat) {
      clean_full    <- sd$clean
      resp_chr_full <- as.character(response_raw_full)
      grp_eligible  <- sapply(grp_levels, function(lvl)
        sum(prep$complete_mask & !is.na(resp_chr_full) & resp_chr_full == lvl, na.rm = TRUE))
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
#' @return data.frame: genotype, overall, stat_g0..stat_gN, response_stat
compute_geno_freq <- function(snp_nm, prep, subpop = FALSE, show_missing = FALSE) {
  sd           <- prep$snp_data[[snp_nm]]
  sm           <- sd$summary_cc
  rtype        <- prep$response_type
  response_cc  <- if (!is.null(prep$response_enc)) prep$response_enc[sd$snp_mask] else NULL

  response_raw_full <- prep$response_raw
  response_raw_cc   <- if (!is.null(response_raw_full)) response_raw_full[sd$snp_mask] else NULL

  geno_obj <- sd$geno_cc
  snp_chr  <- as.character(geno_obj)

  gf <- tryCatch(reorder_geno(sm$genotype.freq, sd$ref, sd$user_levels),
                 error = function(e) sm$genotype.freq)
  if (!is.matrix(gf))
    gf <- matrix(gf, nrow = 1L, dimnames = list("", names(gf)))
  gf <- gf[rownames(gf) != "NA", , drop = FALSE]

  do_strat   <- isTRUE(subpop) && !is.null(response_raw_cc) &&
                rtype %in% c("binary", "categorical")
  grp_levels <- if (do_strat) levels(as.factor(response_raw_full)) else character(0)
  n_grp      <- length(grp_levels)

  resp_chr_cc  <- if (do_strat) as.character(response_raw_cc) else NULL
  strat_totals <- if (do_strat)
    sapply(grp_levels, function(lvl)
      sum(resp_chr_cc == lvl & !is.na(resp_chr_cc) & !is.na(snp_chr), na.rm = TRUE))
  else NULL

  make_row <- function(geno_lbl, overall_str, resp_stat = "", grp_strs = rep("", n_grp)) {
    row <- list(genotype = geno_lbl, overall = overall_str, response_stat = resp_stat)
    for (j in seq_len(n_grp)) row[[paste0("stat_g", j-1L)]] <- grp_strs[j]
    as.data.frame(row, stringsAsFactors = FALSE)
  }

  rows <- lapply(seq_len(nrow(gf)), function(i) {
    geno  <- rownames(gf)[i]
    if (geno == "NA") return(NULL)
    count <- as.integer(gf[i, "Count"])
    prop  <- gf[i, "Proportion"] * 100

    resp_stat <- ""
    if (rtype == "quantitative" && !is.null(response_cc)) {
      mask <- snp_chr == geno & !is.na(snp_chr) & !is.na(response_cc)
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
        idx <- resp_chr_cc == lvl & !is.na(resp_chr_cc)
        n_g <- sum(idx & snp_chr == geno, na.rm = TRUE)
        fmt_cat(n_g, strat_totals[[lvl]])
      }, character(1))
    } else rep("", 0L)

    make_row(geno, fmt_catpct(count, prop), resp_stat, grp_strs)
  })
  rows <- Filter(Negate(is.null), rows)

  if (isTRUE(show_missing) && sd$n_missing > 0) {
    n_elig        <- sum(prep$complete_mask)
    clean_full    <- sd$clean
    resp_chr_full <- if (!is.null(response_raw_full)) as.character(response_raw_full) else NULL

    grp_miss <- if (do_strat) {
      grp_eligible <- sapply(grp_levels, function(lvl)
        sum(prep$complete_mask & !is.na(resp_chr_full) & resp_chr_full == lvl, na.rm = TRUE))
      vapply(grp_levels, function(lvl) {
        n_lv <- sum(is.na(clean_full) & prep$complete_mask &
                    !is.na(resp_chr_full) & resp_chr_full == lvl, na.rm = TRUE)
        fmt_cat(n_lv, grp_eligible[[lvl]])
      }, character(1))
    } else rep("", 0L)

    rows <- c(rows, list(make_row("Missing", fmt_cat(sd$n_missing, n_elig), "", grp_miss)))
  }

  do.call(rbind, rows)
}


#' Compute Hardy-Weinberg test for one SNP, optionally stratified.
#' @return list: $col_labels character(3), $rows data.frame
compute_hwe <- function(snp_nm, prep, subpop = FALSE, show_missing = FALSE) {
  sd           <- prep$snp_data[[snp_nm]]
  geno_obj     <- sd$geno_cc
  response_raw <- if (!is.null(prep$response_raw)) prep$response_raw[sd$snp_mask] else NULL

  hw <- tryCatch(genetics::HWE.exact(geno_obj), error = function(e) NULL)
  if (is.null(hw)) return(NULL)

  get_ordered_counts <- function(go) {
    gf <- tryCatch(reorder_geno(summary(go)$genotype.freq, sd$ref, sd$user_levels),
                   error = function(e) summary(go)$genotype.freq)
    if (!is.matrix(gf))
      gf <- matrix(gf, nrow = 1L, dimnames = list("", names(gf)))
    gf <- gf[rownames(gf) != "NA", , drop = FALSE]
    list(labels = rownames(gf), counts = as.integer(gf[, "Count"]))
  }

  info       <- get_ordered_counts(geno_obj)
  col_labels <- if (length(info$labels) == 3) info$labels else c("AA", "AB", "BB")

  rows     <- list()
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
      if (!is.matrix(gf))
        gf <- matrix(gf, nrow = 1L, dimnames = list("", names(gf)))
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
      snp      = nm, alleles = alleles_lbl,
      group    = if (do_strat) "All" else "",
      n        = res_all$n,
      missing  = if (n_excl > 0L) n_excl else NA_integer_,
      maf      = res_all$maf, geno_counts = res_all$geno_counts,
      hwe_pval = res_all$hwe_pval,
      stringsAsFactors = FALSE))

    if (do_strat) {
      resp_chr      <- as.character(resp_cc)
      resp_chr_full <- as.character(response_raw)
      strat_eligible <- table(factor(
        resp_chr_full[prep$complete_mask & !is.na(resp_chr_full)],
        levels = grp_levels))

      for (lvl in grp_levels) {
        mask <- !is.na(resp_chr) & resp_chr == lvl
        if (sum(mask) == 0) next
        g_sub <- tryCatch(parse_genotype(snp_cc[mask], user_levels), error = function(e) NULL)
        if (is.null(g_sub)) next
        res_s    <- compute_row_stats(g_sub, alt_allele)
        n_excl_s <- max(0L, as.integer(strat_eligible[lvl]) - res_s$n)
        out[[length(out) + 1L]] <- data.frame(
          snp     = nm, alleles = alleles_lbl, group = lvl,
          n       = res_s$n,
          missing = if (n_excl_s > 0L) n_excl_s else NA_integer_,
          maf     = res_s$maf, geno_counts = res_s$geno_counts,
          hwe_pval = res_s$hwe_pval,
          stringsAsFactors = FALSE)
      }
    }
    do.call(rbind, out)
  })

  do.call(rbind, Filter(Negate(is.null), rows))
}


#' Compute covariate descriptive table.
#' @return list: $table data.frame, $grp_levels character, $notes character
compute_cov_desc <- function(prep, subpop = FALSE, mask = NULL) {
  cov_df       <- prep$cov_df
  response_raw <- prep$response_raw
  rtype        <- prep$response_type
  response_var <- prep$response_var

  if (is.null(cov_df)) return(NULL)

  # Whether a Missing row exists is decided from the RAW (pre-mask) column, so
  # the row structure depends only on the data — which .init can read — and not
  # on completeCases/rmSnpMissing. A variable with missing values therefore keeps
  # its Missing row (showing 0 (0.0%)) when a mask removes them, instead of the
  # row appearing and disappearing, which .init could not predict.
  has_miss_resp <- !is.null(response_raw) && anyNA(response_raw)
  has_miss_cov  <- vapply(cov_df, anyNA, logical(1))

  n_total <- prep$n_rows
  if (!is.null(mask) && any(!mask)) {
    if (!is.null(response_raw)) response_raw <- response_raw[mask]
    cov_df <- cov_df[mask, , drop = FALSE]
  }
  n_analyzed <- length(if (!is.null(response_raw)) response_raw else cov_df[[1]])

  do_strat  <- isTRUE(subpop) && !is.null(response_raw) &&
               rtype %in% c("binary", "categorical")
  valid_resp <- if (!is.null(response_raw)) !is.na(response_raw)
                else rep(TRUE, n_analyzed)

  grp_levels <- NULL; mask_list <- NULL; totals <- NULL
  if (do_strat) {
    grp_fac    <- as.factor(response_raw)
    grp_levels <- levels(grp_fac)
    mask_list  <- lapply(grp_levels, function(l) valid_resp & as.character(grp_fac) == l)
    names(mask_list) <- grp_levels
    totals     <- sapply(mask_list, sum)
  }
  n_grp <- length(grp_levels)

  # key identifies the row across runs so .init can pre-create it and .run can
  # setRow into it (see .covdesc_keys in snpStats.b.R). Levels are namespaced
  # under "|lvl|" so a factor level literally named "Missing" cannot collide
  # with the missing-count row's "|miss" key.
  make_row <- function(variable, level, overall,
                       grp_stats = rep("", n_grp), pval = NA_real_, key = NULL) {
    row <- list(key      = as.character(if (is.null(key)) paste0(variable, "|lvl|", level) else key),
                variable = as.character(variable),
                level    = as.character(level),
                overall  = as.character(overall))
    for (j in seq_len(n_grp))
      row[[paste0("stat_g", j - 1L)]] <- as.character(grp_stats[j])
    row$pval <- as.numeric(pval)
    as.data.frame(row, stringsAsFactors = FALSE)
  }

  grp_stats_for <- function(mask, fmt_fn = function(n, tot) fmt_cat(n, tot)) {
    if (!do_strat) return(rep("", 0L))
    vapply(grp_levels, function(lvl)
      fmt_fn(sum(mask & mask_list[[lvl]], na.rm = TRUE), totals[[lvl]]), character(1))
  }
  grp_stats_cont <- function(x_vec) {
    if (!do_strat) return(rep("", 0L))
    vapply(grp_levels, function(lvl) fmt_cont(x_vec[mask_list[[lvl]]]), character(1))
  }

  rows <- list()

  if (!is.null(response_raw) && !is.null(rtype)) {
    if (rtype == "quantitative") {
      rows[[length(rows)+1L]] <- make_row(response_var, "Mean ± SD",
        fmt_cont(response_raw), grp_stats_cont(response_raw),
        key = paste0(response_var, "|mean"))
      mask <- !is.na(response_raw)
      rows[[length(rows)+1L]] <- make_row(response_var, "Valid", as.character(sum(mask)),
        grp_stats_for(mask, function(n, tot) as.character(n)),
        key = paste0(response_var, "|valid"))
      if (has_miss_resp)
        rows[[length(rows)+1L]] <- make_row(response_var, "Missing",
          fmt_cat(sum(!mask), length(response_raw)), grp_stats_for(!mask),
          key = paste0(response_var, "|miss"))
    } else {
      mask <- valid_resp
      rows[[length(rows)+1L]] <- make_row(response_var, "Valid", as.character(sum(mask)),
        grp_stats_for(mask, function(n, tot) as.character(n)),
        key = paste0(response_var, "|valid"))
      if (has_miss_resp)
        rows[[length(rows)+1L]] <- make_row(response_var, "Missing",
          fmt_cat(sum(!valid_resp), length(response_raw)), grp_stats_for(!valid_resp),
          key = paste0(response_var, "|miss"))
    }
  }

  if (ncol(cov_df) > 0) {
    for (v in names(cov_df)) {
      col    <- cov_df[[v]]
      n      <- length(col)
      n_miss <- sum(is.na(col))
      is_cat <- is.factor(col) || is.character(col)
      if (is_cat && !is.factor(col)) col <- factor(col)

      if (is_cat) {
        pval_cat <- if (do_strat) tryCatch({
          ct <- table(col[valid_resp], as.factor(response_raw)[valid_resp])
          suppressWarnings(chisq.test(ct)$p.value)
        }, error = function(e) NA_real_) else NA_real_

        first <- TRUE
        for (lvl in levels(col)) {
          mask <- !is.na(col) & col == lvl
          rows[[length(rows)+1L]] <- make_row(v, lvl, fmt_cat(sum(mask), n),
            grp_stats_for(mask), if (first) pval_cat else NA_real_)
          first <- FALSE
        }
        if (has_miss_cov[[v]])
          rows[[length(rows)+1L]] <- make_row(v, "Missing", fmt_cat(n_miss, n),
            grp_stats_for(is.na(col)), key = paste0(v, "|miss"))
      } else {
        pval_cont <- if (do_strat) tryCatch({
          grps <- split(col[valid_resp], as.factor(response_raw)[valid_resp])
          if (length(grps) == 2) t.test(grps[[1]], grps[[2]])$p.value
          else summary(aov(col ~ as.factor(response_raw)))[[1]][["Pr(>F)"]][1]
        }, error = function(e) NA_real_) else NA_real_

        rows[[length(rows)+1L]] <- make_row(v, "Mean ± SD",
          fmt_cont(col), grp_stats_cont(col), pval_cont, key = paste0(v, "|mean"))
        if (has_miss_cov[[v]])
          rows[[length(rows)+1L]] <- make_row(v, "Missing", fmt_cat(n_miss, n),
            grp_stats_for(is.na(col)), key = paste0(v, "|miss"))
      }
    }
  }

  notes <- if (!is.null(mask) && any(!mask)) {
    n_removed <- n_total - n_analyzed
    paste0("N analyzed: ", n_analyzed,
      if (n_removed > 0) paste0(" (", n_removed, " removed due to missing SNP data)") else "")
  } else character(0)

  list(table = do.call(rbind, rows), grp_levels = grp_levels, notes = notes)
}


# ══════════════════════════════════════════════════════════════════════════════
# 10. Association analysis
# ══════════════════════════════════════════════════════════════════════════════

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

.split_genos <- function(gl) strsplit(gl, "-", fixed = TRUE)[[1]]

.compute_stats <- function(geno_labels, snp_char, response, response_type,
                            response_raw = NULL) {
  split_genos <- .split_genos
  ref_al <- NULL
  for (lbl in geno_labels) {
    for (g in strsplit(lbl, "-", fixed = TRUE)[[1]]) {
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
#' @return list: $snp, $response_type, $col_titles, $notes, $rows data.frame
compute_assoc <- function(snp_nm, prep,
                          models   = c("codominant","dominant","recessive",
                                       "overdominant","logadditive"),
                          ci_width = 95) {
  if (is.null(prep$response_enc) || prep$response_type == "none") return(NULL)
  sd           <- prep$snp_data[[snp_nm]]
  response_enc <- prep$response_enc[sd$snp_mask]
  response_raw <- if (!is.null(prep$response_raw)) prep$response_raw[sd$snp_mask] else NULL
  cov_df_cc    <- if (!is.null(prep$cov_df)) prep$cov_df[sd$snp_mask, , drop=FALSE] else NULL
  rtype        <- prep$response_type
  # Use geno_cc (ref-allele-first orientation) not clean_cc — overdominant counts
  # were wrong when raw data had het alleles in non-ref-first order.
  snp_char     <- as.character(sd$geno_cc)
  ref          <- sd$ref
  user_levels  <- sd$user_levels
  n_miss       <- prep$n_rows - sd$n_typed

  all_genos <- .all_genos_for_snp(user_levels, snp_char, ref)
  is_categorical <- rtype == "categorical"

  col_titles <- list(
    effect = if (rtype %in% c("binary","categorical")) "OR" else "β",
    stat0  = switch(rtype,
      binary      = if (!is.null(response_raw)) levels(as.factor(response_raw))[1] else "Group 0",
      categorical = "N (%)",
      "Mean (SD)"),
    stat1  = switch(rtype, binary = if (!is.null(response_raw)) levels(as.factor(response_raw))[2] else "Group 1", ""))

  notes <- character(0)
  if (!is.null(cov_df_cc) && ncol(cov_df_cc) > 0)
    notes <- c(notes, paste0("Adjusted for: ", paste(names(cov_df_cc), collapse=", ")))
  if (n_miss > 0)
    notes <- c(notes, paste0(n_miss, " observation(s) excluded."))
  if (is_categorical && !is.null(response_raw))
    notes <- c(notes, paste0("Reference category: ‘",
      levels(as.factor(response_raw))[1], "’"))

  all_rows  <- list()
  any_clamp <- FALSE
  fit_diags <- character(0)

  # The covariate-only null model (response ~ covariates) is identical for every
  # genetic model — the SNP is already complete-cased, so all models share one
  # analysis subset. Fit it once here and reuse it for each model's LRT/F test
  # instead of refitting it inside fit_model five times. (Categorical uses
  # multinom and keeps its own null.)
  null_fit <- NULL
  if (rtype %in% c("binary", "quantitative")) {
    ndf <- data.frame(resp = response_enc)
    has_cov <- !is.null(cov_df_cc) && ncol(cov_df_cc) > 0
    if (has_cov) ndf <- cbind(ndf, cov_df_cc)
    ndf   <- ndf[complete.cases(ndf), , drop = FALSE]
    nform <- as.formula(paste("resp ~ 1", if (has_cov) paste("+", safe_rhs(names(cov_df_cc))) else ""))
    null_fit <- tryCatch(
      if (rtype == "binary") glm(nform, data = ndf, family = binomial())
      else                   lm(nform, data = ndf),
      error = function(e) NULL)
  }

  for (mdl in models) {
    snp_enc  <- encode_model(snp_char, ref, mdl, user_levels)
    res_list <- fit_model(snp_enc, response_enc, cov_df_cc, mdl, rtype, ci_width,
                          null_fit = null_fit)
    if (is.null(res_list)) next
    dg <- attr(res_list, "diagnostics")
    if (length(dg) > 0)
      fit_diags <- c(fit_diags, paste0(.MODEL_LABELS[[mdl]], ": ", paste(dg, collapse = "; ")))

    geno_labels <- .geno_labels_for_model(mdl, all_genos, ref)
    aic_val <- { a <- res_list[[1]]$aic; if (!is.null(a) && !is.nan(a)) round(a, 2) else NA_real_ }
    bic_val <- { b <- res_list[[1]]$bic; if (!is.null(b) && !is.nan(b)) round(b, 2) else NA_real_ }

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
            model    = paste0(.MODEL_LABELS[mdl], " — ", cat, " (n=", n_cat, ")"),
            genotype = "Per allele", stat0 = as.character(n_cat), stat1 = "",
            effect   = res$effect, ci_low = res$ci_low, ci_high = res$ci_high,
            pval     = res$pval,   aic    = aic_val,    bic     = bic_val,
            stringsAsFactors = FALSE)
          next
        }
        all_rows[[length(all_rows)+1L]] <- data.frame(
          model    = paste0(.MODEL_LABELS[mdl], " — ", cat, " (n=", n_cat, ")"),
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
            model="", genotype=gl,
            stat0=if ((i+1)<=length(cat_sts)) cat_sts[i+1] else "", stat1="",
            effect=res$effect, ci_low=res$ci_low, ci_high=res$ci_high,
            pval=res$pval, aic="", bic="",
            stringsAsFactors = FALSE)
        }
      }
      next
    }

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
        model="Additive", genotype="Per allele",
        stat0=stat0_val, stat1=stat1_val,
        effect=res$effect, ci_low=res$ci_low, ci_high=res$ci_high,
        # model p-value is the likelihood-ratio (binary/categorical) or F
        # (quantitative) test vs the null, matching the SNPstats reference tool.
        pval=res$global_p, aic=aic_val, bic=bic_val,
        stringsAsFactors = FALSE)
      next
    }
    # Every model reports its LRT/F model p on the reference row (as codominant
    # always did). Codominant additionally shows the per-genotype Wald p on each
    # non-reference row (vs the reference genotype), matching the SNPstats
    # reference tool; the other models leave per-genotype rows blank.
    pval_row1 <- res_list[[1]]$global_p
    all_rows[[length(all_rows)+1L]] <- data.frame(
      model    = .MODEL_LABELS[mdl], genotype = geno_labels[1],
      stat0    = st$s0[1], stat1 = st$s1[1],
      effect   = if (rtype=="binary") 1. else 0.,
      ci_low   = '', ci_high = '',
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
        pval=if (mdl == "codominant") res$pval else NA_real_,
        aic=NA_real_, bic=NA_real_,
        stringsAsFactors = FALSE)
    }
  }

  if (any_clamp)
    notes <- c(notes, "One or more OR/CI suppressed (shown as blank) due to complete or quasi-complete separation.")
  if (length(fit_diags) > 0)
    notes <- c(notes, paste0("Model fit warning(s) — ", paste(fit_diags, collapse = " | "),
                             ". Interpret the affected estimates with caution."))

  list(snp           = snp_nm,
       response_type = rtype,
       col_titles    = col_titles,
       notes         = notes,
       rows          = do.call(rbind, all_rows))
}


# ══════════════════════════════════════════════════════════════════════════════
# 11. Haplotype analysis
# ══════════════════════════════════════════════════════════════════════════════

# Required by haplo.stats::haplo.glm — removes rows with missing response or
# covariates while preserving the geno matrix attributes.
na.geno.keep <- function(m) {
  mf.gindx <- function(m) {
    # data.class() matches haplo.stats internals exactly.
    # ncol==2 check is WRONG for >=2 SNPs (setupGeno produces 2*n_loci columns).
    nvars    <- length(m)
    typevars <- rep("", nvars)
    for (i in seq_len(nvars)) typevars[i] <- data.class(m[[i]])
    gindx <- seq_len(nvars)[typevars == "model.matrix" | typevars == "matrix"]
    if (length(gindx) == 0) stop("No geno matrix found in data frame")
    if (length(gindx) >  1) stop("More than 1 geno matrix in data frame")
    gindx
  }
  gindx    <- mf.gindx(m)
  yxmiss   <- apply(is.na(m[, -gindx, drop=FALSE]), 1, any)
  gmiss    <- apply(is.na(m[,  gindx, drop=FALSE]), 1, all)
  genoAttr <- attributes(m[, gindx])
  allmiss  <- yxmiss | gmiss
  m        <- m[!allmiss, ]
  genoAttr$dim[1] <- genoAttr$dim[1] - sum(allmiss)
  nloc <- ncol(m[, gindx]) / 2
  for (k in seq_len(nloc)) {
    ualleles <- unique(c(m[, gindx][, 2*k-1], m[, gindx][, 2*k]))
    nalleles <- length(genoAttr$unique.alleles[[k]])
    if (length(ualleles) < nalleles)
      genoAttr$unique.alleles[[k]] <-
        genoAttr$unique.alleles[[k]][!is.na(match(seq_len(nalleles), ualleles))]
  }
  for (att in names(genoAttr)) attr(m[, gindx], att) <- genoAttr[[att]]
  attr(m, "yxmiss") <- yxmiss; attr(m, "gmiss") <- gmiss; m
}

decode_haplo_row <- function(codes, label_list) {
  if (all(codes == "*") || any(codes == "*")) return("Rare (combined)")
  parts <- character(length(codes))
  for (i in seq_along(codes)) {
    idx      <- as.numeric(codes[i])
    parts[i] <- if (!is.na(idx)) label_list[[i]][idx] else "?"
  }
  paste(parts, collapse = "-")
}

subset_geno <- function(gs, idx) {
  saved <- attributes(gs); gs2 <- gs[idx, , drop=FALSE]
  for (att in setdiff(names(saved), c("dim","dimnames")))
    attr(gs2, att) <- saved[[att]]
  gs2
}

# haplo.stats::setupGeno() requires equal-length vectors, so every SNP must be
# re-parsed on one shared row mask. haplo.em / haplo.glm handle partially-missing
# genotypes via the EM algorithm, so by default we keep every subject typed at
# >= 1 SNP (union mask) and let missing alleles be estimated — matching the
# SNPstats reference tool. When complete_case is TRUE (the completeCases option)
# we instead keep only subjects typed at ALL SNPs (intersection mask), which for
# assoc/interaction is exactly prep$complete_mask (response/covariates included).
.make_haplo_geno_list <- function(prep, complete_case = FALSE) {
  snp_vars <- prep$snp_vars
  mask <- if (isTRUE(complete_case))
    Reduce(`&`, lapply(snp_vars, function(nm) prep$snp_data[[nm]]$snp_mask))
  else
    Reduce(`|`, lapply(snp_vars, function(nm) !is.na(prep$snp_data[[nm]]$clean)))
  geno_list  <- lapply(snp_vars, function(nm) {
    sd <- prep$snp_data[[nm]]
    parse_genotype(sd$clean[mask], sd$user_levels)
  })
  names(geno_list) <- snp_vars
  geno_list <- Filter(Negate(is.null), geno_list)
  list(geno_list = geno_list, mask = mask)
}

# Rare-haplotype criterion (haploRareCriterion). haplo.freq.min and
# haplo.min.count are mutually exclusive in haplo.glm.control (passing both drops
# the count with a warning), so the user picks one and it governs every haplotype
# analysis. These helpers build the shared control and the matching display
# threshold/label so the three haplo.glm fits and the frequency table agree.
haplo_glm_control <- function(opts) {
  if (identical(opts$haploRareCriterion, "count"))
    haplo.stats::haplo.glm.control(haplo.effect = opts$haploEffect,
                                   haplo.min.count = opts$haploMinCount)
  else
    haplo.stats::haplo.glm.control(haplo.effect = opts$haploEffect,
                                   haplo.freq.min = opts$haploFreqMin)
}

haplo_rare_label <- function(opts) {
  if (identical(opts$haploRareCriterion, "count"))
    paste0("Rare (count<", opts$haploMinCount, ")")
  else
    paste0("Rare (<", opts$haploFreqMin, ")")
}

# Frequency cutoff for the haplo.em frequency table's rare-haplotype pooling.
# For the count criterion, convert the count to an expected frequency using the
# number of subjects (2n chromosomes) in the EM fit, mirroring how haplo.glm
# derives its own freq.min from a count.
haplo_rare_freq_cut <- function(opts, n_subj) {
  if (identical(opts$haploRareCriterion, "count"))
    opts$haploMinCount / (2 * n_subj)
  else
    opts$haploFreqMin
}
