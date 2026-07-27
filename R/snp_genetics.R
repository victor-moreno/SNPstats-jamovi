
# ══════════════════════════════════════════════════════════════════════════════
# snp_genetics.R — biallelic genotype handling, HWE and LD, without `genetics`
#
# WHY THIS FILE EXISTS
#
# The `genetics` package is marked obsolete upstream. jamovi resolves module
# dependencies from a pinned package snapshot, so if a deprecated package stops
# being carried, the module does not degrade — it fails to install. These are
# the only four things SNPstats ever used from it (genotype, allele, HWE.exact,
# LD), so the exposure was never worth the risk.
#
# PROVENANCE
#
# Written from the published statistical definitions, not ported from
# `genetics` (which is GPL). The methods are standard and long predate any
# implementation of them:
#
#   HWE exact test  Wigginton JE, Cutler DJ, Abecasis GR (2005)
#                   A note on exact tests of Hardy-Weinberg equilibrium.
#                   Am J Hum Genet 76:887-893.
#   Haplotype EM    Excoffier L, Slatkin M (1995) Maximum-likelihood estimation
#                   of molecular haplotype frequencies in a diploid population.
#                   Mol Biol Evol 12:921-927.
#   D, D'           Lewontin RC (1964) Genetics 49:49-67.
#   r               Hill WG, Robertson A (1968) Theor Appl Genet 38:226-231.
#
# The output conventions below (allele ordering, table shapes, the error on a
# monomorphic HWE test) were matched to what the module previously observed, so
# result tables and golden values are unchanged. They were established by
# running `genetics` and recording its behaviour, not by reading its code.
# ══════════════════════════════════════════════════════════════════════════════


# ── Genotype objects ─────────────────────────────────────────────────────────

#' Parse "A/B" genotype strings into an internal biallelic genotype object.
#'
#' Alleles are ordered by descending count, ties broken by allele name
#' descending; every genotype string is then written in that order, so a
#' heterozygote is always "<major>/<minor>". Rows with NA stay NA.
#' @return object of class "snpgeno", or stops if nothing is parseable.
snp_genotype <- function(x, sep = "/") {
  x <- as.character(x)
  parts <- strsplit(x, sep, fixed = TRUE)
  ok    <- !is.na(x) & vapply(parts, length, 1L) == 2L

  a <- rep(NA_character_, length(x))
  b <- rep(NA_character_, length(x))
  if (any(ok)) {
    a[ok] <- vapply(parts[ok], `[`, "", 1L)
    b[ok] <- vapply(parts[ok], `[`, "", 2L)
  }

  cnt <- table(c(a[ok], b[ok]))
  if (length(cnt) == 0L) stop("no parseable genotypes")
  # Descending count; ties resolved by descending allele name. sort() on the
  # names first makes the subsequent order() tie-break deterministic.
  nms     <- sort(names(cnt), decreasing = TRUE)
  alleles <- nms[order(-as.integer(cnt[nms]))]

  # Canonicalise each genotype to allele order
  ia <- match(a, alleles); ib <- match(b, alleles)
  swap <- !is.na(ia) & !is.na(ib) & ia > ib
  if (any(swap)) { tmp <- a[swap]; a[swap] <- b[swap]; b[swap] <- tmp[]; }

  geno <- ifelse(is.na(a) | is.na(b), NA_character_, paste0(a, "/", b))

  structure(list(geno = geno, a1 = a, a2 = b, alleles = alleles),
            class = "snpgeno")
}

#' @export
as.character.snpgeno <- function(x, ...) x$geno

#' @export
length.snpgeno <- function(x) length(x$geno)

#' Subset a genotype object. Allele ordering is NOT recomputed: a stratum of a
#' SNP must keep the parent's allele labels or its genotype rows would not line
#' up with the overall table.
#' @export
`[.snpgeno` <- function(x, i, ...) {
  structure(list(geno = x$geno[i], a1 = x$a1[i], a2 = x$a2[i],
                 alleles = x$alleles),
            class = "snpgeno")
}

#' Allele matrix (one row per individual, two columns), for haplo.stats.
snp_allele <- function(g) cbind(g$a1, g$a2)

#' Genotype and allele frequency tables.
#'
#' `genotype.freq` has one row per observed genotype, ordered a1/a1, a1/a2,
#' a2/a2 by the object's allele order; `allele.freq` one row per allele in that
#' order. A trailing "NA" row is added only when something is missing, counting
#' individuals (genotypes) or alleles (2 per missing individual), with
#' Proportion NA. Proportions are over non-missing.
#'
#' Unlike `genetics`, a single observed genotype still returns a 1-row matrix
#' rather than collapsing to a named vector; callers already guard for both.
#' @export
summary.snpgeno <- function(object, ...) {
  geno <- object$geno
  obs  <- !is.na(geno)
  n_typed <- sum(obs)
  n_miss  <- sum(!obs)

  al <- object$alleles
  # Expected genotype labels in allele order: a1/a1, a1/a2, ..., ak/ak
  lab <- character(0)
  for (i in seq_along(al)) for (j in seq(i, length(al)))
    lab <- c(lab, paste0(al[i], "/", al[j]))
  lab <- lab[lab %in% geno[obs]]

  gcount <- as.integer(table(factor(geno[obs], levels = lab)))
  gf <- cbind(Count = gcount,
              Proportion = if (n_typed > 0) gcount / n_typed else rep(NA_real_, length(lab)))
  rownames(gf) <- lab

  acount <- as.integer(table(factor(c(object$a1[obs], object$a2[obs]), levels = al)))
  af <- cbind(Count = acount,
              Proportion = if (n_typed > 0) acount / (2 * n_typed) else rep(NA_real_, length(al)))
  rownames(af) <- al

  if (n_miss > 0) {
    gf <- rbind(gf, "NA" = c(n_miss, NA_real_))
    af <- rbind(af, "NA" = c(2 * n_miss, NA_real_))
  }

  list(genotype.freq = gf, allele.freq = af, n.typed = n_typed)
}


# ── Hardy-Weinberg exact test ────────────────────────────────────────────────

#' Exact test of Hardy-Weinberg equilibrium at a biallelic locus.
#'
#' Conditioning on the observed allele counts, the probability of n12
#' heterozygotes among n individuals is
#'
#'   P(n12 | n, n1) = n! / (n11! n12! n22!) * 2^n12 * n1! n2! / (2n)!
#'
#' and the two-sided p-value is the total probability of every attainable
#' heterozygote count no more probable than the observed one (Wigginton et al.
#' 2005). Evaluated in logs so large samples do not overflow.
#'
#' Stops on a non-biallelic locus, matching what the callers already expect:
#' they wrap this in tryCatch and treat the failure as "no HWE result", which
#' is the right outcome for a monomorphic SNP.
snp_hwe_exact <- function(g) {
  sm <- summary(g)
  af <- sm$allele.freq
  al <- rownames(af)[rownames(af) != "NA"]
  if (length(al) != 2L)
    stop("Exact HWE test can only be computed for 2 markers with 2 alleles")

  geno <- g$geno[!is.na(g$geno)]
  hom1 <- paste0(al[1], "/", al[1])
  hom2 <- paste0(al[2], "/", al[2])
  het  <- paste0(al[1], "/", al[2])

  n11 <- sum(geno == hom1)
  n12 <- sum(geno == het)
  n22 <- sum(geno == hom2)
  if (n11 + n12 + n22 != length(geno))
    stop("Exact HWE test can only be computed for 2 markers with 2 alleles")

  list(p.value = hwe_exact_p(n11, n12, n22),
       observed = c(n11 = n11, n12 = n12, n22 = n22))
}

#' Two-sided HWE exact p-value from the three genotype counts.
hwe_exact_p <- function(n11, n12, n22) {
  n  <- n11 + n12 + n22
  n1 <- 2 * n11 + n12          # copies of allele 1
  n2 <- 2 * n22 + n12
  if (n == 0 || n1 == 0 || n2 == 0) return(NA_real_)

  # n12 shares parity with n1 and cannot exceed either allele's total
  ks  <- seq(n1 %% 2, min(n1, n2), by = 2)
  k11 <- (n1 - ks) / 2
  k22 <- (n2 - ks) / 2

  logp <- lfactorial(n) - lfactorial(k11) - lfactorial(ks) - lfactorial(k22) +
          ks * log(2) + lfactorial(n1) + lfactorial(n2) - lfactorial(2 * n)
  p <- exp(logp - max(logp))
  p <- p / sum(p)

  obs <- p[ks == n12]
  if (length(obs) != 1L) return(NA_real_)
  # 1e-9 slack so a tie in probability is counted, not lost to rounding
  sum(p[p <= obs * (1 + 1e-9)])
}


# ── Pairwise linkage disequilibrium ──────────────────────────────────────────

#' Pairwise LD between two biallelic loci from unphased genotypes.
#'
#' Every two-locus genotype has an unambiguous haplotype composition except the
#' double heterozygote, which is a mixture of the coupling (AB/ab) and repulsion
#' (Ab/aB) phases. EM splits the double heterozygotes between the two phases in
#' proportion to their current probabilities and re-counts until the haplotype
#' frequencies stop moving (Excoffier & Slatkin 1995).
#'
#' With A and B the first (most frequent) allele at each locus:
#'   D  = p(AB) - p(A)p(B)
#'   D' = D / Dmax, Dmax the bound imposed by the allele frequencies
#'   r  = D / sqrt(p(A)p(a)p(B)p(b));  X2 = 2n r^2 on 1 df
#'
#' @return list with D, D', r, R^2, n, X^2 and P-value, or NULL if either locus
#'   is not biallelic or nothing is jointly typed.
snp_ld <- function(g1, g2) {
  keep <- !is.na(g1$geno) & !is.na(g2$geno)
  n    <- sum(keep)
  if (n == 0L) return(NULL)

  a1 <- g1$a1[keep]; a2 <- g1$a2[keep]
  b1 <- g2$a1[keep]; b2 <- g2$a2[keep]

  # Allele order is recomputed on the jointly-typed subset, so the reported D
  # refers to the major allele among the individuals actually used.
  ca <- sort(table(c(a1, a2)), decreasing = TRUE)
  cb <- sort(table(c(b1, b2)), decreasing = TRUE)
  if (length(ca) != 2L || length(cb) != 2L) return(NULL)
  A <- names(ca)[1]; B <- names(cb)[1]

  nA <- (a1 == A) + (a2 == A)          # copies of A, 0/1/2
  nB <- (b1 == B) + (b2 == B)
  pA <- sum(nA) / (2 * n)
  pB <- sum(nB) / (2 * n)

  tab <- table(factor(nA, 0:2), factor(nB, 0:2))

  # Unambiguous haplotype counts; the (1,1) cell is the ambiguous double het
  cAB <- 2 * tab["2","2"] + tab["2","1"] + tab["1","2"]
  cAb <- 2 * tab["2","0"] + tab["2","1"] + tab["1","0"]
  caB <- 2 * tab["0","2"] + tab["0","1"] + tab["1","2"]
  cab <- 2 * tab["0","0"] + tab["0","1"] + tab["1","0"]
  ndh <- tab["1","1"]

  p <- c(pA * pB, pA * (1 - pB), (1 - pA) * pB, (1 - pA) * (1 - pB))
  for (it in seq_len(1000L)) {
    coup <- p[1] * p[4]
    repu <- p[2] * p[3]
    w    <- if (coup + repu > 0) coup / (coup + repu) else 0.5
    cnt  <- c(cAB + ndh * w, cAb + ndh * (1 - w),
              caB + ndh * (1 - w), cab + ndh * w)
    new  <- cnt / sum(cnt)
    if (max(abs(new - p)) < 1e-12) { p <- new; break }
    p <- new
  }

  D    <- p[1] - pA * pB
  Dmax <- if (D > 0) min(pA * (1 - pB), (1 - pA) * pB)
          else       max(-pA * pB, -(1 - pA) * (1 - pB))
  Dprime <- if (Dmax != 0) D / Dmax else NA_real_
  den <- sqrt(pA * (1 - pA) * pB * (1 - pB))
  r   <- if (den > 0) D / den else NA_real_
  X2  <- if (is.na(r)) NA_real_ else 2 * n * r^2

  out <- list(D, Dprime, r, r^2, n, X2,
              if (is.na(X2)) NA_real_ else stats::pchisq(X2, 1, lower.tail = FALSE))
  names(out) <- c("D", "D'", "r", "R^2", "n", "X^2", "P-value")
  out
}
