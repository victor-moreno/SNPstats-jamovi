.libPaths(c("/Users/h501uvma/R/library/4.5",.libPaths()))

library(genetics)
library(haplo.stats)

# ── load your data ─────────────────────────────────────────────
data <- read.delim("../sample1.txt")   # adjust as needed
head(data)

# ── set these to match your actual column names ────────────────
snp_vars       <- c("SNP1", "SNP2", "SNP3")   # your SNP columns
response_var   <- "STATUS"
int_var        <- "SEX"                         # interaction covariate

# ── helpers exactly as in the module ─────────────────────────
parse_genotype <- function(x) {
  tryCatch(genetics::genotype(as.character(x), sep = "/"), error = function(e) NULL)
}

subset_geno <- function(gs, idx) {
  saved_attr <- attributes(gs)
  gs2        <- gs[idx, , drop = FALSE]
  for (att in setdiff(names(saved_attr), c("dim", "dimnames")))
    attr(gs2, att) <- saved_attr[[att]]
  gs2
}

# ── build geno_setup exactly as in the module ─────────────────
geno_list  <- lapply(snp_vars, function(v) parse_genotype(data[[v]]))
names(geno_list) <- snp_vars
allele_mat <- do.call(cbind, lapply(geno_list, function(g) genetics::allele(g)))
geno_setup <- haplo.stats::setupGeno(allele_mat, locus.label = snp_vars)

# ── build model data frame exactly as in the module ───────────
response <- as.integer(as.factor(data[[response_var]])) - 1L
keep     <- !is.na(response) & !is.na(data[[int_var]])

m <- data.frame(y = as.numeric(as.factor(response[keep])) - 1L)
m[[int_var]] <- data[[int_var]][keep]
m$geno <- subset_geno(geno_setup, keep)

m<- m[complete.cases(m), ]

# ── fit ───────────────────────────────────────────────────────
fit <- haplo.stats::haplo.glm(
  as.formula(paste("y ~ geno *", int_var)),
  family    = "binomial",
  data      = m,
  na.action = haplo.stats::na.geno.keep,
  control   = haplo.stats::haplo.glm.control(haplo.freq.min = 0.05)
)

fit0 <- haplo.stats::haplo.glm(
  as.formula(paste("y ~ geno +", int_var)),
  family    = "binomial",
  data      = m,
  na.action = haplo.stats::na.geno.keep,
  control   = haplo.stats::haplo.glm.control(haplo.freq.min = 0.05)
)

fit
fit0

anova(fit0, fit)

─ debug dump ────────────────────────────────────────────────
cat("=== coef rownames ===\n");    print(rownames(summary(fit)$coefficients))
cat("\n=== haplo.common ===\n");   print(fit$haplo.common)
cat("\n=== haplo.base ===\n");     print(fit$haplo.base)
cat("\n=== haplo.names ===\n");    print(fit$haplo.names)
cat("\n=== haplo.unique ===\n");   print(fit$haplo.unique)
cat("\n=== names(fit) ===\n");     print(names(fit))
