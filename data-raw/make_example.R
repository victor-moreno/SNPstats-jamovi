## Generate example dataset matching the SNPstats web app sample format
## Run once to produce SNPstats_example.csv

set.seed(42)
n <- 200

# Helper: simulate bi-allelic SNP under HWE
sim_snp <- function(n, p, alleles = c("A", "G")) {
  a1 <- alleles[1]; a2 <- alleles[2]
  genos <- c(
    paste0(a1, "/", a1),   # freq p^2
    paste0(a1, "/", a2),   # freq 2p(1-p)
    paste0(a2, "/", a2)    # freq (1-p)^2
  )
  probs <- c(p^2, 2*p*(1-p), (1-p)^2)
  g <- sample(genos, n, replace = TRUE, prob = probs)
  # Add ~5% missings
  g[sample(n, floor(n * 0.05))] <- NA
  g
}

status <- sample(c("Case", "Control"), n, replace = TRUE)
sex    <- sample(c("Male", "Female"), n, replace = TRUE)
age    <- round(rnorm(n, mean = 55, sd = 12))

df <- data.frame(
  id     = 1:n,
  snp1   = sim_snp(n, p = 0.3, alleles = c("C", "T")),
  snp2   = sim_snp(n, p = 0.4, alleles = c("A", "G")),
  snp3   = sim_snp(n, p = 0.2, alleles = c("G", "T")),
  snp4   = sim_snp(n, p = 0.5, alleles = c("C", "G")),
  status = status,
  sex    = sex,
  age    = age,
  stringsAsFactors = FALSE
)
dir.create("data", showWarnings = FALSE)
write.csv(df, "data/SNPstats_example.csv", row.names = FALSE, na = "NA")
message("Written data/SNPstats_example.csv with ", n, " rows.")
