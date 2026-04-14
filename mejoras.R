
# check bugs when bad data: zeros, monoallelic, no BB, ...

# replace values > 1000 by Inf or ---

# Genotype frequencies table has error object snp_char not found


# Haplotype frequencies: stratify by response

# Lets revise missing values in response, SNPs and covariates. I want that any missing values are excluded from analysis, but informed.

# SNP summary: add a column with missing values for each SNP. If covariates or response added, Add a note indicating: "Covariates <>and response considered. N observations missing excluded." Only show if missings really exist.

# For SNP associacion and SNP interaction tables, add to Note Model adjusted, an extension "<> observations excluded due to missing values." if any. Include response when counting missings.

# The same note for haplotype frequency, association and interaction tables.For haplotye frequencies, separate missings for SNPs and covariates/response.

# Design
# 	Separate PRS analysis

# 	Consider separating LD + haplotype analysis
	
	
# add reference haplotyype to interaction table




options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2025-05-25"))

install.packages(c("remotes","jmvcore"))
install.packages('jmvtools', repos='https://repo.jamovi.org')
remotes::install_version("node", version="1.2", repos="https://repo.jamovi.org")

install.packages(c("genetics","haplo.stats","arsenal","combinat","gdata","rms","polspline"))

.libPaths(c(
"/Users/h501uvma/Downloads/claude-code/snpstats/jamovi/SNPstats/build/R4.5.0-x64-macos",
 "~/R/library/4.5"
))

.libPaths()
getwd()
jmvtools::install()
