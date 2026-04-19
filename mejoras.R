
# instead of generating the association for the first covariate, populate a drodown with all covariates and let the user choose  

# in this jamovi module, fix the table generated in .fill_cross_class. Now is a modified copy of .fill_strat_by_covariate. The requested fix is that the model terms used for the strata of the covariate, except the reference, needs to be calculated so that the values have a common reference for SNP and covariate. Now, the model used is conditional, but it should be the snp*covariate and calculate the ORs/betas for the levels 2+ of the covariate as beta(snp)+beta(covar)+beta(snp:covar), calculate the appropriate standard error, and exponentiate if response is binary.

# haplotype analysis: 
# show N in Note, in addition to missings
# sort by frequency (freq, assoc & interaction)
# replace NA by ''
# add LRT p value for haplotype association -> first row & note
# add stratified tables for haplotype interaction
# OR are different from web version, check code

# show target response level if binary

# missing values in SNPs are not removed for covariate descriptives, to match 
# the web version. This means that the covariate descriptives are not necessarily 
# on the same subset as the SNP descriptives.
# improvement:
#  - add option to exclude SNP missings in covariate descriptives, and show N in note


# simplify code
# find duplicated code and improve modularity. For example strata labels & missing values for SNP results

# check bugs when bad data: zeros, monoallelic, no BB, ...
# replace values > 1000 by Inf or ---
# Note in association: risk response category
# colorize results
# hide all notes except LRT


# Design
# 	Separate PRS analysis

# 	Consider separating LD + haplotype analysis
	
	




options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2025-05-25"))

install.packages(c("remotes","jmvcore"))
install.packages('jmvtools', repos='https://repo.jamovi.org')
remotes::install_version("node", version="1.2", repos="https://repo.jamovi.org")

install.packages(c("genetics","haplo.stats","arsenal","combinat","gdata","rms","polspline"))

.libPaths(c(
"/Users/h501uvma/Downloads/SNPstats/build/R4.5.0-x64-macos",
 "~/R/library/4.5"
))

.libPaths(c(
"/Applications/jamovi.app/Contents/Resources/modules/base/R",
"/Users/h501uvma/Downloads/SNPstats/build/R4.5.0-arm64-macos"
))

setwd("/Users/h501uvma/Downloads/SNPstats")
source("~/.Rprofile")

.libPaths()
getwd()
jmvtools::install()
