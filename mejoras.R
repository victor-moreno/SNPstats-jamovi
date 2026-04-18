# improve association table as in web version
#  - p value for LRT in ref cat of codominant 

# stratified analysis:
# table labels

# can I have a new function .fill_cross_class that receives the same arguments and generates a table that populates the terms of the interaction model  so that, for a binary response, I get the OR for combinations of SNP and covar levels (assume categorical now). The reference category is the combination reference levels. In rows I get genotypes and in columns the covar categories. for each covar category I get n(%) n(%) OR (95%CI)

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
