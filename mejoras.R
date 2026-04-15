# missings
# Add a check to report SNP missings in descriptive tables as additional row
# Notes: report covariate, SNPs & adjusted. 

# check bugs when bad data: zeros, monoallelic, no BB, ...
# replace values > 1000 by Inf or ---
# Note in association: risk response category
# colorize results
# hide all notes except LRT

# simplify code


# References 
#    Sole et al. 2017 for reporting guidelines
#    tutorial

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
