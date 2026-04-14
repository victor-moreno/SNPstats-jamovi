# missing response/covariate adjust "typed samples" and exclude from all analyses.

# missing notes: combine with adjusted. 
# haplotype frequencies stratified using 2 columns, not additional rows

# add reference haplotyye to interaction table


# check bugs when bad data: zeros, monoallelic, no BB, ...
# replace values > 1000 by Inf or ---


# Design
# 	Separate PRS analysis

# 	Consider separating LD + haplotype analysis
	
	




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
