setwd("/Users/h501uvma/Downloads/claude-code/snpstats-jamovi/SNPstats")


data<-read.delim("data/CRCgenet-SNPs.tsv", header=TRUE, stringsAsFactors = TRUE)
names(data)
#  [1] "phenotype"  "sex"        "age"        "bmi"        "bmi4"       "bmiOMS"     "fit"       
#  [8] "rs12080929" "rs10911251" "rs72647484" "rs6691170"  "rs992157"   "rs11676348" "rs11903757"
# [15] "rs10936599" "rs35360328" "rs812481"   "rs7136702"  "rs3987"     "rs2736100"  "rs647161"  
# [22] "rs4711689"  "rs1321311"  "rs719725"   "rs16892766" "rs6983267"  "rs6469656"  "rs11987193"
# [29] "rs10904849" "rs10795668" "rs704017"   "rs11196172" "rs1035209"  "rs12241008" "rs11190164"
# [36] "rs4246215"  "rs3802842"  "rs174537"   "rs1535"     "rs174550"   "rs3824999"  "rs73208120"
# [43] "rs3217901"  "rs11169552" "rs59336"    "rs10849432" "rs3184504"  "rs11064437" "rs10774214"
# [50] "rs2238126"  "rs3217810"  "rs4444235"  "rs17094983" "rs1957636"  "rs4779584"  "rs9929218" 
# [57] "rs16941835" "rs12603526" "rs4939827"  "rs12970291" "rs3764482"  "rs10411210" "rs2241714" 
# [64] "rs1800469"  "rs2423279"  "rs961253"   "rs6066825"  "rs6061231"  "rs4925386"  "rs4813802" 
# [71] "rs5934683"

# function (data, response = NULL, snps, covariates = NULL, responseType = "auto", 
#     subpop = FALSE, covDesc = FALSE, rmSnpMissing = FALSE, snpSummary = TRUE, 
#     allFreq = FALSE, genoFreq = FALSE, hweTest = FALSE, showMissing = FALSE, 
#     showMissingnessPlot = FALSE, missingnessThreshold = 0.1, 
#     snpAssoc = FALSE, modelCodominant = TRUE, modelDominant = FALSE, 
#     modelRecessive = FALSE, modelOverdominant = FALSE, modelLogAdditive = FALSE, 
#     ciWidth = 95, showAIC = FALSE, snpInteraction = FALSE, interactionType = "multiplicative", 
#     showInteractionTable = FALSE, showInteractionAdjVars = FALSE, 
#     showStratByCovariate = FALSE, showStratByGenotype = FALSE, 
#     showCrossClassTable = FALSE, ldAnalysis = FALSE, ldMatrix = FALSE, 
#     ldMetric = "r2", ldPlot = FALSE, haploFreq = FALSE, haploFreqMin = 0.01, 
#     ldSubpop = FALSE, haploAssoc = FALSE, haploInteraction = FALSE) 

data<-read.delim("data/CRCgenet-SNPs.tsv", header=TRUE, stringsAsFactors = TRUE)

SNPstats::snpStats(data=data, response = "phenotype", snps=names(data)[8:9], covariates = c("sex", "age","bmi","bmi4"), responseType = "auto", 
        covDesc = TRUE, snpSummary = FALSE)

SNPstats::snpStats(data=data, response = "phenotype", snps=names(data)[8:9], covariates = c("sex", "age","bmi","bmi4"), responseType = "auto", 
        covDesc = TRUE, subpop = TRUE, snpSummary = FALSE)



SNPstats::snpStats(data=data, response = "phenotype", snps=names(data)[8:9], covariates = c("sex", "age"), responseType = "auto", 
        snpSummary = TRUE) #, covDesc = TRUE, subpop = TRUE)

SNPstats::snpStats(data=data, response = "phenotype", snps=names(data)[8:9], covariates = c("sex", "age"), 
        snpSummary = TRUE, subpop = TRUE)






SNPstats::snpStats(data=data, response = "phenotype", snps=names(data)[8:9], covariates = c("sex", "age","bmi","bmi4"), responseType = "auto", 
        snpAssoc = TRUE, modelCodominant = TRUE, modelDominant = TRUE, modelRecessive = TRUE, modelOverdominant = TRUE, modelLogAdditive = TRUE) #, covDesc = TRUE, subpop = TRUE)


jmvtools::install()
