

sample1 <- read.delim("../sample1.txt", header = TRUE)
save(sample1, file = "data/sample1.rda", version = 2)

sample2 <- read.delim("../sample2.txt", header = TRUE)
save(sample2, file = "data/sample2.rda", version = 2)

sample3 <- read.delim("../sample3.txt", header = TRUE)
save(sample3, file = "data/sample3.rda", version = 2) 

sample4 <- read.delim("../sample4.txt", header = TRUE)
save(sample4, file = "data/sample4.rda", version = 2)

sample5 <- read.delim("../sample5.txt", header = TRUE)
save(sample5, file = "data/sample5.rda", version = 2)


jmvtools::install()
