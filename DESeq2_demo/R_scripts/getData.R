# getData.R
# script to get data from airway package

# Data: RNA-seq performed on 4 primary human airway smooth muscle cells lines
# treated with 1 micromolar dexamethasone for 18 hours. Each of the four cell 
# lines has a treated and untreated sample 

library(airway)

data(airway)

# get experimental design data
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('utr', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file='~/Desktop/code_demos/DESeq2_demo/Data/sample_info.csv', 
            sep = ',', col.names = T, row.names = T, quote = F)

# get read counts data
countsData <- assay(airway)
write.table(countsData, file='~/Desktop/code_demos/DESeq2_demo/Data/counts_data.csv', 
            sep = ',', col.names = T, row.names = T, quote = F)

