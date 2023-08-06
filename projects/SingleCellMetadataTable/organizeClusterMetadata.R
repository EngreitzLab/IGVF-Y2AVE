#######################################################
## Make a combined cluster metadata table across all clusters
## Jesse Engreitz
## August 1, 2023

library(tidyr)
library(dplyr)

combinedFileList <- read.delim("Y2AVE_SingleCellDatasets.fileList.tsv", stringsAsFactors=F)
datasetTable <- read.delim("Y2AVE_SingleCellDatasets.datasetList.tsv", stringsAsFactors=F)


clusterFiles <- combinedFileList %>% filter(fileType=="cluster metadata table")

readClusterMetadataTable <- function(file) {
	x <- read.delim(paste0("clusterFiles/",file), stringsAsFactors=F)
}