## Metadata organization for Y2AVE Single-Cell Datasets
## Jesse Engreitz
## August 1, 2023

library(tidyr)
library(dplyr)

datasets <- read.delim("Y2AVE_SingleCellDatasets.tsv", stringsAsFactors=F)
datasetsMetadata <- read.delim("Y2AVE_SingleCellDatasets.metadata.txt", stringsAsFactors=F)
for (i in colnames(datasetsMetadata))
    if (is.character(datasetsMetadata[,i]))
        datasetsMetadata[,i] <- gsub("\n",",",datasetsMetadata[,i])


setdiff(datasets$DatasetID, datasetsMetadata$DatasetID)

datasetTable <- datasets %>%
    rename(FilteredDataSynID=synID) %>%
    merge(datasetsMetadata, by="DatasetID") %>%
    select(
        DatasetID,
        DatasetName,
        Assay,
        FilteredDataSynID,
        ProcessedDataSynID,
        Biosample,
        BiosampleOntologyTerm,
        Species,
        ContactNominatedBy,
        ContactDataCurationAndAlignment,
        ContactDataFilteringQC,
        Rationale,
        Seqspec,
        FASTQOriginalLink,
        FASTQSynID,
        Publication,
        JamboreeTeam
        )

stopifnot(nrow(datasets) == nrow(datasetTable))

## Update this file on drive by re-pulling from Synapse, then load it in here:
allFileList <- read.delim("Y2AVE_SingleCellDatasets.allFileTableDump.tsv", stringsAsFactors=F)

########################################################
## Process and organize the filtered data / pseudobulk file list

## This is the list of pseudobulk files in Y2AVE/FilteredData

fileList <- read.delim("Y2AVE_SingleCellDatasets.files.tsv", stringsAsFactors=F)
fileList <- merge(fileList, allFileList %>% rename(FileSynID=id, File=name))

fileList <- fileList %>% 
    rename(FileSynapseName=File,
           FilteredDataSynID=DatasetSynID) %>%
    mutate(FileDownloadName=basename(dataFileKey))

fileList <- fileList %>% 
    mutate(fileType=
        ifelse(grepl("ClusterAssignment.tsv",FileSynapseName), "cluster assignment table",
        ifelse(grepl("ClusterMetadata.tsv",FileSynapseName),   "cluster metadata table",
        ifelse(grepl("sorted.tagAlign.gz.tbi$",FileSynapseName),"cluster pseudobulk ATAC sorted tagAlign Tabix index",
        ifelse(grepl("sorted.tagAlign.gz$",FileSynapseName),    "cluster pseudobulk ATAC sorted tagAlign",
        ifelse(grepl("Cluster.*.atac.filter.fragments.*.tsv.gz.tbi",FileSynapseName), "cluster pseudobulk ATAC sorted tagAlign Tabix index",
        ifelse(grepl("Cluster.*.atac.filter.fragments.*.tsv.gz",FileSynapseName), "cluster pseudobulk ATAC sorted tagAlign",                     
        ifelse(grepl(".atac.filter.cutsites.*.tsv.gz.tbi",FileSynapseName), "cluster pseudobulk ATAC sorted tagAlign Tabix index",
        ifelse(grepl(".atac.filter.cutsites.*.tsv.gz",FileSynapseName), "cluster pseudobulk ATAC sorted tagAlign",                     
        ifelse(grepl(".atac.filter.cutsites.*.tagAlign.gz.tbi",FileSynapseName), "cluster pseudobulk ATAC sorted tagAlign Tabix index",
        ifelse(grepl(".atac.filter.cutsites.*.tagAlign.gz",FileSynapseName), "cluster pseudobulk ATAC sorted tagAlign",                     
        ifelse(grepl("tagAlign_Ma2020_BMMC_SingleDonor_.*tsv.gz.tbi",FileSynapseName), "cluster pseudobulk ATAC sorted tagAlign Tabix index",
        ifelse(grepl("tagAlign_Ma2020_BMMC_SingleDonor_.*tsv.gz",FileSynapseName), "cluster pseudobulk ATAC sorted tagAlign",                     
        ifelse(grepl("tagAlign",FileSynapseName),              "cluster pseudobulk ATAC unsorted tagAlign",
        ifelse(grepl(".combinedFiltered.*.rna.h5",FileSynapseName), "gene expression matrix .h5 combined and filtered",
        ifelse(grepl(".CPM.tsv",FileSynapseName),              "cluster gene expression counts per million",
            "Unrecognized"))))))))))))))))

fileList[fileList$FileDownloadName=="Ma2020_BMMC_SingleDonor_RNA_ClusterMetadata.tsv","fileType"] <- "cluster metadata table alternative"
fileList[fileList$FileDownloadName=="Ma2020_BMMC_SingleDonor_ATAC_ClusterMetadata.tsv","fileType"] <- "cluster metadata table alternative"

## Check if any files are labeled unrecognized — if so, update the logic above
fileList %>% filter(fileType=="Unrecognized")
stopifnot(!any(duplicated(fileList$FileSynapseName)))
stopifnot(!any(duplicated(fileList$FileDownloadName)))
subset(fileList, fileType == "cluster pseudobulk ATAC unsorted tagAlign" & grepl(".tbi",FileSynapseName))
#subset(fileList, fileType == "cluster pseudobulk ATAC unsorted tagAlign")

## Now merge in the "processed data" files (before cell type clustering / pseudobulking)
processedDataList <- allFileList

## Output the unique list of subpools and manually make a subpool -> DatasetID table
write.table(processedDataList %>% pull(subpool) %>% unique(), quote=F, row.names=F, col.names=F, file="SubpoolList.tsv")
subpoolToData <- read.delim("SubpoolToDataset.tsv", stringsAsFactors=F)
setdiff(processedDataList$subpool, subpoolToData$subpool)

processedDataList <- processedDataList %>% 
    filter(type=="file") %>%
    merge(subpoolToData, by="subpool") %>%
    rename(FileSynID=id,
           subpoolOrLaneID=subpool,
           FileSynapseName=name) %>%
    mutate(FileDownloadName=basename(dataFileKey),
           ProcessedDataSynID=parentId
        )

processedDataList <- processedDataList %>% 
    mutate(fileType=
        ifelse(grepl("fragments.*.tsv.gz.tbi",FileSynapseName),         "ATAC fragment file Tabix index",
        ifelse(grepl("fragments.*.tsv.gz",FileSynapseName),             "ATAC fragment file",
        ifelse(grepl("qc.rna.*.barcode.metadata.tsv",FileSynapseName),  "gene expression barcode metrics table",
        ifelse(grepl(".rna.h5",FileSynapseName),                        "gene expression matrix .h5",
        ifelse(grepl(".atac.qc.*.metadata.tsv",FileSynapseName),        "ATAC barcode metrics table",
        ifelse(grepl(".html",FileSynapseName),                          "QC HTML",
        ifelse(grepl(".joint.barcode.metadata.*.csv",FileSynapseName),  "Multiome joint barcode metrics table",
        ifelse(grepl(".sortedByCoord.out.bam",FileSynapseName),         "gene expression sorted BAM",
        ifelse(grepl(".raw.tar.gz",FileSynapseName),                    "gene expression matrix TAR",
            "Unrecognized"))))))))))
processedDataList %>% filter(fileType=="Unrecognized")
stopifnot(!any(duplicated(processedDataList$FileSynapseName)))
stopifnot(!any(duplicated(processedDataList$FileDownloadName)))


combinedFileList <- rbind(
        fileList %>% 
            mutate(ProcessedDataSynID="",subpoolOrLaneID="all") %>% 
            select(DatasetID, FileSynapseName, FileDownloadName, fileType, FileSynID, parentId, projectId, subpoolOrLaneID),
        processedDataList %>% 
            mutate(FilteredDataSynID="") %>%
            select(DatasetID, FileSynapseName, FileDownloadName, fileType, FileSynID, parentId, projectId, subpoolOrLaneID)
    ) %>% 
    arrange(DatasetID, fileType)

write.table(combinedFileList, file="Y2AVE_SingleCellDatasets.fileList.tsv", sep='\t', quote=F, col.names=T, row.names=F)
write.table(datasetTable, file="Y2AVE_SingleCellDatasets.datasetList.tsv", sep='\t', quote=F, col.names=T, row.names=F)


