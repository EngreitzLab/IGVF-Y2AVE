## Metadata organization for Y2AVE Single-Cell Datasets
## Jesse Engreitz
## August 1, 2023

library(tidyr)
library(dplyr)

datasets <- read.delim("Y2AVE_SingleCellDatasets.tsv", stringsAsFactors=F)
datasetsMetadata <- read.delim("Y2AVE_SingleCellDatasets.metadata.txt", stringsAsFactors=F)

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
        ))

stopifnot(nrow(datasets) == nrow(datasetTable))


########################################################
## Process and organize the file list

## This is the list of pseudobulk files in Y2AVE/FilteredData
fileList <- read.delim("Y2AVE_SingleCellDatasets.files.tsv", stringsAsFactors=F)
fileList <- fileList %>% 
    rename(FilteredDataSynID=DatasetSynID) %>%
    mutate(fileType=
        ifelse(grepl("ClusterAssignment.tsv",File), "cluster assignment table",
        ifelse(grepl("ClusterMetadata.tsv",File),   "cluster metadata table",
        ifelse(grepl("sorted.tagAlign.gz.tbi",File),"cluster pseudobulk ATAC sorted tagAlign Tabix index",
        ifelse(grepl("sorted.tagAlign.gz",File),    "cluster pseudobulk ATAC sorted tagAlign",
        ifelse(grepl("tagAlign",File),              "cluster pseudobulk ATAC unsorted tagAlign",
            "Unrecognized"))))))

## Check if any files are labeled unrecognized — if so, update the logic above
fileList %>% filter(fileType=="Unrecognized")
stopifnot(!any(duplicated(fileList$File)))

## Now merge in the "processed data" files (before cell type clustering / pseudobulking)
processedDataList <- read.delim("SYNAPSE_TABLE_QUERY_127648913.csv", sep=',', stringsAsFactors=F)

## Output the unique list of subpools and manually make a subpool -> DatasetID table
write.table(processedDataList %>% pull(subpool) %>% unique(), quote=F, row.names=F, col.names=F, file="SubpoolList.tsv")
subpoolToData <- read.delim("SubpoolToDataset.tsv", stringsAsFactors=F)

processedDataList <- processedDataList %>% 
    filter(type=="file") %>%
    merge(subpoolToData, by="subpool") %>%
    rename(FileSynID=id,
           ProcessedDataSynID=parentId,
           subpoolOrLaneID=subpool,
           File=name)

processedDataList <- processedDataList %>% 
    mutate(fileType=
        ifelse(grepl("fragments.*.tsv.gz.tbi",File),         "ATAC fragment file Tabix index",
        ifelse(grepl("fragments.*.tsv.gz",File),             "ATAC fragment file",
        ifelse(grepl("qc.rna.*.barcode.metadata.tsv",File),  "gene expression barcode metrics table",
        ifelse(grepl(".rna.h5",File),                        "gene expression matrix .h5",
        ifelse(grepl(".atac.qc.*.metadata.tsv",File),        "ATAC barcode metrics table",
        ifelse(grepl(".html",File),                          "QC HTML",
        ifelse(grepl(".joint.barcode.metadata.*.csv",File),  "Multiome joint barcode metrics table",
        ifelse(grepl(".sortedByCoord.out.bam",File),         "gene expression sorted BAM",
        ifelse(grepl(".raw.tar.gz",File),                    "gene expression matrix TAR",
            "Unrecognized"))))))))))
processedDataList %>% filter(fileType=="Unrecognized")
stopifnot(!any(duplicated(processedDataList$File)))

combinedFileList <- rbind(
        fileList %>% 
            mutate(ProcessedDataSynID="",subpoolOrLaneID="all") %>% 
            select(DatasetID, File, fileType, FileSynID, ProcessedDataSynID, FilteredDataSynID, subpoolOrLaneID),
        processedDataList %>% 
            mutate(FilteredDataSynID="") %>%
            select(DatasetID, File, fileType, FileSynID, ProcessedDataSynID, FilteredDataSynID, subpoolOrLaneID)
    ) %>% 
    arrange(DatasetID, fileType)

write.table(combinedFileList, file="Y2AVE_SingleCellDatasets.fileList.tsv", sep='\t', quote=F, col.names=T, row.names=F)
write.table(datasetTable, file="Y2AVE_SingleCellDatasets.datasetList.tsv", sep='\t', quote=F, col.names=T, row.names=F)


