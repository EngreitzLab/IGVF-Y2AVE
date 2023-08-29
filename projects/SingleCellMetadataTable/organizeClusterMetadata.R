#######################################################
## Make a combined cluster metadata table across all clusters
## Jesse Engreitz
## August 1, 2023

library(tidyr)
library(dplyr)

combinedFileList <- read.delim("Y2AVE_SingleCellDatasets.fileList.tsv", stringsAsFactors=F)
datasetTable <- read.delim("Y2AVE_SingleCellDatasets.datasetList.tsv", stringsAsFactors=F)

clusterFiles <- combinedFileList %>% filter(fileType=="cluster metadata table")


findDatasetFile <- function(datasetID, filetype, fileName=NA, id="FileSynID") {
    result <- subset(combinedFileList, fileType == filetype & DatasetID == datasetID)
    if (nrow(result) > 1) {
        print("Found multiple options for ",filetype," for ",datasetID)
        print(result)
        stop()
    } else if (nrow(result) == 0) {
        return("")
    } else {
        return(result[,id])
    }
}

findTagAlignFiles <- function(CellClusterID, datasetID, filetype, fileName=NA, id="FileSynID") {

    ## Handle Eila's tagAlign files, where CellClusterID does not match the filename; in this case, use filename she provided
    if (!is.na(fileName) & fileName != "") 
        result <- subset(combinedFileList, grepl(gsub(".tsv","",fileName), FileSynapseName) & fileType==filetype)

    ## Look for McGinnis dataset files (made by Chloe, but where CellClusterID from Chris is weird)
    else if (grepl("mcginnis_dataset", CellClusterID))
        result <- subset(combinedFileList, 
            fileType==filetype & DatasetID==datasetID & 
            grepl(paste0(gsub("mcginnis_dataset", "CharacterizationMcGinnis_Dataset", CellClusterID),".atac.filter"), FileSynapseName, ignore.case=TRUE))

    ## Look for file by string matching the CellClusterID with a filename
    else
        result <- subset(combinedFileList, 
            fileType==filetype & DatasetID==datasetID & 
            grepl(paste0("tagAlign_",CellClusterID,"\\."), FileSynapseName, ignore.case=T))

    ## Look for file replacing "SeuratCluster" with "Cluster", a la Chloe's files
    if (nrow(result) == 0)
        result <- subset(combinedFileList, 
            fileType==filetype & DatasetID==datasetID & 
            grepl(paste0(gsub("SeuratCluster","Cluster",CellClusterID),"\\."), FileSynapseName, ignore.case=T))

    ## Look for file by string matching the CellClusterID with the Download file name, rather than Synapse file name, just in case
    if (nrow(result) == 0)
        result <- subset(combinedFileList, 
            fileType==filetype & DatasetID==datasetID & 
            grepl(paste0("tagAlign_",CellClusterID,"\\."), FileDownloadName))

    ## Else conclude that file does not exist
    if (nrow(result) == 0)
        result <- ""
    else
        result <- result[,id]

    ## Catch cases where 2 or more files are found - this shouldn't happen
    if (length(result) != 1) {
        print(length(result))
        print(result)
        stopifnot(length(result) == 1)
    }

    return(result)
}



readClusterMetadataTable <- function(metadataFile, datasetID) {
    print(metadataFile)    
    print(datasetID)
    canonicalColumns <- c("CellClusterID", "ManualAnnotationLabel", "nCells", "MeanRNAUMIsPerCell","MeanATACFragmentsPerCell")
    x <- read.delim(paste0("clusterTables/",metadataFile), stringsAsFactors=F)

    if (any(!(canonicalColumns %in% colnames(x)))) {
        missingColumns <- setdiff(canonicalColumns, colnames(x))
        if (missingColumns == "ManualAnnotationLabel") 
            x$ManualAnnotationLabel=""
        else
            stop(paste(
                "Error reading cluster metadata table file:",metadataFile,".\n",
                "Found columns: ", paste(colnames(x), collapse=' '),"\n",
                "Missing columns: ", paste(missingColumns, collapse=' ')))
    }

    extraColumns <- setdiff(colnames(x), canonicalColumns)
    if (!("tagAlignFile" %in% colnames(x))) 
        x$tagAlignFile <- ""

    ## Adjust manual annotation labels for Eila's ENCODE samples to include TissueName
    ## if ("TissueName" %in% colnames(x))
    ##   x <- x %>% mutate(ManualAnnotationLabel=gsub("_", " ", paste0(ManualAnnotationLabel," from ",TissueName)))

    x$Biosample <- subset(datasetTable, DatasetID==datasetID)$Biosample
    x$Assay <- subset(datasetTable, DatasetID==datasetID)$Assay
    x$Species <- subset(datasetTable, DatasetID==datasetID)$Species
    x$ATACtagAlignUnsorted <- mapply(findTagAlignFiles, x$CellClusterID, datasetID, "cluster pseudobulk ATAC unsorted tagAlign", x$tagAlignFile)
    x$ATACtagAlignSorted <- mapply(findTagAlignFiles, x$CellClusterID, datasetID, "cluster pseudobulk ATAC sorted tagAlign", x$tagAlignFile)
    x$ATACtagAlignSortedIndex <- mapply(findTagAlignFiles, x$CellClusterID, datasetID, "cluster pseudobulk ATAC sorted tagAlign Tabix index", x$tagAlignFile)
    x$GeneExpressionCPM <- findDatasetFile(datasetID, "cluster gene expression counts per million")
    x$GeneExpressionFilteredH5 <- findDatasetFile(datasetID, "gene expression matrix .h5 combined and filtered")
    x$PeaksMACS2 <- "In progress"
    x$DatasetID <- datasetID

    tagAlignCols <- c("Biosample","Assay","Species","ATACtagAlignUnsorted","ATACtagAlignSorted","ATACtagAlignSortedIndex","GeneExpressionCPM","GeneExpressionFilteredH5","PeaksMACS2")
    res <- x[,c("DatasetID", canonicalColumns,tagAlignCols)]

    write.table(res, file=paste0("clusterTables/",metadataFile, ".ClusterMetadataClean.tsv"), sep='\t', quote=F, col.names=T, row.names=F)
    return(res)
}

clusterData <- do.call(rbind, mapply(readClusterMetadataTable, clusterFiles$FileDownloadName, clusterFiles$DatasetID, SIMPLIFY=FALSE))
rownames(clusterData) <- NULL
stopifnot(all(!duplicated(clusterData$CellClusterID)))
subset(clusterData, ATACtagAlignSorted == "")

write.table(clusterData, file="Y2AVE_SingleCellDatasets.CellClusterTableRaw.tsv", sep='\t', quote=F, col.names=T, row.names=F)

## Edit script to output metadata tables per experiment, and read them back in — so that we could possibly manually edit and fix things and
##  don't need to go back through the source

#####################################
## Print out list of missing files
expectedFiles <- list(
    `Parse SplitSeq`=c("GeneExpressionCPM","GeneExpressionFilteredH5"),
    `10X snATAC`=c("ATACtagAlignSorted","ATACtagAlignSortedIndex","PeaksMACS2"),
    `10x snATAC`=c("ATACtagAlignSorted","ATACtagAlignSortedIndex","PeaksMACS2"),
    `sciATAC`=c("ATACtagAlignSorted","ATACtagAlignSortedIndex","PeaksMACS2"),
    `MULTI-seq`=c("ATACtagAlignSorted","ATACtagAlignSortedIndex","PeaksMACS2","GeneExpressionCPM","GeneExpressionFilteredH5"),
    `10x Multiome`=c("ATACtagAlignSorted","ATACtagAlignSortedIndex","PeaksMACS2","GeneExpressionCPM","GeneExpressionFilteredH5"),
    `SHARE-seq`=c("ATACtagAlignSorted","ATACtagAlignSortedIndex","PeaksMACS2","GeneExpressionCPM","GeneExpressionFilteredH5")
    )
expectedFileTypes <- unique(unlist(expectedFiles))
allFiles <- clusterData %>% gather(FileType, FileName, one_of(expectedFileTypes)) %>%
    select(DatasetID, CellClusterID, FileType, FileName, Assay) %>% 
    unique() %>%
    mutate(ExpectedFile=FALSE) %>% 
    rowwise() %>%
    mutate(ExpectedFile=FileType %in% expectedFiles[[Assay]]) %>%
    ungroup() %>% as.data.frame()
missingFiles <- allFiles %>% filter(ExpectedFile & FileName == "") %>% select(DatasetID, CellClusterID, FileType)
write.table(missingFiles, file="Y2AVE_SingleCellDatasets.MissingFiles.tsv", sep='\t', quote=F, col.names=T, row.names=F)

#####################################
## Write out Manual Annotation Labels into a separate file
if (FALSE)  ## Only run this if we're sure we're not going to overwrite the file we want
    clusterData %>% 
        select(DatasetID, CellClusterID, ManualAnnotationLabel, Biosample) %>%
        mutate(ManualAnnotationLabel=gsub("_", " ", ManualAnnotationLabel)) %>%
        write.table(file="Y2AVE_SingleCellDatasets.CellClusterManualAnnotationLabels.tsv", sep='\t', quote=F, col.names=T, row.names=F)


#####################################
## Now make edits to the ManualAnnotationLabel column
clusterData <- read.delim("Y2AVE_SingleCellDatasets.CellClusterTableRaw.tsv", stringsAsFactors=F)
newLabels <- read.delim(file="Y2AVE_SingleCellDatasets.CellClusterManualAnnotationLabels.tsv", stringsAsFactors=F)
stopifnot(all(clusterData$CellClusterID %in% newLabels$CellClusterID))
clusterDataNewLabels <- clusterData %>% 
    select(-ManualAnnotationLabel) %>%
    merge(newLabels %>% select(CellClusterID,ManualAnnotationLabel)) %>%
    arrange(DatasetID,CellClusterID) %>%
    select(CellClusterID, DatasetID, ManualAnnotationLabel, everything())

stopifnot(all(!duplicated(clusterDataNewLabels$CellClusterID)))
write.table(clusterDataNewLabels, file="Y2AVE_SingleCellDatasets.CellClusterTable.tsv", sep='\t', quote=F, col.names=T, row.names=F)
