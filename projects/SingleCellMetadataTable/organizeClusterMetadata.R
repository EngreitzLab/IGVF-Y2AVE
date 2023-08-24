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

    x$Assay <- subset(datasetTable, DatasetID==datasetID)$Assay
    x$Species <- subset(datasetTable, DatasetID==datasetID)$Species
    x$ATACtagAlignUnsorted <- mapply(findTagAlignFiles, x$CellClusterID, datasetID, "cluster pseudobulk ATAC unsorted tagAlign", x$tagAlignFile)
    x$ATACtagAlignSorted <- mapply(findTagAlignFiles, x$CellClusterID, datasetID, "cluster pseudobulk ATAC sorted tagAlign", x$tagAlignFile)
    x$ATACtagAlignSortedIndex <- mapply(findTagAlignFiles, x$CellClusterID, datasetID, "cluster pseudobulk ATAC sorted tagAlign Tabix index", x$tagAlignFile)
    x$GeneExpressionCPM <- findDatasetFile(datasetID, "cluster gene expression counts per million")
    x$GeneExpressionFilteredH5 <- findDatasetFile(datasetID, "gene expression matrix .h5 combined and filtered")
    x$PeaksMACS2 <- "In progress"
    x$DatasetID <- datasetID

    tagAlignCols <- c("Assay","Species","ATACtagAlignUnsorted","ATACtagAlignSorted","ATACtagAlignSortedIndex","GeneExpressionCPM","GeneExpressionFilteredH5","PeaksMACS2")
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
## Now make edits to the ManualAnnotationLabel column
#clusterDataRelabeled <- clusterData
#new_str <- gsub('[^[:alnum:] ]','',address_str)

stopifnot(all(!duplicated(clusterData$CellClusterID)))
write.table(clusterData, file="Y2AVE_SingleCellDatasets.CellClusterTable.tsv", sep='\t', quote=F, col.names=T, row.names=F)
