## Read in a submitted Y2AVE variant list
## Parse the header + column names, which are in variable formats depending on the submitters
## Output a file with a single line listing the indices of the desired columns (SPDI, chr, rsID) in that order

library(data.table)

args <- commandArgs(trailingOnly=TRUE)

infile <- args[1]
outfile <- args[2]

requiredColumns <- list(
    SPDI=c("spdi","spdi:::"),
    chr=c("chr","chromosome"))


getColumnNames <- function(filepath, desiredColumns) {
    ## Handle three cases:
    ## 1. File contains header with comment characters ('#'), and row with column names starts with '#'
    ## 2. File contains header with comment characters ('#'), and row with column names does not start with '#'
    ## 3. File does not contain a header, and row with column names does not start with '#'
    ## 4. File does not contain any header, in which case throw a flag but then assume default column order
    
    ## Get lines up to the first line that does not start with '#'
    lines <- c()
    if (grepl(".gz$",filepath))
        incon <- gzfile(filepath, "r")
    else
        incon <- file(filepath, "r")
    while (length(myLine <- scan(incon, what="character",nlines=1,sep=',',quiet=TRUE)) > 0 ) {
        myLine <- paste0(myLine, collapse='')
        lines <- c(lines, myLine)
        if (!grepl('^#', myLine)) 
            break
    }

    colLine <- NULL
    if (grepl("SPDI", lines[length(lines)]))
        colLine <- lines[length(lines)]
    else if (length(lines) > 1)
        if (grepl("SPDI", lines[length(lines)-1]))
            colLine <- lines[length(lines)-1]

    if (is.null(colLine)) {
        print("WARNING: Could not find SPDI column in expected location. Using default file names.")
        print(lines[length(lines)-1])
        colLine <- paste0(c("chrRefSeqID","chr","position","ReferenceAllele","AlternativeAllele","SPDI","rsID"), collapse='\t')
    }
    
    columns <- strsplit(tolower(colLine), "\t")[[1]]

    return(columns)
}

columns <- getColumnNames(infile, desiredColumns)

## Get indices for the desired column n
columnIndices <- list()
for (desiredCol in names(requiredColumns)) {
    index <- which(columns %in% requiredColumns[[desiredCol]])
    if (length(index) == 0) {
        stop(paste("ERROR: Could not find Column ", desiredCol, "\nColumns:", paste(columns,collapse=',')))
    } else if (length(index) > 1) {
        stop(paste("ERROR: Found more than one match for Column ", desiredCol,":",paste(columns[index],collapse=',')))
    } else {
        columnIndices[desiredCol] <- index
    }
}

write.table(paste0(unlist(columnIndices), collapse=','), file=outfile, row.names=F, col.names=F, quote=F)
