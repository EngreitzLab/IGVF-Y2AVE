conda activate bedops   ## includes synapse client

## Curated Y2AVE_SingleCellDatasets.metadata.txt from Google Spreadsheet: https://docs.google.com/spreadsheets/d/1QWa1JUzs7pR02P8uS95MWIqW8-Ldh7BX6D990YXXQbQ/edit#gid=0
dos2unix Y2AVE_SingleCellDatasets.metadata.txt; dos2unix -c mac Y2AVE_SingleCellDatasets.metadata.txt

## Get list of all files in the Y2AVE/FilteredData Synapse folder
synapse list -r syn52167911 | tr -s ' ' '\t' > FilteredDataContents.tsv

## Get list of aligned data files from the Processed Data File View Table
synapse query "SELECT * FROM syn52132754" > Y2AVE_SingleCellDatasets.allFileTableDump.tsv

## Make list of unique datasets
{
    echo -e 'synID\tDatasetID' > Y2AVE_SingleCellDatasets.tsv
    synapse list syn52167911 | grep -v '   ' | grep -v 'ipynb_checkpoints' | grep -v 'try_folder' | sed 's:/::' | tr -s ' ' '\t' >> Y2AVE_SingleCellDatasets.tsv
}

## Make list of files per dataset
{
    echo -e 'DatasetID\tDatasetSynID\tFileSynID\tFile' > Y2AVE_SingleCellDatasets.files.tsv
    sed 1d Y2AVE_SingleCellDatasets.tsv | \
    while read synID DatasetID; do
	synapse list $synID | tr -s ' ' '\t' | grep -v '/' | awk -v ID=$DatasetID -v syn=$synID '{ print ID "\t" syn "\t" $0 }'
    done >> Y2AVE_SingleCellDatasets.files.tsv
}

###############
## Read various tables + metadata tables into R for further processing + column annotation
## Makes final dataset table: Y2AVE_SingleCellDatasets.datasetList.tsv
## Makes final file table: Y2AVE_SingleCellDatasets.fileList.tsv
conda activate EngreitzLab
Rscript organizeMetadata.R

###############
## Get all of the ClusterMetadata tables, merge, and edit a bit for clarity
conda activate synapsetools
mkdir clusterTables
for synID in $(cat Y2AVE_SingleCellDatasets.fileList.tsv | grep "cluster metadata table" | cut -f 4); do
  synapse get $synID --downloadLocation clusterTables/
done

conda activate EngreitzLab
Rscript organizeClusterMetadata.R

##############
## Upload final files to Synapse
conda activate synapsetools
synapse store --id syn52252340 --name Y2AVE_SingleCellDatasets.fileList.tsv --description "List of processed single-cell data files for Y2AVE datasets. One row per file." Y2AVE_SingleCellDatasets.fileList.tsv
synapse store --id syn52252344 --name Y2AVE_SingleCellDatasets.datasetList.tsv --description "List of processed single-cell Y2AVE datasets and associated metadata (including contacts, data provenance, etc.). One row per dataset. Merge with Y2AVE_SingleCellDatasets.fileList.tsv as needed to find files corresponding to each dataset." Y2AVE_SingleCellDatasets.datasetList.tsv
synapse store --id syn52252345 --name Y2AVE_SingleCellDatasets.CellClusterTable.tsv --description "List of cell clusters across all single-cell Y2AVE datasets, one row per cell cluster.\n\nCurrently this includes only datasets with ATAC-seq, but eventually will be edited to include datasets with only scRNA once those data are available. Merge this table with Y2AVE_SingleCellDatasets.datasetList.tsv to get metadata about the datasets corresponding to each file, and with fileList.tsv to get metdata corresponding to the tagAlign files" Y2AVE_SingleCellDatasets.CellClusterTable.tsv
