conda activate synapsetools

## Curated Y2AVE_SingleCellDatasets.metadata.txt from Google Spreadsheet: https://docs.google.com/spreadsheets/d/1QWa1JUzs7pR02P8uS95MWIqW8-Ldh7BX6D990YXXQbQ/edit#gid=0

## Get list of all files in the Y2AVE/FilteredData Synapse folder
synapse list -r syn52167911 | tr -s ' ' '\t' > FilteredDataContents.tsv

## Get list of aligned data files from the Processed Data File View Table
synapse get -q "SELECT * FROM syn52132754"

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

