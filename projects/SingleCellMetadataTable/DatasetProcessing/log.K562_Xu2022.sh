conda activate bedops
synapse get -r syn52333627
echo -e "bc\tCellClusterID" > K562_10XMultiome_Xu2022_CellClusterAssignment.tsv
zcat barcode_cell_types.tsv.gz | sed 1d | awk '{ print $1 "\tK562_10XMultiome_Xu2022_CellCluster0" }' >> K562_10XMultiome_Xu2022_CellClusterAssignment.tsv
gzip K562_10XMultiome_Xu2022_CellClusterAssignment.tsv
synapse store --parentid syn52333627 K562_10XMultiome_Xu2022_CellClusterAssignment.tsv.gz
## Actually no, this is not what is needed. Emailed Tulika instead



conda activate EngreitzLab
SIZES=GRCh38_EBV.chrom.sizes.tsv

## This didn't work because required too much memory
## bedtools sort -i K562_10XMultiome_Xu2022_cluster1.atac.filter.hg38.tagAlign.tsv.gz -g $SIZES | bgzip -c > K562_10XMultiome_Xu2022_CellCluster0.atac.filter.cutsites.hg38.tagAlign.gz

## Instead gotta split by chromosome first
mkdir split
zcat K562_10XMultiome_Xu2022_cluster1.atac.filter.hg38.tagAlign.tsv.gz | awk '{ print $0 > "split/" $1 ".bed" }' &

cat ../$SIZES | 
while read chr size; do 
    bedtools sort -i split/${chr}.bed -g ../$SIZES > split/${chr}.sorted.bed
done

cat ../$SIZES | 
while read chr size; do 
    cat split/${chr}.sorted.bed
done | bgzip -c > K562_10XMultiome_Xu2022_CellCluster0.atac.filter.cutsites.hg38.tagAlign.gz

tabix -p bed K562_10XMultiome_Xu2022_CellCluster0.atac.filter.cutsites.hg38.tagAlign.gz

synapse store --parentid syn52333627 K562_10XMultiome_Xu2022_CellCluster0.atac.filter.cutsites.hg38.tagAlign.gz
synapse store --parentid syn52333627 K562_10XMultiome_Xu2022_CellCluster0.atac.filter.cutsites.hg38.tagAlign.gz.tbi

## Manually edit K562_10XMultiome_Xu2022_CellClusterAssignment.tsv to fix tagAlignFile name
## Then:
synapse store --id syn52333699 K562_10XMultiome_Xu2022_ClusterMetadata.tsv



