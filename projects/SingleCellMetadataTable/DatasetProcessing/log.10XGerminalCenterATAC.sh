## Tabix sort the files for 10X Germinal Center ATAC

## read in cellcluster metadata table

wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv
SIZES=GRCh38_EBV.chrom.sizes.tsv

mkdir GerminalCenter
{
	cd GerminalCenter
	conda activate bedops
	synapse get -r syn52295731
	cd ..
}

conda activate EngreitzLab
sed 1d clusterTables/GerminalCenter10xATAC_combined_ClusterMetadata.tsv |
while read tagAlignFile CellClusterID ManualAnnotationLabel rest; do
	echo $tagAlignFile $CellClusterID
	bedtools sort -i GerminalCenter/$tagAlignFile -g $SIZES | bgzip -c > GerminalCenter/${CellClusterID}.atac.filter.cutsites.hg38.tagAlign.gz && tabix -p bed GerminalCenter/${CellClusterID}.atac.filter.cutsites.hg38.tagAlign.gz
done

## Upload to Synapse
conda activate bedops
sed 1d clusterTables/GerminalCenter10xATAC_combined_ClusterMetadata.tsv |
while read tagAlignFile CellClusterID ManualAnnotationLabel rest; do
	synapse store --parentid syn52295731 GerminalCenter/${CellClusterID}.atac.filter.cutsites.hg38.tagAlign.gz
	synapse store --parentid syn52295731 GerminalCenter/${CellClusterID}.atac.filter.cutsites.hg38.tagAlign.gz.tbi
done

## Still to do 8/21/23: Fix Cluster3 based on tagAlign_GerminalCenter10xATAC_Activated_B.tsv which has a corrupted section
## I Slacked Revathy to ask her to fix this





## Same for GM12878_10xMultiome
mkdir GM12878_10xMultiome
{
	cd GM12878_10xMultiome
	conda activate bedops
	synapse get syn52263983
	cd ..
	conda activate EngreitzLab
	bedtools sort -i GM12878_10xMultiome/tagAlign_GM12878_10XMultiome_GM12878_10XMultiome.tsv -g $SIZES | bgzip -c > GM12878_10xMultiome/GM12878_10XMultiome_Cluster1.atac.filter.cutsites.hg38.tagAlign.gz && tabix -p bed GM12878_10xMultiome/GM12878_10XMultiome_Cluster1.atac.filter.cutsites.hg38.tagAlign.gz
	conda activate bedops
	synapse store --parentid syn52263976 GM12878_10xMultiome/GM12878_10XMultiome_Cluster1.atac.filter.cutsites.hg38.tagAlign.gz
	synapse store --parentid syn52263976 GM12878_10xMultiome/GM12878_10XMultiome_Cluster1.atac.filter.cutsites.hg38.tagAlign.gz.tbi

}
