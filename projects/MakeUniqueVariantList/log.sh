## Jesse Engreitz
## August 1, 2023

#################################################################
## SETUP
git clone git@github.com:mikelove/igvf_spdi_demo.git
conda install -n base conda-libmamba-solver
conda env create -f synapsetools.yml --experimental-solver libmamba  ##needed libmamba speed-up to solve the environment

## Download the Y2AVE Data Table
conda activate synapsetools

################################################################
## Special processing needed for Max Schubach files, which list SPDI files only; need to add the chromosome back in
{
    mkdir syn52271949; synapse get -r syn52271949 --downloadLocation syn52271949
    mkdir syn52271987; synapse get syn52271987 --downloadLocation syn52271987
    outfile=syn52271949/SchubachMPRAVariants.SPDI.tsv
    echo -e '#FileType=Node:Variant
#Contact=Max Schubach (combined by Jesse Engreitz)
#Genome=GRCh38
#Description=Combination of two variant lists from Max Schubach (see syn51227153 and syn30614866) represents sets of variants tested by MPRA (saturation mutagenesis of some regulatory elements, or promoter-proximal variants). Jesse combined lists of SPDIs from Max into proper Y2AVE format
chrRefSeqID\tchr\tposition\tReferenceAllele\tAlternativeAllele\tSPDI' > $outfile

    for file in $(ls -1 syn52271949/*.SPDI.txt syn52271987/*.SPDI.txt); 
    do
cat $file | 
awk '{ 
split($1,x,":");
split(x[1],refseq,".");
split(refseq[1],chrString,"0000");
if (chrString[2] == "23") {
  chromosome = "X"
} else {
  chromosome = sprintf("%d", chrString[2])
}
print x[1] "\t" "chr" chromosome "\t" x[2] "\t" x[3] "\t" x[4] "\t" $1}' >> $outfile
    done
}



#synapse get -q "SELECT * FROM syn51514535" 
## Download table from Synapse manually
# files.syn51514535.tsv

## Download all of the Variants Files
#for synID in $(grep 'Node:Variants' SYNAPSE_TABLE_QUERY_127613705.csv | cut -f 10 -d '"'); do
#  mkdir $synID
#  synapse get $synID --downloadLocation $synID/
#done

FILES=files.syn51514535.tsv

for synID in $(grep 'Node:Variants' $FILES | cut -f 5); do
  mkdir $synID
  destfile=$(synapse get $synID --downloadLocation $synID/ | grep Creating | cut -f 2 -d ' ')
  echo -e "${synID}\t${destfile}" 
done > files.syn51514535.local.tsv

## Now do the combining
conda env create -f bedops.yml --experimental-solver libmamba
conda activate bedops

## Find out where the SPDI and chromosome columns are located
cat files.syn51514535.local.tsv |
while read synID localfile; do
  echo "Starting " $synID " " $localfile
  Rscript getVariantColumnIndices.R $localfile ${localfile}.headers
done

## Extract the SPDI and chromosome columns
cat files.syn51514535.local.tsv |
while read synID localfile; do
  echo "Starting " $synID " " $localfile
  csvtk cut -t -f $(cat ${localfile}.headers) $localfile | 
    ## Make sure that SPDI has been computed 
    grep "^NC_" | 
    ## Change "chr1" to "1"
    sed 's/chr//' | 
    ## Fix rare cases where SPDI position is output in scientific notation
    awk '{ if (index($1,"e+0") != 0) { split($1,a,":"); a[2]=sprintf("%d",a[2]);  $1 = a[1] ":" a[2] ":" a[3] ":" a[4] } print $0 }' > ${localfile}.spdi.tsv
  cat ${localfile} | grep -v "NC" | grep "Error_Ref_Mismatch\|NA" > ${localfile}.failedSpdi.tsv
  #csvtk cut -t -f $(cat ${localfile}.headers) $localfile | grep -n -v "^NC_" > ${localfile}.failedSpdi.tsv
done

cat files.syn51514535.local.tsv | while read synID localfile; do wc -l ${localfile}.failedSpdi.tsv; done > failedSpdiCount.tsv

cat files.syn51514535.local.tsv |
while read synID localfile; do
    cat ${localfile}.spdi.tsv
done > spdiCombined.tsv

{
    LC_ALL=C
    sort -u spdiCombined.tsv > spdiUnique.tsv
}

echo -e '#FileType=Node:Variant
#Contact=Jesse Engreitz (engreitz@stanford.edu)
#Genome=GRCh38
#Description=Unique list of variants submitted for Y2AVE based on the submitted SPDIs. Note that variants that were not submitted with a SPDI were dropped.
chrRefSeqID\tchr\tposition\tReferenceAllele\tAlternativeAllele\tSPDI' > Y2AVECombinedVariants.tsv
cat spdiUnique.tsv | 
awk '{ split($1,x,":"); print "chr" $2 "\t" x[2] "\t" x[2] + 1 "\t" x[3] "\t" x[4] "\t" $1 "\t" x[1]}' | 
sort-bed - |
csvtk cut -t -f 7,1,3,4,5,6 >> Y2AVECombinedVariants.tsv



