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


#####################################################################
## Special processing / combination for the TOPMed files from Kushal
## Download from Synapse into Kushal/
{
    cd Kushal/
    outfile=TOPMed_freeze_8.common_LF.spdi.tsv
    echo -e '#FileType=Node:Variant
#Contact=Kushal Dey
#Genome=GRCh38
#Description=Common and low-frequency variants from TOPMed freeze 8
chrRefSeqID\tchr\tposition\tReferenceAllele\tAlternativeAllele\tSPDI' > $outfile
    for file in chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.TOPMed_freeze_8_spdi.common_LF_modified_updated_converted.vcf; do
	sed 1d $file | awk '{ split($2, x, ":"); print x[1] "\t" $3 "\t" $4 "\t" $6 "\t" $7 "\t" $2 }'
    done >> $outfile
    synapse store --parentid syn52273646 TOPMed_freeze_8.common_LF.spdi.tsv.gz
}



#####################################################################
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

## Convert csv to tsv
cat syn52228735/20230802_IGVFYear2VariantsForAnnotationCAVA.csv | sed 's/,/\t/g' > syn52228735/20230802_IGVFYear2VariantsForAnnotationCAVA.tsv

## Find out where the SPDI and chromosome columns are located
conda activate EngreitzLab
cat files.syn51514535.local.tsv |
while read synID localfile; do
  echo "Starting " $synID " " $localfile
  Rscript IGVF-Y2AVE/projects/MakeUniqueVariantList/getVariantColumnIndices.R $localfile ${localfile}.headers
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

conda activate bedops
cat spdiUnique.tsv | 
awk '{ split($1,x,":"); print "chr" $2 "\t" x[2] "\t" x[2] + 1 "\t" x[3] "\t" x[4] "\t" $1 "\t" x[1]}' | 
sort-bed - > tmp.txt

echo -e '#FileType=Node:Variant
#Contact=Jesse Engreitz (engreitz@stanford.edu)
#Genome=GRCh38
#Description=Unique list of variants submitted for Y2AVE based on the submitted SPDIs. Note that variants that were not submitted with a SPDI were dropped.
chr\tposition\tReferenceAllele\tAlternativeAllele\tSPDI' > Y2AVECombinedVariants.tsv
cat tmp.txt | 
awk '{ n=split($0,a,"\t"); if (n == 7) print $0 }' |
csvtk cut -t -f 1,3,4,5,6 >> Y2AVECombinedVariants.tsv


## Make one more version with all submitted variants except 1000G and TOPMed
{
    cat files.syn51514535.local.tsv |
    grep -v 1000_Genomes |
    grep -v TOPMed_freeze_8 |
    while read synID localfile; do
        cat ${localfile}.spdi.tsv
    done > spdiCombined.noTOPMed1000G.tsv

    {
        LC_ALL=C
        sort -u spdiCombined.noTOPMed1000G.tsv > spdiUnique.noTOPMed1000G.tsv
    }

    conda activate bedops
    cat spdiUnique.noTOPMed1000G.tsv | 
    awk '{ split($1,x,":"); print "chr" $2 "\t" x[2] "\t" x[2] + 1 "\t" x[3] "\t" x[4] "\t" $1 "\t" x[1]}' | 
    sort-bed - > tmp.txt

    echo -e '#FileType=Node:Variant
#Contact=Jesse Engreitz (engreitz@stanford.edu)
#Genome=GRCh38
#Description=Unique list of variants submitted for Y2AVE based on the submitted SPDIs, excluding the large 1000G and TOPMed lists. Note that variants that were not submitted with a SPDI were dropped.
chr\tposition\tReferenceAllele\tAlternativeAllele\tSPDI' > Y2AVECombinedVariants.noTOPMed1000G.tsv
    cat tmp.txt | 
    awk '{ n=split($0,a,"\t"); if (n == 7) print $0 }' |
    csvtk cut -t -f 1,3,4,5,6 >> Y2AVECombinedVariants.noTOPMed1000G.tsv
}


gzip Y2AVECombinedVariants.tsv 
gzip Y2AVECombinedVariants.noTOPMed1000G.tsv
##--parentid syn52279180
synapse store --id syn52279183 --description "Unique list of ~29.1 million variants submitted for Y2AVE based on the submitted SPDIs. Note that variants that were not submitted with a SPDI were dropped." Y2AVECombinedVariants.tsv.gz
synapse store --id syn52279184 --description "Unique list of ~1.79 million variants submitted for Y2AVE based on the submitted SPDIs, excluding the large 1000G and TOPMed lists. Note that variants that were not submitted with a SPDI were dropped." Y2AVECombinedVariants.noTOPMed1000G.tsv.gz

## make BED versions of each
csvtk cut -t -f chr,position,position,SPDI Y2AVECombinedVariants.noTOPMed1000G.tsv.gz | 
sed 1d | awk -v OFS='\t' '{ $2 = $2 - 1; print $0 }' |
gzip > Y2AVECombinedVariants.noTOPMed1000G.bed.gz

csvtk cut -t -f chr,position,position,SPDI Y2AVECombinedVariants.tsv.gz | 
sed 1d | awk -v OFS='\t' '{ $2 = $2 - 1; print $0 }' |
gzip > Y2AVECombinedVariants.bed.gz

synapse store --parentid syn52279180 --description "Unique list of ~29.1 million variants submitted for Y2AVE based on the submitted SPDIs. Note that variants that were not submitted with a SPDI were dropped." Y2AVECombinedVariants.bed.gz
synapse store --parentid syn52279180 --description "Unique list of ~1.79 million variants submitted for Y2AVE based on the submitted SPDIs, excluding the large 1000G and TOPMed lists. Note that variants that were not submitted with a SPDI were dropped." Y2AVECombinedVariants.noTOPMed1000G.bed.gz
