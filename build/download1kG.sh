#!/usr/bin/env bash

echo "Grabbing 1000 genomes files"
echo "Grabbing chr 22"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --quiet
echo "Grabbing chr 22"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --quiet
echo "Grabbing chr 20"
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --quiet
echo "Finished grabbing 1000 genomes files"

#merge the vcfs
echo "Starting merge"
zcat ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep '^#' > merged.vcf
zcat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '^#' >> merged.vcf

echo "Finished merge, cleaning up..."

rm ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
rm ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
rm ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

echo "Finished cleanup, creating gzip"

cat merged.vcf | head -400 | grep "#" | grep "NA" | sed 's/\t/\n/g' | awk -v OFS='\t' '{print $1, $1}' | tail -n +11 > individuals.txt
gzip merged.vcf

echo "Creating reference list"

zcat merged.vcf.gz | grep -v "#"| cut -f 3,4 | uniq >referencelist.txt

echo "Converting to plink format"

# convert to plink format 
plink --vcf merged.vcf.gz --snps-only --biallelic-only strict --thin-count 150000 --make-bed --out subsampled

echo "splitting"

split -l 900 individuals.txt

echo "Grabbing nature dataset"

wget https://media.nature.com/original/nature-assets/nature/journal/v526/n7571/extref/nature15394-s2.zip --quiet

unzip nature15394-s2.zip

echo "Extracted nature dataset, converting from excel to plink"

in2csv nature15394-s2/TABLE_1-Samples_Full_Info.xlsx --sheet B | tail -n +2 > merged.dem

echo "Running plink on reference list"

plink --bfile subsampled --reference-allele referencelist.txt --keep xaa --make-bed --out dset1
plink --bfile subsampled --reference-allele referencelist.txt --keep xab --make-bed --out dset2
plink --bfile subsampled --reference-allele referencelist.txt --keep xac --make-bed --out dset3

split -l 900 merged.dem
mv xaa dset1.ind
mv xab dset2.ind
mv xac dset3.ind

rm merged.dem

echo "finished"
