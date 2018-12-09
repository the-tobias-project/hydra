wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

module load plink/1.90
module load tabix


#merge the vcfs 

tabix ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
zcat ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep '^#' > merged.vcf
zcat ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '^#' >> merged.vcf

rm ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
rm ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

cat merged.vcf.gz | head -400 | grep "#" | grep "NA" | sed 's/\t/\n/g' | awk -v OFS='\t' '{print $1, $1}' | tail -n +11 > individuals.txt
gzip merged.vcf
rm merged.vcf



# convert to plink format 
plink --vcf merged.vcf.gz -thin-count 100000 --make-bed --out subsampled

rm merged.vcf.gz
split -l 900 individuals.txt

wget https://media.nature.com/original/nature-assets/nature/journal/v526/n7571/extref/nature15394-s2.zip

unzip nature15394-s2.zip
in2csv nature15394-s2/TABLE_1-Samples_Full_Info.xlsx --sheet B | tail -n +2 > merged.dem

plink --bfile subsampled --keep xaa --make-bed --out dset1 
plink --bfile subsampled --keep xab --make-bed --out dset2 
plink --bfile subsampled --keep xac --make-bed --out dset3 

split -l 900 merged.dem
mv xaa dset1.ind
mv xab dset2.ind
mv xac dset3.ind

rm merged.dem





