#!/bin/bash
#SBATCH --gres=lscratch:10
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=2:00:00

#Create a project folder, then globus transfer data to the project folder such as eyeGENE_ORF15
#run outside the project folder

TIMESTAMP=$(date "+%Y%m%d-%H%M%S")
ProjectDir=$1
mkdir -p $ProjectDir-$TIMESTAMP/old_bam $ProjectDir-$TIMESTAMP/vcf_nisc/gvcf $ProjectDir-$TIMESTAMP/vcf_nisc/vcf

module load samtools

find $ProjectDir -name "*.bam*" -exec mv {} $ProjectDir-$TIMESTAMP/old_bam \;
find $ProjectDir -name "*.g.vcf.gz*" -exec mv {} $ProjectDir-$TIMESTAMP/vcf_nisc/gvcf \;
find $ProjectDir -name "*.vcf*" -exec mv {} $ProjectDir-$TIMESTAMP/vcf_nisc/vcf \;
for vcf in $ProjectDir-$TIMESTAMP/vcf_nisc/vcf/*.vcf; do bgzip -@ 8 $vcf; tabix -@ 8 $vcf.gz; done
mv $ProjectDir/*.* $ProjectDir-$TIMESTAMP/
rm -rf $ProjectDir
mv $ProjectDir-$TIMESTAMP $ProjectDir
cd $ProjectDir
bash ~/git/NGS_genotype_calling/make_bam_metadatafile.sh


#to be continued