#!/bin/bash
#SBATCH -c8
#SBATCH --mem=64g
#SBATCH --gres=lscratch:50
#SBATCH --time=2:0:0

#Could submit with dependency then start automatically after the previous step completes. 
#sbatch --dependency=afterok:11254323 ~/git/NGS_genotype_calling/dev/bcm_RPGR.porechop.minimap2.alignment.sh OGL
#~20 min per sample, thus adjust the time accordingly. 
#fastq files have to be SampleID-LW and SampleID-MW. If missing either LW or MW, touch an empty file.

set -e
library=$1 #QBplEx or OGL; will be part of @RG
batchname=$(date '+%Y%m%d_%H%M') 
WORK_DIR=/lscratch/$SLURM_JOB_ID

module load porechop minimap2/2.26 sambamba/0.8.2 mosdepth

mkdir -p bam porechop mosdepth
#code below not used
#module load chopper/0.7.0
case "${library^^}" in
	"QBPLEX")
	declare -i maxl=17000
	declare -i minl=200
	;;
	*)
	declare -i maxl=17000
	declare -i minl=12000
	;;
esac

export PORECHOP_ADAPTERS=/data/OGL/resources/porechop/adaptersBCM_E6R.py

if [ -e manifest.csv ];
then
	echo "manifest.csv provided, old bam file used for calling"
	echo "delete manifest.csv to re-align file"
else
	for fastq in fastq/*.fastq.gz; do
		sample=$(basename $fastq | sed "s/_reads.fastq.gz//" | sed "s/.fastq.gz//" | sed "s/-/_/");
		echo -e "$sample" >> manifest.csv
		porechop --discard_unassigned --extra_end_trim 0 -i $fastq -b porechop/$sample
	done
fi

export PORECHOP_ADAPTERS=/data/OGL/resources/porechop/adaptersBCM.py
while read -r sample;
do
porechop --discard_unassigned --extra_end_trim 0 -i porechop/$sample/BCE6R.fastq.gz -b porechop/$sample/bcm
porechop --discard_unassigned --extra_end_trim 0 -i porechop/$sample/BCRPGRr.fastq.gz -b porechop/$sample/RPGR
minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:$library\tID:$sample\tSM:$sample\tPL:ONT" /data/OGL/resources/genomes/GRCh38/OPN1LW/GRCh38Decoy_OPN1LW.mmi porechop/$sample/RPGR/BCRPGRf.fastq.gz \
| sambamba sort -u --compression-level 6 --tmpdir=$WORK_DIR -t 3 -o bam/$sample.RPGR.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin);
minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:$library\tID:$sample\tSM:$sample\tPL:ONT" /data/OGL/resources/genomes/GRCh38/OPN1LW/GRCh38Decoy_OPN1LW.mmi porechop/$sample/bcm/BCLW.fastq.gz \
| sambamba sort -u --compression-level 6 --tmpdir=$WORK_DIR -t 3 -o bam/$sample.LW.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin);
minimap2 -a -x map-ont -Y -t 3 -R "@RG\tLB:$library\tID:$sample\tSM:$sample\tPL:ONT" /data/OGL/resources/genomes/GRCh38/OPN1MW/GRCh38Decoy_OPN1MW.mmi porechop/$sample/bcm/BCMW.fastq.gz \
| sambamba sort -u --compression-level 6 --tmpdir=$WORK_DIR -t 3 -o bam/$sample.MW.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin);
sambamba merge -t 3 -l 6 bam/$sample.bam bam/$sample.RPGR.bam bam/$sample.LW.bam bam/$sample.MW.bam
rm bam/$sample.RPGR.bam* bam/$sample.LW.bam* bam/$sample.MW.bam*
cd mosdepth
mosdepth -t 1 --no-per-base --by /data/OGL/resources/bed/bcmRPGR.bed --use-median --mapq 1 --fast-mode \
$sample.md ../bam/$sample.bam
cd ..
done < manifest.csv

module load R/4.3.0
echo -e "Sample\tRPGR\tOPN1LW\tOPN1MW\tMax_OPN1" > $batchname.coverage.summary.tsv
while read -r sample; do
Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/bcmRPGRcoverageWide.R $sample mosdepth/$sample.md.regions.bed.gz mosdepth/$sample.coverage.wide.tsv
tail -n 1 mosdepth/$sample.coverage.wide.tsv >> $batchname.coverage.summary.tsv
done < manifest.csv

#zcat mosdepth/$sample.md.regions.bed.gz | sed "s/^/$sample\t/" >> mosdepth/coverage.summary.tsv
#test freebayes

module load samtools/1.17 vcflib/1.0.3 vt/0.57721

module load freebayes/1.3.5 deepvariant/1.6.0 whatshap/2.3 annovar/2020-06-08
mkdir -p freebayes deepvariant vcf annotation

#changed -C 5 to 3 to account for low coverage Plasmid Express reads.
while read -r sample;
do
coverage=$(grep $sample $batchname.coverage.summary.tsv | cut -f 5)
echo $coverage
if (( $coverage > 3 )); then
freebayes -f /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --max-complex-gap 80 -p 10 -C 3 -F 0.05 \
--genotype-qualities --strict-vcf --use-mapping-quality --min-base-quality 5 \
--targets /data/OGL/resources/bed/OPN1_e2-e5.bed  \
bam/$sample.bam \
| vcffilter -f "QUAL > 10" \
| bcftools +fill-tags - -Ou -- -t VAF \
| bcftools norm --multiallelics -any --output-type u - \
| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  --output-type u --no-version - \
| bcftools annotate --set-id 'haplo_%CHROM\:%POS%REF\>%ALT' -x ^INFO/AO,^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type u --no-version \
| bcftools norm -d exact --output-type z -o freebayes/$sample.haplo.vcf.gz
tabix -p vcf freebayes/$sample.haplo.vcf.gz
run_deepvariant --model_type ONT_R104 --num_shards 1 \
 --ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
 --regions /data/OGL/resources/bed/OPN1_geneBody_LR-PCR.bed \
 --reads bam/$sample.bam \
 --output_vcf deepvariant/$sample.OPN1.dv.vcf.gz \
 --sample_name $sample \
 --intermediate_results_dir /lscratch/$SLURM_JOB_ID
bcftools norm --multiallelics -any --output-type u deepvariant/$sample.OPN1.dv.vcf.gz \
| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa  --output-type u --no-version - \
| bcftools annotate --set-id 'sml_%CHROM\:%POS%REF\>%ALT' -x ^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type u --no-version \
| bcftools filter --include 'QUAL>0' --output-type u \
| bcftools norm -d exact --output-type z -o deepvariant/$sample.OPN1.dv.filtered.vcf.gz
tabix -p vcf deepvariant/$sample.OPN1.dv.filtered.vcf.gz
whatshap phase --reference /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
 deepvariant/$sample.OPN1.dv.filtered.vcf.gz bam/$sample.bam \
 | bgzip -f > deepvariant/$sample.OPN1.dv.phased.vcf.gz
tabix -p vcf deepvariant/$sample.OPN1.dv.phased.vcf.gz
rm deepvariant/$sample.OPN1.dv.filtered.vcf.gz*
bcftools concat -a --rm-dups none --no-version \
deepvariant/$sample.OPN1.dv.phased.vcf.gz freebayes/$sample.haplo.vcf.gz \
-Oz -o vcf/$sample.OPN1.vcf.gz
tabix -p vcf vcf/$sample.OPN1.vcf.gz
convert2annovar.pl -format vcf4old vcf/$sample.OPN1.vcf.gz -includeinfo --outfile annotation/$sample.avinput
ver=hg38
table_annovar.pl annotation/$sample.avinput \
$ANNOVAR_DATA/$ver \
-buildver $ver \
-remove \
-out annotation/$sample.avinput \
--protocol refGeneWithVer \
-operation g \
--argument '-hgvs' \
--polish -nastring . \
--thread 1 \
--otherinfo
sed -i "1 s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT_FIELDS/" annotation/$sample.avinput."$ver"_multianno.txt
Rscript ~/git/NGS_genotype_calling/NGS_generic_OGL/bcmlocus.R \
/data/OGL/resources/bcmlocus.xlsx \
$sample annotation/$sample.avinput."$ver"_multianno.txt annotation/$sample.bcm.tsv annotation/$sample.bcm.xlsx
rm annotation/$sample.avinput*
else echo -e "####################\n$sample low coverage @ $coverage\n####################";
fi
done < manifest.csv

module load clair3/1.0.10

mkdir -p deepvariant clair3 vcf annotation

while read -r sample; do
coverage=$(grep $sample $batchname.coverage.summary.tsv | cut -f 2)
echo $coverage
if (( $coverage > 3 )); then
run_deepvariant --model_type ONT_R104 --num_shards 1 \
--ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
--regions /data/OGL/resources/bed/RPGR.bed \
--reads bam/$sample.bam \
--output_vcf deepvariant/$sample.RPGR.dv.vcf.gz \
--sample_name $sample \
--intermediate_results_dir $WORK_DIR/dv/$sample
vcffilter -f "QUAL > 0" deepvariant/$sample.RPGR.dv.vcf.gz \
| bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version - \
| bcftools annotate --set-id 'dv_%CHROM\:%POS%REF\>%ALT' -x FORMAT/PL --output-type z -o deepvariant/$sample.RPGR.dv.filtered.vcf.gz
tabix -p vcf deepvariant/$sample.RPGR.dv.filtered.vcf.gz
clair3 --bam_fn=bam/$sample.bam \
--ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
--threads=8 \
--platform=ont \
--bed_fn=/data/OGL/resources/bed/RPGR.bed \
--model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 \
--sample_name=$sample \
--ref_pct_full=1.0 --var_pct_full=0.5 \
--chunk_size=-1 \
--output=$WORK_DIR/clair3/$sample 
###--ctg_name=chrX \ ## --bed_fn=/data/OGL/resources/bed/RPGR_ORF15.bed 
###--var_pct_full=1.0 \ --ref_pct_full=1.0 \
##Will false negative increase if setting minimum read = 1 instead of 2 as in the default??? ##added var_pct_full & ref_pcf_full 8/6/24
cp $WORK_DIR/clair3/$sample/merge_output.vcf.gz clair3/$sample.vcf.gz
cp $WORK_DIR/clair3/$sample/merge_output.vcf.gz.tbi clair3/$sample.vcf.gz.tbi
bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u --no-version clair3/$sample.vcf.gz \
| bcftools annotate --set-id 'clr_%CHROM\:%POS%REF\>%ALT' --output-type z -o clair3/$sample.norm.vcf.gz
tabix -p vcf clair3/$sample.norm.vcf.gz
bcftools concat --threads 1 -a --no-version \
clair3/$sample.norm.vcf.gz deepvariant/$sample.RPGR.dv.filtered.vcf.gz -Oz \
-o vcf/$sample.clair.dv.vcf.gz
tabix -p vcf vcf/$sample.clair.dv.vcf.gz
rm -rf $WORK_DIR/*
convert2annovar.pl -format vcf4old vcf/$sample.clair.dv.vcf.gz -includeinfo --outfile annotation/$sample.avinput
ver=hg38
table_annovar.pl annotation/$sample.avinput \
$ANNOVAR_DATA/$ver \
-buildver $ver \
-remove \
-out annotation/$sample.avinput \
--protocol refGeneWithVer \
-operation g \
--argument '-hgvs' \
--polish -nastring . \
--thread 1 \
--otherinfo
sed -i "1 s/Otherinfo1\tOtherinfo2\tOtherinfo3\tOtherinfo4\tOtherinfo5\tOtherinfo6\tOtherinfo7\tOtherinfo8\tOtherinfo9\tOtherinfo10/CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGT_FIELDS/" annotation/$sample.avinput."$ver"_multianno.txt
Rscript ~/git/NGS_genotype_calling/dev/annovar2excel.R annotation/$sample.avinput."$ver"_multianno.txt annotation/$sample.RPGR.clair.dv.anno.xlsx
rm annotation/$sample.avinput*
else echo -e "####################\n$sample low coverage @ $coverage\n####################";
fi
done < manifest.csv

#clean up directory
#rm -r fastq pod5_pass porechop
rm -rf fastq porechop

##clair3 didn't call some SNP, even with the last line of --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1
# mkdir -p clair3/$sample
# clair3 \
  # --bam_fn=bam/$sample.bam \
  # --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
  # --threads=8 \
  # --platform=ont \
  # --bed_fn=/data/OGL/resources/bed/OPN1_geneBody_LR-PCR.bed \
  # --model_path=/data/OGL/resources/clair3/ont_r1041/r1041_e82_400bps_sup_v420 \
  # --sample_name=$sample \
  # --output=clair3/$sample \
  # --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1
  
##freebayes is not as good as deepvariant
# freebayes -f /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --limit-coverage 1000 --min-alternate-fraction 0.05 -C 3 -p 10 \
# --min-mapping-quality 0 --genotype-qualities --strict-vcf --use-mapping-quality --min-base-quality 5 \
# --targets /data/OGL/resources/bed/OPN1_geneBody_LR-PCR.bed \
# bam/$sample.bam \
# | vcffilter -f "QUAL > 10 & SAF > 0 & SAR > 0 & RPR > 0 & RPL > 0 & AO > 2 & DP > 5" \
# | bcftools +fill-tags - -Ou -- -t VAF \
# | bcftools norm --multiallelics -any --output-type u - \
# | vt decompose_blocksub -p -m -d 2 - \
# | bcftools norm --check-ref s --fasta-ref /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa --output-type u - \
# | bcftools annotate --set-id 'sm_%CHROM\:%POS%REF\>%ALT' -x ^INFO/AO,^FORMAT/GT,FORMAT/DP,FORMAT/VAF --output-type v --no-version \
# | bcftools norm -d exact --output-type z -o freebayes/$sample.small.vcf.gz
# tabix -p vcf freebayes/$sample.small.vcf.gz

#perhaps it's better to use the freebayes raw output to represent the variants, as bcftools norm may shorten the variants.

#rm annotation/$sample.avinput."$ver"_multianno.txt


#find LW/fastq -name "*fastq*" | parallel -j 10 -I% --max-args 1 --tag "minimap2 -a -Y -t 3 /data/OGL/resources/genomes/GRCh38/OPN1LW/GRCh38Decoy_OPN1LW.mmi % | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o %.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)"

#find MW/fastq -name "*fastq*" | parallel -j 10 -I% --max-args 1 --tag "minimap2 -a -Y -t 3 /data/OGL/resources/genomes/GRCh38/OPN1MW/GRCh38Decoy_OPN1MW.mmi % | sambamba sort -u --compression-level 6 --tmpdir=/lscratch/$SLURM_JOB_ID -t 3 -o MW/bam/%.bam <(sambamba view -S -f bam --compression-level 0 -t 3 /dev/stdin)"

#The minimap2 index were generated by: minimap2 -d genome.mmi  /fdb/igenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa 
#The index is now moved to /data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.mmi


