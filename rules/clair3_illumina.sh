sinteractive -c56 --mem=128g --time=16:0:0 --gres=lscratch:200

module load clair3/20250303 whatshap/2.3 samtools/1.21
clair3 --bam_fn sample_bam/D2097_01_BP524859.markDup.bam \
 --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
 --threads=$SLURM_CPUS_PER_TASK --platform=ilmn --gvcf \
 --bed_fn=/data/OGL/resources/bed/xgenV2_OGLv1.500padded.hg38.bed \
 --model_path=/data/OGL/resources/clair3/ilmn \
 --use_whatshap_for_final_output_phasing \
 --sample_name=D2097_01_BP524859 \
 --output=/lscratch/$SLURM_JOB_ID/D2097_01_BP524859
 
#merge_output.gvcf.gz phase?
#merge_output.vcf.gz
#phased_merge_output.vcf.gz

# Used 10G lscrach, 40G mem, 63 max CPU, 40 min (including phasing and haplotagged bam)

AddPairEndAlleleDepth.py --bam_fn sample_bam/D2097_01_BP524859.markDup.bam \
	--clair3_vcf_input /lscratch/$SLURM_JOB_ID/D2097_01_BP524859/phased_merge_output.vcf.gz \
	--vcf_output /lscratch/$SLURM_JOB_ID/D2097_01_BP524859/phased_merge_output.PEAD.vcf.gz \
	--threads $SLURM_CPUS_PER_TASK

#
#PEAD is time consuming. Started at 10:28am;

 
 #phased bam: /lscratch/$SLURM_JOB_ID/D2097_01_BP524859/tmp/phase_output/phase_bam when adding phased bam option
 
 clair3 --bam_fn sample_bam/D1596_1.markDup.bam \
 --ref_fn=/data/OGL/resources/genomes/NCBI/GRCh38Decoy/genome.fa \
 --threads=$SLURM_CPUS_PER_TASK --platform=ilmn --gvcf \
 --model_path=/data/OGL/resources/clair3/ilmn \
 --use_whatshap_for_final_output_phasing \
 --sample_name=D1596_1 \
 --output=/lscratch/$SLURM_JOB_ID/D1596_1

#Genome took about 2 hours including phasing. used 70G lscrach. ~60-78 CPUS, 50 GB memory

 
 
 #https://github.com/HKU-BAL/Clair3/blob/main/docs/split_haplotype_into_haploid_calling.md
 
 --whatshap=STR            Path of whatshap, whatshap >= 1.0 is required.
 
 --cpus-per-task
 $SLURM_CPUS_PER_TASK
 
 
bcftools view --threads 8 -i 'ID~"vg"' D663test.gt3.anno3.vcf.gz -Oz -o D663test.gt3.anno3.dvg.vcf.gz
tabix -@ 8 -f -p vcf D663test.gt3.anno3.dvg.vcf.gz

module load peddy/0.4.8
peddy -p 4 D663test.gt3.anno3.dvg.vcf.gz D663.ped --prefix D663test.gt3_PEDDY --sites hg38
						