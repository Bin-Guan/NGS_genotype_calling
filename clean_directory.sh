#!/bin/bash
#SBATCH --gres=lscratch:10
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=2:00:00

calling_configure_file=$1
variantPrioritization_configure_file=$2

set -e

TIMESTAMP=$(date "+%Y%m%d-%H%M%S")

rm -rf coverage/mean.coverage.done.txt
rm -rf bcmlocus/combine.bcmlocus.done.txt
rm -rf orf15/combine.orf15.done.txt
rm -rf mutserve/haplocheck.done.txt
rm -rf freebayes/freebayes.merge.done.txt
rm -rf prioritization/dv_fb.merge.done.txt
rm -rf .snakemake
rm -rf 00log
rm -rf slurm*.out
rm -rf sample_bam
rm -rf deepvariant/deepvariant_phased_glnexusVcf.merge.done.txt
rm -rf deepvariant/deepvariant.gvcf.merge.done.txt
rm -rf deepvariant/deepvariantVcf.merge.done.txt
rm -rf fastqc/multiqc_report/.snakemake_timestamp
find fastqc -name ".snakemake_timestamp" -exec rm {} \;
rm -rf freebayes/chr_split
rm -rf prioritization/slurm*.out
rm -rf prioritization/.snakemake
rm -rf prioritization/00log
rm -rf prioritization/madeline/madeline.done
rm -rf prioritization/temp
#rm -rf old_bam bam cram
echo "File deletion task done"

module load R/4.3.0
analysis_batch_name=$(grep "^analysis_batch_name:" $1 | head -n 1 | cut -d"'" -f 2)
ngstype=$(grep "^datatype:" $2 | head -n 1 | cut -d"'" -f 2)
ped_file=$(grep "^ped:" $2 | head -n 1 | cut -d"'" -f 2)
#copy sample names to OGL.ped with the last column as analysis_batch_name
#first batch adding column names: awk -F"\t" -v batch="$analysis_batch_name" -v dt="$datatype" 'BEGIN{OFS="\t"} NR==1 {print $0, "batch", "data_type"} NR>1 {print $0, batch, dt}' prioritization/$ped_file >> /data/OGL/resources/OGLsample/OGL.ped
#subsequent batch
awk -F"\t" -v batch="$analysis_batch_name" -v dt="$ngstype" 'BEGIN{OFS="\t"} NR>1 {print $0, batch, dt}' prioritization/$ped_file >> /data/OGL/resources/OGLsample/OGL.ped
chgrp OGL /data/OGL/resources/OGLsample/OGL.ped

#copy relevant files to resources

cp -l -r prioritization/gemini_tsv_filtered /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/temp$TIMESTAMP
chgrp --recursive OGL /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/temp$TIMESTAMP
mv /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/temp$TIMESTAMP/* /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered
rm -r /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/temp$TIMESTAMP

#vcf, SV files to resources

# case "${ngstype^^}" in
	# "PANEL")
		# snakemake -s /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/panel.Snakefile \
		# -pr --local-cores 2 --jobs 1999 \
		# --cluster-config /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/panel.cluster.json \
		# --cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
		# -k --restart-times 1 \
		# --resources res=1 \
		# --configfile $@
		# ;;
	# "AMPLICON"|"AMP")
		# snakemake -s /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/amplicon.Snakefile \
		# -pr --local-cores 2 --jobs 1999 \
		# --cluster-config /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/panel.cluster.json \
		# --cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
		# -k --restart-times 1 \
		# --resources res=1 \
		# --configfile $@
		# ;;
	# "EXOME"|"WES"|"ES")
		# snakemake -s /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/exome.Snakefile \
		# -pr --local-cores 2 --jobs 1999 \
		# --cluster-config /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/exome.cluster.json \
		# --cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
		# -k --restart-times 1 \
		# --resources res=1 \
		# --configfile $@
		# ;;
	# *)
		# snakemake -s /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/Snakefile \
		# -pr --local-cores 2 --jobs 1999 \
		# --cluster-config /home/$USER/git/NGS_genotype_calling/NGS_generic_OGL/cluster.json \
		# --cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
		# -k --restart-times 0 \
		# --resources res=1 \
		# --configfile $@
		# ;;
# esac


#copy index case files to resources and change group ownership to OGL
#make indexSample files

# Rscript ~/git/variant_prioritization/dev/pick_affected_index_sample_from_ped_columnName.R prioritization/$ped_file /data/OGL/resources/OGLsample/OGL.indexSample.$TIMESTAMP.ped $datatype.$analysis_batch_name.indexSample.tsv

# tail -n 2+ /data/OGL/resources/OGLsample/OGL.indexSample.$TIMESTAMP.ped >> /data/OGL/resources/OGLsample/OGL.indexSample.ped

#some scripts are in Z:\resources\OGLsample
#Steps: remove duplicated samples and others "special cases", copy vcf files and SV tsv files to resources, merge samples in another scripts etc. These special cases can be in a file.
#test later.
