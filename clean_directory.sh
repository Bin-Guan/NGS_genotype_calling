#!/bin/bash
#SBATCH --gres=lscratch:10
#SBATCH --cpus-per-task=8
#SBATCH --mem=8g
#SBATCH --time=2:00:00

calling_configure_file=$1
variantPrioritization_configure_file=$2

set -e

TIMESTAMP=$(date "+%Y%m%d-%H%M%S")

rm -rf coverage/mean.coverage.done
rm -rf bcmlocus/combine.bcmlocus.done
rm -rf orf15/combine.orf15.done
rm -rf mutserve/haplocheck.done
rm -rf freebayes/freebayes.merge.done
rm -rf prioritization/dv_fb.merge.done
rm -rf .snakemake
rm -rf 00log
rm -rf slurm*.out
rm -rf sample_bam
rm -rf deepvariant/deepvariant_phased_glnexusVcf.merge.done
rm -rf deepvariant/deepvariant.gvcf.merge.done
rm -rf deepvariant/deepvariantVcf.merge.done
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
analysis_batch_name=$(grep "^analysis_batch_name:" $calling_configure_file | head -n 1 | cut -d"'" -f 2)
ngstype=$(grep "^datatype:" $variantPrioritization_configure_file | head -n 1 | cut -d"'" -f 2)
ped_file=$(grep "^ped:" $variantPrioritization_configure_file | head -n 1 | cut -d"'" -f 2)

export ngstype

#copy sample names to OGL.ped adding datatype and analysis_batch_name
awk -F"\t" -v batch="$analysis_batch_name" -v dt="$ngstype" 'BEGIN{OFS="\t"} NR==1 {print "NGS_type", "batch", $0} NR>1 {print dt, batch, $0}' \
 prioritization/$ped_file | cat /data/OGL/sample_info/master.v3.ped - > /lscratch/$SLURM_JOB_ID/temp.ped
awk '!a[$0]++' /lscratch/$SLURM_JOB_ID/temp.ped > /data/OGL/sample_info/master.v3.ped

#subsequent batch
# awk -F"\t" -v batch="$analysis_batch_name" -v dt="$ngstype" 'BEGIN{OFS="\t"} NR>1 {print $0, batch, dt}' prioritization/$ped_file >> /data/OGL/resources/OGLsample/OGL.ped

chgrp OGL /data/OGL/sample_info/master.ped #what happens for another user?

#copy relevant files to resources

find prioritization/gemini_tsv_filtered/ -name "*.tsv" | parallel -j 8 'gzip -c {} \
 > /data/OGL/resources/GeneSearch/$ngstype/gemini_tsv_filtered/v3_temp/$(echo {/} | cut -d. -f 1).tsv.gz && chgrp OGL /data/OGL/resources/GeneSearch/$ngstype/gemini_tsv_filtered/v3_temp/$(echo {/} | cut -d. -f 1).tsv.gz' 

find prioritization/gemini_tsv_filtered/ -name "*.tsv.gz" | parallel -j 8 'cp {} /data/OGL/resources/GeneSearch/$ngstype/gemini_tsv_filtered/v3_temp/$(echo {/} | cut -d. -f 1).tsv.gz && chgrp OGL /data/OGL/resources/GeneSearch/$ngstype/gemini_tsv_filtered/v3_temp/$(echo {/} | cut -d. -f 1).tsv.gz' 

#chgrp --recursive OGL /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/temp$TIMESTAMP
#mv /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/temp$TIMESTAMP/* /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/v3
#rm -r /data/OGL/resources/GeneSearch/genome/gemini_tsv_filtered/temp$TIMESTAMP

find clinSV/ -name "*.clinsv.SV-CNV.PASS.vcf.gz*" -exec cp {} /data/OGL/resources/clinSV/genome \;
find clinSV/ -name "*.clinsv.SV-CNV.RARE_PASS_GENE.vcf.gz*" -exec cp {} /data/OGL/resources/clinSV/genome \;
find clinSV/ -name "*.clinSV.PASS.annotated.tsv.gz" -exec cp {} /data/OGL/resources/clinSV/genome \;
find clinSV/ -name "*.clinSV.RARE_PASS_GENE.annotated.tsv.gz" -exec cp {} /data/OGL/resources/clinSV/genome \;
find clinSV/ -name "*.clinSV.RARE_PASS_GENE.eG.tsv" -exec cp {} /data/OGL/resources/clinSV/genome \;

cp -a deepvariant/gvcf/* /data/OGL/resources/OGLsample/genome_dv_gvcf
cp -a deepvariant/vcf/*.dv.filtered.vcf.gz* /data/OGL/resources/OGLsample/genome_dv_vcf
cp 


# ( cat $CONTIGFILE | parallel -C "\t" -j 21 --tmpdir $WORK_DIR --eta --halt 2 --line-buffer \
		 	# --tag "whatshap phase --reference $WORK_DIR/$(basename {config[ref_genome]}) \
			# --indels --ignore-read-groups $WORK_DIR/filtered/{{2}}.filtered.vcf.gz $WORK_DIR/{wildcards.sample}.markDup.bam \
			# | bgzip -f --threads $(({threads}-10)) > $WORK_DIR/phased/{{2}}.phased.vcf.gz" \
		# ) && echo "whatshap on chr completed" || exit 7

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
