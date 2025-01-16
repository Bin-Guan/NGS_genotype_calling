#!/bin/bash
#SBATCH --gres=lscratch:10
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --time=2:00:00

rm -rf coverage/mean.coverage.done.txt
rm -rf bcmlocus/combine.bcmlocus.done.txt
rm -rf orf15/combine.orf15.done.txt
rm -rf mutserve/haplocheck.done.txt
rm -rf freebayes/freebayes.merge.done.txt
rm -rf prioritization/dv_fb.merge.done.txt
rm -rf .snakemake
rm -rf 00log
rm -rf slurm*.out
rm -rf sample_bam old_bam
