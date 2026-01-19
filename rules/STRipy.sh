#!/bin/bash
#SBATCH --gres=lscratch:200
#SBATCH --cpus-per-task=28
#SBATCH --mem=64g
#SBATCH --partition=norm
#SBATCH --time=32:0:0

config=$1
sample=$2

#selecting mapped reads using subset_mapped_reads.sh reduced sample processing time to 44 hours.
#selecting reads using subset_chr1-M_reads.sh ???
#The bam file has to be saved in the directory sample_bam and named as *.markDup.bam, unless changing the code below.
set -e
if (( $(module list 2>&1 | grep "stripy-pipeline" | wc -l) == 0 )); then module load stripy-pipeline/2.5; fi
if (( $(module list 2>&1 | grep "bwa/" | wc -l) == 0 )); then module load bwa/0.7.17; fi
ref=$(grep "^ref_genome:" $config | head -n 1 | cut -d"'" -f 2)
mkdir -p STRipy
time python3 $STRIPY_PIPELINE_HOME/stri.py \
 --genome hg38 \
 --reference $ref \
 --output STRipy \
 --locus AFF2,AR,ARX_1,ARX_2,ATN1,ATXN1,ATXN10,ATXN2,ATXN3,ATXN7,ATXN8OS,BEAN1,C9ORF72,CACNA1A,CBL,CNBP,COMP,DAB1,DIP2B,DMD,DMPK,FGF14,FMR1,FOXL2,FXN,GIPC1,GLS,HOXA13_1,HOXA13_2,HOXA13_3,HOXD13,HTT,JPH3,LRP12,MARCHF6,NIPA1,NOP56,NOTCH2NLC,NUTM2B-AS1,PABPN1,PHOX2B,PPP2R2B,PRDM12,RAPGEF2,RFC1,RILPL1,RUNX2,SAMD12,SOX3,STARD7,TBP,TBX1,TCF4,TNRC6A,XYLT1,YEATS2,ZIC2,ZIC3 \
 --config /data/OGL/resources/STRipy/config.json \
 --input sample_bam/$sample.markDup.bam

#Changed threads to $SLURM_CPUS_PER_TASK in the config file, but the program still used 56 CPUs.