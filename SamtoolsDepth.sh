#!/bin/bash

#-----------------------------------------------
#Script to calculate depth for each individual
#from the quality filtered .bam files
#0010-depth
#-----------------------------------------------

#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=16


module load StdEnv/2020 gcc/9.3.0 samtools/1.17

BASEDIR=/home/cnems/projects/rrg-ruzza/cnems/RossCreek/RossCreek_dedups #Where all my quality filtered, overlapped clipped/dedupped/realigned .bam files are (the realigned .bams are also in here)
SAMPLELIST=/home/cnems/projects/rrg-ruzza/cnems/RossCreek/RossCreek_samplemerge.txt #File with individual sample names like 21_RCU_SFO_05_S35 etc and only one sample name for each samples. These are the prefixes to the .bam files I will used here.
REFNAME=SfoGenome #Name of the reference genome being used to add to file names.

for name in `cat $SAMPLELIST`; do
  samtools depth -aa $BASEDIR/${name}'_bt2_'$REFNAME'_minq20_sorted_dedup_overlapclipped_realigned.bam' -q 20 -Q 20 | cut -f 3 | gzip > $BASEDIR/${name}'_bt2_'$REFNAME'_minq20_sorted_dedup_overlapclipped_realigned.bam.depth.gz'
done
#Use samtools depth on the dedup and overlap clipped bam files for each sample where -aa is calculate depth at all positions including unused reference sequences and I think positions with zero depth then;
#-q is the minimum base quality filter and -Q is the minimum mapping quality filter.
#cut -f 3 means select only the 3rd column of the samtools depth .bam file which is the number of reads aligned to the reference at that position, in other words, the depth of coverage for each base indicated.
#zip it using gzip so each sample will have a zipped filed of depth information and have the ending .bam.depth.gz
