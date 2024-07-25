#!/bin/bash

#-----------------------------------------------
#Script to get rid of duplicates and clip
#overlapping reads from the merged bam files
#for each sample
#-----------------------------------------------

#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=20

module load StdEnv/2020
module load nixpkgs/16.09 intel/2018.3
module load bamutil/1.0.14
module load picard


BASEDIR=/home/cnems/projects/def-ruzza/cnems/RossCreek_Upstream #main directory to where my merged files are.
BAMLIST=/home/cnems/projects/def-ruzza/cnems/RossCreek_Upstream/RossCreek_samplemerge.txt #File with individual sample names like 21_RCU_SFO_01 etc and only one same name for each samples
DEDUP=/home/cnems/scratch/RossCreek/RossCreek_dedups #where I want the deduplicated and overlapping clipped merged bams to go in the scratch folder.
REFNAME=SnaGenome1_1 #Name of the reference genome being used.

LINES=$(cat $BAMLIST)

for name in $LINES
do
  #Picard will get rid of duplicate reads and print the dupstat report which tells
java -Xmx40g -jar $EBROOTPICARD/picard.jar MarkDuplicates -I BASEDIR'/RossCreek_mapped/RossCreek_mergedbams/'$name'_merged_paired_bt2_'$REFNAME'_minq20_sorted.bam' -O $DEDUP/$name'_bt2_'$REFNAME'_minq20_sorted_dedup.bam' -M $DEDUP/$name'_bt2_'$REFNAME'_minq20_sorted_dupstat.txt' VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

#Then clip overlapping paired end Reads
bam clipOverlap --in $DEDUP/$name'_bt2_'$REFNAME'_minq20_sorted_dedup.bam' --out $DEDUP/$name'_bt2_'$REFNAME'_minq20_sorted_dedup_overlapclipped.bam' --stats

done
#VALIDATION_STRINGENCY=SILENT this is supposed to improve performance when processing bam files that have variable length data in terms of reads, qualities, and tags.
#M is what I want the stats files to be called for picard.
#--stats for bamutil will print stats on the overlapping reads.
