#!/bin/bash

#----------------------------------------------------------
#Script for running fastp on raw Brook Trout fastq files
#by stream. Naming convention will be 003 as trimmomatic is the third step
#----------------------------------------------------------


#SBATCH --time=00-02:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32

module load StdEnv/2020
module load fastp/0.23.1

#cut_right is the sliding window

BASEDIR=/home/cnems/projects/def-ruzza/cnems/RossCreek_Upstream/RossCreek_trimmed #Cut the poly-g from the adapter trimmed files which are here in trimmed folder
FILE=/home/cnems/projects/def-ruzza/cnems/RossCreek_Upstream/RossCreek_1samplename.txt #file with all samples names including lane number and without the R1_ etc.
POLYTRIM=/home/cnems/projects/def-ruzza/cnems/RossCreek_Upstream/RossCreek_polytrimmed #where I want the trimmed files to go


LINES=$(cat $FILE)

for NAME in $LINES

do
  fastp --threads 16 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 --trim_poly_g -L -A -i $BASEDIR/${NAME}_R1_paired.fastq.gz -I $BASEDIR/${NAME}_R2_paired.fastq.gz -o $POLYTRIM/${NAME}_R1_paired_qual_filtered.fastq.gz -O $POLYTRIM/${NAME}_R2_paired_qual_filtered.fastq.gz -h $POLYTRIM/${NAME}_adapter_clipped_fastp.html
done  |& tee -a $POLYTRIM/RossCreek_polytrim_1sampletest.log
#-L is the length filter and if put -L it disables the length filter which I think we want to do because we are using a sliding window instead
#-A is adapters. Putting -A disables the adapter trimming (already did this)
#Parameters used here are the same as those from the batch effect paper from Therkildsen and Lou
#Two input and output files so need to specify one with -i and -o and the other with -I and -O
#16 is the max number of threads you can ask for in fastp
