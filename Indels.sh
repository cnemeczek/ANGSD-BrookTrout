#!/bin/bash

#-----------------------------------------------
#Script to realign reads from dedup overlapclipped
#bams around indels
#009-indels
#October 21 2023
#-----------------------------------------------

#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c.nemeczek@dal.ca
#SBATCH --account=def-ruzza
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=30

#Run in folder with bams

#Load what is needed for gatk and load gatk
module load nixpkgs/16.09
module load java/1.8.0_121
module load gatk/3.7

BASEDIR=/home/cnems/projects/rrg-ruzza/cnems/RossCreek/RossCreek_dedups #where the bams are
REFERENCE=/home/cnems/projects/rrg-ruzza/cnems/BT_genome/SfoGenome.fna #file name pattern for refenrence files needed.
BAMLIST=/home/cnems/projects/rrg-ruzza/cnems/RossCreek/RossCreek_dedups/RossCreekbam_list_dedup_overlapclipped_sample.list #This is a .list file that shows the entire .bam file path and name for the samples that have been merged, deduped, and overlap clipped.


##Realign around indels across all samples at once

#Make a list of indels
if [ ! -f $BASEDIR/'RossCreekall_samples_for_indel_realigner.intervls' ]; then

java -Xmx120g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REFERENCE -I $BAMLIST -o $BASEDIR/'RossCreekall_samples_for_indel_realigner.intervals' -drf BadMate

fi
## Run the indel realigner tool

java -Xmx120g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REFERENCE -I $BAMLIST -targetIntervals $BASEDIR/'RossCreekall_samples_for_indel_realigner.intervals' --consensusDeterminationModel USE_READS --nWayOut _realigned.bam

#-T is which gatk command you want to run, in this case RealignerTargetCreator and IndelRealigner
#-R is the lake trout reference files needed I think the .fna .fai and .dict files
#-I is the input which are bam alignment files
#-o is where and what I want the output to be called. For RealignerTargetCreator is has to be a .intervals file and -o means I want them all together instead of individual for each sample.
#-drf This makes it so that only reads likely mapped in the right place are used in analysis. If mates map to different contigs, one of them is likely in the wrong place
#If you have a draft genome where the chromosomes are in different contigs then you could have reads mapped correctly but are on different contigs so we want to disbale this.
#For IndelRealigner
#-targetIntervals is the .intervals file made in the previous RealignerTargetCreator step
#--consensusDeterminationModel is USE_READS, which is recommended to balance accuracy and performance.
#--nWayOut is still output but this is specifiying that I want to generate one output for each input
#parameters used are the same as in the therkildsen lcWGS tutorial https://github.com/nt246/lcwgs-guide-tutorial/blob/main/tutorial1_data_processing/markdowns/data_processing.md
