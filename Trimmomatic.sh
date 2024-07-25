#!/bin/bash

#----------------------------------------------------------
#Script for running Trimmomatic on raw Brook Trout fastq files
#by stream. Naming convention will be 002 as trimmomatic is the second step
#----------------------------------------------------------
#SBATCH --time=02-00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32

module load StdEnv/2020
module load trimmomatic/0.39

BASEDIR=/home/cnems/projects/def-ruzza/cnems/raw_fastq/211015_A00987_0289_AHLLJTDSX2 #path all your raw.fastq.gz
ADAPTER=/home/cnems/projects/def-ruzza/cnems/reference/NexteraPE_NT.fa #path to adapters
FILE=/home/cnems/projects/def-ruzza/cnems/RossCreek_Upstream/RossCreekUpstream_SampleNames.txt #file with all samples names including lane number and without the R1_ etc.
TRIMDIR=/home/cnems/projects/def-ruzza/cnems/RossCreek_Upstream/RossCreek_trimmed #where I want the trimmed files to go

LINES=$(cat $FILE)

for name in $LINES

do

java -Xmx120g -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 30 -phredd33 $BASEDIR/${name}_R1_001.fastq.gz $BASEDIR/${name}_R2_001.fastq.gz $TRIMDIR/${name}_R1_paired.fastq.gz $TRIMDIR/${name}_R1_unpaired.fastq.gz $TRIMDIR/${name}_R2_paired.fastq.gz $TRIMDIR/${name}_R2_unpaired.fastq.gz ILLUMINACLIP:$ADAPTER:2:30:10:1:true MINLEN:40
#PE for paired end and phred33 based on the illumina sequencing. In the base dir, get the sample name followed by _R1_001.fastq.gz because that is how the end of the files are named in the basedir and then repeat for the reverse R2 reads. Then put the output in the TRIMDIR using the name which is the sample name followed by R1 unpaired and paired and R2 paired and unpaired.

done  |& tee -a trim.log
