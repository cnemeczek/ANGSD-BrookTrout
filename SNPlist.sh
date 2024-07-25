#!/bin/bash

#-----------------------------------------------
#Script to make SNP list from .mafs.gz
#011_SNPList
#-----------------------------------------------

#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --cpus-per-task=32

module load StdEnv/2020
module load angsd/0.940 #new angsd


STREAM=$1
BASEDIR=/home/cnems/projects/rrg-ruzza/cnems/$STREAM


## Create a SNP list to use in downstream analyses and index is with angsd
gunzip -c $BASEDIR/$STREAM'_SamTools_GL.mafs.gz' | cut -f 1,2,3,4 | tail -n +2 > $BASEDIR/'global_snp_list_'$STREAM'SNPlist.txt'

#This will make a .bin and .idx file with the same prefix as the SNPlist.txt
angsd sites index $BASEDIR/'global_snp_list_'$STREAM'SNPlist.txt'

##Regions format.
cut -f 1,2 $BASEDIR/'global_snp_list_'$STREAM'SNPlist.txt' | sed 's/\t/:/g' > $BASEDIR/'global_snp_list_'$STREAM'SNPlist.regions'

##Get list of chromosomes.
cut -f1 $BASEDIR/'global_snp_list_'$STREAM'SNPlist.txt' | sort | uniq > $BASEDIR/'global_snp_list_'$STREAM'SNPlist.chrs'
