#!/bin/bash

#-----------------------------------------------
#Script to calculate counts for depth filtering
#011-depthANGSD
#-----------------------------------------------

#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --cpus-per-task=16

#This script is to do depth counts on all sites to set depth filters

#Don't add -GL 1 -doMajorMinor 1 -doMaf 1

module load StdEnv/2020 StdEnv/2023
module load angsd/0.940

BASEDIR=/home/cnems/projects/rrg-ruzza/cnems/RossCreek
BAMLIST=$BASEDIR/RossCreek_dedups/RossCreek_realignedbamlist.txt #Shows the paths with names of the bam files I want (realigned.bams) depth.gz??
REFERENCE=/home/cnems/projects/rrg-ruzza/cnems/BT_genome/SfoGenome.fna #want the concatenated SfoGenome fna file.


angsd -nThreads 16 -b $BAMLIST -ref $REFERENCE -out $BASEDIR/RossCreek_dedups/RossCreek_depth/'RossCreekDepthCount' -doDepth 1 -maxDepth 210 -doCounts 1 -dumpCounts 1 -minMapQ 20 -minQ 20 -remove_bads 1 -only_proper_pairs 1 -C 50

#-out path and prefix of what I want the files to be called
#-maxDepth set at 210 (21 indivudals times 10 as an individual max depth cut off= 210 for all)
#-doCounts is count the number of A,C,G,T at all sites
#-dumpCounts 1 is I think generate total sequencing depth for each site creating .pos.gz
#-minMapQ is minimum mapping quality
#-minQ would be individual base quality I believe
#-remove_bads 1 keep those not tagged as bad
#-only_proper_pairs 1 keep only proper pairs
#-C 50 is reduces the effect of reads with excessive mismatching
#This script will make a file .depthGlobal which has the counts for all sites summed across all samples.
