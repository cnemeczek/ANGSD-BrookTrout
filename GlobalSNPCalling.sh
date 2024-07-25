#!/bin/bash

#-----------------------------------------------
#Script to calculate counts for depth filtering
#012_GL_SamTools
#-----------------------------------------------

#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=32

#This script is to generate genotype likelihoods using max and min depth determined through .depthGlobal plot.
#-doGLF 2 so it makes the .beagle files needed for PCA etc.

#Filter bams and do -GL 1 using SamTools
module load StdEnv/2020 StdEnv/2023
module load angsd/0.940 #Use new angsd

BASEDIR=/home/cnems/projects/rrg-ruzza/cnems/Poole #Where I want files to go.
BAMLIST=$BASEDIR/Poole_dedups/Poole_realignedbamlist.txt #Shows names and paths to the bam files for Ross Creek
REFERENCE=/home/cnems/projects/rrg-ruzza/cnems/BT_genome/SfoGenome.fna #want the concatenated SnaGenome fna file.
MINDEP=54
MAXDEP=92
STREAM=Poole

#Filter bams and do -GL 1 using SamTools
angsd -nThreads 30 -b $BAMLIST -ref $REFERENCE -out $BASEDIR/$STREAM'_SamTools_GL' \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minInd 10 -minMapQ 20 -minQ 20 \
-setMinDepth $MINDEP -setMaxDepth $MAXDEP -doCounts 1 -GL 1 -doGLF 2 -doMajorMinor 1 -doIBS 1 -doCov 1 -makeMatrix 1 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05

#-out path and prefix of what I want the files to be called
#-setMinDepth and -setMaxDepth set at these numbers based on plotting and filtering of the .depthGlobal file from ANGSD.
#-doCounts is count the number of A,C,G,T at all sites
#Have to have the -doCounts option for setMinDepth and setMaxDepth
#By giving the reference it looks for both the .fna genome file and the samtools faidx index .fai file
#Filtering options -minInd is remove sites where half the individuals have no data. If set at 96, means only data from at least 96 individuals.
#-minMapQ is minimum mapping quality
#-minQ would be individual base quality I believe
#-remove_bads 1 keep those not tagged as bad
#-only_proper_pairs 1 keep only proper pairs
#-C 50 is reduces the effect of reads with excessive mismatching
#-baq 1 computes base alignment quality which rules out false SNPs close to indels.
#When -doGlf is set to 2 for making the .beagle files this means you also ened to set the -doMajorMinor function
#-doMajorMinor 1 this means infer the major and minor from the genotype likelihoods.

#With -doGlf 2 a beagle.gz file will be made with the likelihoods
