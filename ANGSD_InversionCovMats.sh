#!/bin/bash

#-----------------------------------------------
#Script to use SNP list and get genotype
#likelihoods from those sites by chromosome in beagle format
#Angsd_covmat_bychromo_forPCA.sh
#-----------------------------------------------

#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32

#This script is to generate genotype likelihoods at the SNP list sites generated from the last setMaxDepth
#run this in scratch and its calling on files in projects but will put the new generated files in scratch.

module load StdEnv/2020 StdEnv/2023
module load angsd/0.940 #Use newest angsd.

#do -GL 1 using SamTools
#Run this in the output dir.

OUTPUT=/home/cnems/projects/rrg-ruzza/cnems/InversionCovmats #where I want the cov mats and beagles etc to go.
BAMLIST=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/AllStreams_bamlist.txt #has paths the realiged.bams needed for ANGSD and also the .bai index files of the realigned bams.
REFERENCE=/home/cnems/projects/rrg-ruzza/cnems/BT_genome/SfoGenome.fna #want the concatenated SfoGenome fna file.
SNPLIST=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/Chromosomes/$1 #The SNP list generated from .mafs.gz file. ALso need the index .idx of the SNP list here for the individual chromosomes.
CHRLIST=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/Chromosomes/$2 #The chromosome list also generated from the .mafs.gz file for individual chromosomes.
CHROMO=$3 #name of chromosome to add to output.


angsd -nThreads 30 -b $BAMLIST -ref $REFERENCE -out $OUTPUT/'AllStreams_GLfromSNPList_'$CHROMO'_LDBlock' \
-GL 1 \
-doGlf 2 \
-doMaf 1 \
-doMajorMinor 3 \
-sites $SNPLIST \
-rf $CHRLIST \
-doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1


#When -doGlf is set to 2 for making the .beagle files this means you also ened to set the -doMajorMinor function
#-doMajorMinor is 3 meaning it will use the provided SNP list to calculate allele frequencies at
#-doMaf 1 is calculate allele frequencies based on a fixed major minor and therefor make a .mafs.gz file
#-sites is the SNP list generated from the Global SNP calling script in the previous step.
#-rf is the list of SNPs in chromosomes also generate from the Global SNP calling script.
#make the covariance matrix so I can do a PCA, although the PCA will likely be the same as in the last step.
