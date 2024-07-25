#!/bin/bash

#-----------------------------------------------
#Script to calculate hardy weinberg across
#entire chromos where there are potential inversions
#PotentialInversions_HWE.sh
#-----------------------------------------------

#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=16

module load StdEnv/2020 StdEnv/2023
module load angsd/0.940 #using newest angsd.New angsd only uses 8 threads.


#All my outlier data for inversion testing is on Cedar in the OUTLIERS dir so do this there.

#Arguments to specifcy certain files etc.
BAM=$1 #where this is the name of the bamlist to use whether its left, middle, right group.
GROUP=$2 #where this is the group name I want to add to the output file.
SNP=$3 #where this is the SNP list file, whether inside the LD block or outside the LD block on said chromosome.
CHR=$4 #where this is the chromosome list.chrs file that has the name of the chromosome.
CHROMO=$5 #Name of the chromosome that is used at the beginning of the file name ex. Chromosome9.

#dirs.
BAMLIST=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/$BAM  #I have a bam list for each group of chromosome
SNPLIST=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/Chromosomes/$SNP  #SNP list for entire chromosome.
CHRMLIST=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/Chromosomes/$CHR #chromosome list.
OUTPUT=/home/cnems/projects/rrg-ruzza/cnems/BT_Hobs #where I want the .hwe files to go.


angsd -nThreads 8 \
-doHWE 1 \
-GL 1 \
-minMapQ 20 \
-minQ 20 \
-remove_bads 1 \
-doMajorMinor 3 \
-sites $SNPLIST  \
-b $BAMLIST \
-rf $CHRMLIST \
-out $OUTPUT/$CHROMO'_EntireChromo_'$GROUP'_hwe'
