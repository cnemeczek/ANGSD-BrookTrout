#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=16


#Load modules
module load StdEnv/2020 StdEnv/2023
module load angsd/0.940 #use new angsd

REF=/home/cnems/projects/rrg-ruzza/cnems/BT_genome/SfoGenome.fna #want the concatenated SnaGenome fna file.
BAMLIST=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/AllStreams_bamlist.txt #has paths the realiged.bams needed for ANGSD and also the .bai index files of the realigned bams.
PRUNEDSNPLIST=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/AllStreams_getbeagle/LD/AllChromosomes_Pruned_SNP_list.txt #List of all chromosomes pruned SNPs (those that are unlinked)
OUT=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/AllStreams_getbeagle/LD #Where I want the output files to go

#First need to index the pruned SNP list. The index command will make a .bin and .idx file that ANGSD will need in the following commmand. So run this script where the files will be put
angsd sites index $PRUNEDSNPLIST

#Then run ANGSD supplying the SNP list and ANGSD will use the .idx file generated in the previous index command.
angsd -nThreads 8 -bam $BAMLIST -anc $REF -out $OUT/'AllStreams_SamTools_LDpruned' -sites $PRUNEDSNPLIST \
-GL 1 -doGlf 2 -doMajorMinor 3 -doMaf 1 -doPost 1 -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1

#when -doGlf is set to 2 for making the .beagle files this means you also need to set the -doMajorMinor function
#-doMajorMinor is 3 meaning it will use the provided SNP list to calculate allele frequencies at
#-doMaf 1 is calculate allele frequencies based on a fixed major minor and therefor make a .mafs.gz file
#-doPost 1 is caculate posterior probabilties using frequency as a prior.
#Have to have -doCounts 1 in order to have the IBS and Cov mats produced. Also have to have -makeMatrix 1
