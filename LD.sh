#!/bin/bash

#-----------------------------------------------
#Script to calculate LD for individual chromosome
#beagle files
#AllStreams_LDByChromosome.sh
#-----------------------------------------------

#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=16

#This script is designed to calculate LD from the beagle files produced for each chromosome in AllStreams_BeagleByChromosome.sh
#It will make a matrix which I can use for the LD blocks plotting as a heatmap.

module load gcc gsl

LG_FILE=/home/cnems/projects/rrg-ruzza/cnems/ChromosomeList.txt #A list of all 42 (excluding mtDNA) brook trout chromosomes which will will use as input for generating beagle by chromosome
BEAGLEDIR=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/Chromosomes #This is where all the beagles by chromosome are located.
LD_DIR=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/Chromosomes/LD #where I want the new .ld files to go.
N_IND=192 #number of brookies
CHROMOSOME=$1 #Arugment to supply the name of the chromosome so I can run them in parallel

#Take the individual beagle and pos files produced in the previous step AllStreams_BeagleByChromosome.sh and we will be making the .ld files for each chromosome but this is what we want them to be called.
BEAGLE_FILE=$BEAGLEDIR/'AllStreams_SamTools_GLfromSNPList_'$CHROMOSOME'_formatforLD.beagle.gz' #Where these are the beagles per chromosome. i is the name of the chromosome from the LG_FILE
POS_FILE=$BEAGLEDIR/'AllStreams_SamTools_GLfromSNPList_'$CHROMOSOME'.pos.gz' #The position file from each individual chromosome maf file.
LDOUT_FILE=$LD_DIR/'AllStreams_SamTools_GLfromSNPList_'$CHROMOSOME'.ld' #What I want the .ld files to be named.
#calculate the number of sites has to be zcat beacsue the pos file is a .gz
N_SITES=$(zcat $POS_FILE | tail -n +2 | wc -l | cut -d " " -f 1) #Have to have the tail -n +2 so it will skip the header. If not the number of sites wont be correct.
echo "calculating LD on " $CHROMOSOME

#run ngsLD
ngsLD --n_threads 8 --geno $BEAGLE_FILE \
--posH $POS_FILE \
--probs \
--n_ind $N_IND \
--n_sites $N_SITES \
--max_kb_dist 0 \
--max_snp_dist 0 \
--rnd_sample 0.5 \
--out $LDOUT_FILE


#zip the .ld file that is produced.
#gzip $LDOUT_FILE


#When ngsLD is run the number of sites will be grabbed from the N_SITES where it counts the lines in each position file.
#It is going to calculate ld with no max_kb_dist and snp_dist so it will measure between every position
#posH beacuse the position file has a header.
#Randomly sample half of the data as there are a lot of SNPs and it will likely make massive matrices.
