#!/bin/bash

#-----------------------------------------------
#Script to compute LD
#AllStreamsByPopulation_017_LD_OnSamTools_SNPList_AfterDecay.sh
#-----------------------------------------------

#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=16

module load gcc gsl

#Can run this in the basedir so no need to have cd BASEDIR in the script.

#This script is going to calculate LD using .beagle.gz formatted genotype likelihood files from ANGSD. I also need the .mafs.gz files from the same run.
#I want the files from the SNP calling with the SNP list using SamTools


BASEDIR=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/AllStreams_getbeagle #This is where the .beagle.gz and .mafs.gz files from samtools SNP list run are located
OUTPUT=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/AllStreams_getbeagle/LD #Where I want the output .tsv files to go.

ngsLD --n_threads 16 --geno $BASEDIR/'AllStreams_SamTools_GLfromSNPList_minus4ColforLD.beagle.gz' \
--posH $BASEDIR/'AllStreams_SamTools_GLfromSNPList.pos.gz' \
--probs \
--n_ind 192 \
--n_sites 2620282 \
--max_kb_dist 15 \
--out $OUTPUT/'AllStreams_maxkb15_SamToolsSNPlist_afterdecay.ld'

#I am supplying genotype likelihoods so I use --geno which is the beagle.gz
#--posH is the pos.gz file with chromosome and position which was cut from the .mafs.gz. It has a header so H included in --posH
#--n_ind is number of brookies
#--n_sites is number of SNPs which you get from counting lines in the .mafs.gz
#--max_kb_dist is set to 10 so that it assumes SNPs in less than 10kb windows are linked.
#--ouput prefix for output files which will end in .tsv.
#As output ngsLD makes a .tsv file with LD results for all pairs of sites for which LD was calculated.
