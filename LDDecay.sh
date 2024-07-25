#!/bin/bash

#-----------------------------------------------
#Script to compute LD decay on LD file
#AllStreams_LDDecay.sh
#-----------------------------------------------

#SBATCH --time=03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --account=
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --cpus-per-task=16

#Run R through Compute Canada
module load StdEnv/2020 r/4.1.2
module load gcc gsl


#Run in the basedir where the .ld file is
#Prior to running this script. need to compress the .ld file.

BASEDIR=/home/cnems/projects/rrg-ruzza/cnems/AllStreams/AllStreams_getbeagle/LD

NGSLD=/home/cnems/ngsLD/scripts #this is where the script for LD decay is located


#read the ld.gz file which is in BASEDIR and then bash the decay script

ls $BASEDIR/'AllStreams_maxkb100_fordecay.ld.gz' | Rscript --vanilla --slave $NGSLD/'fit_LDdecay.R' \
--ld r2 \
--n_ind 192 \
--fit_level 100 \
--fit_boot 100 \
--max_kb_dist 200 \
--fit_bin_size 100 \
--plot_data \
--plot_scale 3 \
--out $BASEDIR/'100kbld_200kbdecay_fitbinsize1000_Decayplot.pdf'

#Where the ld files is a zipped version of the .ld file to use for caculating decay.
#specifying number of individuals is relevant for fitting the r squared decay model which is the first of two models the decay does.
#fit_level 1 is fit using Nelder-Mead, 2 is Nelder-Mead and BFGS, 3 is Nelder-Mead, BFGS, L-BFGS-B. Setting at 100 means run through the 3 algorithms 100 times and pick the best
#--fit_boot is do 100 bootstrap replicates
#max_kb_dist is max distance between SNPs in KB
#--fit_bin_size is bin data into fixed-sized windows for fitting where the default is 250bp
#--plot_data means put the data points on
#--plot_scale means increase the size of the plot so everything fits nicely.
