#!/bin/bash

#SBATCH --time=3-00:00
#SBATCH --job-name=ld_pruning
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --nodes=1
#SBATCH --array=1-43
#SBATCH --ntasks=1
#SBATCH --mem=10000M
#SBATCH --tmp=120G

## Define some variables
INPUT_PATH=/home/cnems/scratch/All_Streams/All_Streams_bams/LD/ #Where the .ld input file is located.
LD=AllStreams_SamToolsSNPlist_AfterDecay.ld #This is the .ld file generated after determining that the max_kb_dist for brook trout should be 20kb. (second run of LD after first testing with 100-500kb distance)
MAXDIST=20 #kb distance between SNPs to say there is no decay. This value was chosen from the ld decay curve in step 018 where decay for brook trout dropped substantially at about 20kb.
MINWEIGHT=0.4 #Setting this at 0.4 because based on my decay graph, this is the r squared or LD rate at which the decay rate drops and it seems like other papers are choosing values that match the rate of decay. --min_weight is the minimum weight of an edge (linkage levels between them) to assume SNPs are connected
PRUNE_GRAPH=/home/cnems/ngsLD/scripts/prune_graph.pl #this is where the script for LD decay is located
LG_LIST=/home/cnems/scratch/All_Streams/All_Streams_bams/Chromosome_List.txt #This chrom list was made from the .pos.gz of the .mafs.gz from step 2 get_beagle script (_012)


##################################################

## Keep a record of the Job ID
echo $SLURM_JOB_ID
## Create and move to working directory for job
WORKDIR=${SLURM_TMPDIR}/
cd $WORKDIR
## Transfer the input files
cp $INPUT_PATH$LD $WORKDIR
## Select the LG from the input file. LG is the linkage group? ie the chromosome.
LG=`head $LG_LIST -n $SLURM_ARRAY_TASK_ID | tail -n 1`
grep ^${LG}: $WORKDIR$LD > $WORKDIR${LD%%.*}_${LG}.ld
## Define the output name
OUT=${LD%%.*}_unlinked_${LG}.id
## Run the perl script
perl $PRUNE_GRAPH \
--in_file $WORKDIR${LD%%.*}_${LG}.ld \
--max_kb_dist $MAXDIST \
--min_weight $MINWEIGHT \
--out $WORKDIR$OUT
## Move output files back
cp $WORKDIR$OUT $INPUT_PATH
