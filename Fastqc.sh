

BASEDIR=/home/cnems/projects/def-ruzza/cnems
SAMPLELIST=$BASEDIR/sample_lists/sample_listfinal.txt #sample list has the full names of all my files

#cat $SAMPLELIST | sed -n '1,800p'| while read -r SAMPLE #Run the first 800 samples from samplelist.
cat $SAMPLELIST | sed -n '801,$p'| while read -r SAMPLE #$ is the last line. Print from line 801 to the end
do #This chunk is reading through the samplelist, then only using the first 800 lines, then while it is reading, for each SAMPLE it will then make a script
#with the below echo command

echo '#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32
#SBATCH -o '$BASEDIR'/error_reports/fastqc_'$SAMPLE'.o

#Load the modules

module load StdEnv/2020
module load nixpkgs/16.09
module load fastqc/0.11.9


fastqc '$BASEDIR'/raw_fastq/211015_A00987_0289_AHLLJTDSX2/'$SAMPLE' -o '$BASEDIR'/fastqc/

scontrol show job ${SLURM_JOB_ID}' > $BASEDIR/shell/fastqc/fastqc_$SAMPLE.sh #scontrol is the details from the job itself, run time etc and it puts it in error report folder. ' is the end of the echo command, writing a file. Then
#it is putting script files into the shell/fastqc folder. fastqc_$SAMPLE is how I wanted it to be named.


sbatch $BASEDIR/shell/fastqc/fastqc_$SAMPLE.sh

done
