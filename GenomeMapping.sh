#!/bin/bash

#----------------------------------------------------------
#Script for mapping paired reads to SfoGenome using bowtie2
#with end to end (--very-sensitive)
#005-genomemapping
#----------------------------------------------------------

#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --account=
#SBATCH --mail-user=
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32

module load StdEnv/2020
module load gcc/9.3.0
module load bowtie2/2.4.4
module load samtools/1.17


BASEDIR=/home/cnems/projects/rrg-ruzza/cnems/RossCreek
FILE=$BASEDIR/RossCreekUpstream_SampleNames.txt
FASTQFILT=/home/cnems/projects/def-ruzza/cnems/RossCreek_Upstream/RossCreek_polytrimmed #Where the adaptor trimmed and polyg trimmed files are. Only want the paired end data for R1 and R2.
REFERENCE=/home/cnems/projects/rrg-ruzza/cnems/BT_genome/SfoGenome #want the index files .bt2 Bowtie will know to use the bt2 files even though there are other files types here so just put the name of the genome which is at the front of all file names in the folder
REFNAME=SfoGenome #The name of my reference genome
MAPPEDBAM=$BASEDIR/RossCreek_mapped #where I want the mapped and filtered reads to go.

LINES=$(cat $FILE)

for name in $LINES
do
 bowtie2 -q --phred33 --very-sensitive -p 30 -I 0 -X 1500 --fr --rg-id $name --rg SM:$name --rg LB:$name -x $REFERENCE -1 $FASTQFILT/${name}_R1_paired_qual_filtered.fastq.gz -2 $FASTQFILT/${name}_R2_paired_qual_filtered.fastq.gz -S $MAPPEDBAM/${name}'_paired_bt2_'$REFNAME'.sam'

#Convert to bam file for storage (including the mapped reads)
 samtools view -@ 30 -bS -F 4 $MAPPEDBAM/${name}'_paired_bt2_'$REFNAME'.sam' > $MAPPEDBAM/${name}'_paired_bt2_'$REFNAME'.bam' #make the .sam file to a .bam file
 rm -f $MAPPEDBAM/${name}'_paired_bt2_'$REFNAME'.sam' #remove the sam files.
 #Filter the mapped reads to only retain reads that have quality of 20
 #Filter the bam files to get rid of non-unique mappings and those with quality less than 20.
 samtools view -@ 30 -h -q 20 $MAPPEDBAM/${name}'_paired_bt2_'$REFNAME'.bam' | samtools view -@ 30 -buS - | samtools sort -m 2G -@ 30 -o $MAPPEDBAM/${name}'_paired_bt2_'$REFNAME'_minq20_sorted.bam'
done

#-fr is the option for expected relative orientation of the mates in sets.
#-I and -X is the expected range of inter-mate distances. This if for concordant pairs- Reads are concordant if a pair aligns with the mate orientation (-I and -X) and with the expected range of distances.
#-I is the minimum fragment length for valid paired-end alignments. I think setting this at 0 means there can't be gaps between alignments?
#when you set -I and -X bowtie scans that size window to figure out if there is concordant alignments.
#-q indicates that the input reads are fasta files
#p is threads
#-fr if there is a candidate paired-end alignment where 1 mate is upstream of the reverse complement of mate 2 and -I and -X are met, then the alignment is valid.
#--rg-id makes the read id text that will attach to the SAM output. *It makes the header for the SAM file. SM is Sample, LB is library, PU is platform unit (eg flow cell lane). unique ID, PL is platform used to make the reads
#I think this is essentially telling what the SAM file should hold.
#-x is the name of the reference genome.
#-1 is the R1 paired adapter/polyg trimmed fastq files
#2 is the R2 paired adapter/polyg trimmed fastq files
#-S is what I want the SAM files to be named and where to put them
#samtool view -bS -b is bam and -S is automatically detects the type of file? -F is do not output alignment with any bits set in FLAG present in the FLAG field. WHat is 4?Is this only the reads that mapped successfully?
#samtools -h means include the header -q skip alignments with MAPQ smaller than 20 -buS is output in bam that are -u uncompressed and automatically detect sample file (needed if files are SAM)
