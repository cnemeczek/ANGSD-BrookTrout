#Script to generate the pruned SNP list for each chromosome using the SNPlist.mafs.gz and LD pruned files from the ngsLD perl prune_graph.pl script

basedir="/home/cnems/projects/rrg-ruzza/cnems/AllStreams/Chromosomes/" #This is where the individual chromosome.mafs.gz, Chromosome_List.txt
chromdir="/home/cnems/projects/rrg-ruzza/cnems/" #where the chromosome list is
LDdir="/home/cnems/projects/rrg-ruzza/cnems/AllStreams/AllStreams_getbeagle/LD/" #This is where the .id pruned files are per chromosome.

#chromosomes that not have 0 unlinked
chromosome_list <- read.delim(paste0(chromdir, "ChromosomeList.txt"), header = F)$V1
typeof(chromosome_list)

for (name in chromosome_list) {

  pruned_position <- as.integer(gsub(".*:","", readLines(paste0(LDdir,"AllStreams_maxkb15_SamToolsSNPlist_afterdecay_unlinked_",name,".id")))) #the files of unlinked SNPs where name is the name of the chromosome

  snp_list <- read.table(paste0(basedir,"AllStreams_SamTools_GLfromSNPList_",name,".mafs.gz", sep=""), stringsAsFactors = F, header = F) [,1:4] #This is the .mafs.gz produced from step 2 of ANGSD pipeline (012)

  pruned_snp_list <- snp_list[snp_list$V2 %in% pruned_position, ]

  write.table(pruned_snp_list, paste0(LDdir,"AllStreams_maxkb15_SamTools_LDpruned_SNP_list_",name, sep = ""), col.names = F, row.names = F, quote = F, sep =  "\t")

}

#Read through the chromosome list and for each chromosome name
#for pruned position, take each individual chromosome unlinked. id file which starts with the name of the .ld file that was used for the previous pruning step
#for snp_list it is reading the .mafs.gz from 012_getbeagle SNP list which gives the positions of all the SNPs where column 1:4 are chromosome, position, major and minor alleles.
#pruned_snp_list is I think taking the position ($V2) from the original SNP list and finding it in the pruned_position dataframe.
#paste0 means put together name and the AllStreams together to make the new file name.
#gsub is replace .*: in the unlinked.id files with just a no space, so get rid of them.
