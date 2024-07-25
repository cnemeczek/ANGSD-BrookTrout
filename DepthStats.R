

##Ross Creek####
bamdir <- "/home/cnems/projects/rrg-ruzza/cnems/RossCreek/RossCreek_dedups" #Where all my quality filtered, overlapped clipped and dedupped .bam files are (the realigned .bams are also in here)
basedir <- "/home/cnems/projects/rrg-ruzza/cnems/RossCreek" #Where the sample list and merged sample lists are
bam_list <- data.frame(Indiv=readLines(paste0(basedir, "/RossCreek_samplemerge.txt"))) #List of the names of samples like 21_RCU_SFO_05_S22 for example.
bam_list$Loc <- "RossCreek"

depthout <- data.frame(matrix(nrow = length(bam_list$Indiv),ncol = 3))
colnames(depthout) <- c("Loc","Indiv","MeanDepthNonZero")

i <- 1
for (i in 1:length(bam_list$Indiv)){ #For each sample in the list do the following
  print(i)
  bamfile = bam_list[i,1] #make an object call bamfile which refers to the bam list.
  # Compute depth stats
  depth <- readLines(paste0(bamdir, "/", bamfile, "_bt2_SfoGenome_minq20_sorted_dedup_overlapclipped_realigned.bam.depth.gz"))
  print("File is loaded!")
  depth <- as.numeric(depth)

  print("Calculating mean!")
  #adding indiv name and mean depth to output file
  depthout[i,1] <- bam_list[i,2]
  depthout[i,2] <- bamfile
  depthout[i,3] <- mean(depth[depth > 0])
}

write.csv(depthout, "/home/cnems/projects/rrg-ruzza/cnems/RossCreek/RossCreek_dedups/RossCreekDepthSfoGenome.csv") #Export as a .csv
