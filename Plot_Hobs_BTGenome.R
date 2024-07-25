#This script is for plotting the Hobs results of hwe files

library(RColorBrewer) #1.1.3
library(viridis) #0.6.5
library(ggplot2) #3.5.0
library(patchwork)


####Chromosome 12####
#CM055694.1

######Left Group#####

#read in the data

hwe_chromo12_left <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome12_EntireChromo_leftgroup_hwe.hwe.slidingwindows', header = T)


#Now add on the windscanr results.

Win_chromo12_left <- ggplot(hwe_chromo12_left, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 12", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 12 non-inverted homokaryotypes", breaks =  seq(from =0, to = 65000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo12_left

ggsave("BTChromo12_AverageHobsEntireChromosome_LeftGroup.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)



######Middle Group#####

#read in the data

hwe_chromo12_middle <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome12_EntireChromo_middlegroup_hwe.hwe.slidingwindows', header = T)

#Now add on the windscanr results.

Win_chromo12_middle <- ggplot(hwe_chromo12_middle, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 12", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 12 heterokaryotypes", breaks =  seq(from =0, to = 65000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo12_middle

ggsave("BTChromo12_AverageHobsEntireChromosome_MiddleGroup.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)


#####Middle with Outliers####
#Middle with outliers refers to the addition of a few individuals that looked like they were between the heterkaryotypes and homokarytypes
#in the AllStreamsPCA_InversionGroupsCircled_March13_2024 PCA.png.

#read in the data

hwe_chromo12_middleout <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome12_EntireChromo_middlegroupoutliers_hwe.hwe.slidingwindows', header = T)

#Now add on the windscanr results.

Win_chromo12_middleout <- ggplot(hwe_chromo12_middleout, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 12", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 12 outlier heterokaryotypes", breaks =  seq(from =0, to = 65000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo12_middleout

ggsave("BTChromo12_AverageHobsEntireChromosome_MiddleGroupOutliers.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

#Stick these two middle group plots together.

patchwork = Win_chromo12_middle + Win_chromo12_middleout + plot_layout(guides = "collect", ncol = 3)
patchwork

ggsave("Chrom12AverageHobs_middlegroupvsoutliers.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 14, height = 8, units = c("in"), dpi = 300)



######Right Group#####

hwe_chromo12_right <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome12_EntireChromo_rightgroup_hwe.hwe.slidingwindows', header = T)


#Now add on the windscanr results.

Win_chromo12_right <- ggplot(hwe_chromo12_right, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 12", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 12 potential inverted homokaryotypes", breaks =  seq(from =0, to = 65000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  
Win_chromo12_right

ggsave("BTChromo12_AverageHobsEntireChromosome_RightGroup.png", path = "C:/Users/cneme/Documents/Thesis/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

######Raw Right group data####
#I think I need a dotplot of the raw data, not the windowed version so that I can pinpoint if my inversion area was correct
#I need the .Hobs file before it was run through windowscanr for this. 

chromo12_right_hobs <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome12_EntireChromo_rightgroup_hwe.hwe.Hobs', header = T)


chromo12_right_hobs <- ggplot(chromo12_right_hobs, aes(Position, Hobs)) +
  geom_point() +
  #xlim(0, 40000000)
  scale_x_continuous(breaks = seq(0,40000000, 2000000)) 
chromo12_right_hobs

#From this plot can see that the SNP list should be adjusted to be from 14-32 million bp. Original based on LD heatmap was 14-35.

####Chromosome 19####
#CM055701.1

######Left Group#####

#read in the data

hwe_chromo19_left <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome19_EntireChromo_leftgroup_hwe.hwe.slidingwindows', header = T)


#Now add on the windscanr results.

Win_chromo19_left <- ggplot(hwe_chromo19_left, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 19", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 19 non-inverted homokaryotypes", breaks =  seq(from =0, to = 60000000, by = 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo19_left

ggsave("BTChromo19_AverageHobsEntireChromosome_LeftGroup.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)



######Middle Group#####

#read in the data

hwe_chromo19_middle <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome19_EntireChromo_middlegroup_hwe.hwe.slidingwindows', header = T)

#Now add on the windscanr results.

Win_chromo19_middle <- ggplot(hwe_chromo19_middle, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 19", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 19 heterokaryotypes", breaks =  seq(from =0, to = 60000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo19_middle

ggsave("BTChromo12_AverageHobsEntireChromosome_MiddleGroup.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)




#Stick these two middle group plots together.

patchwork = Win_chromo12_middle + Win_chromo12_middleout + plot_layout(guides = "collect", ncol = 3)
patchwork

ggsave("Chrom12AverageHobs_middlegroupvsoutliers.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 14, height = 8, units = c("in"), dpi = 300)


######Right Group#####

hwe_chromo19_right <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome19_EntireChromo_rightgroup_hwe.hwe.slidingwindows', header = T)


#Now add on the windscanr results.

Win_chromo19_right <- ggplot(hwe_chromo19_right, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 19", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 19 potential inverted homokaryotypes", breaks =  seq(from =0, to = 60000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  
Win_chromo19_right
#THis one is interesting, there are at least 2 dips in Hobs

ggsave("BTChromo19_AverageHobsEntireChromosome_RightGroup.png", path = "C:/Users/cneme/Documents/Thesis/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

######Raw Right group data####
#I think I need a dotplot of the raw data, not the windowed version so that I can pinpoint if my inversion area was correct
#I need the .Hobs file before it was run through windowscanr for this. 

chromo19_right_hobs <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome19_EntireChromo_rightgroup_hwe.hwe.Hobs', header = T)


chromo19_right_hobs <- ggplot(chromo19_right_hobs, aes(Position, Hobs)) +
  geom_point() +
  #xlim(0, 40000000)
  scale_x_continuous(breaks = seq(0,60000000, 2000000)) 
chromo19_right_hobs

#This one is interesting, there's almost 3 sections of dips in Hobs
# but I am going to go with the one that corresponds to the darkest
#blob of LD heatmap which would be between 21-33mp.




####Chromosome 27####
#CM055709.1

######Left Group#####

#read in the data

hwe_chromo27_left <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome27_EntireChromo_leftgroup_hwe.hwe.slidingwindows', header = T)


#Now add on the windscanr results.

Win_chromo27_left <- ggplot(hwe_chromo27_left, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 27", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 27 potential inverted homokaryotypes", breaks =  seq(from =0, to = 55000000, by = 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo27_left

ggsave("BTChromo27_AverageHobsEntireChromosome_LeftGroup.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)



######Middle Group#####

#read in the data

hwe_chromo27_middle <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome27_EntireChromo_middlegroup_hwe.hwe.slidingwindows', header = T)

#Now add on the windscanr results.

Win_chromo27_middle <- ggplot(hwe_chromo27_middle, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 27", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 27 heterokaryotypes", breaks =  seq(from =0, to = 55000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo27_middle

ggsave("BTChromo27_AverageHobsEntireChromosome_MiddleGroup.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)


#####Middle with Outliers####
#Middle with outliers refers to the addition of a few individuals that looked like they were between the heterkaryotypes and homokarytypes
#in the AllStreamsPCA_InversionGroupsCircled_March13_2024 PCA.png.

#read in the data

hwe_chromo27_middleout <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome27_EntireChromo_middlegroupoutliers_hwe.hwe.slidingwindows', header = T)

#Now add on the windscanr results.

Win_chromo27_middleout <- ggplot(hwe_chromo27_middleout, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 27", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 27 outlier heterokaryotypes", breaks =  seq(from =0, to = 55000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo27_middleout

ggsave("BTChromo27_AverageHobsEntireChromosome_MiddleGroupOutliers.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

#Stick these two middle group plots together.

patchwork2 = Win_chromo27_middle + Win_chromo27_middleout + plot_layout(guides = "collect", ncol = 3)
patchwork2

ggsave("Chrom27AverageHobs_middlegroupvsoutliers.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 14, height = 8, units = c("in"), dpi = 300)



######Right Group#####

hwe_chromo27_right <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome27_EntireChromo_rightgroup_hwe.hwe.slidingwindows', header = T)


#Now add on the windscanr results.

Win_chromo27_right <- ggplot(hwe_chromo27_right, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 27", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 27 non-inverted homokaryotypes", breaks =  seq(from =0, to = 55000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  
Win_chromo27_right
#THis one is interesting, there are at least 2 dips in Hobs

ggsave("BTChromo27_AverageHobsEntireChromosome_RightGroup.png", path = "C:/Users/cneme/Documents/Thesis/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

######Raw left group data####
#I think I need a dotplot of the raw data, not the windowed version so that I can pinpoint if my inversion area was correct
#I need the .Hobs file before it was run through windowscanr for this. 

chromo27_left_hobs <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome27_EntireChromo_leftgroup_hwe.hwe.Hobs', header = T)


chromo27_left_hobs <- ggplot(chromo27_left_hobs, aes(Position, Hobs)) +
  geom_point() +
  #xlim(0, 40000000)
  scale_x_continuous(breaks = seq(0,50000000, 2000000)) 
chromo27_left_hobs

#Inversion chunk between 21-33mb.



####Chromosome 31####
#CM055713.1

######Left Group#####

#read in the data

hwe_chromo31_left <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome31_EntireChromo_leftgroup_hwe.hwe.slidingwindows', header = T)


#Now add on the windscanr results.

Win_chromo31_left <- ggplot(hwe_chromo31_left, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 31", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 31 potential inverted homokaryotypes", breaks =  seq(from =0, to = 60000000, by = 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo31_left

ggsave("BTChromo31_AverageHobsEntireChromosome_LeftGroup.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)



######Middle Group#####

#read in the data

hwe_chromo31_middle <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome31_EntireChromo_middlegroup_hwe.hwe.slidingwindows', header = T)

#Now add on the windscanr results.

Win_chromo31_middle <- ggplot(hwe_chromo31_middle, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 31", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 31 heterokaryotypes", breaks =  seq(from =0, to = 60000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo31_middle

ggsave("BTChromo31_AverageHobsEntireChromosome_MiddleGroup.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

#####Middle with Outliers####
#Middle with outliers refers to the addition of a few individuals that looked like they were between the heterkaryotypes and homokarytypes
#in the AllStreamsPCA_InversionGroupsCircled_March13_2024 PCA.png.

#read in the data

hwe_chromo31_middleout <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome31_EntireChromo_middlegroupoutliers_hwe.hwe.slidingwindows', header = T)

#Now add on the windscanr results.

Win_chromo31_middleout <- ggplot(hwe_chromo31_middleout, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 31", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 31 outlier heterokaryotypes", breaks =  seq(from =0, to = 60000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  

Win_chromo31_middleout

ggsave("BTChromo31_AverageHobsEntireChromosome_MiddleGroupOutliers.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

#Stick these two middle group plots together.

patchwork3 = Win_chromo31_middle + Win_chromo31_middleout + plot_layout(guides = "collect", ncol = 3)
patchwork3

ggsave("Chrom31AverageHobs_middlegroupvsoutliers.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 14, height = 8, units = c("in"), dpi = 300)




######Right Group#####

hwe_chromo31_right <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome31_EntireChromo_rightgroup_hwe.hwe.slidingwindows', header = T)


#Now add on the windscanr results.

Win_chromo31_right <- ggplot(hwe_chromo31_right, aes(win_mid, Hobs_mean)) + geom_point(alpha = 0.5, colour = "grey") +
  geom_point() +
  expand_limits(x=0,y=0) +
  ylim(0,1) +
  labs(x= "Chromosome 31", y= "Average Hobs") +
  scale_x_continuous(name = "Chromosome 31 non-inverted homokaryotypes", breaks =  seq(from =0, to = 60000000, by= 5000000)) + #to is the length of the entire chromosome in bp.
  theme(text = element_text(size = 12), axis.text.x = element_blank())  
Win_chromo31_right
#THis one is interesting, there are at least 2 dips in Hobs

ggsave("BTChromo31_AverageHobsEntireChromosome_RightGroup.png", path = "C:/Users/cneme/Documents/Thesis/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

######Raw left group data####
#I think I need a dotplot of the raw data, not the windowed version so that I can pinpoint if my inversion area was correct
#I need the .Hobs file before it was run through windowscanr for this. 

chromo31_left_hobs <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome31_EntireChromo_leftgroup_hwe.hwe.Hobs', header = T)


chromo31_left_hobs <- ggplot(chromo31_left_hobs, aes(Position, Hobs)) +
  geom_point() +
  #xlim(0, 40000000)
  scale_x_continuous(breaks = seq(0,60000000, 2000000)) 
chromo31_left_hobs

#Inversion chunk between 14-40mb.

##Plot of Hobs accross entire 4 Inversion Chromos####

#This will make it so all inverted groups are on the right side of the plot and non-inverted on the left

patchwork = Win_chromo12_left + Win_chromo12_middle + Win_chromo12_right + 
  Win_chromo19_left + Win_chromo19_middle + Win_chromo19_right+
  Win_chromo27_right+ Win_chromo27_middle + Win_chromo27_left +
  Win_chromo31_right + Win_chromo31_middle + Win_chromo31_left + plot_layout(guides = "collect", ncol = 3)
patchwork


#This one is good regarding text size but I need to try to get rid of duplicate x axis labels

patchwork[[2]] =patchwork[[2]] + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank() )

patchwork[[3]] =patchwork[[3]] + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank() )

patchwork[[5]] =patchwork[[5]] + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank() )

patchwork[[6]] =patchwork[[6]] + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank() )

patchwork[[8]] =patchwork[[8]] + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank() )

patchwork[[9]] =patchwork[[9]] + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank() )

patchwork[[11]] =patchwork[[11]] + theme(axis.text.y = element_blank(),
                                           axis.ticks.y = element_blank(),
                                           axis.title.y = element_blank() )

patchwork[[12]] =patchwork[[12]] + theme(axis.text.y = element_blank(),
                                           axis.ticks.y = element_blank(),
                                           axis.title.y = element_blank() )
patchwork

ggsave("4Chromos_AverageHobs_noxaxislabelsInverted_NotInverted.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 14, height = 8, units = c("in"), dpi = 300)


###Boxplots for Hobs####

#For my thesis I re-ran Hobs with adjusted SNP lists based on the above scatter plots because you can see more precisely where the inversions
#start and stop and then did the boxplots based on those SNP lists.

#The data was produced using a 10kb window with 10kb slide. For my thesis I also ended up doing it just with a window the size of the inversions region
#Below data was produced using a window the size of each inversion with no slide.

######Chromosome 12####
#Read in the left group Hobs data that went through windowscanr

chromo12_left <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome12_inversionregion_rightgroup.hwe.Hobs.windowrolling', header = T)

#add a column called group to indicate which group from the PCA this is to combine data from left,middle,right
#where for this chromo left is inverted homos, middle is heteros, and right is non-inverted homos.

chromo12_left$group <- 'Non-inv.'

#read in the middle group Hobs data that went through windowscanr.

chromo12_middle <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome12_inversionregion_middlegroup.hwe.Hobs.windowrolling', header = T)

chromo12_middle$group <- 'Heterokaryotypes'


#read in the right group Hobs data that went through windowscanr.

chromo12_right <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome12_inversionregion_rightgroup.hwe.Hobs.windowrolling', header = T)

chromo12_right$group <- 'Inv.'

#Bind the 3 dataframes.

Chromo12 <- rbind(chromo12_left, chromo12_middle, chromo12_right)

#Then make the group a factor so it will plot in the order I want along the x axis.

Chromo12$group <- factor(Chromo12$group, levels = c("Inv.", "Heterokaryotypes", "Non-inv."))


#Now try to make a box plot with the 3 different groups
#Note that this will be the mean of means if I am using the windowscanr results.

box_chromo12 <- ggplot(Chromo12, aes(x= group, y = Hobs_mean, fill = group)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(limits = c(0,0.6)) +
  #scale_colour_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  labs(x="Group", y= "Average Hobs", title = "a) Chromosome 12") +
  theme(axis.text.x = element_text(size=19), axis.text.y = element_text(size = 19), axis.title = element_text(size = 19),
        plot.title = element_text(size = 19), axis.title.x = element_blank(),
        panel.background = element_blank())
box_chromo12

#This works but I have a bunch of Hobs of 1 and close to 1 which is because of the window....

ggsave("BTChromo12_AverageHobs_Boxplot.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)


######Chromosome 19####
#Read in the left group Hobs data that went through windowscanr

chromo19_left <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome19_inversionregion_leftgroup.hwe.Hobs.windowrolling', header = T)

#add a column called group to indicate which group from the PCA this is to combine data from left,middle,right
#where for this chromo left is inverted homos, middle is heteros, and right is non-inverted homos.

chromo19_left$group <- 'Non-inv.'

#read in the middle group Hobs data that went through windowscanr.

chromo19_middle <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome19_inversionregion_middlegroup.hwe.Hobs.windowrolling', header = T)

chromo19_middle$group <- 'Heterokaryotypes'

#read in the right group Hobs data that went through windowscanr.

chromo19_right <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome19_inversionregion_rightgroup.hwe.Hobs.windowrolling', header = T)

chromo19_right$group <- 'Inv.'

#Bind the 3 dataframes.

Chromo19 <- rbind(chromo19_left, chromo19_middle, chromo19_right)

#Then make the group a factor so it will plot in the order I want along the x axis.

Chromo19$group <- factor(Chromo19$group, levels = c("Inv.", "Heterokaryotypes", "Non-inv."))


#Now try to make a box plot with the 3 different groups
#Note that this will be the mean of means if I am using the windowscanr results.

box_chromo19 <- ggplot(Chromo19, aes(x= group, y = Hobs_mean, fill = group)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(limits = c(0,0.6)) +
  #scale_colour_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  labs(x="Group", y= "Average Hobs", title = "b) Chromosome 19") +
  theme(axis.text.x = element_text(size=19), axis.text.y = element_text(size = 19), axis.title = element_text(size = 19),
        plot.title = element_text(size = 19), axis.title.x = element_blank(),
        panel.background = element_blank())
box_chromo19

#This works but I have a bunch of Hobs of 1 and close to 1 which is because of the window....

ggsave("BTChromo12_AverageHobs_Boxplot.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)



######Chromosome 27####
#Read in the left group Hobs data that went through windowscanr

chromo27_left <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome27_inversionregion_leftgroup.hwe.Hobs.windowrolling', header = T)

#add a column called group to indicate which group from the PCA this is to combine data from left,middle,right
#where for this chromo left is inverted homos, middle is heteros, and right is non-inverted homos.

chromo27_left$group <- 'Inv.'

#read in the middle group Hobs data that went through windowscanr.

chromo27_middle <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome27_inversionregion_middlegroup.hwe.Hobs.windowrolling', header = T)

chromo27_middle$group <- 'Heterokaryotypes'

#read in the right group Hobs data that went through windowscanr.

chromo27_right <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome27_inversionregion_rightgroup.hwe.Hobs.windowrolling', header = T)

chromo27_right$group <- 'Non-inv.'

#Bind the 3 dataframes.

Chromo27 <- rbind(chromo27_left, chromo27_middle, chromo27_right)

#Then make the group a factor so it will plot in the order I want along the x axis.

Chromo27$group <- factor(Chromo27$group, levels = c("Inv.", "Heterokaryotypes", "Non-inv."))


#Now try to make a box plot with the 3 different groups
#Note that this will be the mean of means if I am using the windowscanr results.

box_chromo27 <- ggplot(Chromo27, aes(x= group, y = Hobs_mean, fill = group)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(limits = c(0,0.6)) +
  #scale_colour_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  labs(x="Group", y= "Average Hobs", title = "c) Chromosome 27") +
  theme(axis.text.x = element_text(size=19), axis.text.y = element_text(size = 19), axis.title = element_text(size = 19),
        plot.title = element_text(size = 19), axis.title.x = element_blank(),
        panel.background = element_blank())
box_chromo27

#This works but I have a bunch of Hobs of 1 and close to 1 which is because of the window....

ggsave("BTChromo12_AverageHobs_Boxplot.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)



######Chromosome 31####
#Read in the left group Hobs data that went through windowscanr

chromo31_left <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome31_inversionregion_leftgroup.hwe.Hobs.windowrolling', header = T)

#add a column called group to indicate which group from the PCA this is to combine data from left,middle,right
#where for this chromo left is inverted homos, middle is heteros, and right is non-inverted homos.

chromo31_left$group <- 'Inv.'

#read in the middle group Hobs data that went through windowscanr.

chromo31_middle <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome31_inversionregion_middlegroup.hwe.Hobs.windowrolling', header = T)

chromo31_middle$group <- 'Heterokaryotypes'

#read in the right group Hobs data that went through windowscanr.

chromo31_right <- read.table(file='C:/Users/cneme/Documents/Hobs/Chromosome31_inversionregion_rightgroup.hwe.Hobs.windowrolling', header = T)

chromo31_right$group <- 'Non-inv.'

#Bind the 3 dataframes.

Chromo31 <- rbind(chromo31_left, chromo31_middle, chromo31_right)

#Then make the group a factor so it will plot in the order I want along the x axis.

Chromo31$group <- factor(Chromo31$group, levels = c("Inv.", "Heterokaryotypes", "Non-inv."))


#Now try to make a box plot with the 3 different groups
#Note that this will be the mean of means if I am using the windowscanr results.

box_chromo31 <- ggplot(Chromo31, aes(x= group, y = Hobs_mean, fill = group)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(limits = c(0,0.6)) +
  #scale_colour_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  labs(x="Group", y= "Average Hobs", title = "d) Chromosome 31") +
  theme(axis.text.x = element_text(size=19), axis.text.y = element_text(size = 19), axis.title = element_text(size = 19),
        plot.title = element_text(size = 19), axis.title.x = element_blank(),
        panel.background = element_blank())
box_chromo31

#This works but I have a bunch of Hobs of 1 and close to 1 which is because of the window....

ggsave("BTChromo12_AverageHobs_Boxplot.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)


###Facet 4 boxplots together####

box_chromo12 + box_chromo19 + box_chromo27 + box_chromo31 + plot_layout(guides = "collect", ncol = 2)
ggsave("InverChromos_AverageHobs_Boxplot_biggerfont.png", path = "C:/Users/cneme/Documents/Hobs/Plots", bg = "transparent", width = 18, height = 8, units = c("in"), dpi = 300)



####Chromosome 9####
#checking hobs for chromosome 9 as it has a little LD block and the PCA is convincing. Plot raw hobs.

#####Raw right group####
#I think I need a dotplot of the raw data, not the windowed version so that I can pinpoint if my inversion area was correct
#I need the .Hobs file before it was run through windowscanr for this. 

chromo9_right_hobs <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome9_EntireChromo_rightgroup_hwe.hwe.Hobs', header = T)


chromo9_right_hobs <- ggplot(chromo9_right_hobs, aes(Position, Hobs)) +
  geom_point() +
  #xlim(0, 40000000)
  scale_x_continuous(breaks = seq(0,70000000, 2000000)) 
chromo9_right_hobs

#Not convincing that there is a potential inversion
#I do not see a distinct dip in hobs like we do for the other 4 major chromos.



####Chromosome 25####

#####Raw left group####
#I think I need a dotplot of the raw data, not the windowed version so that I can pinpoint if my inversion area was correct
#I need the .Hobs file before it was run through windowscanr for this. 

chromo25left_hobs <- read.table(file='C:/Users/cneme/Documents/Hobs/BTChromosome25_EntireChromo_leftgroup_hwe.hwe.Hobs', header = T)


chromo25left_hobs <- ggplot(chromo25left_hobs, aes(Position, Hobs)) +
  geom_point() +
  #xlim(0, 40000000)
  scale_x_continuous(breaks = seq(0,45000000, 2000000)) 
chromo25left_hobs

#This one is also not convincing...I do not see any real dip in Hobs. 

