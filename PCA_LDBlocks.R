###PCAs for chromosomes with LD blocks###

library(tidyverse) #version 1.3.1
library(ggplot2) #version 3.3.6
library(dplyr) #version 1.0.8
library(RColorBrewer) #version 1.1.3
library(ggrepel) #
library(ggpubr) #version 0.4.0 to do the facet wraps.
library(patchwork) #version 1.1.1 also for facetting

#set the working directory

setwd("C:/Users/cneme/Documents/PCA_LDheatmaps")

#now read in the .covmats

####CH12 CM055694.1 All individuals####
#This .covmat is from running angsd on SNPs within the LD block on chromosome 12

#Read in the data
CH12_cov <- as.matrix(read.table(paste0("AllStreams_GLfromSNPList_CM055694.1_LDBlock.covMat"), header= F))

#Read in file of sample names to add to the plot

#I need to bring in data so I can have a population column and individual ID

#Bring in all the sample merge texts from compute canada and merge those to start

BH.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/BlackHole_samplemerge.txt", header = F)
CV.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/ChurchVault_samplemerge.txt", header = F)
HE.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Healeys_samplemerge.txt", header = F)
PL.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Poole_samplemerge.txt", header = F)
RB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Robinson_samplemerge.txt", header = F)
RC.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/RossCreek_samplemerge.txt", header = F)
SB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SaundersBrookDownstream_samplemerge.txt", header = F)
SS.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SheepShearer_samplemerge.txt", header = F)
WW.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/WoodworthDownstream_samplemerge.txt", header = F)

#Add a column to each sample table to specify the population or stream name
BH.names$population <- "BlackHole"
CV.names$population <- "ChurchVault"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Robinson"
RC.names$population <- "RossCreek"
SB.names$population <- "Saunders"
SS.names$population <- "SheepShearer"
WW.names$population <- "Woodworth"



#bind all tables
AllStreamsNames <- rbind(BH.names,CV.names,HE.names,PL.names,RB.names,RC.names,SB.names,SS.names,WW.names)


#add column title for indivudal id by specifying column 1 in the data frame. 

colnames(AllStreamsNames)[1] <- paste("individual")

#calculate the eigen vectors and values from the covariance matrix.
CH12 <- eigen(CH12_cov)


#extract eigenvectors which are the principal components.
eigenvectors = CH12$vectors 

#This a data frame where each column is labelled X1, X2 and so on. These will be the principle components. So X1 is PC1 etc.
pca.vectors = as_tibble(data.frame(eigenvectors)) 

#This will bind the sample name table adn the eigen vectors together so I can put on sample labels. 
pca.vectors2 = as_tibble(cbind(AllStreamsNames, data.frame(eigenvectors)))

#This will produce objects that give the percent variance for each PC.
#Make these and then just change the axis title to show these percents. 


pca.eigenval.sum = sum(CH12$values) #sum of eigenvalues where the eigen values are the squares of the standard deviation (variance)
varPC1 <- (CH12$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1 =19.5%
varPC2 <- (CH12$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2 =6.0%
varPC3 <- (CH12$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3 =4.3%
varPC4 <- (CH12$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4 =3.4%

#order with factors so can plot from east to west

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("RossCreek", "Woodworth", "BlackHole", "ChurchVault", "Saunders", "Robinson", "SheepShearer", "Healeys", "Poole"))


#Where X1 and X2 are PC1 and PC2
pca1 = ggplot(data = pca.vectors2, aes(x=X1, y=X2)) + 
  geom_point(aes(colour = population), size=1.8, show.legend = FALSE) + #Have to have this legend is false for when individually labelling the points.) +
  geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Poole",]) + #This will make only points with pc1 less than this value show up
  #scale_color_manual(values = c("#00AFBB", "#00AFBB", "#00AFBB","#00AFBB","#00AFBB","#00AFBB","#7F00FF", "#E7B800", "#FC4E07"))+  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Saunders",]) + #This will make only points with pc1 less than this value show up
  #labs(title= "Chromosome 12", x = "PC1 (19.6% explained variance)", y= "PC2 (6.0% explained variance)") +
  theme_gray()
pca1


ggsave("PCA_LDBlock_Chromosome12_Allindividuals.png", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)


####CH19 CM055701.1 All individuals####

#Read in the data
CH19_cov <- as.matrix(read.table(paste0("AllStreams_GLfromSNPList_CM055701.1_LDBlock.covMat"), header= F))

#Read in file of sample names to add to the plot

#I need to bring in data so I can have a population column and individual ID

#Bring in all the sample merge texts from compute canada and merge those to start

BH.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/BlackHole_samplemerge.txt", header = F)
CV.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/ChurchVault_samplemerge.txt", header = F)
HE.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Healeys_samplemerge.txt", header = F)
PL.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Poole_samplemerge.txt", header = F)
RB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Robinson_samplemerge.txt", header = F)
RC.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/RossCreek_samplemerge.txt", header = F)
SB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SaundersBrookDownstream_samplemerge.txt", header = F)
SS.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SheepShearer_samplemerge.txt", header = F)
WW.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/WoodworthDownstream_samplemerge.txt", header = F)

#Add a column to each sample table to specify the population or stream name
BH.names$population <- "BlackHole"
CV.names$population <- "ChurchVault"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Robinson"
RC.names$population <- "RossCreek"
SB.names$population <- "Saunders"
SS.names$population <- "SheepShearer"
WW.names$population <- "Woodworth"



#bind all tables
AllStreamsNames <- rbind(BH.names,CV.names,HE.names,PL.names,RB.names,RC.names,SB.names,SS.names,WW.names)


#add column title for indivudal id by specifying column 1 in the data frame. 

colnames(AllStreamsNames)[1] <- paste("individual")

#calculate the eigen vectors and values from the covariance matrix.
CH19 <- eigen(CH19_cov)


#extract eigenvectors which are the principal components.
eigenvectors = CH19$vectors 

#This a data frame where each column is labelled X1, X2 and so on. These will be the principle components. So X1 is PC1 etc.
pca.vectors = as_tibble(data.frame(eigenvectors)) 

#This will bind the sample name table adn the eigen vectors together so I can put on sample labels. 
pca.vectors2 = as_tibble(cbind(AllStreamsNames, data.frame(eigenvectors)))

#This will produce objects that give the percent variance for each PC.
#Make these and then just change the axis title to show these percents. 


pca.eigenval.sum = sum(CH19$values) #sum of eigenvalues where the eigen values are the squares of the standard deviation (variance)
varPC1 <- (CH19$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1 =18.5%
varPC2 <- (CH19$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2 =6.5%
varPC3 <- (CH19$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3 =4.9%
varPC4 <- (CH19$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4 =3.0%

#order with factors so can plot from east to west

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("RossCreek", "Woodworth", "BlackHole", "ChurchVault", "Saunders", "Robinson", "SheepShearer", "Healeys", "Poole"))


#Where X1 and X2 are PC1 and PC2
pca2 = ggplot(data = pca.vectors2, aes(x=X1, y=X2)) + 
  geom_point(aes(colour = population), size=1.8, show.legend = FALSE) + #Have to have this legend is false for when individually labelling the points.) +
  geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 30, data = pca.vectors2[pca.vectors2$population %in% "SheepShearer",]) + #This will make only points with pc1 less than this value show up
  #scale_color_manual(values = c("#00AFBB", "#00AFBB", "#00AFBB","#00AFBB","#00AFBB","#00AFBB","#7F00FF", "#E7B800", "#FC4E07"))+  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Saunders",]) + #This will make only points with pc1 less than this value show up
  labs(title= "Chromosome 19", x = "PC1 (18.5% explained variance)", y= "PC2 (6.5% explained variance)") +
  theme_gray()
pca2


ggsave("PCA_LDBlock_Chromosome19_Allindividuals.png", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)




####CH27 CM055709.1 All individuals####

#Read in the data
CH27_cov <- as.matrix(read.table(paste0("AllStreams_GLfromSNPList_CM055709.1_LDBlock.covMat"), header= F))

#Read in file of sample names to add to the plot

#I need to bring in data so I can have a population column and individual ID

#Bring in all the sample merge texts from compute canada and merge those to start

BH.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/BlackHole_samplemerge.txt", header = F)
CV.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/ChurchVault_samplemerge.txt", header = F)
HE.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Healeys_samplemerge.txt", header = F)
PL.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Poole_samplemerge.txt", header = F)
RB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Robinson_samplemerge.txt", header = F)
RC.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/RossCreek_samplemerge.txt", header = F)
SB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SaundersBrookDownstream_samplemerge.txt", header = F)
SS.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SheepShearer_samplemerge.txt", header = F)
WW.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/WoodworthDownstream_samplemerge.txt", header = F)

#Add a column to each sample table to specify the population or stream name
BH.names$population <- "BlackHole"
CV.names$population <- "ChurchVault"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Robinson"
RC.names$population <- "RossCreek"
SB.names$population <- "Saunders"
SS.names$population <- "SheepShearer"
WW.names$population <- "Woodworth"



#bind all tables
AllStreamsNames <- rbind(BH.names,CV.names,HE.names,PL.names,RB.names,RC.names,SB.names,SS.names,WW.names)


#add column title for indivudal id by specifying column 1 in the data frame. 

colnames(AllStreamsNames)[1] <- paste("individual")

#calculate the eigen vectors and values from the covariance matrix.
CH27 <- eigen(CH27_cov)


#extract eigenvectors which are the principal components.
eigenvectors = CH27$vectors 

#This a data frame where each column is labelled X1, X2 and so on. These will be the principle components. So X1 is PC1 etc.
pca.vectors = as_tibble(data.frame(eigenvectors)) 

#This will bind the sample name table adn the eigen vectors together so I can put on sample labels. 
pca.vectors2 = as_tibble(cbind(AllStreamsNames, data.frame(eigenvectors)))

#This will produce objects that give the percent variance for each PC.
#Make these and then just change the axis title to show these percents. 


pca.eigenval.sum = sum(CH27$values) #sum of eigenvalues where the eigen values are the squares of the standard deviation (variance)
varPC1 <- (CH27$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1 =20.9%
varPC2 <- (CH27$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2 =6.2%
varPC3 <- (CH27$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3 =5.4%
varPC4 <- (CH27$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4 =4.0%

#order with factors so can plot from east to west

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("RossCreek", "Woodworth", "BlackHole", "ChurchVault", "Saunders", "Robinson", "SheepShearer", "Healeys", "Poole"))


#Where X1 and X2 are PC1 and PC2
pca3 = ggplot(data = pca.vectors2, aes(x=X1, y=X2)) + 
  geom_point(aes(colour = population), size=1.8, show.legend = FALSE) + #Have to have this legend is false for when individually labelling the points.) +
  geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 30, data = pca.vectors2[pca.vectors2$population %in% "Poole",]) + #This will make only points with pc1 less than this value show up
  #scale_color_manual(values = c("#00AFBB", "#00AFBB", "#00AFBB","#00AFBB","#00AFBB","#00AFBB","#7F00FF", "#E7B800", "#FC4E07"))+  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Saunders",]) + #This will make only points with pc1 less than this value show up
  labs(title= "Chromosome 27", x = "PC1 (20.9% explained variance)", y= "PC2 (6.2% explained variance)") +
  theme_gray()
pca3


ggsave("PCA_LDBlock_Chromosome27_Allindividuals.png", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)



####CH31 CM055713.1 All individuals####

#Read in the data
CH31_cov <- as.matrix(read.table(paste0("AllStreams_GLfromSNPList_CM055713.1_LDBlock.covMat"), header= F))

#Read in file of sample names to add to the plot

#I need to bring in data so I can have a population column and individual ID

#Bring in all the sample merge texts from compute canada and merge those to start

BH.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/BlackHole_samplemerge.txt", header = F)
CV.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/ChurchVault_samplemerge.txt", header = F)
HE.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Healeys_samplemerge.txt", header = F)
PL.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Poole_samplemerge.txt", header = F)
RB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Robinson_samplemerge.txt", header = F)
RC.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/RossCreek_samplemerge.txt", header = F)
SB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SaundersBrookDownstream_samplemerge.txt", header = F)
SS.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SheepShearer_samplemerge.txt", header = F)
WW.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/WoodworthDownstream_samplemerge.txt", header = F)

#Add a column to each sample table to specify the population or stream name
BH.names$population <- "BlackHole"
CV.names$population <- "ChurchVault"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Robinson"
RC.names$population <- "RossCreek"
SB.names$population <- "Saunders"
SS.names$population <- "SheepShearer"
WW.names$population <- "Woodworth"



#bind all tables
AllStreamsNames <- rbind(BH.names,CV.names,HE.names,PL.names,RB.names,RC.names,SB.names,SS.names,WW.names)


#add column title for indivudal id by specifying column 1 in the data frame. 

colnames(AllStreamsNames)[1] <- paste("individual")

#calculate the eigen vectors and values from the covariance matrix.
CH31 <- eigen(CH31_cov)


#extract eigenvectors which are the principal components.
eigenvectors = CH31$vectors 

#This a data frame where each column is labelled X1, X2 and so on. These will be the principle components. So X1 is PC1 etc.
pca.vectors = as_tibble(data.frame(eigenvectors)) 

#This will bind the sample name table adn the eigen vectors together so I can put on sample labels. 
pca.vectors2 = as_tibble(cbind(AllStreamsNames, data.frame(eigenvectors)))

#This will produce objects that give the percent variance for each PC.
#Make these and then just change the axis title to show these percents. 


pca.eigenval.sum = sum(CH31$values) #sum of eigenvalues where the eigen values are the squares of the standard deviation (variance)
varPC1 <- (CH31$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1 =17.3%
varPC2 <- (CH31$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2 =6.4%
varPC3 <- (CH31$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3 =4.5%
varPC4 <- (CH31$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4 =3.8%

#order with factors so can plot from east to west

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("RossCreek", "Woodworth", "BlackHole", "ChurchVault", "Saunders", "Robinson", "SheepShearer", "Healeys", "Poole"))


#Where X1 and X2 are PC1 and PC2
pca4 = ggplot(data = pca.vectors2, aes(x=X1, y=X2)) + 
  geom_point(aes(colour = population), size=1.8, show.legend = FALSE) + #Have to have this legend is false for when individually labelling the points.) +
  geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 30, data = pca.vectors2[pca.vectors2$population %in% "Woodworth",]) + #This will make only points with pc1 less than this value show up
  #scale_color_manual(values = c("#00AFBB", "#00AFBB", "#00AFBB","#00AFBB","#00AFBB","#00AFBB","#7F00FF", "#E7B800", "#FC4E07"))+  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Saunders",]) + #This will make only points with pc1 less than this value show up
  labs(title= "Chromosome 31", x = "PC1 (17.3% explained variance)", y= "PC2 (6.4% explained variance)") +
  theme_gray()
pca4


ggsave("PCA_LDBlock_Chromosome31_Allindividuals.png", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

#####Group PCA####

figure <- ggarrange(pca1, pca2, pca3, pca4,
                    labels = c("A", "B", "C", "D", "E"),
                    common.legend = TRUE, legend = "right")
figure

ggsave("Facet_AllStreams_PCA_LDblockSectionsOFChromoms_March11_2024.png", path = "C:/Users/cneme/Documents/PCA_LDheatmaps", bg = "white", width = 12, height = 8, units = c("in"), dpi = 300)


####Same plot with adjustments to size of symbols and background/legend.####


#now read in the .covmats

####CH12 CM055694.1 All individuals####
#This .covmat is from running angsd on SNPs within the LD block on chromosome 12

#Read in the data
CH12_cov <- as.matrix(read.table(paste0("AllStreams_GLfromSNPList_CM055694.1_LDBlock.covMat"), header= F))

#Read in file of sample names to add to the plot

#I need to bring in data so I can have a population column and individual ID

#Bring in all the sample merge texts from compute canada and merge those to start

BH.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/BlackHole_samplemerge.txt", header = F)
CV.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/ChurchVault_samplemerge.txt", header = F)
HE.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Healeys_samplemerge.txt", header = F)
PL.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Poole_samplemerge.txt", header = F)
RB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Robinson_samplemerge.txt", header = F)
RC.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/RossCreek_samplemerge.txt", header = F)
SB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SaundersBrookDownstream_samplemerge.txt", header = F)
SS.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SheepShearer_samplemerge.txt", header = F)
WW.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/WoodworthDownstream_samplemerge.txt", header = F)

#Add a column to each sample table to specify the population or stream name
BH.names$population <- "Easternmost Populations"
CV.names$population <- "Easternmost Populations"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Easternmost Populations"
RC.names$population <- "Easternmost Populations"
SB.names$population <- "Easternmost Populations"
SS.names$population <- "Sheep Shearer"
WW.names$population <- "Easternmost Populations"



#bind all tables
AllStreamsNames <- rbind(BH.names,CV.names,HE.names,PL.names,RB.names,RC.names,SB.names,SS.names,WW.names)


#add column title for indivudal id by specifying column 1 in the data frame. 

colnames(AllStreamsNames)[1] <- paste("individual")

#calculate the eigen vectors and values from the covariance matrix.
CH12 <- eigen(CH12_cov)


#extract eigenvectors which are the principal components.
eigenvectors = CH12$vectors 

#This a data frame where each column is labelled X1, X2 and so on. These will be the principle components. So X1 is PC1 etc.
pca.vectors = as_tibble(data.frame(eigenvectors)) 

#This will bind the sample name table adn the eigen vectors together so I can put on sample labels. 
pca.vectors2 = as_tibble(cbind(AllStreamsNames, data.frame(eigenvectors)))

#This will produce objects that give the percent variance for each PC.
#Make these and then just change the axis title to show these percents. 


pca.eigenval.sum = sum(CH12$values) #sum of eigenvalues where the eigen values are the squares of the standard deviation (variance)
varPC1 <- (CH12$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1 =19.5%
varPC2 <- (CH12$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2 =6.0%
varPC3 <- (CH12$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3 =4.3%
varPC4 <- (CH12$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4 =3.4%

#order with factors so can plot from east to west

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("Easternmost Populations", "Sheep Shearer", "Healeys", "Poole"))


#Where X1 and X2 are PC1 and PC2
pca1 = ggplot(data = pca.vectors2, aes(x=X1, y=X2)) + 
  geom_point(aes(colour = population), size=3) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_point(aes(colour = population), size=2.3, show.legend = FALSE) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Poole",]) + #This will make only points with pc1 less than this value show up
  scale_color_manual(values = c("#00AFBB","#7F00FF", "#E7B800", "#FC4E07"))+  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Saunders",]) + #This will make only points with pc1 less than this value show up
  labs(x = "PC1", y= "PC2") +
  ggtitle("a) Chromosome 12") +
  #labs(title= "Chromosome 12", x = "PC1 (19.6% explained variance)", y= "PC2 (6.0% explained variance)") +
  theme(panel.background = element_blank(), plot.title = element_text(size = 16), #hjust 1 is left 0 is right
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        axis.text = element_text(size = 14))
pca1


ggsave("PCA_LDBlock_Chromosome12_Allindividuals.png", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)


####CH19 CM055701.1 All individuals####

#Read in the data
CH19_cov <- as.matrix(read.table(paste0("AllStreams_GLfromSNPList_CM055701.1_LDBlock.covMat"), header= F))

#Read in file of sample names to add to the plot

#I need to bring in data so I can have a population column and individual ID

#Bring in all the sample merge texts from compute canada and merge those to start

BH.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/BlackHole_samplemerge.txt", header = F)
CV.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/ChurchVault_samplemerge.txt", header = F)
HE.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Healeys_samplemerge.txt", header = F)
PL.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Poole_samplemerge.txt", header = F)
RB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Robinson_samplemerge.txt", header = F)
RC.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/RossCreek_samplemerge.txt", header = F)
SB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SaundersBrookDownstream_samplemerge.txt", header = F)
SS.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SheepShearer_samplemerge.txt", header = F)
WW.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/WoodworthDownstream_samplemerge.txt", header = F)

#Add a column to each sample table to specify the population or stream name
BH.names$population <- "Easternmost Populations"
CV.names$population <- "Easternmost Populations"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Easternmost Populations"
RC.names$population <- "Easternmost Populations"
SB.names$population <- "Easternmost Populations"
SS.names$population <- "Sheep Shearer"
WW.names$population <- "Easternmost Populations"



#bind all tables
AllStreamsNames <- rbind(BH.names,CV.names,HE.names,PL.names,RB.names,RC.names,SB.names,SS.names,WW.names)


#add column title for indivudal id by specifying column 1 in the data frame. 

colnames(AllStreamsNames)[1] <- paste("individual")

#calculate the eigen vectors and values from the covariance matrix.
CH19 <- eigen(CH19_cov)


#extract eigenvectors which are the principal components.
eigenvectors = CH19$vectors 

#This a data frame where each column is labelled X1, X2 and so on. These will be the principle components. So X1 is PC1 etc.
pca.vectors = as_tibble(data.frame(eigenvectors)) 

#This will bind the sample name table adn the eigen vectors together so I can put on sample labels. 
pca.vectors2 = as_tibble(cbind(AllStreamsNames, data.frame(eigenvectors)))

#This will produce objects that give the percent variance for each PC.
#Make these and then just change the axis title to show these percents. 


pca.eigenval.sum = sum(CH19$values) #sum of eigenvalues where the eigen values are the squares of the standard deviation (variance)
varPC1 <- (CH19$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1 =18.5%
varPC2 <- (CH19$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2 =6.5%
varPC3 <- (CH19$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3 =4.9%
varPC4 <- (CH19$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4 =3.0%

#order with factors so can plot from east to west

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("Easternmost Populations", "Sheep Shearer", "Healeys", "Poole"))


#Where X1 and X2 are PC1 and PC2
pca2 = ggplot(data = pca.vectors2, aes(x=X1, y=X2)) + 
  geom_point(aes(colour = population), size=3) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_point(aes(colour = population), size=1.8, show.legend = FALSE) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 30, data = pca.vectors2[pca.vectors2$population %in% "SheepShearer",]) + #This will make only points with pc1 less than this value show up
  scale_color_manual(values = c("#00AFBB","#7F00FF", "#E7B800", "#FC4E07"))+  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Saunders",]) + #This will make only points with pc1 less than this value show up
  labs(title= "b) Chromosome 19", x = "PC1", y= "PC2") +
  theme(panel.background = element_blank(), plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        axis.text = element_text(size = 14))
pca2


ggsave("PCA_LDBlock_Chromosome19_Allindividuals.png", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)




####CH27 CM055709.1 All individuals####

#Read in the data
CH27_cov <- as.matrix(read.table(paste0("AllStreams_GLfromSNPList_CM055709.1_LDBlock.covMat"), header= F))

#Read in file of sample names to add to the plot

#I need to bring in data so I can have a population column and individual ID

#Bring in all the sample merge texts from compute canada and merge those to start

BH.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/BlackHole_samplemerge.txt", header = F)
CV.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/ChurchVault_samplemerge.txt", header = F)
HE.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Healeys_samplemerge.txt", header = F)
PL.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Poole_samplemerge.txt", header = F)
RB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Robinson_samplemerge.txt", header = F)
RC.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/RossCreek_samplemerge.txt", header = F)
SB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SaundersBrookDownstream_samplemerge.txt", header = F)
SS.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SheepShearer_samplemerge.txt", header = F)
WW.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/WoodworthDownstream_samplemerge.txt", header = F)

#Add a column to each sample table to specify the population or stream name
BH.names$population <- "Easternmost Populations"
CV.names$population <- "Easternmost Populations"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Easternmost Populations"
RC.names$population <- "Easternmost Populations"
SB.names$population <- "Easternmost Populations"
SS.names$population <- "Sheep Shearer"
WW.names$population <- "Easternmost Populations"



#bind all tables
AllStreamsNames <- rbind(BH.names,CV.names,HE.names,PL.names,RB.names,RC.names,SB.names,SS.names,WW.names)


#add column title for indivudal id by specifying column 1 in the data frame. 

colnames(AllStreamsNames)[1] <- paste("individual")

#calculate the eigen vectors and values from the covariance matrix.
CH27 <- eigen(CH27_cov)


#extract eigenvectors which are the principal components.
eigenvectors = CH27$vectors 

#This a data frame where each column is labelled X1, X2 and so on. These will be the principle components. So X1 is PC1 etc.
pca.vectors = as_tibble(data.frame(eigenvectors)) 

#This will bind the sample name table adn the eigen vectors together so I can put on sample labels. 
pca.vectors2 = as_tibble(cbind(AllStreamsNames, data.frame(eigenvectors)))

#This will produce objects that give the percent variance for each PC.
#Make these and then just change the axis title to show these percents. 


pca.eigenval.sum = sum(CH27$values) #sum of eigenvalues where the eigen values are the squares of the standard deviation (variance)
varPC1 <- (CH27$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1 =20.9%
varPC2 <- (CH27$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2 =6.2%
varPC3 <- (CH27$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3 =5.4%
varPC4 <- (CH27$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4 =4.0%

#order with factors so can plot from east to west

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("Easternmost Populations", "Sheep Shearer", "Healeys", "Poole"))


#Where X1 and X2 are PC1 and PC2
pca3 = ggplot(data = pca.vectors2, aes(x=X1, y=X2)) + 
  geom_point(aes(colour = population), size=3) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_point(aes(colour = population), size=1.8, show.legend = FALSE) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 30, data = pca.vectors2[pca.vectors2$population %in% "Poole",]) + #This will make only points with pc1 less than this value show up
  scale_color_manual(values = c("#00AFBB","#7F00FF", "#E7B800", "#FC4E07"))+  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Saunders",]) + #This will make only points with pc1 less than this value show up
  labs(title= "c) Chromosome 27", x = "PC1", y= "PC2") +
  theme(panel.background = element_blank(), plot.title = element_text(size = 16),
             axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
             legend.text = element_text(size = 16),
             legend.title = element_blank(),
             axis.text = element_text(size = 14))
pca3


ggsave("PCA_LDBlock_Chromosome27_Allindividuals.png", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)



####CH31 CM055713.1 All individuals####

#Read in the data
CH31_cov <- as.matrix(read.table(paste0("AllStreams_GLfromSNPList_CM055713.1_LDBlock.covMat"), header= F))

#Read in file of sample names to add to the plot

#I need to bring in data so I can have a population column and individual ID

#Bring in all the sample merge texts from compute canada and merge those to start

BH.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/BlackHole_samplemerge.txt", header = F)
CV.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/ChurchVault_samplemerge.txt", header = F)
HE.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Healeys_samplemerge.txt", header = F)
PL.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Poole_samplemerge.txt", header = F)
RB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Robinson_samplemerge.txt", header = F)
RC.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/RossCreek_samplemerge.txt", header = F)
SB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SaundersBrookDownstream_samplemerge.txt", header = F)
SS.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SheepShearer_samplemerge.txt", header = F)
WW.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/WoodworthDownstream_samplemerge.txt", header = F)

#Add a column to each sample table to specify the population or stream name
BH.names$population <- "Easternmost Populations"
CV.names$population <- "Easternmost Populations"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Easternmost Populations"
RC.names$population <- "Easternmost Populations"
SB.names$population <- "Easternmost Populations"
SS.names$population <- "Sheep Shearer"
WW.names$population <- "Easternmost Populations"



#bind all tables
AllStreamsNames <- rbind(BH.names,CV.names,HE.names,PL.names,RB.names,RC.names,SB.names,SS.names,WW.names)


#add column title for indivudal id by specifying column 1 in the data frame. 

colnames(AllStreamsNames)[1] <- paste("individual")

#calculate the eigen vectors and values from the covariance matrix.
CH31 <- eigen(CH31_cov)


#extract eigenvectors which are the principal components.
eigenvectors = CH31$vectors 

#This a data frame where each column is labelled X1, X2 and so on. These will be the principle components. So X1 is PC1 etc.
pca.vectors = as_tibble(data.frame(eigenvectors)) 

#This will bind the sample name table adn the eigen vectors together so I can put on sample labels. 
pca.vectors2 = as_tibble(cbind(AllStreamsNames, data.frame(eigenvectors)))

#This will produce objects that give the percent variance for each PC.
#Make these and then just change the axis title to show these percents. 


pca.eigenval.sum = sum(CH31$values) #sum of eigenvalues where the eigen values are the squares of the standard deviation (variance)
varPC1 <- (CH31$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1 =17.3%
varPC2 <- (CH31$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2 =6.4%
varPC3 <- (CH31$values[3]/pca.eigenval.sum)*100 #Variance explained by PC3 =4.5%
varPC4 <- (CH31$values[4]/pca.eigenval.sum)*100 #Variance explained by PC4 =3.8%

#order with factors so can plot from east to west

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("Easternmost Populations", "Sheep Shearer", "Healeys", "Poole"))


#Where X1 and X2 are PC1 and PC2
pca4 = ggplot(data = pca.vectors2, aes(x=X1, y=X2)) + 
  geom_point(aes(colour = population), size=3) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_point(aes(colour = population), size=1.8, show.legend = FALSE) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 30, data = pca.vectors2[pca.vectors2$population %in% "Woodworth",]) + #This will make only points with pc1 less than this value show up
  scale_color_manual(values = c("#00AFBB","#7F00FF", "#E7B800", "#FC4E07"))+  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Saunders",]) + #This will make only points with pc1 less than this value show up
  labs(title= "d) Chromosome 31", x = "PC1", y= "PC2") +
  theme(panel.background = element_blank(), plot.title = element_text(size = 16),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        axis.text = element_text(size = 14))
pca4


ggsave("PCA_LDBlock_Chromosome31_Allindividuals.png", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

#####Group PCA####

figure <- ggarrange(pca1, pca2, pca3, pca4,
                    common.legend = TRUE, legend = "bottom")
figure

ggsave("Facet_AllStreams_PCA_LDblockSectionsOFChromoms_April18_2024.png", path = "C:/Users/cneme/Documents/PCA_LDheatmaps", bg = "white", width = 12, height = 8, units = c("in"), dpi = 300)


##Testing chromosome 9 and 25 for 3 potential groups based on small LD blocks####

##Chromosome 9 CM055691.1####
#This .covmat is from running angsd on SNPs within the LD block on chromosome 12

#Read in the data
CH9_cov <- as.matrix(read.table(paste0("AllStreams_GLfromSNPList_CM055691.1_LDBlock.covMat"), header= F))

#Read in file of sample names to add to the plot


#Read in file of sample names to add to the plot

#I need to bring in data so I can have a population column and individual ID

#Bring in all the sample merge texts from compute canada and merge those to start

BH.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/BlackHole_samplemerge.txt", header = F)
CV.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/ChurchVault_samplemerge.txt", header = F)
HE.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Healeys_samplemerge.txt", header = F)
PL.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Poole_samplemerge.txt", header = F)
RB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Robinson_samplemerge.txt", header = F)
RC.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/RossCreek_samplemerge.txt", header = F)
SB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SaundersBrookDownstream_samplemerge.txt", header = F)
SS.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SheepShearer_samplemerge.txt", header = F)
WW.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/WoodworthDownstream_samplemerge.txt", header = F)

#Add a column to each sample table to specify the population or stream name
BH.names$population <- "Easternmost Populations"
CV.names$population <- "Easternmost Populations"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Easternmost Populations"
RC.names$population <- "Easternmost Populations"
SB.names$population <- "Easternmost Populations"
SS.names$population <- "Sheep Shearer"
WW.names$population <- "Easternmost Populations"



#bind all tables
AllStreamsNames <- rbind(BH.names,CV.names,HE.names,PL.names,RB.names,RC.names,SB.names,SS.names,WW.names)


#add column title for indivudal id by specifying column 1 in the data frame. 

colnames(AllStreamsNames)[1] <- paste("individual")

#calculate the eigen vectors and values from the covariance matrix.
CH9 <- eigen(CH9_cov)


#extract eigenvectors which are the principal components.
eigenvectors = CH9$vectors 

#This a data frame where each column is labelled X1, X2 and so on. These will be the principle components. So X1 is PC1 etc.
pca.vectors = as_tibble(data.frame(eigenvectors)) 

#This will bind the sample name table adn the eigen vectors together so I can put on sample labels. 
pca.vectors2 = as_tibble(cbind(AllStreamsNames, data.frame(eigenvectors)))

#This will produce objects that give the percent variance for each PC.
#Make these and then just change the axis title to show these percents. 


pca.eigenval.sum = sum(CH9$values) #sum of eigenvalues where the eigen values are the squares of the standard deviation (variance)
varPC1 <- (CH9$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1 =21.9%
varPC2 <- (CH9$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2 =7.2%


#order with factors so can plot from east to west

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("Easternmost Populations", "Sheep Shearer", "Healeys", "Poole"))


#Where X1 and X2 are PC1 and PC2
pca9 = ggplot(data = pca.vectors2, aes(x=X1, y=X2)) + 
  geom_point(aes(colour = population), size=3) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_point(aes(colour = population), size=2.3, show.legend = FALSE) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Sheep Shearer",]) + #This will make only points with pc1 less than this value show up
  scale_color_manual(values = c("#00AFBB","#7F00FF", "#E7B800", "#FC4E07"))+  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Saunders",]) + #This will make only points with pc1 less than this value show up
  labs(x = "PC1", y= "PC2") +
  ggtitle("a) Chromosome 9") +
  #labs(title= "Chromosome 12", x = "PC1 (19.6% explained variance)", y= "PC2 (6.0% explained variance)") +
  theme(panel.background = element_blank(), plot.title = element_text(size = 16), #hjust 1 is left 0 is right
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        axis.text = element_text(size = 14))
pca9

ggsave("PCA_LDBlock_Chromosome9_Allindividuals.png", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)


##Chromosome 25 CM055707.1####
#This .covmat is from running angsd on SNPs within the LD block on chromosome 12

#Read in the data
CH25_cov <- as.matrix(read.table(paste0("AllStreams_GLfromSNPList_CM055707.1_LDBlock.covMat"), header= F))

#Read in file of sample names to add to the plot

#I need to bring in data so I can have a population column and individual ID

#Bring in all the sample merge texts from compute canada and merge those to start

BH.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/BlackHole_samplemerge.txt", header = F)
CV.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/ChurchVault_samplemerge.txt", header = F)
HE.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Healeys_samplemerge.txt", header = F)
PL.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Poole_samplemerge.txt", header = F)
RB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/Robinson_samplemerge.txt", header = F)
RC.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/RossCreek_samplemerge.txt", header = F)
SB.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SaundersBrookDownstream_samplemerge.txt", header = F)
SS.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/SheepShearer_samplemerge.txt", header = F)
WW.names <- read.table("C:/Users/cneme/Documents/PCA_LDheatmaps/WoodworthDownstream_samplemerge.txt", header = F)

#Add a column to each sample table to specify the population or stream name
BH.names$population <- "Easternmost Populations"
CV.names$population <- "Easternmost Populations"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Easternmost Populations"
RC.names$population <- "Easternmost Populations"
SB.names$population <- "Easternmost Populations"
SS.names$population <- "Sheep Shearer"
WW.names$population <- "Easternmost Populations"

#to figure out which individual is the outlier in the right group
BH.names$population <- "BlackHole"
CV.names$population <- "ChurchVault"
HE.names$population <- "Healeys"
PL.names$population <- "Poole"
RB.names$population <- "Robinson"
RC.names$population <- "RossCreek"
SB.names$population <- "Saunders"
SS.names$population <- "SheepShearer"
WW.names$population <- "Woodworth"



#bind all tables
AllStreamsNames <- rbind(BH.names,CV.names,HE.names,PL.names,RB.names,RC.names,SB.names,SS.names,WW.names)


#add column title for indivudal id by specifying column 1 in the data frame. 

colnames(AllStreamsNames)[1] <- paste("individual")

#calculate the eigen vectors and values from the covariance matrix.
CH25 <- eigen(CH25_cov)


#extract eigenvectors which are the principal components.
eigenvectors = CH25$vectors 

#This a data frame where each column is labelled X1, X2 and so on. These will be the principle components. So X1 is PC1 etc.
pca.vectors = as_tibble(data.frame(eigenvectors)) 

#This will bind the sample name table adn the eigen vectors together so I can put on sample labels. 
pca.vectors2 = as_tibble(cbind(AllStreamsNames, data.frame(eigenvectors)))

#This will produce objects that give the percent variance for each PC.
#Make these and then just change the axis title to show these percents. 


pca.eigenval.sum = sum(CH25$values) #sum of eigenvalues where the eigen values are the squares of the standard deviation (variance)
varPC1 <- (CH25$values[1]/pca.eigenval.sum)*100 #Variance explained by PC1 =19.7%
varPC2 <- (CH25$values[2]/pca.eigenval.sum)*100 #Variance explained by PC2 =8.4%


#order with factors so can plot from east to west

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("Easternmost Populations", "Sheep Shearer", "Healeys", "Poole"))


#for outlier

pca.vectors2$population <- factor(pca.vectors2$population, levels = c("RossCreek", "Woodworth", "BlackHole", "ChurchVault", "Saunders", "Robinson", "SheepShearer", "Healeys", "Poole"))


#Where X1 and X2 are PC1 and PC2
pca10 = ggplot(data = pca.vectors2, aes(x=X1, y=X2)) + 
  geom_point(aes(colour = population), size=3) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_point(aes(colour = population), size=2.3, show.legend = FALSE) + #Have to have this legend is false for when individually labelling the points.) +
  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Healeys",]) + #This will make only points with pc1 less than this value show up
  scale_color_manual(values = c("#00AFBB","#7F00FF", "#E7B800", "#FC4E07"))+  #geom_text_repel(aes(x= X1, y= X2, label = individual), size= 3, max.overlaps = 21, data = pca.vectors2[pca.vectors2$population %in% "Saunders",]) + #This will make only points with pc1 less than this value show up
  labs(x = "PC1", y= "PC2") +
  ggtitle("a) Chromosome 25") +
  #labs(title= "Chromosome 12", x = "PC1 (19.6% explained variance)", y= "PC2 (6.0% explained variance)") +
  theme(panel.background = element_blank(), plot.title = element_text(size = 16), #hjust 1 is left 0 is right
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        axis.text = element_text(size = 14))
pca10

ggsave("PCA_LDBlock_Chromosome25_Allindividuals.png", bg = "transparent", width = 12, height = 8, units = c("in"), dpi = 300)

