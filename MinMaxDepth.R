####Plotting the .depthGlobal files ####
#This is to then filter and figure out what to set the min and max depth to for further ananalyses in ANGSD.

#load libraries

library(tidyr)
library(readr)
library(data.table)
library(ggplot2)
library(na.tools)
library(dplyr)



#Read in the .depthGlobal file produced by ANGSD.#

RC <- read_tsv("/home/danielruzzante/BrookTroutlcWGS_2021/ComputeCanadaData/BTGenome/RossCreekDepthCount.depthGlobal", col_names = F)


#The way this reads in is that there is 212 columns representing the max depth I set in the ANGSD script. So first column is 0 depth, 2nd is 1x and so on
#I need this to be one column instead of many. The transpose function will  work for this.

test <- transpose(RC)

#Change column V1 to number of sites

colnames(test)[colnames(test) == "V1"] <- "NumberOfSites"


#Now it is all in one column though I might now need a column for depth indicating depth 0-212

test3 <- as.data.frame(seq(0,211))
colnames(test3)[colnames(test3) == "seq(0, 211)"] <- "depth"


#Now bind it to test

RCbind <- cbind(test, test3)


#Plot
#After plotting find the mean and +1 -1 sd. Don't want the ones around 1 where this is a tone of reads.

plot <- ggplot(data = RCbind, aes(x=depth, y=NumberOfSites)) +
  scale_x_continuous(breaks = seq(0,200, by=10)) +
  geom_freqpoly(stat = "identity")
plot


#Now that I have the plot. Have to run some filters and then get mean and standard deviation so they can be put on the plot
#The filter that is set below depends on the original plot. In this case I will set at about 35 and 90 to capture the normal distribution.
#Don't want the ones with really high depth because they are duplicates and the ones at the beginning have 0 or near 0 depth.
#Calculate mean and standard deviation of the RCbind data frame

#Filter
#First get rid of the na for 212 that has no data.
RCbind2 <- na.omit(RCbind)

#For the filter function- since my data is all numeric, have to supply the column that I want to filter by and use logical filter conditions like > <=.

RCbind2 <- filter(RCbind2, depth > 35, depth <= 90)


#Mean and SD

meanRC <- mean(RCbind2$depth)
sdRC <- sqrt(var(RCbind2$depth))

#Now +1 and -1 standard deviation which will be maximum and minimum filters for ANGSD.

max_depth_RC <- meanRC + sdRC
min_depth_RC <- meanRC - sdRC

#Now plot with the means and sd +1 -1

plot2 <- ggplot(data = RCbind2, aes(x=depth, y=NumberOfSites)) +
  geom_freqpoly(stat = "identity") +
  geom_vline(xintercept = c(min_depth_RC, max_depth_RC), color = "red")
plot2

#Focusing on the normal distribution and putting the red lines for the max and min depth that I should use
#where max and min are now numbers in the values section under the global environment on the right
#For Ross Creek min depth should be 46 and max should be 79.
