######Chromosome 9#####

#Read in the gene table for chromosome 9 inverted region of 22Mb-42mb.

CH9_Genes <- read.table(file = 'C:/Users/cneme/Downloads/CHROMO_9_Genes', header = F)

#Add headers manually based on the describe table schema from table brower of USSC genome browser.

colnames(CH9_Genes) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes', 'chromStarts',
                         'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames', 'geneName', 'geneName2')

Chrom9_GeneList <- sort(unique(CH9_Genes$geneName2))

#Turn list into dataframe

List <- as.data.frame(Chrom9_GeneList)

#Maybe filter out the loci that aren't known genes. 
#Easiest way for now is to look which rows have LOC and just select those out. 

List2 <- as.data.frame(List[-(31:286),])

write.table(List2, file="Chromo9_Genes.txt",quote=FALSE, row.names=FALSE, sep='\t')



######Chromosome 14####

#Read in the gene table for chromosome 9 inverted region of 22Mb-42mb.

CH14_Genes <- read.table(file = 'C:/Users/cneme/Downloads/CHROMO_14_Genes', header = F)

#Add headers manually based on the describe table schema from table brower of USSC genome browser.

colnames(CH14_Genes) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes', 'chromStarts',
                         'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames', 'geneName', 'geneName2')

Chrom14_GeneList <- sort(unique(CH14_Genes$geneName2))

#Turn list into dataframe

Chromo14_List <- as.data.frame(Chrom14_GeneList)

#Maybe filter out the loci that aren't known genes. 
#Easiest way for now is to look which rows have LOC and just select those out. 

Chromo14_Filter <- as.data.frame(Chromo14_List[-(32:479),])


write.table(Chromo14_Filter, file="Chromo14_Genes.txt",quote=FALSE, row.names=FALSE, sep='\t')


#####Chromosome 19####

#Read in the gene table for chromosome 9 inverted region of 22Mb-42mb.

CH19_Genes <- read.table(file = 'C:/Users/cneme/Downloads/CHROMO_19_Genes', header = F)

#Add headers manually based on the describe table schema from table brower of USSC genome browser.

colnames(CH19_Genes) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes', 'chromStarts',
                          'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames', 'geneName', 'geneName2')

Chrom19_GeneList <- sort(unique(CH19_Genes$geneName2))

#Turn list into dataframe

Chromo19_List <- as.data.frame(Chrom19_GeneList)

#Maybe filter out the loci that aren't known genes. 
#Easiest way for now is to look which rows have LOC and just select those out. 

Chromo19_Filter <- as.data.frame(Chromo19_List[-(37:303),])

write.table(Chromo19_Filter, file="Chromo19_Genes.txt",quote=FALSE, row.names=FALSE, sep='\t')


#####Chromosome 24####

#Read in the gene table for chromosome 9 inverted region of 22Mb-42mb.

CH24_Genes <- read.table(file = 'C:/Users/cneme/Downloads/CHROMO_24_Genes', header = F)

#Add headers manually based on the describe table schema from table brower of USSC genome browser.

colnames(CH24_Genes) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes', 'chromStarts',
                          'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames', 'geneName', 'geneName2')

Chrom24_GeneList <- sort(unique(CH24_Genes$geneName2))

#Turn list into dataframe

Chromo24_List <- as.data.frame(Chrom24_GeneList)

#Maybe filter out the loci that aren't known genes. 
#Easiest way for now is to look which rows have LOC and just select those out. 

Chromo24_Filter <- as.data.frame(Chromo24_List[-(11:181),])

write.table(Chromo24_Filter, file="Chromo24_Genes.txt",quote=FALSE, row.names=FALSE, sep='\t')
