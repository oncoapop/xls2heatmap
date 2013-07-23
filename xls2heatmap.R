# R-script to automatedly plot the bionominal output of Deep Seq data
# Change the variables between the ####### to make the script run correctly

# source("http://bioconductor.org/biocLite.R")
# BiocLite()
# source("http://www.bioconductor.org/biocLite.R"); biocLite("VariantAnnotation")
# biocLite(c("GenomicFeatures", "AnnotationDbi"))
# install.packages("XLConnect") 
# install.packages("rJava")

library("VariantAnnotation")
library("IRanges")
library("GenomicRanges")
library(foreign)
library(lattice)
library(XLConnect)
require(Heatplus)
library(limma)


basedir <- "C:/Users/dyap_000/Documents/R"
dirname <- paste(basedir, "BioNom_QC", sep="/")

setwd(dirname)
# This is where the output of this R-script will be
outpath <- dirname

SA="SA029_MS130708"

#File names matching a pattern
pat1=paste("*",SA,sep="")
pat=paste(pat1, "*.*xls", sep="_")
file_names=list.files(pattern=pat)
file_names
#########################################################################


############# DO NOT CHANGE ANYTHING BELOW THIS LINE #################

file = file_names[1]

#Load workbook one by one
wb <- loadWorkbook ( file , create = FALSE )

#Load the second worksheet which is variant allelic freq
sheets <- c(getSheets(wb)[2],getSheets(wb)[1])
sheets

dflist <- readWorksheet(wb, sheet = sheets, header = TRUE, startRow = c(1,1), startCol = c(1,1), endCol = c(49,49), endRow = c(289,289))

########################################################
# dflist[2] is the reads (depth)
# Col 1 is the ID, the data starts from col 2
ef <- data.matrix(dflist[[2]][2:ncol(dflist[[1]])])

# Put the ID back as rownames
rownames(ef) <- dflist[[2]]$target

# Selecting only those position with at least 100 reads in the normal
# Normal sample
# Uncommet this section out if normal controls failed for some reason
sample=colnames(ef)[1]
ff <- ef[ef[, sample] > 100,]
gf <- ff[complete.cases(ff),]
# Select single cells positions (at least one read out in a sample)
sel <- subset(ef[rowSums(!is.na(ef)) > 1,])

reads=paste(SA,"reads.csv",sep="-")
write.csv(sel, file=reads)

hmcols<-colorRampPalette(c("dark green","red"))(1000)
readpdf=paste(SA,"reads.pdf",sep="-")
pdf(readpdf, width=6, height=6)
title=paste(SA, "reads", sep=" ")
levelplot(sel, cex.axis=0.8, col.regions = hmcols, xlab="SNV positions", ylab="Samples", main=title, aspect="fill", cuts=1000)

dev.off()

##########################################################
# dflist[1] is the allele Freq sheet
# Col 1 is the ID, the data starts from col 2
jf <- data.matrix(dflist[[1]][2:ncol(dflist[[2]])])

# Put the ID back as rownames
rownames(jf) <- dflist[[1]]$target

# Selecting only those position with at least MAF > 0.2 in the normal
# ff <- jf[jf[, sample] > 0.2,]
# kf <- jf[complete.cases(jf),]
self <- subset(jf[rowSums(!is.na(ef)) > 2,])

freq=paste(SA,"freq.csv",sep="-")
write.csv(self, file=freq)

hmcols<-colorRampPalette(c("dark green","red"))(1000)
freqpdf=paste(SA,"freq.pdf",sep="-")
pdf(freqpdf, width=6, height=6)
title=paste(SA, "Alt. Allele Freq", sep=" ")
levelplot(self, cex.axis=0.8, col.regions = hmcols, xlab="SNV positions", ylab="Samples", main=title, aspect="fill", cuts=1000, at=c(seq(0.2,1, by=0.1)))

dev.off()

freqpdf=paste(SA,"freq-clus.pdf",sep="-")
pdf(freqpdf, width=6, height=6)
reg2 = regHeatmap(self)
plot(reg2)
dev.off()



#####################################################################################################

