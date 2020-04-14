#Calling copy-number using sWGS or SNP arrays: Identify copy-number changes in sWGS 
#Relative copy-number calling using shallow whole-genome sequencing

#!/usr/bin/Rscript --slave

#source("https://bioconductor.org/biocLite.R")
#biocLite("QDNAseq")

argv <- commandArgs(TRUE)

library("methods")
library('matrixStats')
library("QDNAseq")

#Get 30kb bin annotations for hg19 genome
bins <- getBinAnnotations(binSize=30)

#Plot the readcounts with filtered reads highlighted.
readCounts <- binReadCounts(bins, bamfiles=argv[1])
png("plot_read_counts.png", width=1000)
plot(readCounts, logTransform=FALSE)
highlightFilters(readCounts, logTransform=TRUE, residual=TRUE, blacklist=TRUE)
dev.off()

#Apply QDNAseq filters
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
isobarPlot(readCountsFiltered)

#Estimate GC correction and noise plot
readCountsFiltered <- estimateCorrection(readCountsFiltered)
noisePlot(readCountsFiltered)

#Normalise and smooth outliers.
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
png("plot_copynumber.png", width=1000)
plot(copyNumbersSmooth)

#Segment the copy-number profile.
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented)

#Call copy-number
copyNumbersCalled <- callBins(copyNumbersSegmented)
#Plot final profile.
plot(copyNumbersCalled)

dev.off()
exportBins(copyNumbersSmooth, file=argv[2])

