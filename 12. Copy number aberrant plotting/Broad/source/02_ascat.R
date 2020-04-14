source("/home/audrey/Script_Blanpainlab/R_library/ascat_2.5.2.R")

#If the only tumor is given, further step is not working so you should put files for two times. See line 15.
ascat.bc = ascat.loadData("logR","BAF", "logR","BAF", 
                          chrs = c(1:22,"X","Y"), gender = NULL, sexchromosomes = c("X","Y"))
ascat.plotRawData(ascat.bc)

#ascat.gg corresponds to germline genotypes
ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=NULL)
ascat.plotSegmentedData(ascat.bc)

#Calculating the allele-specific copy numbers
ascat.output = ascat.runAscat(ascat.bc)


#Error log
#ascat.bc = ascat.loadData("logR_SCS18P5R1.fastq.txt.bz2", "BAF_SCS18P5R1.fastq.txt.bz2", chrs = c(1:22,"X","Y"), gender = NULL, sexchromosomes = c("X","Y"))
#ascat.plotRawData(ascat.bc)
#ascat.bc = ascat.aspcf(ascat.bc, ascat.gg=NULL, selectsamples = 1)
#Error in 1:dim(hom)[2] : argument of length 0
