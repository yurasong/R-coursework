# FOCAL

sample="Esophagus"
setwd('results_Esophagus')
getwd()
genome <- read.table("../genome.hg19.easy.tsv", header=T, sep="\t")
genomeC <- read.table("../genome.hg19.tsv", header=T, sep="\t")

#setwd("~/Documents/Basic_workflow/Plot_Gistic")
# Input file are compressed just like the once send to GISTIC (log2(CN) - 1)
# loc <- "~/Documents/Basic_workflow/Plot_Gistic/Compressed/"
# postfix <- ".compressed"

# There are 3 important columns: 1, 5 and 6; Sample Name, group and Ploidy
# They should be in this order! other columns can be anything
# it also assumes this file has a header!!
#files <- read.table("samples.tsv", header=T, sep="\t")

# This file represents the genome and the cytobands
# and looks like this
# chrom  chromStart  chromEnd   type
# 1      0           121236957  short_arm


#nbrFiles=dim(files)[1]

# COL<-colorRampPalette(c('red','white','darkBlue'), bias=1, interpolate = c("spline"))
COL<-colorRampPalette(c('Navy','royalBlue','lightyellow','#AC0A00','darkred'), bias=1, interpolate = c("spline"))
# COL<-colorRampPalette(c('Navy','lightyellow','darkred'), bias=1, interpolate = c("spline"))

COLgroup <- c("#EC1E27","#2FAB49","#2F3592","#F57616","#83D6E4","#BD9E61")



### human 
cumChrLength = c(0,247199719,489950867,689397693,880660755,1061498620,1232395611,1391217034,1537491859,1677766684,1813141420,1947593803,	
                 2079880336,2194008315,2300368899,2400707813,2489530066,2568184807,	2644301959,2708108609,2770544572,2817488894,2867080325)

## ,3021664561)  X chromosome

#mouse mm10
#cumChrLength = c(0, 195471971, 377585195, 537624875, 694132991, 845967675, 995704221, 1141145680, 1270546893, 1395142003, 1525836996, 1647919539, 1768048561, 1888470200, 2013372444, 2117416129, 2215623897, 2310611168, 2401313807, 2462745373)
#2633776672, 2725521370, 2725537669)  mm10 X, Y, MT

# change scale here !!!!!!!!!!
# 15 is 10e-15
# 8 will be 10e-8
# adapt the thick marks also - line 171
plotrangeA=5
plotrangeD=5
barrange=sum(plotrangeA + plotrangeD)/20  ### 5% of sum plotrangeA + plotrange D

chrRange <- c(cumChrLength[1], sum(-cumChrLength[length(cumChrLength)]-120000000))  ## de 12m is toegevoegd om de labels te kunnen plaatsen  ### tsjek dit bij andere figuren
sampleRange <- c(-plotrangeD,plotrangeA)



#quartz()
#png("Figure2B.png", width=2000, height=3000)
pdf(paste0(sample,".scale_focal.pdf"), width=6, height=10)


############### pdf("Figure2B.pdf FLT N=38.pdf", width=10, height=6)
par(mar = c(0,0,0,0), oma=c(0,0,0,0), cex=1.5)
plot(sampleRange,chrRange,type = "n",ylab="",main="", xaxt = "n", yaxt="n", frame.plot=F)

#### Dit is de Chromosome bar! de breedte wordt in rect onderaan aangegeven

#for(j in 1:dim(genome)[1]) {
for(j in 1:22) {
  if(genome[j,1]=="X") {
    chrOffset=cumChrLength[23]
    col="black"
    colplot="#EBEBEB"
  } else {
    chr <- as.numeric(as.character(genome[j,1]))
    chrOffset=cumChrLength[chr]
    if(chr%%2==0) {
      col="#EBEBEB"
      colplot="black"	
    } else {
      col="black"
      colplot="#EBEBEB"
    }	
  }
  begin=as.numeric(genome[j,2])+chrOffset
  end=as.numeric(genome[j,3])+chrOffset
  begin=-begin
  end=-end
  
  type=as.character(genome[j,4])
  if(type=="centromere" || type=="heterochromatin") {
    col="darkgrey"
  }
  rect(-barrange,begin,barrange,end,col=col, border=col)
  
  if (col=="black") {		### deze if staat er omdat elk chromosome een label geven te groot is voor de mensen)
    text(0,((begin+end)/2), label= j, cex=0.8,col=colplot)
  }
  
  
  if (col=="black") {
    
    rect(-plotrangeD,begin,-barrange,end,col=colplot,border=colplot)
    rect(plotrangeA,begin,barrange,end,col=colplot,border=colplot)
  }
}


##########   hier begint FOCAAL Specifiek intermezzo

### in deze file is de q steeds verwijderd ## voorlopig nog manueel  ### voor human linken aan genome file...

focal<-read.delim("scores.gistic", header=T)

for(j in 1:dim(focal)[1]) {
  if(focal[j,1]=="Amp") { 
    # on the y-axis
    begin=sum (cumChrLength[as.numeric(as.character(focal[j,2]))],as.numeric(as.character(focal[j,3])))
    end=sum (cumChrLength[as.numeric(as.character(focal[j,2]))],as.numeric(as.character(focal[j,4])))
    
    begin=-begin
    end=-end
    
    cnv=as.numeric(focal[j,5])
    cnv=sum(cnv,barrange) # omdat er op 1 begonnnen wordt
    col="firebrick2"		
    segments(barrange,begin,cnv,begin,col=col, lwd=2)
    segments(cnv,begin,cnv,end,col=col, lwd=2)
    segments(barrange,end,cnv,end,col=col, lwd=2)
    # if (cnv>plotrangeA){
    #   newy=begin+30000000
    #   text(plotrangeA/8*6,newy,labels= "q-value<10e-40" , col="gold", cex=0.5)}
    
  }
  if(focal[j,1]=="Del") { 
    # on the y-axis
    begin=sum (cumChrLength[as.numeric(as.character(focal[j,2]))],as.numeric(as.character(focal[j,3])))
    end=sum (cumChrLength[as.numeric(as.character(focal[j,2]))],as.numeric(as.character(focal[j,4])))
    
    
    begin=-begin
    end=-end
    cnv=as.numeric(focal[j,5]) 
    ### hier hebben we reeds de te plotten waarde
    cnv=sum(-cnv,-barrange)		# omdat er op 1 begonnnen wordt
    col="dodgerblue"		
    segments(-barrange,begin,cnv,begin,col=col, lwd=2)
    segments(cnv,begin,cnv,end,col=col, lwd=2)
    segments(-barrange,end,cnv,end,col=col, lwd=2)
    
  }
}

## Ad threshold of -log10 (Q-value)  = 0.25
threshold=-log(0.25,10)
thresholdamp=sum(barrange,threshold)
thresholddel=sum(-barrange,-threshold)	
col="darkgreen"
segments(thresholdamp,min(cumChrLength),thresholdamp,sum(-max(cumChrLength),-70000000),col=col, lwd=1, lty ="dashed")
segments(thresholddel,min(cumChrLength),thresholddel,sum(-max(cumChrLength),-70000000),col=col, lwd=1, lty="dashed")
text(thresholdamp,-max(cumChrLength)-130000000,labels= "0.25" , col="darkgreen", cex=0.65)
text(thresholddel,-max(cumChrLength)-130000000,labels= "0.25" , col="darkgreen", cex=0.65)

###  Kader
rect(-plotrangeD,min(cumChrLength), plotrangeA,-max(cumChrLength))

### X-as  .. Add thicks     ###### Ad a IF ..thick function that adds thicks as long as the thicks are in the plot range,

for (P in c(1,2,3,4)){
  Q=10^-P
  
  X=-log(Q,10)
  X=sum(X,barrange)	
  
  
  if(X< plotrangeA ) { 
    segments(X,-max(cumChrLength),X,sum(-max(cumChrLength),-20000000),col="darkgreen", lwd=1)
    text(X,(-max(cumChrLength)-60000000), labels= Q , cex=0.5, col="darkgreen") 
  }
  ## part2 deletion
  if(X< plotrangeD ) { 	
    X=-log(Q,10)
    X=sum(-X,-barrange)	
    segments(X,-max(cumChrLength),X,sum(-max(cumChrLength),-20000000),col="darkgreen", lwd=1)
    text(X,(-max(cumChrLength)-60000000), labels= Q , cex=0.5, col="darkgreen") 
    
    
  }}



# line for the "1" bar
segments(barrange,-max(cumChrLength),barrange,sum(-max(cumChrLength),-20000000),col="darkgreen", lwd=1)
text(barrange,(-max(cumChrLength)-60000000), labels= "1" , cex=0.5, col="darkgreen") 
segments(-barrange,-max(cumChrLength),-barrange,sum(-max(cumChrLength),-20000000),col="darkgreen", lwd=1)
text(-barrange,(-max(cumChrLength)-60000000), labels= "1" , cex=0.5, col="darkgreen") 

# range <- 1e3
# for(i in 1:range) {
# 	rect(nbrFiles/range*(i-1),sum(chrRange)*1.01,nbrFiles/range*i,sum(chrRange)*1.03, col=COL(range)[i], border=NA)
# }

#title
mtext(paste0("Shallow seq - ",sample," - focal amplifications (red) and deletions (blue)"),cex.main=0.65, side = 3, line = -0.8, outer = TRUE)

dev.off()


