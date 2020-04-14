# BROAD

sample="Esophagus"
setwd('results_Esophagus')
getwd()
genome <- read.table("../genome.hg19.easy.tsv", header=T, sep="\t")
genomeC <- read.table("../genome.hg19.tsv", header=T, sep="\t")

COL<-colorRampPalette(c('Navy','royalBlue','lightyellow','#AC0A00','darkred'), bias=1, interpolate = c("spline"))
COLgroup <- c("#EC1E27","#2FAB49","#2F3592","#F57616","#83D6E4","#BD9E61")

### human 
cumChrLength = c(0,247199719,489950867,689397693,880660755,1061498620,1232395611,1391217034,1537491859,1677766684,1813141420,1947593803,
                 2079880336,2194008315,2300368899,2400707813,2489530066,2568184807,	2644301959,2708108609,2770544572,2817488894,2867080325)

## ,3021664561)  X chromosome
#mouse mm10
#cumChrLength = c(0, 195471971, 377585195, 537624875, 694132991, 845967675, 995704221, 1141145680, 1270546893, 1395142003, 1525836996, 1647919539, 1768048561, 1888470200, 2013372444, 2117416129, 2215623897, 2310611168, 2401313807, 2462745373)
#2633776672, 2725521370, 2725537669)  mm10 X, Y, MT

# the grey background seperated for each chromosome

# change scale here !!!!!!!!!!
# 15 is 10e-15
# 8 will be 10e-8
# adapt the thick marks also - line 298
plotrangeA=2
plotrangeD=2
barrange=sum(plotrangeA + plotrangeD)/20  ### 5% of sum plotrangeA + plotrange D - chr bar

chrRange <- c(cumChrLength[1], sum(-cumChrLength[length(cumChrLength)]-120000000))  ## de 12m is toegevoegd om de labels te kunnen plaatsen  ### tsjek dit bij andere figuren
sampleRange <- c(-plotrangeD,plotrangeA)

# to create a pdf
pdf(paste0(sample,".scale_broad.pdf"), width=6, height=10)


#par(mfrow=c(1,2))
par(mar = c(0,0,0,0), oma=c(0,0,0,0), cex=1.5)
plot(sampleRange,chrRange, type = "n",ylab="",xaxt = "n", yaxt="n", frame.plot=F)

#the chromosome bar - with 22 chromosomes in black background
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
  begin=as.numeric(genome[j,2])+chrOffset-10000
  end=as.numeric(genome[j,3])+chrOffset
  begin=-begin
  end=-end
  
  type=as.character(genome[j,4])
  if(type=="centromere" || type=="heterochromatin") {
    col="darkgrey"
  }
  rect(-barrange,begin,barrange,end,col=col, border=col)
  
  if (col=="black") {		### deze if staat er omdat elk chromosome een label geven te groot is voor de mensen)
    text(0,((begin+end)/2), label= j, cex=0.65, col=colplot)
  }
  
  
  if (col=="black") {
    
    rect(-plotrangeD,begin,-barrange,end,col=colplot,border=colplot)
    rect(plotrangeA,begin,barrange,end,col=colplot,border=colplot)
  }
}

##########   hier begint FOCAAL Specifiek intermezzo

# broad and focal amplifications/deletions of a chromosome region containing multiple genes
# export them from GISTIC results - broad significance results an scores gistic
broad<-read.delim("broad_significance_results.txt", header=T)
#focal<-read.delim(paste0(sample,"_GIS"), header=T)

### Aplifications, red
for (k in 1:22){	
  ##	Per chrom laden we P an Q in 
  for(l in 1:dim(genomeC)[1])	{
    
    if((genomeC[l,1])==k) {					
      if 	( genomeC[l,4]=="long_arm") { 
        print.noquote("long")
        Qstart=genomeC[l,2]
        Qstop=genomeC[l,3]
      }  		
      
      if 	( genomeC[l,4]=="short_arm") { 
        Pstart=genomeC[l,2]
        Pstop=genomeC[l,3]
        print.noquote("short")
      }  
    }
  }
  
  ### hier weet hij start en stop voor P en Q voor chrom K
  
  ## nu zoeken we de CNV
  # see in broad_sign_results obtained from GISTIC if lines 6 and 10 actually correspond 
  # to amplification q-values and delection q-values
  for(j in 1:dim(broad)[1]) {
    parm=paste(k,"p",sep="")
    
    if (parm == "13p") {
      pACNV = 999
    }
    if (parm == "14p") {
      pACNV = 999
    }
    if (parm == "15p") {
      pACNV = 999
    }
    if (parm == "21p") {
      pACNV = 999
    }
    if (parm == "22p") {
      pACNV = 999
    }
    if (broad[j,1] == parm) {
      pACNV <-(broad[j,6])  
    } 
    
    
    qarm=paste(k,"q",sep="")
    if (broad[j,1] == qarm) {
      qACNV <-(broad[j,6])  
    } 
    
  }
  
  for(j in 1:dim(broad)[1]) {
    parm=paste(k,"p",sep="")
    
    if (parm == "13p") {
      pDCNV = 999
    }
    if (parm == "14p") {
      pDCNV = 999
    }
    if (parm == "15p") {
      pDCNV = 999
    }
    if (parm == "21p") {
      pDCNV = 999
    }
    if (parm == "22p") {
      pDCNV = 999
    }
    if (broad[j,1] == parm) {
      pDCNV <-(broad[j,10]) 
    }
    
    qarm=paste(k,"q",sep="")
    if (broad[j,1] == qarm) {
      qDCNV <-(broad[j,10]) 
    } 
  }
  
  
  # blue (del) and red (amp) rectangles 
  beginQ=sum(cumChrLength[k],Qstart)
  endQ=sum(cumChrLength[k],Qstop)
  
  beginP=sum(cumChrLength[k],Pstart)
  endP=sum(cumChrLength[k],Pstop)
  
  if (pACNV <0.00000000000000000000000000000000000000000000000000000000001) {
    pAcnv = plotrangeA
  } else {
    pAcnv=-log(pACNV,10)
  }
  if ( parm == "13p") {
    pAcnv = 0 }
  if ( parm == "14p") {
    pAcnv = 0 }
  if ( parm == "15p") {
    pAcnv = 0 }
  if ( parm == "21p") {
    pAcnv = 0 }
  if ( parm == "22p") {
    pAcnv = 0 }
  if (qACNV <0.00000000000000000000000000000000000000000000000000000000001){
    qAcnv = plotrangeA
  } else {
    qAcnv=-log(qACNV,10)
  }
  
  pAcnv=sum(pAcnv,barrange)
  qAcnv=sum(qAcnv,barrange)
  
  col="firebrick2"	
  rect(barrange,-beginQ,qAcnv,-endQ,col=col, border=col)
  rect(barrange,-beginP,pAcnv,-endP,col=col, border=col)
  
  # add an "*" if q-values = 0
  if (qACNV<0.00000000000000000000000000000000000000000000000000000000001){
    newyQ=-(beginQ+endQ)/2
    text(plotrangeA/8*6,newyQ,labels= "q-value<10e-40" , col="gold", cex=0.5)}
  # if(qACNV< 10^-13 & qACNV > 10^-15){
  #   qACNV =10^-13
  #   newyQ=-(beginQ+endQ)/2
  #   text(plotrangeA/8*6,newyQ,labels= "q-value < 10E-14" , col="gold", cex=0.5)}
  if (pACNV <0.00000000000000000000000000000000000000000000000000000000001){
    newyP=-(beginP+endP)/2
    text(plotrangeA/8*6,newyP,labels= "q-value<10e-40" , col="gold", cex=0.5)}
  
  col="dodgerblue"
  if (pDCNV <0.00000000000000000000000000000000000000000000000000000000001) {
    pDcnv = plotrangeD
  } else {
    pDcnv=-log(pDCNV,10)
  }
  if ( parm == "13p") {
    pDcnv = 0 }
  if ( parm == "14p") {
    pDcnv = 0 }
  if ( parm == "15p") {
    pDcnv = 0 }
  if ( parm == "21p") {
    pDcnv = 0 }
  if ( parm == "22p") {
    pDcnv = 0 }
  
  if (qDCNV <0.00000000000000000000000000000000000000000000000000000000001){
    qDcnv = plotrangeD
  } else {
    qDcnv=-log(qDCNV,10)
  }
  
  print("bar")
  print(barrange)
  print(parm)
  print("pACNV")
  print(pACNV)
  print("pAcnv")
  print(pAcnv)
  print("pDCNV")
  print(pDCNV)
  print("pDcnv")
  print(pDcnv)
  
  #####
  #############
  pDcnv=sum(pDcnv,barrange)
  qDcnv=sum(qDcnv,barrange)
  rect(-barrange,-beginQ,-qDcnv,-endQ,col=col, border=col)
  rect(-barrange,-beginP,-pDcnv,-endP,col=col, border=col)
  if (qDCNV <0.00000000000000000000000000000000000000000000000000000000001){
    newyQ=-(beginQ+endQ)/2
    text(-plotrangeD/8*6,newyQ,labels= "q-value<10e-40" , col="gold", cex=0.5)
  }
  if (pDCNV<0.00000000000000000000000000000000000000000000000000000000001){
    newyP=-(beginP+endP)/2
    text(-plotrangeD/8*6,newyP,labels= "q-value<10e-40" , col="gold", cex=0.5)
  }
  # print(parm)
  # print("pACNV")
  # print(pACNV)
  # print("pAcnv")
  # print(pAcnv)
  # print("pDCNV")
  # print(pDCNV)
  # print("pDcnv")
  # print(pDcnv)}
}

#dev.off()

## Ad threshold of -log10 (Q-value)  = 0.25
threshold=-log(0.25,10)
thresholdamp=sum(barrange,threshold)
thresholddel=sum(-barrange,-threshold)	
col="darkgreen"
segments(thresholdamp,min(cumChrLength),thresholdamp,sum(-max(cumChrLength),-70000000),col=col, lwd=1, lty ="dashed")
segments(thresholddel,min(cumChrLength),thresholddel,sum(-max(cumChrLength),-70000000),col=col, lwd=1, lty="dashed")

text(thresholdamp,-max(cumChrLength)-110000000,labels= "0.25" , col="darkgreen", cex=0.65)
text(thresholddel,-max(cumChrLength)-110000000,labels= "0.25" , col="darkgreen", cex=0.65)

###  black rectangle around the plot
rect(-plotrangeD,min(cumChrLength), plotrangeA,-max(cumChrLength))

#add thick - the black bars at the bottom of the plot
# Ad a IF..thick function that adds thicks as long as the thicks are in the plot range

for (P in c(1,3,6,9,12)){
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

# add a segment and text for the "1" y axis
segments(barrange,-max(cumChrLength),barrange,sum(-max(cumChrLength),-20000000),col="darkgreen", lwd=1)
text(barrange,-max(cumChrLength)-60000000, labels= "1" , cex=0.5, col="darkgreen") 
segments(-barrange,-max(cumChrLength),-barrange,sum(-max(cumChrLength),-20000000),col="darkgreen", lwd=1)
text(-barrange,-max(cumChrLength)-60000000, labels= "1" , cex=0.5, col="darkgreen") 
# range <- 1e3
# for(i in 1:range) {
# 	rect(nbrFiles/range*(i-1),sum(chrRange)*1.01,nbrFiles/range*i,sum(chrRange)*1.03, col=COL(range)[i], border=NA)
# }

### 2nd info in the same graph

pct<-read.delim(paste0("pct_",sample,".txt"), header=T)

### Aplifications, red
for (k in 1:22){	
  ##	Per chrom laden we P an Q in 
  for(l in 1:dim(genomeC)[1])	{
    
    
    if((genomeC[l,1])==k) {					
      if 	( genomeC[l,4]=="long_arm") { 
        print.noquote("long")
        Qstart=genomeC[l,2]
        Qstop=genomeC[l,3]
      }  		
      
      if 	( genomeC[l,4]=="short_arm") { 
        Pstart=genomeC[l,2]
        Pstop=genomeC[l,3]
        print.noquote("short")
      }  
    }
  }
  
  ### hier weet hij start en stop voor P en Q voor chrom K
  
  ## nu zoeken we de CNV
  # see in broad_sign_results obtained from GISTIC if lines 6 and 10 actually correspond 
  # to amplification q-values and delection q-values
  
  #  for(j in 1:dim(broad)[1]) {
  #      if(broad[j,2:dim(broad)]>0.1) {
  #        #text(plotrangeA,chrOffset,labels=broad[0,2:dim(broad)], cex=0.4)
  #        }
  #    else{
  #    }
  # }
  #    }
  #    
  for(j in 1:dim(pct)[1]) {
    parm=paste(k,"p",sep="")
    if (pct[j,1] == parm) 
    {pAPCT <-(pct[j,3])}
    
    qarm=paste(k,"q",sep="")
    if (pct[j,1] == qarm) 
    {qAPCT <-(pct[j,3])}
  }
  
  for(j in 1:dim(pct)[1]) {
    parm=paste(k,"p",sep="")
    if (pct[j,1] == parm) 
    {pDPCT <-(pct[j,5])}
    
    qarm=paste(k,"q",sep="")
    if (pct[j,1] == qarm) 
    {qDPCT <-(pct[j,5])}
  }
  
  #
  # rectangles (del) and (amp)  
  # y axis
  beginQ=sum(cumChrLength[k],Qstart)
  endQ=sum(cumChrLength[k],Qstop)
  beginP=sum(cumChrLength[k],Pstart)
  endP=sum(cumChrLength[k],Pstop)
  # x axis
  pApct=(plotrangeA-barrange)/1*pAPCT
  #-log(pAPCT,10)
  pApct=sum(pApct,barrange)
  qApct=(plotrangeA-barrange)/1*qAPCT
  #-log(qAPCT,10)
  qApct=sum(qApct,barrange)
  
  col="darkslateblue"	
  
  # rect(barrange,-beginQ,qApct,-endQ,col=NA, border="black")
  # rect(barrange,-beginP,pApct,-endP,col=NA, border="black")
  # lines.default(qApct,sum(-beginQ,-endQ),col="gold")
  # lines.default(pApct,sum(-beginP,-endP),col="gold")
  # abline(v=qApct,col="black")
  # abline(v=pApct,col="black")
  
  # p arm is first
  # q arm is second
  # segments(x0=qApct,x1=qApct,y0=-beginQ,y1=-endQ,col="darkslateblue")
  # segments(x0=pApct,x1=pApct,y0=-beginP,y1=-endP,col="darkslateblue")
  # segments(x0=barrange,x1=qApct,y0=-beginQ,y1=-beginQ,col="darkslateblue", lty="dotted")
  # segments(x0=barrange,x1=pApct,y0=-endP,y1=-endP,col="darkslateblue", lty="dotted")
  # segments(x0=barrange,x1=qApct,y0=-endQ,y1=-endQ,col="darkslateblue", lty="dotted")
  # segments(x0=barrange,x1=pApct,y0=-beginP,y1=-beginP,col="darkslateblue", lty="dotted")
  # 
  
  segments(x0=qApct,x1=qApct,y0=-beginQ,y1=-endQ,col="darkslateblue")
  segments(x0=pApct,x1=pApct,y0=-beginP,y1=-endP,col="darkslateblue")
  segments(x0=barrange,x1=qApct,y0=-beginQ,y1=-beginQ,col="darkslateblue", lwd=0.5)
  segments(x0=barrange,x1=pApct,y0=-endP,y1=-endP,col="darkslateblue", lwd=0.5)
  segments(x0=barrange,x1=qApct,y0=-endQ,y1=-endQ,col="darkslateblue", lwd=0.5)
  segments(x0=barrange,x1=pApct,y0=-beginP,y1=-beginP,col="darkslateblue", lwd=0.5)
  
  
  pDpct=(plotrangeD-barrange)/1*(pDPCT)
  #-log(pDPCT,10)
  pDpct=sum(pDpct,barrange)
  qDpct=(plotrangeD-barrange)/1*(qDPCT)
  #-log(qDPCT,10)
  qDpct=sum(qDpct,barrange)
  # rect(-barrange,-beginQ,-qDpct,-endQ,col=NA, border="black")
  # rect(-barrange,-beginP,-pDpct,-endP,col=NA, border="black")
  # lines.default(qDpct,sum(-beginQ,-endQ),col="gold")
  # lines.default(pDpct,sum(-beginP,-endP),col="gold")
  # abline(v=qDpct,col="black")
  # abline(v=pDpct,col="black")
  # segments(x0=-qDpct,x1=-qDpct,y0=-beginQ,y1=-endQ,col="darkslateblue")
  # segments(x0=-pDpct,x1=-pDpct,y0=-beginP,y1=-endP,col="darkslateblue")
  # segments(x0=-barrange,x1=-qDpct,y0=-beginQ,y1=-beginQ,col="darkslateblue", lty="dotted")
  # segments(x0=-barrange,x1=-pDpct,y0=-endP,y1=-endP,col="darkslateblue", lty="dotted")
  # segments(x0=-barrange,x1=-qDpct,y0=-endQ,y1=-endQ,col="darkslateblue", lty="dotted")
  # segments(x0=-barrange,x1=-pDpct,y0=-beginP,y1=-beginP,col="darkslateblue", lty="dotted")
  # 
  segments(x0=-qDpct,x1=-qDpct,y0=-beginQ,y1=-endQ,col="darkslateblue")
  segments(x0=-pDpct,x1=-pDpct,y0=-beginP,y1=-endP,col="darkslateblue")
  segments(x0=-barrange,x1=-qDpct,y0=-beginQ,y1=-beginQ,col="darkslateblue", lwd=0.5)
  segments(x0=-barrange,x1=-pDpct,y0=-endP,y1=-endP,col="darkslateblue", lwd=0.5)
  segments(x0=-barrange,x1=-qDpct,y0=-endQ,y1=-endQ,col="darkslateblue", lwd=0.5)
  segments(x0=-barrange,x1=-pDpct,y0=-beginP,y1=-beginP,col="darkslateblue", lwd=0.5)
  
}

## Add a 10% mark
threshold=(plotrangeA-barrange)/10
thresholdamp=sum(barrange,threshold)
thresholddel=sum(-barrange,-threshold)	
col="darkslateblue"
segments(thresholdamp,min(cumChrLength),thresholdamp,min(cumChrLength)-20000000,col=col, lwd=1)
segments(thresholddel,min(cumChrLength),thresholddel,min(cumChrLength)-20000000,col=col, lwd=1)
text(thresholdamp,min(cumChrLength)+20000000,labels= "10%" , col=col, cex=0.5)
text(thresholddel,min(cumChrLength)+20000000,labels= "10%" , col=col, cex=0.5)

# add 50% mark
threshold=(plotrangeA-barrange)/2
#-log(0.5,10)
thresholdamp=sum(barrange,threshold)
thresholddel=sum(-barrange,-threshold)	
col="darkslateblue"
segments(thresholdamp,min(cumChrLength),thresholdamp,min(cumChrLength)-20000000,col=col, lwd=1)
segments(thresholddel,min(cumChrLength),thresholddel,min(cumChrLength)-20000000,col=col, lwd=1)
text(thresholdamp,min(cumChrLength)+20000000,labels= "50%" , cex=0.5, col=col)
text(thresholddel,min(cumChrLength)+20000000,labels= "50%" , cex=0.5, col=col)
# add a segment and text for the "1" y axis
# segments(barrange,-max(cumChrLength),barrange,sum(-max(cumChrLength),-20000000),col="black", lwd=1)
# text(barrange,(-max(cumChrLength)-60000000), labels= "1" , cex=0.5) 
# segments(-barrange,-max(cumChrLength),-barrange,sum(-max(cumChrLength),-20000000),col="black", lwd=1)
# text(-barrange,(-max(cumChrLength)-60000000), labels= "1" , cex=0.5) 

# add a segment and text for the "100%" y axis
col="darkslateblue"
segments(plotrangeA,-max(cumChrLength),plotrangeA,-max(cumChrLength),col=col, lwd=1)
info=plotrangeA
text(info,min(cumChrLength)+20000000, labels= "100%" , cex=0.5, col=col) 
info=-plotrangeD
segments(-plotrangeD,-max(cumChrLength),-plotrangeD,-max(cumChrLength),col=col, lwd=1)
text(info,min(cumChrLength)+20000000, labels= "100%" , cex=0.5, col=col)

mtext(paste0("Shallow seq - ",sample," - amplifications (red) and deletions (blue)"),cex.main=0.65, side = 3, line = -0.8, outer = TRUE)
# I added a title to the plot
# mtext("green: q-values",cex=0.7, side = 1, col="darkgreen", at=-10, line=-1.8)
# mtext("purple: % sample concerned",cex=0.7, side = 1, col="darkslateblue", at=-10,line=-1.3)

mtext("green: q-values",cex=0.7, side = 1, col="darkgreen", at=0, line=-1.5)
mtext("purple: % sample concerned",cex=0.7, side = 1, col="darkslateblue", at=0,line=-1)


#mtext("yellow: q-values almost 0",cex=0.7, side = 1, col="darkgoldenrod2", at =10,line=-1.3)

#dev.off()
