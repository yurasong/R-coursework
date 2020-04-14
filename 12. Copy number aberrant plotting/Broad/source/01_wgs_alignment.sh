#!/bin/bash

rougefonce='\e[0;31m'
blanc='\e[1;37m'
bleuclair='\e[1;34m'
vertfonce='\e[0;32m'
rose='\e[1;35m'

##################################################################################################
#INSTRUCTIONS
echo -e "${rougefonce}INSTRUCTIONS TO USE THIS SCRIPT :"
echo -e "The sample name should be in one word without special characters. For example P1et2blood not P1-2-blood"
echo -e "To have downloaded and install : BWA, R, picardtools, samtools, QDNA seq, Ascat"
echo -e "For analysis, please input hg19 as reference genome."
echo -e "Parameters of the script should be in this order genome_name location_of_the_genome_indexes filetumor_R1.fastq location_of_the_fastq_files outputdirectory directory_with_homozygous.bed_file" 
echo -e "${vertfonce}To interupt this script at any time you can just type Ctrl+C"
#######################################################################################

##############################################################################
# Verifying the arguments of the script
if [[ $# -ne 6 ]]
    then 
    echo "Bad number of argument"
    exit
fi
if [[ "$1" != mm10 && "$1" != mm9 && "$1" != hg19 ]]
    then
    echo "Bad genome : please choose mm10 or mm9 or hg19"
    exit
fi
if [[ "$3" != ${3%.*}.fastq ]]
    then
    echo "Bad tumor R1 file name. Should be a fastq file"
    exit
fi
if [[ ! -d "$2" ]]
    then 
    echo "directory with genome indexes doesn't exist"
    exit
fi
if [[ ! -f "$2/$1.amb" || ! -f "$2/$1.ann" || ! -f "$2/$1.pac" || ! -f "$2/$1.sa" ]]
    then 
    echo "index files are not computed. You need to run bwa index first"
    exit
fi
if [[ ! -d "$4" ]]
    then 
    echo "directory with fastq files doesn't exist"
    exit
fi
if [[ ! -f "$4/$3" ]]
    then 
    echo "tumor R1 file doesn't exist"
    exit
fi
if [[ ! -d "$5" ]]
    then 
    echo "output directory doesn't exist"
    exit
fi
if [[ ! -w "$5" ]]
    then 
    echo "You don't have the permission to write in $5"
    exit
fi

if [[ ! -d "$6" ]]
    then 
    echo "directory with homozygous.bed file doesn't exist"
    exit
fi
if [[ ! -f "$6/homozygous.bed" ]]
    then 
    echo "homozygous.bed file doesn't exist"
    exit
fi
#############################################################################

#############################################################################
# The file hg38 of the 8th of february 2017 correspond to primary dna of hg38 version 87 genome downloaded from ensembl ftp.

echo -e "${vertfonce}Read groups${bleuclair}"
read l1 < $4/$3
flowcell=$(echo $l1 | cut -d ':' -f 3)
lane=$(echo $l1 | cut -d ':' -f 4)
index=$(echo $l1 | cut -d ':' -f 10)
samplename=$(echo $3 | cut -d '-' -f 2 | cut -d '_' -f 1)
RG=`echo "@RG\tID:$flowcell.$lane.$index\tSM:$samplename\tPL:illumina\tLB:lib1\tPU:$flowcell.$lane.$index" `

#Fastqc
#echo -e "${vertfonce}Quality controls :${bleuclair}"
#fastqc $4/$3

echo -e "${vertfonce}Remove adapter sequences :${bleuclair}"
TrimmomaticSE -threads 4 -phred33 $4/$3 $5/trim_R1_$samplename.fq ILLUMINACLIP:/home/audrey/reference_genomes/GATK/adapters_exome.txt:2:30:10 HEADCROP:5 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Remove mouse contamination
echo -e "${vertfonce}Remove mouse contamination :${bleuclair}"
bbsplit.sh ambiguous=best ambiguous2=toss maxindel=900000 refstats=stat_alignment_$samplename.txt threads=4 qtrim=f untrim=f minratio=0.65 in1=$5/trim_R1_$samplename.fq ref=$2/$1.fa,$2/mm9.fa basename=decontamin_$samplename.%_#.fq outu1=unmapped1_$samplename.fq -Xmx50g
##-Xmx option should be at least 50g, if not, OutOfMemoryError: Java heap space.

# Realigning files on human genome
echo -e "${vertfonce}Realigning human specific reads with an other alignment software${bleuclair}"
bwa mem -t 5 $2/$1 $5/decontamin_$samplename.hg19_1.fq > $5/alignment_$samplename.sam

# Converting sam to bam
echo -e "${vertfonce}Converting sam to bam ${bleuclair}"
samtools view -bS "$5/alignment_$samplename.sam" > "$5/alignment_$samplename.bam"

# Replace read groups in a BAM file.
##This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.
echo -e "${vertfonce}Replacing read groups in a BAM file.${bleuclair}"
picard.jar AddOrReplaceReadGroups I=$5/alignment_$samplename.bam O=$5/alignment_replace_$samplename.bam SORT_ORDER=coordinate RGPL=Illumina RGSM=$samplename RGLB=LB4 RGPU=Lane VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

# Sorting bam file
echo -e "${vertfonce}Sorting bam${bleuclair}"
picard.jar SortSam I=$5/alignment_replace_$samplename.bam O=$5/alignment_$samplename.sort.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT

#picard metrics
repavant=`echo "$5/Metrics_beforeMarkDup_$samplename"`
mkdir $repavant
picard.jar CollectMultipleMetrics I=$5/alignment_$samplename.sort.bam O=$repavant"/multiple_metrics" R=$2/$1".fa" PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=LENIENT


#Marking duplicates
echo -e "${vertfonce}Mark duplicates${bleuclair}"
repapres=`echo "$5/Metrics_afterMarkDup_$samplename"`
mkdir $repapres
picard.jar MarkDuplicates I=$5/alignment_$samplename.sort.bam O=$5/alignment_$samplename.sort.dedup.bam M=$repapres/marked_dup_metrics_$samplename.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true

#Indexing bam file
echo -e "${vertfonce}index${bleuclair}"
samtools index -b $5/alignment_$samplename.sort.dedup.bam $5/alignment_$samplename.sort.dedup.bai

echo -e "${vertfonce}Alignment has been finished!${bleuclair}"

##qDNAseq 
echo -e "${vertfonce}QDNA Rscript${bleuclair}"
Rscript /home/audrey/Script_Blanpainlab/R_library/qDNASeq.R $5/alignment_$samplename.sort.dedup.bam $5/alignment_$samplename.sort.dedup.txt
homozygous=`echo "$6/homozygous.bed"`

##Detecting and plotting Copy number variation
# Create ASCAT input QDNAseq
echo -e "Making BAF and logR files"

awk 'BEGIN{FS="\t";OFS="\t"} FNR==1{INDEX++} INDEX==1{gsub("chr","",$2); if(FNR==1){STR="Name" OFS "Chr" OFS "Position"} else {STR=$1 OFS $2 OFS ($3+$4-1)/2}; for(x=5;x<=NF;x+=5){if(FNR==1){VAL=$x} else{VAL=0.5} STR=STR OFS VAL}  print STR; oNF=NF} INDEX==2{gsub("chr","",$2); STR=$1 OFS $2 OFS $3; for(x=5;x<=oNF;x+=5) {STR=STR OFS 0} print STR}' <(paste alignment_$samplename.sort.dedup.txt) $homozygous | tr '\r' ',' | sed 's/,\t/\t/g' | bzip2 > BAF_${samplename}.txt.bz2

awk 'BEGIN{FS="\t";OFS="\t"} FNR==1{INDEX++} INDEX==1{gsub("chr","",$2); if(FNR==1){STR="Name" OFS "Chr" OFS "Position"} else {STR=$1 OFS $2 OFS ($3+$4-1)/2}; for(x=5;x<=NF;x+=5){STR=STR OFS $x}  print STR; oNF=NF} INDEX==2{gsub("chr","",$2); STR=$1 OFS $2 OFS $3; for(x=5;x<=oNF;x+=5) {STR=STR OFS 0} print STR}' <(paste alignment_$samplename.sort.dedup.txt) $homozygous | tr '\r' ',' | sed 's/,\t/\t/g' | bzip2 > logR_${samplename}.txt.bz2


# Run Ascat
echo -e "${Vertfonce}Ascat${bleuclair}"

mkdir $5/Ascat_output
mkdir $5/Plot

echo `date`" - Alignment and qDNAseq have been Finished!"
