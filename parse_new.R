#!/usr/bin/env Rscript

########################
## Load the libraries ##
########################
wd=getwd()
suppressMessages(library("optparse"))
suppressMessages(library("openxlsx"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyr"))
suppressMessages(library("gtools"))
suppressMessages(library("knitr"))
suppressMessages(library("kableExtra"))
suppressMessages(library("yaml"))
suppressMessages(library("rmarkdown"))
suppressMessages(library("recombinator",lib.loc=wd))
suppressMessages(library("xlsx"))

source(paste0(wd,RepeatsFunctions.R))

##############################
## Parsing the arguments    ## 
## passed from command line ##
############################## 
option_list <- list(
    make_option(c("-f", "--file"), action="store", default=NULL,type="character",
        help="Input file to identify Repeats. Should be an output from RepeatMasker."),
    make_option(c("-c","--chrfile"),action="store",default="true",
        type="character", help="Print output files for each chromosome separately. Yes will output the files for each chromosome along with the merged file. No will output just the merged file"),
    make_option(c("-r", "--reptype"), action="store", default="FullLengthLTR",
        type="character", help="Repeat type to identify. Accepted arguments: FullLengthLTR, soloLTR, FullLengthLTR_soloLTR,LINE, SINE. Default is FullLengthLTR"),
    make_option(c("-t", "--RMtype"), type="character", default=NULL, action="store",
        help="Output from RepeatMasker. Genome file or Variations file. Accepted arguments: genome, variation"),
    make_option(c("-bp", "--basepairs"), type="integer", default=25, action="store",
        help="Number of basepairs in the start and end. Default is 25"),    
    make_option(c("-d", "--diffbp"), type="integer", default=1000, action="store",
        help="Difference between the lengths of LTR begin and LTR end in full length LTR. Default is 1000.")
)

opt <- parse_args(OptionParser(usage="%prog [options] file",option_list=option_list),positional_arguments=TRUE)


#`%notin%` <- Negate(`%in%`)

print(opt$options)

if(is.null(opt$options$RMtype)&is.null(opt$options$file)){
print("please enter RMtype and an input file or --help for options")
}else if(is.null(opt$options$RMtype)&!is.null(opt$options$file)){
print("please enter RMtype or --help for options")
}else if(!is.null(opt$options$RMtype)&is.null(opt$options$file)){
print("please enter an input file or --help for options")
}


#if(opt$options$RMtype == "genome")
#{
#  Genome <- na.omit(as.data.frame(read.delim2(paste0(getwd(),"/",file),sep="\t",stringsAsFactors=default.stringsAsFactors())))
#  attach(Genome)
#  Genome$GeneType <- unlist(lapply(strsplit(as.character(Genome$genoName),"_"), function(x) if(length(x) > 1) "AlternateGene" else "Gene"))
#  Genome <- as.data.frame(Genome %>% filter(GeneType == "Gene" & genoName != "chrM"))[,-18]
#  Genome$repClass <- gsub(pattern="[?]",replacement='',FullLength.Genome$repClass)
#  Genome$query_sequence <- paste0(FullLength.Genome$genoName,":",FullLength.Genome$genoStart,"-",FullLength.Genome$genoEnd)
#  if(opt$options$reptype=="FullLengthLTR"){
#   FullLengthERV <- c()
#     Frequency <- c()
#     Summary <- c()
#     system(paste0("mkdir ",pwd,"/OutputFiles_FullLengthLTR"))
#      for (i in unique(Genome$genoName)) {
#       chr <- Genome %>% filter(genoName == i & !repClass %in% c("Simple_repeat","Low_complexity","Other","Unknown")) %>% arrange(genoStart,genoEnd)
#       Genome.FullLength <- genome.FullLength.ERV(chr,opt$options$basepairs,opt$options$diffbp)
#       if(!is.null(Genome.FullLength)){
#          FullLengthERV <- rbind(FullLengthERV,Genome.FullLength[[1]])
#          Frequency <- rbind(Frequency,cbind(Genome.FullLength[[2]],rep(paste0(unique(chr$genoName)),nrow(Genome.FullLength[[2]]))))
#          Summary <- rbind(Summary,cbind(Genome.FullLength[[3]],rep(paste0(unique(chr$genoName)),nrow(Genome.FullLength[[3]]))))
#          if((opt$options$chrfile=="True") | (opt$options$chrfile=="T") | (opt$options$chrfile=="true") | (opt$options$chrfile=="TRUE") | (opt$options$chrfile=="t")){
#          write.table(Genome.FullLength[[1]],paste0(pwd,"/OutputFiles_FullLengthLTR/",unique(chr$genoName),"_FullLengthERVs.txt"),row.names=FALSE,quote=FALSE,sep="\t")
#          write.table(Genome.FullLength[[2]],paste0(pwd,"/OutputFiles_FullLengthLTR/",unique(chr$genoName),"_Frequency.txt"),row.names=FALSE,quote=FALSE,sep="\t")
#          write.table(Genome.FullLength[[3]],paste0(pwd,"/OutputFiles_FullLengthLTR/",unique(chr$genoName),"_Summary.txt"),row.names=FALSE,quote=FALSE,sep="\t")
#          }
#        }else{
#          next
#        }
#      }
#       levels(Frequency$Chromosome) <- mixedsort(levels(Frequency$Chromosome))
#       bed <- Frequency[,c(3,4,5,1)]
#       colnames(bed) <- c("chromStart","chromEnd","chrom","name")
#       colnames(Frequency)[ncol(Frequency)] <- "Chromosome"
#       colnames(Summary)[ncol(Summary)] <- "Chromosome"
#       write.table(FullLengthERV,paste0(pwd,"/OutputFiles_FullLengthLTR/FullLengthERVs.txt"),row.names=FALSE,quote=FALSE,sep="\t")
#       write.table(Frequency,paste0(pwd,"/OutputFiles_FullLengthLTR/Frequency.txt"),row.names=FALSE,quote=FALSE,sep="\t")
#       write.table(Summary,paste0(pwd,"/OutputFiles_FullLengthLTR/Summary.txt"),row.names=FALSE,quote=FALSE,sep="\t")
#       write.table(bed,paste0(pwd,"/OutputFiles_FullLengthLTR/FullLengthERV.bed"),row.names=FALSE,quote=FALSE,sep="\t")
#       GenomeWideFreq <- table(Frequency$Matching_Repeats,Frequency$Chromosome)
#       GenomeWideFreq <- tibble::rownames_to_column(as.data.frame.matrix(GenomeWideFreq),"Repeats")
#       GenomeWideFreq$TotalRepeats <- rowSums(GenomeWideFreq[,-1])
#       write.table(GenomeWideFreq,paste0(pwd,"/OutputFiles_FullLengthLTR/GenomeWideFrequencies.txt"),row.names=FALSE,quote=FALSE,sep="\t")
#    }else if(opt$options$reptype=="soloLTR"){
#       
#    }
#}






