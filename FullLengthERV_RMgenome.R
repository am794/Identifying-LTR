#!/usr/bin/env Rscript

########################
## Load the libraries ##
########################
wd=getwd()
library("openxlsx")
library("dplyr")
library("tidyr")
library("gtools")
library("knitr")
library("kableExtra")
library("yaml")
library("rmarkdown")
library("recombinator",lib.loc=paste0(wd))
library("xlsx")
library("gtools")

######################
## Read the files   ##
## Data Preparation ##
######################
file <- commandArgs()
pwd <- getwd()
mm10.genome <- na.omit(as.data.frame(read.delim2(paste0(pwd,"/InputFiles/",file[6]),sep="\t",stringsAsFactors=default.stringsAsFactors())))
attach(mm10.genome)
## mm10.genome <- na.omit(as.data.frame(read.delim2(paste0(pwd,"/InputFiles/repeatmaskermm10wholegenome13Nov2020"),sep="\t",stringsAsFactors=default.stringsAsFactors())))
#mm10.genome$Chr <- unlist(lapply(strsplit(as.character(mm10.genome$genoName),"_"), function(x) x[1]))
#mm10.genome <- mm10.genome %>% filter(!repClass %in% c("Simple_repeat","Low_complexity","Other","Unknown")) %>% arrange(genoStart,genoEnd)
mm10.genome$GeneType <- unlist(lapply(strsplit(as.character(mm10.genome$genoName),"_"), function(x) if(length(x) > 1) "AlternateGene" else "Gene"))
mm10.genome <- as.data.frame(mm10.genome %>% filter(GeneType == "Gene" & genoName != "chrM"))[,-18]
mm10.genome$repClass <- gsub(pattern="[?]",replacement='',mm10.genome$repClass)
#mm10.genome$query_sequence <- paste0(mm10.genome$genoName,":",mm10.genome$genoStart,"-",mm10.genome$genoEnd)
mm10.genome$query_sequence <- paste0(mm10.genome$genoStart,"-",mm10.genome$genoEnd)

########################################
## Function to identify FullLengthERV ##
########################################
genome.FullLength.ERV <- function(genome_table,cbp_window,LTR_diff){
   final_list <- fullLen.ERV <- NULL
   val=0
   for(j in 1:nrow(genome_table)){
     if(j>val & genome_table[j,'repClass']=="LTR"){
       if(abs(genome_table[j,14]) < cbp_window & abs(genome_table[j,16]) < cbp_window){
         vals <- which((genome_table[,11] == genome_table[j,11]) & (abs(genome_table[,15]-genome_table[j,15]) < LTR_diff))
         if(vals[1] != j){
           index=which(vals==j)
           vals=vals[-(1:index-1)]
         }
         for (k in 1:length(vals)) {
           if((j < vals[k]) & k>1){
             semi <- genome_table[j:vals[k],]
             semi$repName <- factor(semi$repName)
             if(length(unique(semi$strand)) == 1 & all(semi$repClass == "LTR")){
               if(k!=length(vals)){
              next
               }
               else if(k==length(vals)){
                 fullLen.ERV <- rbind(fullLen.ERV,semi)
                 val=vals[k]
                 #semi_start <- as.numeric(strsplit(semi$query_sequence,"[_:.-]")[[1]][2])+semi[1,6]
                 #semi_end <- as.numeric(strsplit(semi$query_sequence,"[_:.-]")[[1]][2])+semi[nrow(semi),7]
                 semi_start <- semi[1,7]
                 semi_end <- semi[nrow(semi),8]
                 temp <- paste0(data.frame(table(semi$repName))[,1],":",data.frame(table(semi$repName))[,2],collapse=" ; ")
                 final_list = as.data.frame(rbind(final_list,cbind(paste(unique(semi$repName),collapse = ","),temp,semi_start,semi_end)))
                 break
               }
             }
             else if (k>2){
                semi <- genome_table[j:vals[k-1],]
                semi$repName <- factor(semi$repName)
                fullLen.ERV <- rbind(fullLen.ERV,semi)
                #semi_start <- as.numeric(strsplit(semi$query_sequence,"[_:.-]")[[1]][2])+semi[1,6]
                #semi_end <- as.numeric(strsplit(semi$query_sequence,"[_:.-]")[[1]][2])+semi[nrow(semi),7]
                semi_start <- semi[1,7]
                semi_end <- semi[nrow(semi),8]
                temp <- paste0(data.frame(table(semi$repName))[,1],":",data.frame(table(semi$repName))[,2],collapse=" ; ")
                final_list = as.data.frame(rbind(final_list,cbind(paste(unique(semi$repName),collapse = ","),temp,semi_start,semi_end)))
                val=vals[k-1]
                break
             }
             else{
               break
             }
           }}}}}
if(!is.null(final_list)){
  colnames(final_list) <- c("Matching_Repeats","Repeat_frequency","ERV_Start","ERV_End")
  ERV_Freq <- as.data.frame(table(final_list[,1]))
  colnames(ERV_Freq)=c("ERV_Name","Freq")
  ERV_Freq <- ERV_Freq %>% arrange(desc(Freq))
  return(list(fullLen.ERV,final_list,ERV_Freq))
#return(final_list)
}
else{
 return(NULL)
}
}

#######################################
## Function for identifying solo LTR ##
#######################################
solo.LTR <- function(rm_data, cbp_window, Full.length.ERV){
solo.LTR_Filtered <- anti_join(rm_data,Full.length.ERV) 
solo.LTR_list <- solo.final_list <- NULL
#sltr <- table(ltr$query_sequence)[which(table(ltr$query_sequence) == 1)]
sltr <- unique(solo.LTR_Filtered$query_sequence)
for(m in sltr){
  sub <- solo.LTR_Filtered %>% filter(query_sequence == m)
  if(nrow(sub) == 1 && sub[1,12] == "LTR"){
    if(abs(sub[1,14]) < cbp_window & abs(sub[1,16]) < cbp_window){
          solo.LTR_list <- rbind(solo.LTR_list,sub)
          #semi_start <- as.numeric(strsplit(sub$query_sequence,"[_:.-]")[[1]][2])+sub[1,7]
          #semi_end <- as.numeric(strsplit(sub$query_sequence,"[_:.-]")[[1]][2])+sub[nrow(sub),8]
          #temp <- paste0(data.frame(table(sub$matching_repeat))[,1],":",data.frame(table(sub$matching_repeat))[,2],collapse=" ; ")
          semi_start <- sub[1,7]
          semi_end <- sub[nrow(sub),8]
          solo.final_list = as.data.frame(rbind(solo.final_list,cbind(as.character(sub[1,6]),m,as.character(sub[1,11]),semi_start,semi_end)))
    }
  }
  else if(nrow(sub) > 1){
    for(n in 1:nrow(sub)){
      #if(sub[n,16] == "LTR"){
      if((abs(sub[n,14]) < cbp_window && abs(sub[n,16]) < cbp_window)&& (sub[n,12] == "LTR")){
        solo.LTR_list <- rbind(solo.LTR_list,sub[n,])
        #semi_start <- as.numeric(strsplit(sub[n,18],"[_:.-]")[[1]][2])+sub[n,7]
        #semi_end <- as.numeric(strsplit(sub[n,18],"[_:.-]")[[1]][2])+sub[n,8]
        #temp <- paste0(data.frame(table(sub[n,10]))[,1],":",data.frame(table(sub[n,10]))[,2],collapse=" ; ")
        semi_start <- sub[n,7]
        semi_end <- sub[n,8]
        solo.final_list = as.data.frame(rbind(solo.final_list,cbind(as.character(sub[1,6]),m,as.character(sub[n,11]),semi_start,semi_end)))
      #}
       # break #Use this break to exit the query as soon as it identifies a solo LTR
    }
    }
  }
}
  colnames(solo.final_list) <- c("Chromosome","QuerySequence","Matching_Repeats","LTR_Start","LTR_End")
  LTR_Freq <- as.data.frame(table(solo.final_list[,2]))
  colnames(LTR_Freq)=c("LTR_Name","Freq")
  LTR_Freq <- LTR_Freq %>% arrange(desc(Freq))
  return(list(solo.LTR_list,solo.final_list,LTR_Freq))
}

##################################################
## Using the above functions for each chromosome ##
##################################################
FullLengthERV <- c()
Frequency <- c()
Summary <- c()
system(paste0("mkdir ",pwd,"/OutputFiles"))
for (i in unique(mm10.genome$genoName)) {
chr <- mm10.genome %>% filter(genoName == i & !repClass %in% c("Simple_repeat","Low_complexity","Other","Unknown")) %>% arrange(genoStart,genoEnd)
mm10Genome.FullLength <- genome.FullLength.ERV(chr,25,1000)
if(!is.null(mm10Genome.FullLength)){
  FullLengthERV <- rbind(FullLengthERV,mm10Genome.FullLength[[1]])
  Frequency <- rbind(Frequency,cbind(mm10Genome.FullLength[[2]],rep(paste0(unique(chr$genoName)),nrow(mm10Genome.FullLength[[2]]))))
  Summary <- rbind(Summary,cbind(mm10Genome.FullLength[[3]],rep(paste0(unique(chr$genoName)),nrow(mm10Genome.FullLength[[3]]))))
  #write.xlsx(mm10Genome.FullLength[[1]],paste0(pwd,"/",unique(chr$Chr),"_FullLengthERVs.xlsx"),row.names=FALSE)
  write.table(mm10Genome.FullLength[[1]],paste0(pwd,"/OutputFiles/",unique(chr$genoName),"_FullLengthERVs.txt"),row.names=FALSE,quote=FALSE,sep="\t")
  write.table(mm10Genome.FullLength[[2]],paste0(pwd,"/OutputFiles/",unique(chr$genoName),"_Frequency.txt"),row.names=FALSE,quote=FALSE,sep="\t")
  write.table(mm10Genome.FullLength[[3]],paste0(pwd,"/OutputFiles/",unique(chr$genoName),"_Summary.txt"),row.names=FALSE,quote=FALSE,sep="\t")
}else{
next
}}
levels(Frequency$Chromosome) <- mixedsort(levels(Frequency$Chromosome))
bed <- Frequency[,c(5,3,4,1)]
colnames(bed) <- c("chrom","chromStart","chromEnd","name")
colnames(Frequency)[ncol(Frequency)] <- "Chromosome"
colnames(Summary)[ncol(Summary)] <- "Chromosome"
write.table(FullLengthERV,paste0(pwd,"/OutputFiles/FullLengthERVs.txt"),row.names=FALSE,quote=FALSE,sep="\t")
write.table(Frequency,paste0(pwd,"/OutputFiles/Frequency.txt"),row.names=FALSE,quote=FALSE,sep="\t")
write.table(Summary,paste0(pwd,"/OutputFiles/Summary.txt"),row.names=FALSE,quote=FALSE,sep="\t")
write.table(bed,paste0(pwd,"/OutputFiles/FullLengthERV.bed"),row.names=FALSE,quote=FALSE,sep="\t")
GenomeWideFreq <- table(Frequency$Matching_Repeats,Frequency$Chromosome)
GenomeWideFreq <- tibble::rownames_to_column(as.data.frame.matrix(GenomeWideFreq),"Repeats")
GenomeWideFreq$TotalRepeats <- rowSums(GenomeWideFreq[,-1])
write.table(GenomeWideFreq,paste0(pwd,"/OutputFiles/GenomeWideFrequencies.txt"),row.names=FALSE,quote=FALSE,sep="\t")


################################################
## Using soloLTR function for each chromosome ##
################################################

##system(paste0("if [! -d ",pwd,"/OutputFiles ]"))
#system(paste0("mkdir ",pwd,"/OutputFiles_soloLTR"))
#soloLTR <- c()
#Frequency <- c()
#for (i in unique(mm10.genome$genoName)){
#chr <- mm10.genome %>% filter(genoName == i & repClass %in% c("LTR")) %>% arrange(genoStart,genoEnd)
#fullLen <- read.table(paste0(pwd,"/OutputFiles/",i,"_FullLengthERVs.txt"),sep="\t",header=TRUE)
#fullLen$query_sequence <- factor(unlist(lapply(fullLen$query_sequence,function(x){unlist(strsplit(as.character(x),":"))[2]})))
#mm10Genome.soloLTR <- solo.LTR(chr,25,fullLen)
#if(!is.null(mm10Genome.soloLTR)){
# soloLTR <- rbind(soloLTR,mm10Genome.soloLTR[[1]])
# Frequency <- rbing(Frequency,mm10Genome.soloLTR[[2]])
# write.table(mm10Genome.soloLTR[[1]],paste0(pwd,"/OutputFiles_soloLTR/",unique(chr$genoName),"_soloLTR.txt"),row.names=FALSE,quote=FALSE,sep="\t")
# write.table(mm10Genome.soloLTR[[1]],paste0(pwd,"/OutputFiles_soloLTR/",unique(chr$genoName),"_Frequency.txt"),row.names=FALSE,quote=FALSE,sep="\t")
#}else{
#next
#}}
#levels(Frequency$Chromosome) <- mixedsort(levels(Frequency$Chromosome))
#bed <- Frequency[,c(1,4,5,3)]
#colnames(bed) <- c("chrom","chromStart","chromEnd","name")
#write.table(soloLTR,paste0(pwd,"/OutputFiles_soloLTR/soloLTR.txt"),row.names=FALSE,quote=FALSE,sep="\t")
#write.table(Frequency,paste0(pwd,"/OutputFiles_soloLTR/Frequency.txt"),row.names=FALSE,quote=FALSE,sep="\t")
#write.table(bed,paste0(pwd,"/OutputFiles_soloLTR/soloLTR.bed"),row.names=FALSE,quote=FALSE,sep="\t")
#GenomeWideFreq <- table(Frequency$Matching_Repeats,Frequency$Chromosome)
#GenomeWideFreq <- tibble::rownames_to_column(as.data.frame.matrix(GenomeWideFreq),"Repeats")
#GenomeWideFreq$TotalRepeats <- rowSums(GenomeWideFreq[,-1])
#write.table(GenomeWideFreq,paste0(pwd,"/OutputFiles_soloLTR/GenomeWideFrequencies.txt"),row.names=FALSE,quote=FALSE,sep="\t")





