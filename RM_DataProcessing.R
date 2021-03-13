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
suppressMessages(library("recombinator"))
suppressMessages(library("xlsx"))

##############################
## Parsing the arguments    ## 
## passed from command line ##
############################## 
option_list = list(
    make_option(c("-f","--file"), type='character', action="store_true",
    default=NA, help="Output from RepeatMasker"),
    #make_option(c("-f","--fill"), type="character", default=NULL,
    #help="RepeatMasker file or Variation file", metavar="character")
    ##make_option(c("-t","--fill"), type='character', action="store_true",
    ##default=NULL, help="Output from Masker"),
    make_option(c("-r","--rep"), type='character', action="store_true",
    default=NULL, help="Repetition Name to Identify Use FullLengthERV soloLTR SINE or LINE"))
    ##make_option(c("-f","--filetype"),type='character',action="store",default=NULL,
    ##help="RepeatMasker Genome file or Variation file",metavar="character"))

opt_parser=OptionParser(option_list=option_list)
opt=parse_args(opt_parser,args=c("--file","--rep"))

if(is.null(opt$file)){
  print_help(opt_parser)
  stop("No arguments entered. Please enter arguments to continue. Enter --help for the list of arguments", call.=FALSE)
}

#print(sum(unlist(opt)))
#print(opt_parser)
#####################################
## Functions to identify rep class ##
#####################################

# FullLengthERV from RepeatMasker Genome file
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

# FullLengthERV for RepeatMasker Variation file
Full.length.ERV <- function(rm_data, cbp_window, LTR_diff){
final_list <- fullLen.ERV <- NULL
chr_coord <- table(rm_data$query_sequence)[which(table(rm_data$query_sequence) > 1)]
for(i in names(chr_coord)){
   sub_filt <- rm_data  %>% filter(query_sequence == i)
   val=0
   for(j in 1:nrow(sub_filt)){
     if(j>val & sub_filt[j,'class']=="LTR"){
       if(abs(sub_filt[j,12]) < cbp_window & abs(sub_filt[j,14]) < cbp_window){
         vals <- which((sub_filt[,10] == sub_filt[j,10]) & (abs(sub_filt[,13]-sub_filt[j,13]) < LTR_diff))
         for (k in 1:length(vals)) {
           if((j < vals[k]) & k>1){
             semi <- sub_filt[j:vals[k],]
             if(length(unique(semi$s_s)) == 1 & all(semi$class == "LTR")){
               if(k!=length(vals)){
              next
               }
               else if(k==length(vals)){
                 fullLen.ERV <- rbind(fullLen.ERV,semi)
                 val=vals[k]
                 semi_start <- as.numeric(strsplit(semi$query_sequence,"[_:.-]")[[1]][2])+semi[1,6]
                 semi_end <- as.numeric(strsplit(semi$query_sequence,"[_:.-]")[[1]][2])+semi[nrow(semi),7]
                 temp <- paste0(data.frame(table(semi$matching_repeat))[,1],":",data.frame(table(semi$matching_repeat))[,2],collapse=" ; ")
                 final_list = as.data.frame(rbind(final_list,cbind(i, paste(unique(semi$matching_repeat),collapse = ","),temp,semi_start,semi_end)))
                 break
               }
             }
             else if (k>2){
                semi <- sub_filt[j:vals[k-1],]
                fullLen.ERV <- rbind(fullLen.ERV,sub_filt[j:vals[k-1],])
                semi_start <- as.numeric(strsplit(semi$query_sequence,"[_:.-]")[[1]][2])+semi[1,6]
                semi_end <- as.numeric(strsplit(semi$query_sequence,"[_:.-]")[[1]][2])+semi[nrow(semi),7]
                temp <- paste0(data.frame(table(semi$matching_repeat))[,1],":",data.frame(table(semi$matching_repeat))[,2],collapse=" ; ")
                final_list = as.data.frame(rbind(final_list,cbind(i, paste(unique(semi$matching_repeat),collapse = ","),temp,semi_start,semi_end)))
                val=vals[k]
                break
          }
          }}}}}}
colnames(final_list) <- c("QuerySequence","Matching_Repeats","Repeat_frequency","ERV_Start","ERV_End")
ERV_Freq <- as.data.frame(table(final_list[,2]))
colnames(ERV_Freq)=c("ERV_Name","Freq")
ERV_Freq <- ERV_Freq %>% arrange(desc(Freq))
return(list(fullLen.ERV,final_list,ERV_Freq))
}


# soloLTR
solo.LTR <- function(rm_data, cbp_window, Full.length.ERV){
solo.LTR_Filtered <- anti_join(rm_data,Full.length.ERV) 
solo.LTR_list <- solo.final_list <- NULL
#sltr <- table(ltr$query_sequence)[which(table(ltr$query_sequence) == 1)]
sltr <- unique(solo.LTR_Filtered$query_sequence)
for(m in sltr){
  sub <- solo.LTR_Filtered %>% filter(query_sequence == m)
  if(nrow(sub) == 1 && sub[1,16] == "LTR"){
    if(abs(sub[1,12]) < cbp_window & abs(sub[1,14]) < cbp_window){
          solo.LTR_list <- rbind(solo.LTR_list,sub)
          semi_start <- as.numeric(strsplit(sub$query_sequence,"[_:.-]")[[1]][2])+sub[1,6]
          semi_end <- as.numeric(strsplit(sub$query_sequence,"[_:.-]")[[1]][2])+sub[nrow(sub),7]
          #temp <- paste0(data.frame(table(sub$matching_repeat))[,1],":",data.frame(table(sub$matching_repeat))[,2],collapse=" ; ")
          solo.final_list = as.data.frame(rbind(solo.final_list,cbind(m, sub[1,10],semi_start,semi_end)))
    }
  }
  else if(nrow(sub) > 1){
    for(n in 1:nrow(sub)){
      #if(sub[n,16] == "LTR"){
      if((abs(sub[n,12]) < cbp_window && abs(sub[n,14]) < cbp_window)&& (sub[n,16] == "LTR")){
        solo.LTR_list <- rbind(solo.LTR_list,sub[n,])
        semi_start <- as.numeric(strsplit(sub[n,5],"[_:.-]")[[1]][2])+sub[n,6]
        semi_end <- as.numeric(strsplit(sub[n,5],"[_:.-]")[[1]][2])+sub[n,7]
        #temp <- paste0(data.frame(table(sub[n,10]))[,1],":",data.frame(table(sub[n,10]))[,2],collapse=" ; ")
        solo.final_list = as.data.frame(rbind(solo.final_list,cbind(m,sub[n,10],semi_start,semi_end)))
      #}
       # break #Use this break to exit the query as soon as it identifies a solo LTR
    }
    }
  }
}
  colnames(solo.final_list) <- c("QuerySequence","Matching_Repeats","LTR_Start","LTR_End")
  LTR_Freq <- as.data.frame(table(solo.final_list[,2]))
  colnames(LTR_Freq)=c("LTR_Name","Freq")
  LTR_Freq <- LTR_Freq %>% arrange(desc(Freq))
  return(list(solo.LTR_list,solo.final_list,LTR_Freq))
}

##################################
## Data Preprocessing Functions ##
##################################
file.out <- function(file0){
file <- strsplit(noquote(file0[2:nrow(file0),]),"[  ]+")
file <- lapply(file,function(x){if(x[1]=='') c(x[-1], x[1]) else x})
file <- as.data.frame(do.call(rbind,lapply(file,as.vector)))[,c(1:15)]
colnames(file) <- c("SW_score","perc_div.","perc_del.","perc_ins.","query_sequence","position_begin_query","in_end_query","query_(left)","s_s","matching_repeat","repeat_class/family" ,"position_begin","in_end","repeat_(left)","ID" )
#file[,c(8,12,14)] <- as.data.frame(apply(file[,c(8,12,14)],2,function(x){gsub("[()]","",x)}))
#file[,c(1:4,6:8,12:15)] <- as.data.frame(apply(file[,c(1:4,6:8,12:15)], 2, function(x){x= as.numeric(gsub("[()]","",x))}))
file[,c(1:4,6:8,12:15)] <- as.data.frame(apply(file[,c(1:4,6:8,12:15)], 2, function(x){x=gsub("[)]","",x); x=as.numeric(gsub("[()]","-",x))}))
file$class <- unlist(lapply(strsplit(file$`repeat_class/family`,"/"),'[[', 1))
file <- file[mixedorder(file$query_sequence),]
return(file)
}

file.excel <- function(file0){
colnames(file0) <- paste(file0[1,],file0[2,],sep="_")
colnames(file0)[6:7] <- paste0(colnames(file0)[6:7],"_query")
file0 <- file0[-c(1,2),]
file0[,c(1:4,6:8,12:15)] <- as.data.frame(apply(file0[,c(1:4,6:8,12:15)], 2, function(x){as.numeric(x)}))
file0$class <- unlist(lapply(strsplit(file0$`repeat_class/family`,"/"),'[[', 1))
file <- file0[mixedorder(file0$query_sequence),-16]
return(file)
}




 
