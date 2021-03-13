#!/usr/bin/env Rscript

#######################################################
# Script prints mean, median and standard dev of      #
# millidiv column for each repname in each chromosome #
#######################################################
wd=getwd()
library(dplyr)
solo <- read.table(paste0(wd,"OutputFiles_soloLTR/soloLTR_2.txt"),header=TRUE)
#solo_summary <- solo %>% group_by(genoName,repName) %>% summarise(mean = mean(milliDiv), n = n(),std_div=sd(milliDiv),median_val=median(milliDiv))
solo_summary <- solo %>% group_by(repName) %>% 
  summarise(Mean = mean(milliDiv),Std.dev=sd(milliDiv),
            Median=median(milliDiv),Range_Min=min(milliDiv),
            Range_Max=max(milliDiv), Frequency = n()) %>%
  arrange(desc(Frequency))

my_list <- c("IAPLTR2a", "IAPLTR1_Mm", "MTA_Mm", "RLTR10", "IAPLTR1a_Mm",
             "IAPLTR2b","IAPLTR2a2_Mm","IAPLTR2_Mm","ERVB7_1-LTR_MM","MT2_Mm")

my_plots <- function(x){
  df <- solo %>% filter(solo$repName == x)
  ggplot(df, aes(x=milliDiv))+geom_histogram(binwidth = 5)+ggtitle(paste0(x))+
    theme(axis.text.x = element_text(face="bold",  size=18),
          axis.text.y = element_text(face="bold",  size=18),
          axis.title=element_text(size=18,face="bold"),
          plot.title = element_text(size = 20, face = "bold"))
  ggsave(paste0(wd,"OutputFiles_soloLTR/Histogram_",x,".jpeg"))
  ggsave(paste0(wd,"OutputFiles_soloLTR/Histogram_",x,".png"))
  ggsave(paste0(wd,"OutputFiles_soloLTR/Histogram_",x,".pdf"))
} 
lapply(my_list, my_plots)

write.table(as.data.frame(solo_summary),paste0(wd,"OutputFiles_soloLTR/soloLTR_summaryStats.txt"),row.names=FALSE,quote=FALSE,sep="\t")

