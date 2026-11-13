####################################################################################################################
###                                               Start_Up                                                      ####
####################################################################################################################
#module purge; module load r/4.4.3-k6zf4xz; R

rm(list=ls())

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])

if (length(script_path) == 0) {
  stop("Unable to determine script path")
}

script_path <- normalizePath(script_path, mustWork = TRUE)
script_dir <- dirname(script_path)
parent_dir <- dirname(script_dir)

source(file.path(parent_dir, "parameters.R"))

###   Directory creation  ### 
dir <- PROJECT_DIR
setwd(dir)

print(paste0("Directory Set to ", dir))

###      Packages         ### 

if(!require(tidyverse)){
    install.packages("tidyverse")
    library(tidyverse)
}
if(!require(data.table)){
    install.packages("data.table")
    library(data.table)
}
if(!require(ggraph)){
    install.packages("ggraph")
    library(ggraph)
}
print("Modules Loaded")

####################################################################################################################
###                                              ggplot theme		                                            ####
####################################################################################################################

ggplot_theme_basic <- theme(
					panel.grid = element_blank(),
					plot.title = element_text(size=30, face="bold",hjust = 0.5),
					axis.text.x = element_text(size=24, margin = margin(t = 10, r = 0, b = 0, l = 0)),
					axis.text.y  = element_text(size=24, margin = margin(t = 10, r = 0, b = 0, l = 0),hjust = 0),
					axis.title.y = element_text(size=28, margin = margin(t = 0, r = 20, b = 0, l = 0)),
					axis.title.x = element_text(size=28, margin = margin(t = 20, r = 0, b = 0, l = 0)),
					plot.margin = margin(t = 20, r = 20, b = 10, l = 10, unit = "pt"),
					axis.ticks=element_blank(),
					panel.spacing = unit(1, "lines"),
					strip.text = element_text(angle=0, size=28),
					plot.background = element_rect(fill = 'white', colour = 'white'),
					panel.border = element_rect(fill=NA, colour = "black"),
					panel.background = element_rect(fill = 'white', colour = 'white',linetype="solid"),
					strip.background = element_rect(fill="white",colour="white",linetype="solid"),
					legend.background = element_rect(fill = "white"),
					legend.title = element_text(color = "black", size = 30, margin = margin(t = 0, r = 0, b = 0, l = 0)),
					legend.text = element_text(color = "black", size = 28),
					#legend.position = "none" ,
					legend.key = element_rect(fill = "white",colour="white"),
					legend.margin = margin(t = 0, r = 20, b = 0, l = 20),
					axis.line = element_blank()
				)

####################################################################################################################
###                                               Loadin Data	                                                ####
####################################################################################################################

fqd <- fread(file = paste0(T_META_DIR,"/Sequenced_Samples_Info.tsv"), sep = '\t')
sample_info <- fqd %>% filter(Capture == "Genomic" & Sequence_Platform == "Illumina") %>% select(Individual_ID,Father_ID,Mother_ID,Sex,Affection_Status,Project_Name,Date_sequenced) %>% unique()

fastqscreen <- fread(file = paste0(T_SEQ_DIR,"/hg38/QualityChecks/FastQScreen/Genomic/MultiQC_Total/multiqc_data/multiqc_fastq_screen.txt"), sep = '\t')
fastqscreen <- fastqscreen %>% filter(!is.na(`Oral counts`)) %>% mutate(Individual_ID = gsub("_.*","",Sample))
sample_contam <- fastqscreen %>% select(Individual_ID,total_reads,"Human counts","Oral counts","No hits counts") %>% group_by(Individual_ID) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
	mutate(estimated_bacteria_percentage=(`Oral counts`+`No hits counts`)/total_reads) %>% mutate(oral_percentage=`Oral counts`/total_reads) %>% mutate(human_percentage=`Human counts`/total_reads) %>%
	mutate(Unmapped_percentage=`No hits counts`/total_reads)
  
sample_full <- left_join(sample_contam,sample_info,by="Individual_ID")

write.table(sample_full, file=paste0("All_Samples_Summary.tsv"), quote=FALSE, sep='\t', row.names=F)

####################################################################################################################
###                                               Loadin Data	                                                ####
####################################################################################################################

quit()

