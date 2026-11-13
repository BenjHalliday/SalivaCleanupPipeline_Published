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
topdir <- paste0(PROJECT_DIR,"/Analyses/")
dir <- paste0(topdir,"Compute_Comparison")
Dirs <- c("graphs","data_clean","data_raw")
ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir)),TRUE)
setwd(dir)
for(subDir in Dirs) {
    ifelse(!dir.exists(file.path(dir, subDir)), dir.create(file.path(dir, subDir)), TRUE)
}

print(paste0("Directory Set to ", dir))

###      Packages         ### 

if (!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)

if (!require(data.table)) install.packages("data.table")
library(data.table)

if (!require(ggraph)) install.packages("ggraph")
library(ggraph)

print("Modules Loaded")

####################################################################################################################
###                                              ggplot theme		                                            ####
####################################################################################################################

ggplot_theme_basic <- theme(
	panel.grid = element_blank(),
	plot.title = element_text(size=30, face="bold",hjust = 0.5),
	axis.text.x = element_text(size=24, margin = margin(t = 10, r = 0, b = 0, l = 0)),
	axis.text.y  = element_text(size=24, margin = margin(t = 0, r = 10, b = 0, l = 0),hjust = 0),
	axis.title.x = element_text(size=28, margin = margin(t = 20, r = 0, b = 0, l = 0)),
	axis.title.y = element_text(size=28, margin = margin(t = 0, r = 20, b = 0, l = 0)),
	plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt"),
	axis.ticks=element_blank(),
	panel.spacing = unit(1, "lines"),
	strip.text = element_text(angle=0, size=28),
	plot.background = element_rect(fill = 'white', colour = 'white'),
	panel.border = element_rect(fill=NA, colour = "black", linewidth = 1.5),
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
###                                              Saliva Data Read In                                            ####
####################################################################################################################

samplemap <- paste0(PROJECT_DIR,"/MetaData/Sample_Map.tsv")

stats <- grep("sacct",list.files(path=paste0(dir,"/data_raw"),recursive=F,full.names=F), value = TRUE)
remove <- c("")

negative_intersect <- stats[!(stats %in% remove)]

for( stat in negative_intersect ) {
	file <- paste0(dir,"/data_raw/",stat)
	name <- gsub("_sacct_output.tsv","",stat)
	if(file.exists(file)) {
		sacct <- fread(file,fill=TRUE) 
	
		sample_id <- name %>% gsub("_.*","",.)
		tag <- name %>% gsub(paste0(sample_id,"_NoFastP_"),"",.) %>% gsub("_k.*","",.)
		seedlength <- name %>% gsub(paste0(sample_id,"_"),"",.) %>% gsub(".*_k","",.)

		new_row <- sacct %>% filter(str_detect(JobName,"Aligning")) %>% 
			mutate(Alignment=name) %>%
			mutate(Fastq = gsub("_Aligning_.*","",JobName)) %>%
			mutate(Chunk = gsub(".*_Aligning_","",JobName)) %>%
			mutate(Sample_ID=sample_id) %>%
			mutate(Tag=tag) %>%
			mutate(SeedLength=seedlength)
			
  		assign(paste0("data_",name),new_row)
	} else {
		print(paste0(file," file doesn't exist"))
	}
}
listOfDataFrames <- ls(pattern=paste0("^data_*"))
Combined <- bind_rows(mget(listOfDataFrames)) %>% mutate(Group=ifelse(Tag == "NoKraken","SeedLength",ifelse(Tag == "Krakentaxid","Kraken",ifelse(Tag == "Kraken_trained_taxid","KrakenTrained",Tag)))) %>%
	mutate(SampleType=ifelse(substr(Sample_ID, nchar(Sample_ID), nchar(Sample_ID)) == "S","Saliva",ifelse(substr(Sample_ID, nchar(Sample_ID), nchar(Sample_ID)) == "B","Blood","ERROR"))) %>%
	mutate(Individual_ID = str_sub(Sample_ID, 1, str_length(Sample_ID) - 1))

SampleMap <- fread(file = samplemap, sep = '\t', colClasses = "character") 
Combined <- Combined %>% mutate(Individual_ID_Recoded = recode(Individual_ID, !!!setNames(SampleMap$Paper, SampleMap$Original)))

write.table(Combined, file = "data_clean/Alignment_Speed_Raw.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#note should be State not MaxRSS, columns shifted incorrextly
Combined <- Combined %>% filter(MaxRSS == "COMPLETED") %>%
  mutate(TotalCPU_Seconds = sapply(TotalCPU, function(x) {
    sum(as.numeric(unlist(strsplit(x, ":"))) * c(3600, 60, 1))
  }))

write.table(Combined, file = "data_clean/Alignment_Speed.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################################################################
###                                               Compute Data Read In Seed_Length                              ####
####################################################################################################################

Combined_SeedLength <- Combined %>% filter(Group == "SeedLength") %>% mutate(SeedLength = as.numeric(SeedLength))

Combined_Plot_Time <- Combined_SeedLength %>% group_by(Alignment,SeedLength,Sample_ID,Individual_ID,SampleType,Individual_ID_Recoded) %>% summarise(TotalCPU_Seconds_Alignment = sum(TotalCPU_Seconds))
Combined_Plot_Count <- Combined_SeedLength %>% group_by(Alignment,SeedLength,Sample_ID,Individual_ID,SampleType,Individual_ID_Recoded) %>% summarize(row_count = n(), .groups = "drop")

Combined_Plot <- full_join(Combined_Plot_Time,Combined_Plot_Count)

Combined_Plot <- Combined_Plot %>% mutate(TotalCPU_Hours_Alignment = TotalCPU_Seconds_Alignment / 3600)

Combined_Plot %>% filter(SeedLength == 19) %>% ungroup() %>% 
	summarise( mean = mean(TotalCPU_Hours_Alignment, na.rm = TRUE), sd = sd(TotalCPU_Hours_Alignment, na.rm = TRUE), median = median(TotalCPU_Hours_Alignment, na.rm = TRUE),n = n()) %>%
	write.table("data_clean/Alignment_Speed_SeedLength_Summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

model <- lm(TotalCPU_Hours_Alignment ~ Sample_ID - 1 , data = Combined_Plot)
intercept <- model$coefficients[1]
Combined_Plot$Normalized_TotalCPU_Hours_Alignment <- residuals(model) + intercept

write.table(Combined_Plot, file = "data_clean/Alignment_Speed_SeedLength_Summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

summary(lm(TotalCPU_Hours_Alignment ~ Individual_ID -1 + SeedLength + SampleType,data = Combined_Plot))$coefficients %>% as.data.frame() %>% tibble::rownames_to_column("Term") %>%
	write_tsv(., "data_clean/Alignment_Speed_SeedLength_LinearModel.tsv")
summary(lm(TotalCPU_Hours_Alignment ~ Individual_ID -1 + SeedLength + SampleType,data = Combined_Plot %>% filter(Individual_ID_Recoded != "Sample_3")))$coefficients %>% as.data.frame() %>% tibble::rownames_to_column("Term") %>%
	write_tsv(., "data_clean/Alignment_Speed_SeedLength_LinearModel_HeavyContamDrop.tsv")

plot <-  ggplot(Combined_Plot, aes(x = SeedLength, y = TotalCPU_Hours_Alignment, group = Sample_ID, color = Sample_ID)) +
	geom_line() +
	labs(x = "BWA Alignment Seed Length", y = "Total CPU Time for Alignment (Hours)") +
	scale_y_continuous() +
	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/SeedLength_Compute_Split.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

plot <-  ggplot(Combined_Plot, aes(x = SeedLength, y = Normalized_TotalCPU_Hours_Alignment, group = Sample_ID, color = Sample_ID)) +
	geom_line() +
	labs(x = "BWA Alignment Seed Length", y = "Normalized Total CPU Time for Alignment (Hours)") +
	scale_y_continuous() +
	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/SeedLength_Compute_Split_Normalised.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

plot <- ggplot(Combined_Plot %>% mutate(SeedLength = as.character(SeedLength)), aes(x = SeedLength, y = Normalized_TotalCPU_Hours_Alignment, group = interaction(SeedLength, SampleType), fill = SampleType, colour = SampleType)) +
    geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, position = position_dodge(width = 0.75)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), size = 2, alpha = 0.7) +
    labs(x = "BWA Alignment Seed Length", y = "Normalised Total CPU Time for Alignment (Hours)") +
	scale_colour_manual(values = c("Blood" = "#E63946", "Saliva" = "#1D8CE0")) +
	scale_fill_manual(values = c("Blood" = "#E63946", "Saliva" = "#1D8CE0")) +
    guides(alpha = "none") +
    ggplot_theme_basic
ggsave(plot = plot, filename = "graphs/SeedLength_Compute_Split_Normalised_Boxplot.pdf", device = "pdf", width = 40, height = 30, units = "cm", dpi = 640)

####################################################################################################################
###                                               Fold Change		     	                                    ####
####################################################################################################################

# Choose a reference length 
reference_length <- min(Combined_Plot$SeedLength)

fold_changes <- Combined_Plot %>% group_by(Sample_ID) %>% mutate(
    	Reference_Time = Normalized_TotalCPU_Hours_Alignment[SeedLength == reference_length],
    	Fold_Change = Normalized_TotalCPU_Hours_Alignment / Reference_Time ,
    	Fold_Difference = abs(Fold_Change - 1) + 1
  	) %>%
  	ungroup()

write.table(fold_changes, file = "data_clean/Alignment_Speed_SeedLength_Summary_FoldChange.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

plot <- ggplot(fold_changes, aes(x = SeedLength, y = Fold_Difference, group = Sample_ID, color = Sample_ID)) +
	geom_line() +
	geom_point() +
	labs(x = "BWA Alignment Seed Length", y = "Sample Normalised Fold Reduction in Total CPU Time for Alignment (Hours)") +
	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/SeedLength_Compute_Split_Normalised_FoldChange.pdf"),device="pdf", width = 60, height = 40, units = "cm",dpi = 640)

plot <- ggplot(fold_changes %>% mutate(SeedLength = as.character(SeedLength)), aes(x = SeedLength, y = Fold_Difference, group = SeedLength)) +
	geom_boxplot() +
	#geom_point() +
	labs(x = "BWA Alignment Seed Length", y = "Sample Normalised Fold Reduction in Total CPU Time for Alignment (Hours)") +
	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/SeedLength_Compute_Split_Normalised_FoldChange_Boxplot.pdf"),device="pdf", width = 60, height = 40, units = "cm",dpi = 640)

####################################################################################################################
###                                               Compute Data Read In Cross-Comp                               ####
####################################################################################################################

# Remove Trained Kraken as they were aligned with 16G per chunk instead of 8G per chunk
Combined_Cross <- Combined %>% filter(SeedLength %in% c("19","23","24")) %>% filter(Group != "KrakenTrained") %>% mutate(Group = factor(Group, levels = c("SeedLength", "Decoy", "Kraken")))

Combined_Plot_Time <- Combined_Cross %>% group_by(Alignment,SeedLength,Sample_ID,Individual_ID,Group,SampleType,Individual_ID_Recoded) %>% summarise(TotalCPU_Seconds_Alignment = sum(TotalCPU_Seconds))
Combined_Plot_Count <- Combined_Cross %>% group_by(Alignment,SeedLength,Sample_ID,Individual_ID,Group,SampleType,Individual_ID_Recoded) %>% summarize(row_count = n(), .groups = "drop")

Combined_Plot <- full_join(Combined_Plot_Time,Combined_Plot_Count)

Combined_Plot <- Combined_Plot %>% mutate(TotalCPU_Hours_Alignment = TotalCPU_Seconds_Alignment / 3600)
model <- lm(TotalCPU_Hours_Alignment ~ Sample_ID - 1 , data = Combined_Plot)
intercept <- model$coefficients[1]
Combined_Plot$Normalized_TotalCPU_Hours_Alignment <- residuals(model) + intercept

write.table(Combined_Plot, file = "data_clean/Alignment_Speed_Cross_Summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

summary(lm(TotalCPU_Hours_Alignment ~ Individual_ID -1 + SeedLength + Group + SampleType,data = Combined_Plot))$coefficients %>% as.data.frame() %>% tibble::rownames_to_column("Term") %>%
	write_tsv(., "data_clean/Alignment_Speed_Cross_LinearModel.tsv")
summary(lm(TotalCPU_Hours_Alignment ~ Individual_ID -1 + SeedLength + Group + SampleType,data = Combined_Plot %>% filter(Individual_ID_Recoded != "Sample_3")))$coefficients %>% as.data.frame() %>% tibble::rownames_to_column("Term") %>%
	write_tsv(., "data_clean/Alignment_Speed_Cross_LinearModel_HeavyContamDrop.tsv")
summary(lm(TotalCPU_Hours_Alignment ~ Individual_ID -1 + SeedLength + Group ,data = Combined_Plot %>% filter(SampleType == "Saliva")))$coefficients %>% as.data.frame() %>% tibble::rownames_to_column("Term") %>%
	write_tsv(., "data_clean/Alignment_Speed_Cross_LinearModel_SalivaOnly.tsv")

plot <-  ggplot(Combined_Plot, aes(x = SeedLength, y = TotalCPU_Hours_Alignment, group = interaction(Sample_ID, Group), colour=Group)) +
	geom_line() +
	labs(x = "BWA Alignment Seed Length", y = "Total CPU Time for Alignment (Hours)") +
	scale_y_continuous() +
	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Cross_Compute_Split.pdf"),device="pdf", width = 60, height = 40, units = "cm",dpi = 640)

plot <-  ggplot(Combined_Plot, aes(x = SeedLength, y = Normalized_TotalCPU_Hours_Alignment, group = interaction(Sample_ID, Group), colour=Group)) +
	geom_line() +
	labs(x = "BWA Alignment Seed Length", y = "Normalized Total CPU Time for Alignment (Hours)") +
	scale_y_continuous() +
	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Cross_Compute_Split_Normalised.pdf"),device="pdf", width = 60, height = 40, units = "cm",dpi = 640)

plot <-  ggplot(Combined_Plot, aes(x = SeedLength, y = Normalized_TotalCPU_Hours_Alignment)) +
	geom_boxplot() +
	facet_wrap(~ Group) +
	labs(x = "BWA Alignment Seed Length", y = "Normalized Total CPU Time for Alignment (Hours)") +
	scale_y_continuous() +
	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Cross_Compute_Split_Normalised_Boxplot_SeedLength.pdf"),device="pdf", width = 60, height = 40, units = "cm",dpi = 640)

plot <-  ggplot(Combined_Plot, aes(x = Group , y = Normalized_TotalCPU_Hours_Alignment)) +
	geom_boxplot() +
	facet_wrap(~ SeedLength) +
	labs(x = "Cleanup Method", y = "Normalized Total CPU Time for Alignment (Hours)") +
	scale_y_continuous() +
	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Cross_Compute_Split_Normalised_Boxplot_Group.pdf"),device="pdf", width = 60, height = 40, units = "cm",dpi = 640)

####################################################################################################################
###                                           		    Kraken Burden 			                                ####
####################################################################################################################

for( stat in negative_intersect %>% keep(~ str_detect(., "Krakentaxid")) ) {
	file <- paste0(dir,"/data_raw/",stat)
	name <- gsub("_sacct_output.tsv","",stat)
	if(file.exists(file)) {
		sacct <- fread(file,fill=TRUE) 
	
		sample_id <- name %>% gsub("_.*","",.)
		tag <- name %>% gsub(paste0(sample_id,"_NoFastP_"),"",.) %>% gsub("_k.*","",.)
		seedlength <- name %>% gsub(paste0(sample_id,"_"),"",.) %>% gsub(".*_k","",.)

		new_row <- sacct %>% filter(str_detect(JobName,"Kraken2")) %>% 
			mutate(Alignment=name) %>%
			mutate(Fastq = gsub("_Kraken2_.*","",JobName)) %>%
			mutate(Chunk = gsub(".*-","",JobName)) %>%
			mutate(Sample_ID=sample_id) %>%
			mutate(Tag=tag) %>%
			mutate(SeedLength=seedlength)
			
  		assign(paste0("kraken_",name),new_row)
	} else {
		print(paste0(file," file doesn't exist"))
	}
}
listOfDataFrames <- ls(pattern=paste0("^kraken_*"))
Combined_Kraken <- bind_rows(mget(listOfDataFrames)) %>% mutate(Group=ifelse(Tag == "NoKraken","SeedLength",ifelse(Tag == "Krakentaxid","Kraken",ifelse(Tag == "Kraken_trained_taxid","KrakenTrained",Tag)))) %>%
	mutate(SampleType=ifelse(substr(Sample_ID, nchar(Sample_ID), nchar(Sample_ID)) == "S","Saliva",ifelse(substr(Sample_ID, nchar(Sample_ID), nchar(Sample_ID)) == "B","Blood","ERROR"))) %>%
	mutate(Individual_ID = str_sub(Sample_ID, 1, str_length(Sample_ID) - 1))

SampleMap <- fread(file = samplemap, sep = '\t', colClasses = "character") 
Combined_Kraken <- Combined_Kraken %>% mutate(Individual_ID_Recoded = recode(Individual_ID, !!!setNames(SampleMap$Paper, SampleMap$Original)))

write.table(Combined_Kraken, file = "data_clean/Alignment_Speed_Kraken.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

Combined_Kraken <- Combined_Kraken %>% filter(MaxRSS == "COMPLETED") %>%
  mutate(
    TotalCPU_Seconds = sapply(TotalCPU, function(x) {
      parts <- as.numeric(unlist(strsplit(x, ":")))
      if (length(parts) == 3) {
        # H:M:S
        parts[1] * 3600 + parts[2] * 60 + parts[3]
      } else if (length(parts) == 2) {
        # M:S
        parts[1] * 60 + parts[2]
      } else {
        parts[1]
      }
    })
  )

Combined_Kraken_Plot_Time <- Combined_Kraken %>% group_by(Alignment,SeedLength,Sample_ID,Individual_ID,SampleType,Individual_ID_Recoded) %>% summarise(TotalCPU_Seconds_Alignment = sum(TotalCPU_Seconds))
Combined_Kraken_Plot_Count <- Combined_Kraken %>% group_by(Alignment,SeedLength,Sample_ID,Individual_ID,SampleType,Individual_ID_Recoded) %>% summarize(row_count = n(), .groups = "drop")

Combined_Kraken_Plot <- full_join(Combined_Kraken_Plot_Time,Combined_Kraken_Plot_Count)

Combined_Kraken_Plot <- Combined_Kraken_Plot %>% mutate(TotalCPU_Hours_Alignment = TotalCPU_Seconds_Alignment / 3600) %>%
	mutate(SampleType=ifelse(substr(Sample_ID, nchar(Sample_ID), nchar(Sample_ID)) == "S","Saliva",ifelse(substr(Sample_ID, nchar(Sample_ID), nchar(Sample_ID)) == "B","Blood","ERROR"))) %>%
	mutate(Individual_ID = str_sub(Sample_ID, 1, str_length(Sample_ID) - 1))

write.table(Combined_Kraken_Plot, file = "data_clean/Alignment_Speed_Kraken_Summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

summary_all <- Combined_Kraken_Plot %>%
  ungroup() %>%
  summarise(
    Group = "All",
    mean = mean(TotalCPU_Hours_Alignment, na.rm = TRUE),
    sd = sd(TotalCPU_Hours_Alignment, na.rm = TRUE),
    median = median(TotalCPU_Hours_Alignment, na.rm = TRUE),
    n = n()
  )
summary_by_subgroups <- Combined_Kraken_Plot %>%
  group_by(SampleType) %>%
  summarise(
    mean = mean(TotalCPU_Hours_Alignment, na.rm = TRUE),
    sd = sd(TotalCPU_Hours_Alignment, na.rm = TRUE),
    median = median(TotalCPU_Hours_Alignment, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

bind_rows(summary_all, summary_by_subgroups) %>%
  write.table("data_clean/Alignment_Speed_Kraken_Summary.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)

####################################################################################################################
###                                              	 fin			                                            ####
####################################################################################################################

quit()
