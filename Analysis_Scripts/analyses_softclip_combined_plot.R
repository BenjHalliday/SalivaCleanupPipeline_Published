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
dir <- paste0(PROJECT_DIR,"/Analyses/SoftClip_Comparison_Analysis")
Dirs <- c("graphs","data_clean","data_raw")
ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir)),TRUE)
setwd(dir)
for(subDir in Dirs) {
    ifelse(!dir.exists(file.path(dir, subDir)), dir.create(file.path(dir, subDir)), TRUE)
}

print(paste0("Directory Set to ", dir))

###      Packages         ### 

if(!require(tidyverse)){
    install.packages("tidyverse")
}
library(tidyverse)
if(!require(data.table)){
    install.packages("data.table")
}
library(data.table)
if(!require(ggraph)){
    install.packages("ggraph")
}
library(ggraph)
if(!require(ggpubr)){
    install.packages("ggpubr")
}
library(ggpubr)

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
###                                               Loadin Data	                                                ####
####################################################################################################################

files <- list.files(paste0(PROJECT_DIR,"/Alignments/"), pattern = "_DoubleSoftClipped_Repeat_Filtered_Mappability_FreqTable\\.tsv$", recursive = TRUE, full.names = TRUE)
file_tags <- files %>% gsub(".*/CollectCigarMetrics/","",.) %>% gsub("_DoubleSoftClipped_.*","",.)

for(file in files) {
	file_tag <- file_tags[which(files == file)]
	file.copy(from = file, to = paste0("data_raw/",file_tag,"_DoubleSoftClipped_Repeat_Filtered_Mappability_FreqTable.tsv"), overwrite=T)
	
	if(file.exists(file)){
	    if(file.info(file)$size > 0){
	        data <- fread(file = file, sep = ' ',col.names=c("Frequency","Chromosome","Start","End","Length")) %>% filter(Frequency >= 2)

	        id <- file_tag %>% gsub("_.*","",.)
			group <- file_tag %>% gsub(paste0(id,"_NoFastP_"),"",.) %>% gsub("_k.*","",.)
			seedlength <- file_tag %>% gsub(paste0(id,"_"),"",.) %>% gsub(".*_k","",.)
			
			data$Sample_ID <- id
			data$Tag <- group
			data$SeedLength <- seedlength
						
	        assign(paste0("ft_",file_tag),data)
	    } else {
	        print(paste0(file," is empty!!!"))
	    } 
	} else {
	    print(paste0(file," is missing!!!"))
	}
}
ls(pattern=paste0("^ft_*")) -> listOfDataFrames
bind_rows(mget(listOfDataFrames)) -> Combined

Combined <- Combined %>% mutate(Group=ifelse(Tag == "NoKraken","SeedLength",ifelse(Tag == "Krakentaxid","Kraken",ifelse(Tag == "Kraken_trained_taxid","KrakenTrained",Tag)))) %>% 
	mutate(SampleType=ifelse(substr(Sample_ID, nchar(Sample_ID), nchar(Sample_ID)) == "S","Saliva",ifelse(substr(Sample_ID, nchar(Sample_ID), nchar(Sample_ID)) == "B","Blood","ERROR"))) %>%
	mutate(SeedLength=as.numeric(SeedLength))

write.table(Combined, file=paste0("data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_FreqTable.tsv"), quote=FALSE, sep='\t', row.names=F)

####################################################################################################################
###                                                 Create Mask Region                                          ####
####################################################################################################################

Combined_Filtered <- Combined %>% filter(Group != "Minimap2") #%>% filter(Frequency >= 5)
Blood_Mask <- Combined_Filtered %>% filter(SampleType == "Blood") %>% select(Chromosome,Start,End,Sample_ID) %>% unique() %>% group_by(Chromosome,Start,End) %>% count()

write.table(Blood_Mask, file=paste0("data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_BloodMask.tsv"), quote=FALSE, sep='\t', row.names=F)

Saliva_Masked <- anti_join(Combined_Filtered,Blood_Mask) %>% filter(SampleType == "Saliva") 

write.table(Saliva_Masked, file=paste0("data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Masked_Saliva.tsv"), quote=FALSE, sep='\t', row.names=F)

Saliva_Masked_24 <- Saliva_Masked %>% filter(SeedLength == 24)

write.table(Saliva_Masked_24, file=paste0("data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Masked_Saliva_k24.tsv"), quote=FALSE, sep='\t', row.names=F)

Saliva_Masked_Summary <- Saliva_Masked %>% filter(Frequency >= 4 & Group == "SeedLength" & SeedLength == "19" & Length >= 28) %>%
  count(Sample_ID) %>%
  summarise(
    median = median(n, na.rm = TRUE),
    Q1     = quantile(n, 0.25, na.rm = TRUE),
    Q3     = quantile(n, 0.75, na.rm = TRUE),
    IQR    = IQR(n, na.rm = TRUE),
    min    = min(n, na.rm = TRUE),
    max    = max(n, na.rm = TRUE)
  )
write.table(Saliva_Masked_Summary, file = "data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Masked_Saliva_Large_Summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################################################################
###                                                 	Seed Length                                             ####
####################################################################################################################

Seed_Length_Comp <- Saliva_Masked %>% filter(Group == "SeedLength") %>% filter(Frequency >= 4)
Seed_Length_Comp_Count <- Seed_Length_Comp %>% select(Sample_ID,SeedLength,Chromosome,Start,End) %>% group_by(Sample_ID,SeedLength) %>% count()

write.table(Seed_Length_Comp_Count, file=paste0("data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_SeedLengthComparison.tsv"), quote=FALSE, sep='\t', row.names=F)

Seed_Length_Comp_Count <- Seed_Length_Comp_Count %>% mutate(n_k = n / 1000)

plot <- ggplot(Seed_Length_Comp_Count, aes(x = SeedLength, y = n_k, colour = Sample_ID)) +
		  	geom_line() +
		  	labs(x = "BWA Alignmnent Seed Length", y = "Genomic Regions with Misalignmnet (1000's)") +
		  	scale_y_continuous() +
		  	scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_SeedLength_Comparison.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

Seed_Length_Comp_Count_Percentage <- Seed_Length_Comp_Count %>%
  group_by(Sample_ID) %>%
  mutate(n_percentage = n / first(n) * 100) %>%
  ungroup() 
  
write.table(Seed_Length_Comp_Count_Percentage, file=paste0("data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Masked_Saliva_SeedLength_Percentage.tsv"), quote=FALSE, sep='\t', row.names=F)
  
plot <- ggplot(Seed_Length_Comp_Count_Percentage, aes(x = SeedLength, y = n_percentage, colour = Sample_ID)) +
		  	geom_line() +
		  	labs(x = "BWA Alignmnent Seed Length", y = "Genomic Regions with Misalignmnet (%)") +
		  	scale_y_continuous() +
		  	scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_SeedLength_Comparison_Percentage.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

plot <- ggplot(Seed_Length_Comp_Count_Percentage, aes(x = SeedLength, y = n_percentage, group=SeedLength)) +
    		geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#1B9E77", fill = "#1B9E77") +  # hide outliers since we show all points
    		geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#1B9E77", fill = "#1B9E77") +            # show all points
		  	labs(x = "BWA Alignmnent Seed Length", y = "Genomic Regions with Misalignmnet (%)") +
		  	scale_y_continuous() +
		  	scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  	ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_SeedLength_Comparison_Percentage_Boxplot.pdf"),device="pdf", width = 27.5, height = 30, units = "cm",dpi = 640)

write.table(Seed_Length_Comp_Count_Percentage, file = "data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Masked_Saliva_SeedLength_Percentage_PlotData.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
Seed_Length_Comp_Count_Percentage_Summary <- Seed_Length_Comp_Count_Percentage %>% group_by(SeedLength) %>% summarise(n_percentage_median=median(n_percentage),n_percentage_sd = sd(n_percentage))
write.table(Seed_Length_Comp_Count_Percentage_Summary, file = "data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Masked_Saliva_SeedLength_Percentage_Summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################################################################
###                                                 	Complete Comp                                           ####
####################################################################################################################

groups <- c("Decoy", "Kraken", "KrakenTrained","SeedLength")
seedlengths <- c(19,23,24)
seedlengths_num <- as.integer(seedlengths)
plus1_map <- tibble(
    SeedLength = as.character(seedlengths_num),
    SeedLengthPlus1 = as.character(seedlengths_num + 1)
)

Combined_Filtered_Plus1 <- anti_join(Combined_Filtered, Blood_Mask) %>%
	filter(SampleType == "Saliva") %>% 
	filter(Group == "SeedLength") %>%
	filter(SeedLength %in% plus1_map$SeedLengthPlus1) %>% 
    filter(Frequency >= 4) %>%
    mutate( Group = "SeedLength + 1",SeedLength = as.integer(SeedLength) - 1 ) %>% filter(SeedLength %in% seedlengths)
     
Complete_Comp <- anti_join(Combined_Filtered, Blood_Mask) %>%
    filter(SampleType == "Saliva") %>%
    filter(Group %in% c("Decoy", "Kraken", "KrakenTrained", "SeedLength")) %>%
    filter(SeedLength %in% seedlengths) %>%
    filter(Frequency >= 4) %>%
    bind_rows(Combined_Filtered_Plus1 %>% filter(SeedLength %in% seedlengths))

Complete_Comp_Count <- Complete_Comp %>%
    select(Sample_ID, Group, SeedLength, Chromosome, Start, End) %>%
    group_by(Sample_ID, SeedLength, Group) %>%
    count()

Complete_Comp_Count_Normalized <- Complete_Comp_Count %>%
    group_by(Sample_ID, SeedLength) %>%
    mutate(n_normalized = (n / n[Group == "SeedLength"]) * 100) %>%
    filter(Group != "SeedLength") %>% filter(Group != "SeedLength + 1") %>%
    mutate(Group_Plot = recode(Group,
        "Decoy"          = "Decoy\nAlignment",
        "Kraken"         = "Kraken2\nStandard",
        "KrakenTrained"  = "Kraken2\nTrained"
    ))
write.table(Complete_Comp_Count_Normalized, file=paste0("data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_CompleteComparison.tsv"), quote=FALSE, sep='\t', row.names=F)

plot <- ggplot((Complete_Comp_Count_Normalized), aes(x = Group_Plot, y = n_normalized, group=Group_Plot, colour=Group_Plot, fill=Group_Plot, alpha=0.3)) +
		  geom_boxplot(outlier.size = 5) +
		  labs(x = "Targeted Cleanup Method", y = "Misalignment Regions Remaining after Targeted Cleanup Method\nCompared with Matched Seed Length Alone (%)") +
		  scale_y_continuous(limits=c(30,70)) +
  		  scale_color_brewer(palette = "Dark2") +
  		  scale_fill_brewer(palette = "Dark2") +
		  stat_compare_means( method = "t.test", comparisons = list(c("Decoy\nAlignment", "Kraken2\nStandard"), c("Decoy\nAlignment", "Kraken2\nTrained")),size = 12, label = "p.signif" , paired = TRUE) +
		  facet_wrap(~ SeedLength, nrow = 1, scales = "free_x", labeller = labeller(SeedLength = function(x) paste("Alignment Seed\nLength:", x))) +
		  guides(alpha = "none") +
		  ggplot_theme_basic +
		  theme(
			legend.position = "none",
			strip.text = element_text(margin = margin(t = 5, b = 10))
		  )
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Total_Comparison_Normalised_Boxplot.pdf"),device="pdf", width = 45, height = 40, units = "cm",dpi = 640)

paired_t_by_seed <- function(df, group_a, group_b) {
  df %>%
    filter(Group %in% c(group_a, group_b)) %>%
    select(SeedLength, Sample_ID, Group, n_normalized) %>%
    pivot_wider(names_from = Group, values_from = n_normalized) %>%
    group_by(SeedLength) %>%
    summarise(
      test     = list(t.test(.data[[group_a]], .data[[group_b]], paired = TRUE)),
      p_value  = test[[1]]$p.value,
      t_stat   = unname(test[[1]]$statistic),
      df       = unname(test[[1]]$parameter),
      n_pairs  = df + 1,
      est_diff = unname(test[[1]]$estimate),     # mean of (group_a - group_b)
      ci_low   = test[[1]]$conf.int[1],
      ci_high  = test[[1]]$conf.int[2],
      .groups = "drop"
    ) %>%
    select(-test)
}

t_test_trained <- paired_t_by_seed(Complete_Comp_Count_Normalized, "KrakenTrained", "Decoy")
print(t_test_trained)
write.table(t_test_trained, file = "data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Comparison_Normalised_t-test-Trained.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

t_test_kraken <- paired_t_by_seed(Complete_Comp_Count_Normalized, "Kraken", "Decoy")
print(t_test_kraken)
write.table(t_test_kraken, file = "data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Comparison_Normalised_t-test-Untrained.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

####################################################################################################################
###                                               Contamination		                                            ####
####################################################################################################################

for(id in unique(Complete_Comp_Count_Normalized$Sample_ID)) {
        rm(list = ls(pattern="^kraken_*"))
        path <- paste0(PROJECT_DIR,"/Alignments/",id,"_NoFastP_Kraken_trained_taxid_k19/Kraken2/")
        files <- list.files(path, pattern = "_Report\\.tsv$")
        for(file in files) {
                num <- which(file == files)
                data <- fread(file = paste0(path,file), sep = '\t', quote = "", stringsAsFactors = FALSE)
                human <- data %>% filter(V6 == "Homo sapiens") %>% pull(V3)
                total <- data %>% select(V3) %>% sum()

                assign( paste0("kraken_",num), data.frame(num,human,total))

        }

        bind_rows(mget(ls(pattern="^kraken_*"))) -> combined_kraken

        tot_human <- sum(combined_kraken$human)
        tot_tot <- sum(combined_kraken$total)

        Human_Percentage <- tot_human / tot_tot
        Sample_ID <- id

        assign( paste0("percentage_",id), data.frame(Sample_ID,Human_Percentage))
}
bind_rows(mget(ls(pattern="^percentage_*")))  %>% mutate( Microbiome_Percentage = (1 - Human_Percentage) * 100 ) -> combined_human_percentages

Complete_Comp_Count_contam <- Complete_Comp_Count %>% filter(SeedLength == 19) %>% filter(Group == "SeedLength") %>% left_join(combined_human_percentages, by = "Sample_ID")
write.table(Complete_Comp_Count_contam, file=paste0("data_clean/Sample_Contamination_Estimates.tsv"), quote=FALSE, sep='\t', row.names=F)

cor(Complete_Comp_Count_contam$n,Complete_Comp_Count_contam$Microbiome_Percentage)
test <- cor.test(Complete_Comp_Count_contam$n,Complete_Comp_Count_contam$Microbiome_Percentage)
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA,
  p_value = format(test$p.value, scientific = FALSE, digits = 6), estimate = as.numeric(test$estimate), conf_low = test$conf.int[1], conf_high = test$conf.int[2], alternative = test$alternative
)
readr::write_tsv(test_out, "data_clean/Sample_Contamination_Estimatess_ContaminationCorrelation.tsv")

####################################################################################################################
###                                                 	Fin			                                            ####
####################################################################################################################

quit()
