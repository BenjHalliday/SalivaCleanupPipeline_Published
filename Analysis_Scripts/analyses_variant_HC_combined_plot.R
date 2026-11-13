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
dir <- paste0(PROJECT_DIR,"/Analyses/Variant_Comparison_Analysis/")
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

if (!require(ggpubr)) install.packages("ggpubr")
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
	legend.title = element_text(color = "black", size = 26, margin = margin(t = 0, r = 0, b = 20, l = 0)),
	legend.text = element_text(color = "black", size = 24),
	#legend.position = "none" ,
	legend.key = element_rect(fill = "white",colour="white"),
	legend.margin = margin(t = 0, r = 20, b = 0, l = 20),
	axis.line = element_blank()
)

####################################################################################################################
###                                               Loadin Data	                                                ####
####################################################################################################################

files <- list.files("data_raw/", pattern = "_out\\.tsv$")
file_tags <- files %>% gsub("_SoftClip.*","",.)

for(file in files) {
	file_tag <- file_tags[which(files == file)]
	if(file.exists(paste0("data_raw/",file))) {
	    if(file.info(paste0("data_raw/",file))$size > 0){
	        data <- fread(file = paste0("data_raw/",file), sep = '\t', colClasses = "character") 

	        id <- file_tag %>% gsub("_.*","",.)
			group <- file_tag %>% gsub(paste0(id,"_NoFastP_"),"",.) %>% gsub("_k.*","",.)
			seedlength <- file_tag %>% gsub(paste0(id,"_"),"",.) %>% gsub(".*_k","",.)
			
			colnames(data) <- colnames(data) %>% gsub(".*]","",.) %>% gsub(id,"Individual",.)
			
			data$Sample_ID <- id
			data$Tag <- group
			data$SeedLength <- seedlength
					
	        assign(paste0("HC_",file_tag),data)
	    } else {
	        print(paste0(file," is empty!!!"))
	    } 
	} else {
	    print(paste0(file," is missing!!!"))
	}
}
ls(pattern=paste0("^HC_*")) -> listOfDataFrames
bind_rows(mget(listOfDataFrames)) -> Combined

Combined <- Combined %>% mutate(Group=ifelse(Tag == "NoKraken","SeedLength",ifelse(Tag == "Krakentaxid","Kraken",ifelse(Tag == "Kraken_trained_taxid","KrakenTrained",Tag)))) %>% 
	mutate(SampleType=ifelse(substr(Sample_ID, nchar(Sample_ID), nchar(Sample_ID)) == "S","Saliva",ifelse(substr(Sample_ID, nchar(Sample_ID), nchar(Sample_ID)) == "B","Blood","ERROR"))) %>%
	mutate(SeedLength=as.numeric(SeedLength)) %>% filter(CHROM %in% paste0("chr",seq(1,22)))

write.table(Combined, file=paste0("data_clean/Combined_DoubleSoftClippedRegions_Variants_HC.tsv"), quote=FALSE, sep='\t', row.names=F)

####################################################################################################################
###                                              Anti_Join Saliva-Blood                                         ####
####################################################################################################################

Combined %>% filter(SampleType == "Saliva") %>% mutate( Individual_ID = substr(Sample_ID, 1, nchar(Sample_ID) - 1) ) -> Combined_Saliva
Combined %>% filter(SampleType == "Blood") %>% mutate( Individual_ID = substr(Sample_ID, 1, nchar(Sample_ID) - 1) ) -> Combined_Blood

anti_join(Combined_Saliva,Combined_Blood,by=c("Individual_ID","Tag","SeedLength","Group","CHROM","POS","REF","ALT")) -> Combined_Saliva_Unique
anti_join(Combined_Blood,Combined_Saliva,by=c("Individual_ID","Tag","SeedLength","Group","CHROM","POS","REF","ALT")) -> Combined_Blood_Unique
semi_join(Combined_Blood,Combined_Saliva,by=c("Individual_ID","Tag","SeedLength","Group","CHROM","POS","REF","ALT")) -> Combined_Union

write.table(Combined_Blood_Unique, file=paste0("data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_BloodUnique.tsv"), quote=FALSE, sep='\t', row.names=F)
write.table(Combined_Union, file=paste0("data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_Union.tsv"), quote=FALSE, sep='\t', row.names=F)

####################################################################################################################
###                                                 	Seed Length                                             ####
####################################################################################################################

Seed_Length_Comp <- Combined_Saliva_Unique %>% filter(Group == "SeedLength") %>% group_by(Sample_ID,SeedLength) %>% count()

write.table(Seed_Length_Comp, file=paste0("data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison.tsv"), quote=FALSE, sep='\t', row.names=F)

plot <- ggplot(Seed_Length_Comp, aes(x = SeedLength, y = n, colour = Sample_ID)) +
		  geom_line() +
		  labs(x = "BWA Alignmnent Seed Length", y = "Number of Unique Saliva Variants in Double Soft Clipped Regions") +
		  scale_y_continuous() +
		  scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

Seed_Length_Comp_Percentage <- Seed_Length_Comp %>%
  group_by(Sample_ID) %>%
  mutate(n_percentage = n / first(n) * 100) %>%
  ungroup() 
    
plot <- ggplot(Seed_Length_Comp_Percentage, aes(x = SeedLength, y = n_percentage, colour = Sample_ID)) +
		  geom_line() +
		  labs(x = "BWA Alignmnent Seed Length", y = "Unique Saliva Variants in Double Soft Clipped Regions (%)") +
		  scale_y_continuous() +
		  scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Percentage.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

plot <- ggplot(Seed_Length_Comp_Percentage, aes(x = SeedLength, y = n_percentage, group=SeedLength)) +
    		geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#1B9E77", fill = "#1B9E77") +  # hide outliers since we show all points
    		geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#1B9E77", fill = "#1B9E77") +            # show all points
		  labs(x = "BWA Alignmnent Seed Length", y = "Unique Saliva Variants in Misalignment Regions (%)") +
		  scale_y_continuous() +
		  scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Percentage_Boxplot.pdf"),device="pdf", width = 27.5, height = 30, units = "cm",dpi = 640)

write.table(Seed_Length_Comp_Percentage, file = "data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Percentage_PlotData.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################################################################
###                                                 	Seed Length  Blood                                      ####
####################################################################################################################

Seed_Length_Comp_Union <- Combined_Union %>% filter(Group == "SeedLength") %>% group_by(Sample_ID,SeedLength) %>% count()

write.table(Seed_Length_Comp_Union, file=paste0("data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Union.tsv"), quote=FALSE, sep='\t', row.names=F)

plot <- ggplot(Seed_Length_Comp_Union, aes(x = SeedLength, y = n, colour = Sample_ID)) +
		  geom_line() +
		  labs(x = "BWA Alignmnent Seed Length", y = "Number of Variants at Double Soft Clipped Regions") +
		  scale_y_continuous() +
		  scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Union.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

Seed_Length_Comp_Percentage_Union <- Seed_Length_Comp_Union %>%
  group_by(Sample_ID) %>%
  mutate(n_percentage = n / first(n) * 100) %>%
  ungroup() 
  
plot <- ggplot(Seed_Length_Comp_Percentage_Union, aes(x = SeedLength, y = n_percentage, colour = Sample_ID)) +
		  geom_line() +
		  labs(x = "BWA Alignmnent Seed Length", y = "Number of Variants at Double Soft Clipped Regions") +
		  scale_y_continuous() +
		  scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Percentage_Union.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

plot <- ggplot(Seed_Length_Comp_Percentage_Union, aes(x = SeedLength, y = n_percentage, group=SeedLength)) +
		  geom_boxplot() +
		  labs(x = "BWA Alignmnent Seed Length", y = "Number of Variants at Double Soft Clipped Regions") +
		  scale_y_continuous() +
		  scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Percentage_Boxplot_Union.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

####################################################################################################################
###                                                 	Seed Length  Union                                      ####
####################################################################################################################

Seed_Length_Comp_Blood <- Combined_Blood_Unique %>% filter(Group == "SeedLength") %>% group_by(Sample_ID,SeedLength) %>% count()

write.table(Seed_Length_Comp_Blood, file=paste0("data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Blood.tsv"), quote=FALSE, sep='\t', row.names=F)

plot <- ggplot(Seed_Length_Comp_Blood, aes(x = SeedLength, y = n, colour = Sample_ID)) +
		  geom_line() +
		  labs(x = "BWA Alignmnent Seed Length", y = "Number of Variants at Double Soft Clipped Regions") +
		  scale_y_continuous() +
		  scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Blood.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

Seed_Length_Comp_Percentage_Blood <- Seed_Length_Comp_Blood %>%
  group_by(Sample_ID) %>%
  mutate(n_percentage = n / first(n) * 100) %>%
  ungroup() 
  
plot <- ggplot(Seed_Length_Comp_Percentage_Blood, aes(x = SeedLength, y = n_percentage, colour = Sample_ID)) +
		  geom_line() +
		  labs(x = "BWA Alignmnent Seed Length", y = "Number of Variants at Double Soft Clipped Regions") +
		  scale_y_continuous() +
		  scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Percentage_Blood.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

plot <- ggplot(Seed_Length_Comp_Percentage_Blood, aes(x = SeedLength, y = n_percentage, group=SeedLength)) +
		  geom_boxplot() +
		  labs(x = "BWA Alignmnent Seed Length", y = "Number of Variants at Double Soft Clipped Regions") +
		  scale_y_continuous() +
		  scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Percentage_Boxplot_Blood.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

####################################################################################################################
###                                                 	Complete Comp                                           ####
####################################################################################################################

groups <- c("Decoy", "Kraken", "KrakenTrained","SeedLength")
seedlengths <- c("19", "23", "24")

Complete_Comp <- Combined_Saliva_Unique %>% filter(Group %in% groups) %>% filter(SeedLength %in% seedlengths)
Complete_Comp_Count <- Complete_Comp %>% group_by(Sample_ID,SeedLength,Group) %>% count()

Complete_Comp_Count_Normalized <- Complete_Comp_Count %>%
  group_by(Sample_ID, SeedLength) %>%
  mutate(n_normalized = n / n[Group == "SeedLength"]) %>% filter(Group != "SeedLength")

write.table(Complete_Comp_Count_Normalized, file=paste0("data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_Complete_Comparison.tsv"), quote=FALSE, sep='\t', row.names=F)

plot <- ggplot(Complete_Comp_Count_Normalized, aes(x = Group, y = n_normalized, group=Group)) +
		  geom_boxplot() +
		  labs(x = "Targetted Cleanup Method", y = "Percentage of Variants at Double Soft Clipped\nRegions Remaining After Cleanup Method") +
		  scale_y_continuous(limits=c(0.1,1.1)) +
		  stat_compare_means( method = "t.test", comparisons = list(c("Decoy", "Kraken"), c("Decoy", "KrakenTrained")), label = "p.signif") +
		  facet_wrap(~ SeedLength, nrow = 1, scales = "free_x", labeller = labeller(SeedLength = function(x) paste("Alignment Seed\nLength:", x))) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_Complete_Comparison.pdf"),device="pdf", width = 60, height = 45, units = "cm",dpi = 640)

plot <- ggplot((Complete_Comp_Count_Normalized), aes(x = Group, y = n_normalized, group=Group, colour=Group,fill=Group,alpha=0.3)) +
		  geom_boxplot(outlier.size = 5) +
		  labs(x = "Targeted Cleanup Method", y = "Misalignment Remaining after Targeted Cleanup Method\nRelative to No Cleanup (%)") +
		  scale_y_continuous(limits=c(0.4,1.1)) +
  		  scale_color_brewer(palette = "Dark2") +
  		  scale_fill_brewer(palette = "Dark2") +
		  stat_compare_means( method = "t.test", comparisons = list(c("Decoy", "Kraken"), c("Decoy", "KrakenTrained")),size = 12, label = "p.signif" , paired = TRUE) +
		  facet_wrap(~ SeedLength, nrow = 1, scales = "free_x", labeller = labeller(SeedLength = function(x) paste("Alignment Seed\nLength:", x))) +
		  guides(alpha = "none") +
		  ggplot_theme_basic +
		  theme(
		  	axis.text.x = element_text(size=24, margin = margin(t = 10, r = 0, b = 0, l = 0)),
			axis.text.y  = element_text(size=30, margin = margin(t = 10, r = 0, b = 0, l = 0),hjust = 0),
			axis.title.y = element_text(size=36, margin = margin(t = 0, r = 20, b = 0, l = 0)),
			axis.title.x = element_text(size=36, margin = margin(t = 20, r = 0, b = 0, l = 0)),
			strip.text = element_text(angle=0, size=38),
			legend.title = element_text(color = "black", size = 36, margin = margin(t = 0, r = 0, b = 0, l = 0)),
			legend.text = element_text(color = "black", size = 32)
		  )
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_Complete_Comparison_Normalised_Boxplot.pdf"),device="pdf", width = 60, height = 45, units = "cm",dpi = 640)


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
write.table(t_test_trained, file = "data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_Complete_Comparison_Normalised_t-test-Trained.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

t_test_trained <- paired_t_by_seed(Complete_Comp_Count_Normalized, "Kraken", "Decoy")
print(t_test_trained)
write.table(t_test_trained, file = "data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_Complete_Comparison_Normalised_t-test-Untrained.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

####################################################################################################################
###                                                 	Variability	                                            ####
####################################################################################################################

Combined_Saliva_Unique_VariantConsistency <- inner_join(rbind(Combined_Saliva_Unique,Combined_Blood_Unique),rbind(Combined_Saliva_Unique,Combined_Blood_Unique) %>% filter(Group == "SeedLength") %>% group_by(Sample_ID,CHROM,POS,REF,ALT) %>% count(), by=c("Sample_ID","CHROM","POS","REF","ALT")) %>% 
	ungroup() %>% distinct(Sample_ID,CHROM,POS,REF,ALT,SeedLength,.keep_all = TRUE) 
write.table(Combined_Saliva_Unique_VariantConsistency, file=paste0("data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_VariantConsistency.tsv"), quote=FALSE, sep='\t', row.names=F)

Combined_Saliva_Unique_VariantConsistency_Count <- Combined_Saliva_Unique_VariantConsistency %>% group_by(Sample_ID,SeedLength,n) %>% count(name = "n_count")

dark2_9 <- c(RColorBrewer::brewer.pal(8, "Dark2"), "#4D4D4D")

for( sample in Combined_Saliva_Unique_VariantConsistency_Count %>% ungroup() %>% pull(Sample_ID) %>% unique()){
	plot <- ggplot(filter(Combined_Saliva_Unique_VariantConsistency_Count, Sample_ID == sample), aes(x = factor(SeedLength), y = n_count, fill = factor(n), colour = factor(n))) +
		geom_bar(width = 0.6, alpha = 0.3, stat = "identity") +
  		labs(x = "BWA Alignmnent Seed Length", y = "Number of Unique Saliva Variants\nin Misalignment Regions", fill = "Variant Occurence\nAcross Conditions", colour = "Variant Occurence\nAcross Conditions") +
		scale_colour_manual(values = dark2_9) +
		scale_fill_manual(values = dark2_9) +
		scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
		guides(
    		fill = guide_legend(ncol = 3, byrow = TRUE, keywidth = unit(1.5, "cm"), keyheight = unit(1.5, "cm")),
    		colour = guide_legend(ncol = 3, byrow = TRUE, keywidth = unit(1.5, "cm"), keyheight = unit(1.5, "cm"))
		) +
		theme(legend.spacing.y = unit(0.5, "cm")) +
		ggplot_theme_basic 
	ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_BarPlot_",sample,".pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)
}

plot <- ggplot(Seed_Length_Comp_Percentage, aes(x = SeedLength, y = n_percentage, group=SeedLength)) +
    		geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#1B9E77", fill = "#1B9E77") +  # hide outliers since we show all points
    		geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#1B9E77", fill = "#1B9E77") +            # show all points
		  labs(x = "BWA Alignmnent Seed Length", y = "Unique Saliva Variants in Double Soft Clipped Regions (%)") +
			scale_y_continuous() +
		  scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		  ggplot_theme_basic
ggsave(plot=plot,filename=paste0("graphs/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Percentage_Boxplot.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)

####################################################################################################################
###                                                 	Fin			                                            ####
####################################################################################################################

quit()
