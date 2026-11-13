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
topdir <- PROJECT_DIR
dir <- paste0(topdir,"Analyses/Kmer_Comparison_Analysis")
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
if(!require(patchwork)){
    install.packages("patchwork")
}
library(patchwork)

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
###                                              Meta Read In                                                   ####
####################################################################################################################

	Meta <- paste0(PROJECT_DIR,"/MetaData/Comparison_Sequenced_Samples_Cleanup.tsv")

### Saliva kmer Graph ###

	dataframe_Normalised_Percentage <- read.table( file = "data_clean/SeedLength_Count_kmer_Uniq_Saliva_Boxplot_Normalised_Percentage_PlotData.tsv", sep = "\t", header = TRUE, row.names = NULL, quote = "" )  
    
    plot_A <- ggplot(dataframe_Normalised_Percentage, aes(x = SeedLength_Num, y = Uniq_Saliva_Count_normalized_percentage, group = SeedLength_Num)) +
    	geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#1B9E77", fill = "#1B9E77") +  # hide outliers since we show all points
    	geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#1B9E77", fill = "#1B9E77") +            # show all points
        labs(x = "BWA Alignment Seed Length", y = paste0("Normalised Unique Saliva k-mers (%)")) +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
 		scale_y_continuous() +
        ggplot_theme_basic

### Common kmer Graph ###

 	Analyses_SeedLength_Adjusted <- read.table( file = "data_clean/SeedLength_Count_kmer_Common_Boxplot_Normalised_PlotData.tsv", sep = "\t", header = TRUE, row.names = NULL, quote = "" )  
   
    plot_B <- ggplot(Analyses_SeedLength_Adjusted, aes(x = SeedLength_Num, y = Common_Count_adjusted, group = SeedLength_Num)) +
    	geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#D95F02", fill = "#D95F02") +  # hide outliers since we show all points
    	geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#D95F02", fill = "#D95F02") +            # show all points
        labs(x = "BWA Alignment Seed Length", y = paste0("Adjusted Common k-mers (Millions)")) +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		scale_y_continuous(labels = function(x) formatC(x/1e6, format = "f", digits = 1), expand = expansion(mult = c(0.1, 0.07))) +
        ggplot_theme_basic
  
### Saliva alignment Graph ###

	Seed_Length_Comp_Count_Percentage <- read.table( file = paste0(topdir, "/Analyses/SoftClip_Comparison_Analysis/data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Masked_Saliva_SeedLength_Percentage_PlotData.tsv"), sep = "\t", header = TRUE, row.names = NULL, quote = "" )  

	plot_C <- ggplot(Seed_Length_Comp_Count_Percentage, aes(x = SeedLength, y = n_percentage, group=SeedLength)) +
    	geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#1B9E77", fill = "#1B9E77") +  # hide outliers since we show all points
    	geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#1B9E77", fill = "#1B9E77") +            # show all points
		labs(x = "BWA Alignmnent Seed Length", y = "Genomic Regions with Misalignmnet (%)") +
		scale_y_continuous() +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		ggplot_theme_basic

### Saliva variant Graph ###


	Seed_Length_Comp_Percentage <- read.table( file = paste0(topdir, "/Analyses/Variant_Comparison_Analysis/data_clean/Combined_DoubleSoftClippedRegions_Variants_HC_SeedLength_Comparison_Percentage_PlotData.tsv"), sep = "\t", header = TRUE, row.names = NULL, quote = "" )  

	plot_D <- ggplot(Seed_Length_Comp_Percentage, aes(x = SeedLength, y = n_percentage, group=SeedLength)) +
    	geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#1B9E77", fill = "#1B9E77") +  # hide outliers since we show all points
    	geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#1B9E77", fill = "#1B9E77") +            # show all points
		labs(x = "BWA Alignmnent Seed Length", y = "Unique Variants in Misaligned Regions (%)") +
		scale_y_continuous() +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		ggplot_theme_basic
  
### MultiPlot ###
 
plot_A <- plot_A + labs(tag = "A)") + theme(plot.tag = element_text(size = 28, face = "bold"))
plot_B <- plot_B + labs(tag = "B)") + theme(plot.tag = element_text(size = 28, face = "bold"))
plot_C <- plot_C + labs(tag = "C)") + theme(plot.tag = element_text(size = 28, face = "bold"))
plot_D <- plot_D + labs(tag = "D)") + theme(plot.tag = element_text(size = 28, face = "bold"))
plot <- plot_A + plot_B + plot_C + plot_D + plot_layout(ncol = 2, widths = c(0.9, 1))
ggsave(plot=plot, filename="graphs/Joint_MultiPlot.pdf", device="pdf", width = 56, height = 56, units = "cm", dpi = 640)

####################################################################################################################
###                                              	 fin			                                            ####
####################################################################################################################

quit()
