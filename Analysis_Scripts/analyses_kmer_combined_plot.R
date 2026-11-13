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

kmin=5 #4
kmax=50 #149

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

individuals <- fread(file = Meta, sep = '\t', header = T ) %>% select(Individual) %>% mutate(Individual = sub("[SB]$", "", Individual)) %>% distinct(.,Individual,.keep_all=T) %>% pull(Individual)
conditions <- c("Decoy", "Kraken", "KrakenTrained")
seedlengths <- c("k19", "k23", "k24")

Analyses1 <- expand.grid(Individual = individuals, Condition = conditions, SeedLength = seedlengths)
Analyses1 <- Analyses1 %>% mutate( Key=paste0("kmer_",Individual,"_NoFastP_",ifelse(Condition == "Decoy","Decoy", ifelse(Condition == "Kraken","Krakentaxid", ifelse(Condition == "KrakenTrained","Kraken_trained_taxid", "NoKraken" ) ) ) , "_" , SeedLength) )

individuals <- fread(file = Meta, sep = '\t', header = T ) %>% select(Individual) %>% mutate(Individual = sub("[SB]$", "", Individual)) %>% distinct(.,Individual,.keep_all=T) %>% pull(Individual)
conditions <- c("SeedLength")
seedlengths <- c("k19", "k20", "k21", "k22", "k23", "k24", "k25", "k26", "k27")

Analyses2 <- expand.grid(Individual = individuals, Condition = conditions, SeedLength = seedlengths)
Analyses2 <- Analyses2 %>% mutate( Key=paste0("kmer_",Individual,"_NoFastP_",ifelse(Condition == "Decoy","Decoy", ifelse(Condition == "Kraken","Krakentaxid", ifelse(Condition == "KrakenTrained","Kraken_trained_taxid", "NoKraken" ) ) ) , "_" , SeedLength) )

Analyses <- rbind(Analyses1,Analyses2)

####################################################################################################################
###                                              Contamination Estimation                                       ####
####################################################################################################################

for(id in unique(individuals)) {
	for (suffix in c("B", "S")) {
        rm(list = ls(pattern="^kraken_*"))
        sample_id <- paste0(id, suffix)
        message("Processing: ", sample_id)
        
        path <- paste0(PROJECT_DIR,"/Alignments/",sample_id,"_NoFastP_Krakentaxid_k19/Kraken2/")
        files <- list.files(path, pattern = "_Report\\.tsv$") %>% discard(~ startsWith(.x, "V"))
        for(file in files) {
                num <- which(file == files)
                data <- fread(file = paste0(path,file), sep = '\t', quote = "", stringsAsFactors = FALSE)
                human <- data %>% filter(V6 == "Homo sapiens") %>% pull(V2)
                bacteria <- data %>% filter(V6 == "Bacteria") %>% pull(V2)
                unclassified <- data %>% filter(V6 == "unclassified") %>% pull(V2)
                total <- data %>% select(V3) %>% sum()
                assign( paste0("kraken_",num), data.frame(num,human,total,bacteria,unclassified))
        }
        bind_rows(mget(ls(pattern="^kraken_*"))) -> combined_kraken

        tot_human <- sum(combined_kraken$human)
        tot_unclassified <- sum(combined_kraken$unclassified)
        tot_bacteria <- sum(combined_kraken$bacteria)
        tot_tot <- sum(combined_kraken$total)
        Human_Percentage <- (tot_human / tot_tot) * 100
        Unclassified_Percentage <- (tot_unclassified / tot_tot) * 100
        Bacteria_Percentage <- (tot_bacteria / tot_tot) * 100
        Individual <- sample_id
        Type <- ifelse(suffix == "B", "Blood", "Saliva")
        assign( paste0("percentage_",sample_id), data.frame(Individual,Type,Unclassified_Percentage,Bacteria_Percentage,Human_Percentage))
    }
}
bind_rows(mget(ls(pattern="^percentage_*"))) -> combined_human_percentages

combined_diffvcontam <- combined_human_percentages %>% mutate( NonHuman_Percentage = (100 - Human_Percentage))
readr::write_tsv(combined_diffvcontam, "data_clean/Contamination_Estimates.tsv")

combined_diffvcontam_summary <- combined_diffvcontam %>% group_by(Type) %>%
    summarise(
        n        = n(),
        median   = median(NonHuman_Percentage, na.rm = TRUE),
        q25      = quantile(NonHuman_Percentage, 0.25, na.rm = TRUE),
        q75      = quantile(NonHuman_Percentage, 0.75, na.rm = TRUE),
        iqr      = IQR(NonHuman_Percentage, na.rm = TRUE),
        min      = min(NonHuman_Percentage, na.rm = TRUE),
        max      = max(NonHuman_Percentage, na.rm = TRUE),
        .groups  = "drop"
    )
readr::write_tsv(combined_diffvcontam_summary, "data_clean/Contamination_Estimates_Summary.tsv")

####################################################################################################################
###                                              Data Read In                                            	    ####
####################################################################################################################

# Saliva Unique

for( key in Analyses$Key ) {
	condition <- Analyses %>% filter(Key==key) %>% pull(Condition) %>% unique()
	individual <- Analyses %>% filter(Key==key) %>% pull(Individual) %>% unique()
	saliva <- gsub(individual,paste0(individual,"S"),key) %>% gsub("kmer_","",.)
	blood <- gsub(individual,paste0(individual,"B"),key) %>% gsub("kmer_","",.)
	
	file_path <- paste0(dir,"/data_raw/",condition,"/",key,"/data_clean/unique_to_",saliva,"_dataset_kmc_k-mer_Counts.csv")
	if(file.exists(file_path)) {
		uniq_saliva <- fread(file_path) 
		uniq_saliva_count <- uniq_saliva %>% filter(count <= kmax & count >= kmin) %>% sum()
		Analyses$Uniq_Saliva_Count[Analyses$Key == key] <- uniq_saliva_count
	} else {
		print(paste0(key," file doesn't exist"))
	}
}

for( key in Analyses$Key ) {
	condition <- Analyses %>% filter(Key==key) %>% pull(Condition) %>% unique()
	individual <- Analyses %>% filter(Key==key) %>% pull(Individual) %>% unique()
	saliva <- gsub(individual,paste0(individual,"S"),key) %>% gsub("kmer_","",.)
	blood <- gsub(individual,paste0(individual,"B"),key) %>% gsub("kmer_","",.)

	file_path <- paste0(dir,"/data_raw/",condition,"/",key,"/data_clean/unique_to_",blood,"_dataset_kmc_k-mer_Counts.csv")
	if(file.exists(file_path)) {
		uniq_blood <- fread(file_path) 
		uniq_blood_count <- uniq_blood %>% filter(count <= kmax & count >= kmin) %>% sum()
		Analyses$Uniq_Blood_Count[Analyses$Key == key] <- uniq_blood_count
	} else {
		print(paste0(key," file doesn't exist"))
	}
}

for( key in Analyses$Key ) {
	condition <- Analyses %>% filter(Key==key) %>% pull(Condition) %>% unique()
	individual <- Analyses %>% filter(Key==key) %>% pull(Individual) %>% unique()
	saliva <- gsub(individual,paste0(individual,"S"),key) %>% gsub("kmer_","",.)
	blood <- gsub(individual,paste0(individual,"B"),key) %>% gsub("kmer_","",.)

	file_path <- paste0(dir,"/data_raw/",condition,"/",key,"/data_clean/common_kmers_k-mer_Counts.csv")
	if(file.exists(file_path)) {
		common <- fread(file_path) 
		common_count <- common %>% filter(count <= kmax & count >= kmin) %>% sum()
		Analyses$Common_Count[Analyses$Key == key] <- common_count
	} else {
		print(paste0(key," file doesn't exist"))
	}
}
write.table(Analyses, file = "data_clean/Kmer_Combined_Counts.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

Analyses <- Analyses %>% mutate(SeedLength_Num = SeedLength %>% gsub("k","",.) %>% as.numeric())

####################################################################################################################
###                                   	    Seedlength Graph			                                         ###
####################################################################################################################

summary_fn <- function(x) {
    x <- as.numeric(x)
    tibble(
        n        = sum(!is.na(x)),
        mean     = mean(x, na.rm = TRUE),
        sd       = sd(x, na.rm = TRUE),
        median   = median(x, na.rm = TRUE),
        q25      = quantile(x, 0.25, na.rm = TRUE),
        q75      = quantile(x, 0.75, na.rm = TRUE),
        min      = min(x, na.rm = TRUE),
        max      = max(x, na.rm = TRUE),
    )
}

Analyses_SeedLength <- Analyses %>% filter(Condition == "SeedLength")
cor(Analyses_SeedLength$Common_Count,Analyses_SeedLength$SeedLength_Num)

lm_fit <- lm(Common_Count ~ SeedLength_Num + Individual, data = Analyses_SeedLength)
m <- round(as.numeric(coef(lm_fit)[2]), 2)  # Slope (m)
b <- round(as.numeric(coef(lm_fit)[1]), 2)  # Intercept (b)
r2 <- round(summary(lm_fit)$r.squared, 3)  # R² value
m_pct <- round((m / median(Analyses_SeedLength$Common_Count)) * 100, 10)
message("Slope as % of median: ", m_pct, "% per unit seed length")

individual_slopes <- Analyses_SeedLength %>%
    group_by(Individual) %>%
    summarise(
        slope = coef(lm(Common_Count ~ SeedLength_Num))[2],
        .groups = "drop"
    )

median_slope <- round(median(individual_slopes$slope, na.rm = TRUE))
sd_slope     <- round(sd(individual_slopes$slope, na.rm = TRUE))

message("Median ± SD slope: ", median_slope, " ± ", sd_slope)

### Summaries ###

Analysis <- c("Uniq_Saliva_Count","Uniq_Blood_Count","Common_Count") 
for( i in 1:length(Analysis)){
    analysis <- Analysis[i]
        
	summary <- Analyses_SeedLength %>% group_by(SeedLength_Num) %>% summarise(across(.data[[analysis]], summary_fn, .unpack = TRUE), .groups = "drop")
    write.table(summary, file = paste0("data_clean/Kmer_SeedLength_Summary_",analysis,".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
}

### Saliva Graphs ###

    plot <- ggplot(Analyses_SeedLength, aes(x = SeedLength_Num, y = Uniq_Saliva_Count, group = Individual, colour = Individual)) +
        geom_line(linewidth = 3) +
        labs(x = "BWA Alignment Seed Length", y = paste0("Unique Saliva k-mers (Millions)")) +
        scale_y_continuous(labels = scales::comma_format(scale = 1/1e6)) +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
        ggplot_theme_basic
    ggsave(plot = plot, filename = paste0("graphs/SeedLength_Count_kmer_Uniq_Saliva_Lineplot.pdf"), device = "pdf", width = 45, height = 40, units = "cm", dpi = 640)

    plot <- ggplot(Analyses_SeedLength, aes(x = SeedLength_Num, y = Uniq_Saliva_Count, group = SeedLength_Num)) +
    	geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#1B9E77", fill = "#1B9E77") +  # hide outliers since we show all points
    	geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#1B9E77", fill = "#1B9E77") +            # show all points
        labs(x = "BWA Alignment Seed Length", y = paste0("Unique Saliva k-mers (Millions)")) +
        scale_y_continuous(labels = scales::comma_format(scale = 1/1e6)) +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
        ggplot_theme_basic
    ggsave(plot = plot, filename = paste0("graphs/SeedLength_Count_kmer_Uniq_Saliva_Boxplot.pdf"), device = "pdf", width = 45, height = 40, units = "cm", dpi = 640)

	median_last <- Analyses_SeedLength %>% group_by(Individual) %>%
    	summarise(last_val = last(Uniq_Saliva_Count), .groups = "drop") %>%
    	summarise(median_last = median(last_val)) %>%
    	pull(median_last)
    
    dataframe_Normalised_Percentage <- Analyses_SeedLength %>%
    	group_by(Individual) %>%
    	mutate(
        	Uniq_Saliva_Count_normalized = Uniq_Saliva_Count - last(Uniq_Saliva_Count),
        	Uniq_Saliva_Count_normalized_percentage = Uniq_Saliva_Count_normalized / first(Uniq_Saliva_Count_normalized) * 100
    	) %>%
    	ungroup()

	write.table(dataframe_Normalised_Percentage, file = "data_clean/SeedLength_Count_kmer_Uniq_Saliva_Boxplot_Normalised_Percentage_PlotData.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  
    plot_A <- ggplot(dataframe_Normalised_Percentage, aes(x = SeedLength_Num, y = Uniq_Saliva_Count_normalized_percentage, group = SeedLength_Num)) +
    	geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#1B9E77", fill = "#1B9E77") +  # hide outliers since we show all points
    	geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#1B9E77", fill = "#1B9E77") +            # show all points
        labs(x = "BWA Alignment Seed Length", y = paste0("Normalised Unique Saliva k-mers (%)")) +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
 		scale_y_continuous() +
        ggplot_theme_basic
    ggsave(plot = plot_A, filename = paste0("graphs/SeedLength_Count_kmer_Uniq_Saliva_Boxplot_Normalised_Percentage.pdf"), device = "pdf", width = 40, height = 40, units = "cm", dpi = 640)

### Blood Graphs ###

    plot <- ggplot(Analyses_SeedLength, aes(x = SeedLength_Num, y = Uniq_Blood_Count, group = Individual, colour = Individual)) +
        geom_line(linewidth = 3) +
        labs(x = "BWA Alignment Seed Length", y = paste0("Unique Blood k-mers (Millions)")) +
        scale_y_continuous(labels = scales::comma_format(scale = 1/1e6)) +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
        ggplot_theme_basic
    ggsave(plot = plot, filename = paste0("graphs/SeedLength_Count_kmer_Uniq_Blood_Lineplot.pdf"), device = "pdf", width = 45, height = 40, units = "cm", dpi = 640)

	lm_fit <- lm(Uniq_Blood_Count ~ 0 + Individual, data = Analyses_SeedLength)
	coefs <- coef(lm_fit)
	individual_coefs <- tibble(Individual = gsub("Individual", "", names(coefs[grep("Individual", names(coefs))])),coef_adjustment = coefs[grep("Individual", names(coefs))])
    median_raw <- median(Analyses_SeedLength$Uniq_Blood_Count, na.rm = TRUE)
    Analyses_SeedLength_Adjusted <- Analyses_SeedLength %>% left_join(individual_coefs, by = "Individual") %>% mutate(Uniq_Blood_Count_adjusted = Uniq_Blood_Count - coef_adjustment + median_raw)
    
    plot <- ggplot(Analyses_SeedLength_Adjusted, aes(x = SeedLength_Num, y = Uniq_Blood_Count_adjusted, group = SeedLength_Num)) +
    	geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#1B9E77", fill = "#1B9E77") +  # hide outliers since we show all points
    	geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#1B9E77", fill = "#1B9E77") +            # show all points
        labs(x = "BWA Alignment Seed Length", y = paste0("Normalised Unique Blood k-mers (Millions)")) +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
        scale_y_continuous(labels = scales::comma_format(scale = 1/1e6)) +
        ggplot_theme_basic
    ggsave(plot = plot, filename = paste0("graphs/SeedLength_Count_kmer_Uniq_Blood_Boxplot_Normalised.pdf"), device = "pdf", width = 45, height = 40, units = "cm", dpi = 640)

### Common Graphs ###

    plot <- ggplot(Analyses_SeedLength, aes(x = SeedLength_Num, y = Common_Count, group = Individual, colour = Individual)) +
        geom_line(linewidth = 3) +
        labs(x = "BWA Alignment Seed Length", y = paste0("Common k-mers (Millions)")) +
        scale_y_continuous(labels = scales::comma_format(scale = 1/1e6)) +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
        ggplot_theme_basic
    ggsave(plot = plot, filename = paste0("graphs/SeedLength_Count_kmer_Common_Lineplot.pdf"), device = "pdf", width = 45, height = 40, units = "cm", dpi = 640)

	lm_fit <- lm(Common_Count ~ 0 + Individual, data = Analyses_SeedLength)
	coefs <- coef(lm_fit)
	individual_coefs <- tibble(Individual = gsub("Individual", "", names(coefs[grep("Individual", names(coefs))])),coef_adjustment = coefs[grep("Individual", names(coefs))])
    median_raw <- median(Analyses_SeedLength$Common_Count, na.rm = TRUE)
    Analyses_SeedLength_Adjusted <- Analyses_SeedLength %>% left_join(individual_coefs, by = "Individual") %>% mutate(Common_Count_adjusted = Common_Count - coef_adjustment + median_raw)
    
	write.table(Analyses_SeedLength_Adjusted, file = "data_clean/SeedLength_Count_kmer_Common_Boxplot_Normalised_PlotData.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    
    plot_B <- ggplot(Analyses_SeedLength_Adjusted, aes(x = SeedLength_Num, y = Common_Count_adjusted, group = SeedLength_Num)) +
    	geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA, color = "#D95F02", fill = "#D95F02") +  # hide outliers since we show all points
    	geom_jitter(width = 0.01, size = 2, alpha = 0.7, color = "#D95F02", fill = "#D95F02") +            # show all points
        labs(x = "BWA Alignment Seed Length", y = paste0("Adjusted Common k-mers (Millions)")) +
		scale_x_continuous(breaks = seq(19, 27), labels = seq(19, 27)) +
		scale_y_continuous(labels = function(x) formatC(x/1e6, format = "f", digits = 1), expand = expansion(mult = c(0.1, 0.07))) +
        ggplot_theme_basic
    ggsave(plot = plot_B, filename = paste0("graphs/SeedLength_Count_kmer_Common_Boxplot_Normalised.pdf"), device = "pdf", width = 45, height = 40, units = "cm", dpi = 640)
  
### MultiPlot ###
 
plot_A <- plot_A + labs(tag = "A)") + theme(plot.tag = element_text(size = 28, face = "bold"))
plot_B <- plot_B + labs(tag = "B)") + theme(plot.tag = element_text(size = 28, face = "bold"))
plot <- plot_A + plot_B  + plot_layout(ncol = 2, widths = c(0.95, 1))
ggsave(plot=plot, filename="graphs/SeedLength_Count_kmer_MultiPlot.pdf", device="pdf", width = 55, height = 30, units = "cm", dpi = 640)

plot_A <- plot_A + labs(tag = "A)") + theme(plot.tag = element_text(size = 28, face = "bold"), axis.title.y = element_text(size=28, margin = margin(t = 0, r = 30, b = 0, l = 20)))
plot_B <- plot_B + labs(tag = "B)") + theme(plot.tag = element_text(size = 28, face = "bold"), axis.title.y = element_text(size=28, margin = margin(t = 0, r = 30, b = 0, l = 20)))
plot <- plot_A + plot_B  + plot_layout(ncol = 1, axes = "collect", widths = c(1, 1))
ggsave(plot=plot, filename="graphs/SeedLength_Count_kmer_MultiPlot_Horizontal.pdf", device="pdf", width = 40, height = 40, units = "cm", dpi = 640) 

####################################################################################################################
###                                              	 fin			                                            ####
####################################################################################################################

quit()
