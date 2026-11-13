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
if(!require(ggpubr)){
    install.packages("ggpubr")
    library(ggpubr)
}
if(!require(patchwork)){
    install.packages("patchwork")
    library(patchwork)
}
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
###                                               Load in                                                       ####
####################################################################################################################

# GTEx
GTEx <- read.table( file = paste0(GENRES_DIR, "/Resources/GTEx/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count_avg.gct"))
colnames(GTEx) <- c("Ensembl_Transcript_ID","Ensembl_Gene_ID","Avg_Transcript_Count")
sub("\\..*","\\",GTEx$Ensembl_Transcript_ID) -> GTEx$Ensembl_Transcript_ID

# gnomAD Constraint
gnomAD <- read.table( file = paste0(GENRES_DIR, "/Resources/gnomAD/4.1/gnomad.v4.1.constraint_metrics.tsv"), header = TRUE,sep = "\t")

####################################################################################################################
###                                               Loadin Data Het                                               ####
####################################################################################################################

files <- list.files("data_raw/", pattern = "h.*_out\\.tsv$") #%>% keep(~ str_detect(.x, "trained_taxid_k24")) 
file_tags <- files %>% gsub("_SoftClip.*","",.)
ids <- files %>% gsub("_SoftClip.*","",.) %>% gsub("_rare_coding.*","",.) %>% gsub(".*Norm_ann_","",.)

rm(list = ls(pattern="^df_*"))
Combined=NULL
for(file in files) {    
    lines <- readLines(paste0("data_raw/",file), n = 2)

    if (length(lines) == 0) {
    	print(paste0(file," File Missing"))
	} else if (length(lines) == 1) {
    	print(paste0(file," File is Empty"))	
	} else {
        data <- fread(file = paste0("data_raw/",file), sep = '\t', header = TRUE, colClasses = c("character"), quote = "", stringsAsFactors = FALSE) %>% mutate(across(everything(), ~ na_if(.x, ""))) %>% select(where(~ !all(is.na(.))))
 		id <- ids[file == files]
      	colnames(data) <- colnames(data) %>% gsub(".*]","",.) %>% gsub(":","_",.) 

		data$Individual <- id  
		data$Set <- file %>% gsub(".*_coding_","",.) %>% gsub("_out.tsv","",.)
		data$Tag <- file_tags[file == files] %>% gsub("_Merged_Norm_ann.*","",.)
		names <- colnames(data)
		fields <- names[str_detect(names, id)] 
		for(field in fields) {    
			new_field <- field %>% gsub(id,"Individual",.)
			data <- data %>% mutate(!!new_field := !!sym(field))
		}
		data_processed <- data %>%
 			separate(Individual_AD, into = c("First", "Second"), sep = ",", convert = TRUE) %>%
  			mutate(Individual_VAF := Second / (First + Second)) %>%
  			mutate(Individual_AD := paste0(First,",",Second)) %>%
  			select(-First, -Second) %>%
  			filter(Individual_VAF >= 0.2)

        if (dim(data_processed)[1] != 0) {
	        assign( paste0("df_",file), data_processed )
        }
    }
}
ls(pattern="^df_*") -> listOfDataFrames
bind_rows(mget(listOfDataFrames)) -> Combined

write.table(Combined, file="data_clean/Genome_Singleton_Coding.tsv", quote=FALSE, sep='\t', row.names=F)

##########################################################################################################
###########################           Data Cleaning                 ######################################
##########################################################################################################

#############################
#### Split arrays ####
#############################

colnames(Combined) -> names
names[str_detect(names, "DecipherG2P")] -> DecipherG2Ps
Combined %>%  separate_rows(CSQ,sep=",") %>%  separate_rows(DecipherG2Ps,sep=",")  -> Combined_Split

###############################
### Split CSQ field  ####
###############################

column_names <- sprintf("CSQ%d", 1:53)
Combined_Split  %>% tidyr::separate(col=CSQ, sep ='\\|',into=column_names, convert = TRUE, remove = FALSE, fill = "right") -> Combined_Split_CSQ

apply(Combined_Split_CSQ,2,as.character) %>% as.data.frame -> Combined_Split_CSQ
Combined_Split_CSQ$gnomADv4_1_AF_grpmax_joint <- as.numeric(as.character(Combined_Split_CSQ$gnomADv4_1_AF_grpmax_joint))
Combined_Split_CSQ$gnomADv4_1_AF_grpmax_joint[is.na(Combined_Split_CSQ$gnomADv4_1_AF_grpmax_joint)] <- 0

Combined_Split_CSQ %>% filter(CSQ33 == "Ensembl" & CSQ8 == "protein_coding" & CSQ22 != "Clone_based_ensembl_gene") %>% filter(CSQ21 == "") %>%
     filter(gnomADv4_1_AF_grpmax_joint <= 0.001) -> Combined_Split_CSQ_Filtered
Combined_Split_CSQ_Filtered$DecipherG2P_gene_symbol <- as.character(Combined_Split_CSQ_Filtered$DecipherG2P_gene_symbol)
Combined_Split_CSQ_Filtered$CSQ4 <- as.character(Combined_Split_CSQ_Filtered$CSQ4)

write.table(Combined_Split_CSQ_Filtered, file="data_clean/Genome_Singleton_Coding_Filtered.tsv", quote=FALSE, sep='\t', row.names=F)

###########################################################################
#########    Filter Transcripts based on gnomAD Canonical        ##########
###########################################################################

Combined_Split_CSQ_Filtered$gnomADv4_1_AC_joint <- as.numeric(as.character(Combined_Split_CSQ_Filtered$gnomADv4_1_AC_joint))
Combined_Split_CSQ_Filtered$gnomADv4_1_AC_joint[is.na(Combined_Split_CSQ_Filtered$gnomADv4_1_AC_joint)] <- 0
Combined_Split_CSQ_Filtered_Canonical <- Combined_Split_CSQ_Filtered

factor(Combined_Split_CSQ_Filtered_Canonical$CSQ3,levels=c("HIGH","MODERATE","LOW","MODIFIER")) -> Combined_Split_CSQ_Filtered_Canonical$CSQ3

gnomAD$mane_select[match(Combined_Split_CSQ_Filtered_Canonical$CSQ7,gnomAD$transcript)] -> Combined_Split_CSQ_Filtered_Canonical$MANE_Transcript_Canonical
gnomAD$lof.pLI[match(Combined_Split_CSQ_Filtered_Canonical$CSQ7,gnomAD$transcript)] -> Combined_Split_CSQ_Filtered_Canonical$gnomADv4_Transcript_pLI
gnomAD$lof.z_score[match(Combined_Split_CSQ_Filtered_Canonical$CSQ7,gnomAD$transcript)] -> Combined_Split_CSQ_Filtered_Canonical$gnomADv4_Transcript_lof_z
gnomAD$mis.z_score[match(Combined_Split_CSQ_Filtered_Canonical$CSQ7,gnomAD$transcript)] -> Combined_Split_CSQ_Filtered_Canonical$gnomADv4_Transcript_mis_z
gnomAD$syn.z_score[match(Combined_Split_CSQ_Filtered_Canonical$CSQ7,gnomAD$transcript)] -> Combined_Split_CSQ_Filtered_Canonical$gnomADv4_Transcript_syn_z

Combined_Split_CSQ_Filtered_Canonical_output <- Combined_Split_CSQ_Filtered_Canonical %>% arrange(desc(MANE_Transcript_Canonical)) %>% distinct(CHROM, POS,Tag,Individual, .keep_all = TRUE) %>% filter(CSQ3 == "MODERATE" | CSQ3 == "HIGH")

write.table(Combined_Split_CSQ_Filtered_Canonical_output, file="data_clean/Genome_Singleton_Coding_Filtered_Canonical.tsv", quote=FALSE, sep='\t', row.names=F)

###########################################################################
#########    					Graphing			 	         ##########
###########################################################################

Group_Tally <- Combined_Split_CSQ_Filtered_Canonical_output %>% filter(Individual_GT == "0/1") %>% group_by(Individual,Tag) %>% tally() %>% mutate(Group = ifelse(Tag == "NoFastP_Kraken_trained_taxid_k24", "Cleanup","No Cleanup"))
readr::write_tsv(Group_Tally, "data_clean/Comparison_HetCandidates.tsv")
Group_Tally_Summary <- Group_Tally %>% group_by(Group) %>% summarise( count = n(), mean = mean(n), median = median(n), sd = sd(n), min = min(n), max = max(n), iqr = IQR(n), q25 = quantile(n, 0.25), q75 = quantile(n, 0.75), .groups = "drop")
readr::write_tsv(Group_Tally_Summary, "data_clean/Comparison_HetCandidates_Summary.tsv")

#### Stats Test ####

test <- Group_Tally %>% group_by(Group) %>% summarise(method = shapiro.test(n)$method, p_value = shapiro.test(n)$p.value, W = shapiro.test(n)$statistic)
readr::write_tsv(test, "data_clean/Comparison_HetCandidates_Shapiro.tsv")
	# Not normal - use wilcoxon
#test <- wilcox.test(n ~ Group, data = Group_Tally, paired = TRUE,conf.int = TRUE)
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA,
  p_value = format(test$p.value, scientific = FALSE, digits = 6), estimate = as.numeric(test$estimate), conf_low = test$conf.int[1], conf_high = test$conf.int[2], alternative = test$alternative
)
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_Wilcoxon.tsv")

#### Graph ####

Group_Tally <- Group_Tally %>% mutate(Group = factor(Group,levels=c("No Cleanup","Cleanup")))

plot <- ggplot(Group_Tally, aes(x = Group, y = n, group = Individual)) +
  geom_line(aes(color = Individual), alpha = 0.5) +
  geom_point() +
	guides(alpha = "none",fill = "none",colour = "none") +
  ggplot_theme_basic +
  labs(y = "Number of Candidate Variants",x=NULL) 
ggsave(plot=plot,filename=paste0("graphs/Validation_CleanUp_Comp-Individual.pdf"),device="pdf", width = 20, height = 20, units = "cm",dpi = 640)

Group_Tally <- Group_Tally %>% mutate(Group_Plot = ifelse(Group == "Cleanup", "Trained Kraken\nSeed Length 24", "No Kraken\n Seed Length 19"))

plot_group_boxplot <- ggplot(Group_Tally, aes(y = n, x = Group_Plot, color = Group_Plot, fill = Group_Plot)) +
    geom_boxplot(width = 0.6, alpha = 0.3, outlier.shape = NA) +  # hide outliers since we show all points
    geom_jitter(width = 0.15, size = 2, alpha = 0.7) +            # show all points
    labs(x = "Cleanup Method (n=21)", y = "Number of Candidate Variants") +
    stat_compare_means(method = "wilcox.test", comparisons = list(c("Trained Kraken\nSeed Length 24", "No Kraken\n Seed Length 19")), label = "p.signif", size = 7) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    guides(alpha = "none", fill = "none", colour = "none") +
    ggplot_theme_basic
ggsave(plot=plot_group_boxplot,filename=paste0("graphs/Validation_CleanUp_Comp-Group.pdf"),device="pdf", width = 25, height = 30, units = "cm",dpi = 640)

###########################################################################
#########    			Variant Comparisons			 	         ##########
###########################################################################

Combined_Split_CSQ_Filtered_Canonical_output_Post <- Combined_Split_CSQ_Filtered_Canonical_output %>% filter(Individual_GT == "0/1") %>% filter(Tag == "NoFastP_Kraken_trained_taxid_k24")
Combined_Split_CSQ_Filtered_Canonical_output_Pre <- Combined_Split_CSQ_Filtered_Canonical_output %>% filter(Individual_GT == "0/1") %>% filter(Tag == "NoFastP_NoKraken_k19")

Combined_Split_CSQ_Filtered_Canonical_output_Post_new <- anti_join(Combined_Split_CSQ_Filtered_Canonical_output_Post,Combined_Split_CSQ_Filtered_Canonical_output_Pre, by=c("Individual","Individual_GT","CHROM","POS","REF","ALT"))
Combined_Split_CSQ_Filtered_Canonical_output_Pre_new <- anti_join(Combined_Split_CSQ_Filtered_Canonical_output_Pre,Combined_Split_CSQ_Filtered_Canonical_output_Post, by=c("Individual","Individual_GT","CHROM","POS","REF","ALT"))
Combined_Split_CSQ_Filtered_Canonical_output_Overlap <- inner_join(Combined_Split_CSQ_Filtered_Canonical_output_Post, Combined_Split_CSQ_Filtered_Canonical_output_Pre,by = c("Individual", "Individual_GT", "CHROM", "POS", "REF", "ALT", "CSQ4", "CSQ3","gnomADv4_1_AC_joint","AC","DecipherG2P_disease_name","gnomADLOFmetricsbyGene_pLI"))

write.table(Combined_Split_CSQ_Filtered_Canonical_output_Post_new, file="data_clean/Genome_Singleton_Coding_Filtered_Canonical_PostUnique.tsv", quote=FALSE, sep='\t', row.names=F)
write.table(Combined_Split_CSQ_Filtered_Canonical_output_Pre_new, file="data_clean/Genome_Singleton_Coding_Filtered_Canonical_PreUnique.tsv", quote=FALSE, sep='\t', row.names=F)
write.table(Combined_Split_CSQ_Filtered_Canonical_output_Overlap, file="data_clean/Genome_Singleton_Coding_Filtered_Canonical_Overlap.tsv", quote=FALSE, sep='\t', row.names=F)

Combined_Split_CSQ_Filtered_Canonical_output_Post_new$Set <- "New"
Combined_Split_CSQ_Filtered_Canonical_output_Pre_new$Set <- "Removed"
Combined_Split_CSQ_Filtered_Canonical_output_Overlap$Set <- "Both"

Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined <- bind_rows(
  Combined_Split_CSQ_Filtered_Canonical_output_Post_new,
  Combined_Split_CSQ_Filtered_Canonical_output_Pre_new,
  Combined_Split_CSQ_Filtered_Canonical_output_Overlap
) %>% mutate(CSQ3 = as.character(CSQ3))

summary_table <- Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined %>% group_by(Set) %>% summarise(count = n(), .groups = "drop")
readr::write_tsv(summary_table, "data_clean/Comparison_HetCandidates_SummaryNumbers.tsv")

## HIGH Impact

summary_table <- Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined %>% group_by(Set, CSQ3) %>% summarise(count = n(), .groups = "drop") %>% group_by(Set) %>% mutate(proportion = count / sum(count))
contingency <- table(Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined$Set, Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined$CSQ3)

test <- chisq.test(contingency[c("Both", "New"), ])
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA, p_value = format(test$p.value, scientific = FALSE, digits = 6))
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_Impact_ChiSquared_Added.tsv")

test <- chisq.test(contingency[c("Both", "Removed"), ])
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA, p_value = format(test$p.value, scientific = FALSE, digits = 6))
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_Impact_ChiSquared_Removed.tsv")

## Novelty

Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined <- Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined %>% mutate(GnomAD_Novel = gnomADv4_1_AC_joint == 0)

summary_table <- Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined %>% group_by(Set, GnomAD_Novel) %>% summarise(count = n(), .groups = "drop") %>% group_by(Set) %>% mutate(proportion = count / sum(count))
contingency <- table(Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined$Set, Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined$GnomAD_Novel)

test <- chisq.test(contingency[c("Both", "New"), ])
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA, p_value = format(test$p.value, scientific = FALSE, digits = 6))
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_Novelty_ChiSquared_Added.tsv")

test <- chisq.test(contingency[c("Both", "Removed"), ])
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA, p_value = format(test$p.value, scientific = FALSE, digits = 6))
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_Novelty_ChiSquared_Removed.tsv")

## Shared

Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined <- Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined %>% mutate(Shared = AC > 1)

summary_table <- Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined %>% group_by(Set, Shared) %>% summarise(count = n(), .groups = "drop") %>% group_by(Set) %>% mutate(proportion = count / sum(count))
contingency <- table(Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined$Set, Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined$Shared)

test <- chisq.test(contingency[c("Both", "New"), ])
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA, p_value = format(test$p.value, scientific = FALSE, digits = 6))
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_Shared_ChiSquared_Added.tsv")

test <- chisq.test(contingency[c("Both", "Removed"), ])
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA, p_value = format(test$p.value, scientific = FALSE, digits = 6))
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_Shared_ChiSquared_Removed.tsv")

## Disease Gene

Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined <- Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined %>% mutate(DiseaseGene = DecipherG2P_disease_name != ".")

summary_table <- Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined %>% group_by(Set, DiseaseGene) %>% summarise(count = n(), .groups = "drop") %>% group_by(Set) %>% mutate(proportion = count / sum(count))
contingency <- table(Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined$Set, Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined$DiseaseGene)

test <- chisq.test(contingency[c("Both", "New"), ])
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA, p_value = format(test$p.value, scientific = FALSE, digits = 6))
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_DiseaseGene_ChiSquared_Added.tsv")

test <- chisq.test(contingency[c("Both", "Removed"), ])
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA, p_value = format(test$p.value, scientific = FALSE, digits = 6))
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_DiseaseGene_ChiSquared_Removed.tsv")

## Constraint

Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined <- Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined %>% mutate(ConstraintGene = str_split(gnomADLOFmetricsbyGene_pLI, ",") %>% map_lgl(~ any(suppressWarnings(as.numeric(.x)) > 0.9, na.rm = TRUE)) )

summary_table <- Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined %>% group_by(Set, ConstraintGene) %>% summarise(count = n(), .groups = "drop") %>% group_by(Set) %>% mutate(proportion = count / sum(count))
contingency <- table(Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined$Set, Combined_Split_CSQ_Filtered_Canonical_output_Post_Combined$ConstraintGene)

test <- chisq.test(contingency[c("Both", "New"), ])
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA, p_value = format(test$p.value, scientific = FALSE, digits = 6))
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_ConstraintGene_ChiSquared_Added.tsv")

test <- chisq.test(contingency[c("Both", "Removed"), ])
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA, p_value = format(test$p.value, scientific = FALSE, digits = 6))
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_ConstraintGene_ChiSquared_Removed.tsv")

###########################################################################
#########    				Recurrent FP Genes		 	         ##########
###########################################################################

RFP <- read.table( file = paste0(PROJECT_DIR,"/Validation/Resources/Field_2019_RFP.tsv"), header=T )

Combined_Split_CSQ_Filtered_Canonical_output_Post_new_noRFP <- Combined_Split_CSQ_Filtered_Canonical_output_Post_new %>% filter( !str_split(CSQ4, ",") %>% map_lgl(~ any(.x %in% RFP$Gene) ) )
Combined_Split_CSQ_Filtered_Canonical_output_Pre_new_noRFP <- Combined_Split_CSQ_Filtered_Canonical_output_Pre_new %>% filter( !str_split(CSQ4, ",") %>% map_lgl(~ any(.x %in% RFP$Gene) ) )
Combined_Split_CSQ_Filtered_Canonical_output_Overlap_noRFP <- Combined_Split_CSQ_Filtered_Canonical_output_Overlap %>% filter( !str_split(CSQ4, ",") %>% map_lgl(~ any(.x %in% RFP$Gene) ) )

Post <- data.frame(Set = "Added", No_Filt = nrow(Combined_Split_CSQ_Filtered_Canonical_output_Post_new), RFP_Filt = nrow(Combined_Split_CSQ_Filtered_Canonical_output_Post_new_noRFP))
Pre <- data.frame(Set = "Removed", No_Filt = nrow(Combined_Split_CSQ_Filtered_Canonical_output_Pre_new), RFP_Filt = nrow(Combined_Split_CSQ_Filtered_Canonical_output_Pre_new_noRFP))
Both <- data.frame(Set = "Both", No_Filt = nrow(Combined_Split_CSQ_Filtered_Canonical_output_Overlap), RFP_Filt = nrow(Combined_Split_CSQ_Filtered_Canonical_output_Overlap_noRFP))

combined_counts <- rbind(Post, Pre, Both)

combined_counts %>% mutate(PropRemaining = RFP_Filt / No_Filt)
readr::write_tsv(combined_counts, "data_clean/Comparison_HetCandidates_RFP.tsv")
write.table(Combined_Split_CSQ_Filtered_Canonical_output_Post_new_noRFP, file="data_clean/Genome_Singleton_Coding_Filtered_Canonical_PostUnique_RFP.tsv", quote=FALSE, sep='\t', row.names=F)
write.table(Combined_Split_CSQ_Filtered_Canonical_output_Pre_new_noRFP, file="data_clean/Genome_Singleton_Coding_Filtered_Canonical_PreUnique_RFP.tsv", quote=FALSE, sep='\t', row.names=F)

###########################################################################
#########    					Contam v Diff		 	         ##########
###########################################################################

for(id in unique(ids)) {
        rm(list = ls(pattern="^kraken_*"))
        path <- paste0(PROJECT_DIR,"/Validation/Alignments/CMM_Alignments/",id,"_NoFastP_Kraken_trained_taxid_k24/Kraken2/")
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
        Individual <- id

        assign( paste0("percentage_",id), data.frame(Individual,Human_Percentage))
}
bind_rows(mget(ls(pattern="^percentage_*"))) -> combined_human_percentages

combined_diff <- Group_Tally %>%
  select(Individual, Group, n) %>%
  pivot_wider(names_from = Group, values_from = n) %>%
  mutate(Diff = `No Cleanup` - Cleanup )
  
combined_diffvcontam <- left_join(combined_diff,combined_human_percentages,by="Individual") %>% mutate( Microbiome_Percentage = (1 - Human_Percentage) * 100 )
readr::write_tsv(combined_diffvcontam, "data_clean/Comparison_HetCandidates_Contamination.tsv")

lm_fit <- lm(Diff ~ Microbiome_Percentage, data = combined_diffvcontam)
m <- round(as.numeric(coef(lm_fit)[2]), 2)  # Slope (m)
b <- round(as.numeric(coef(lm_fit)[1]), 2)  # Intercept (b)
r2 <- round(summary(lm_fit)$r.squared, 3)  # R² value
equation <- if (b < 0) {
  as.character(as.expression(
    substitute(italic(y) == m ~ italic(x) ~ "-" ~ abs_b * "," ~~ italic(R)^2 ~ "=" ~ r2, 
               list(m = m, abs_b = abs(b), r2 = r2))))
} else {
  as.character(as.expression(
    substitute(italic(y) == m ~ italic(x) ~ "+" ~ b * "," ~~ italic(R)^2 ~ "=" ~ r2, 
               list(m = m, b = b, r2 = r2))))
}

test <- cor.test(combined_diffvcontam$Microbiome_Percentage, combined_diffvcontam$Diff, method = "pearson")
test_out <- tibble::tibble( method = test$method, statistic = as.numeric(test$statistic), df = if (!is.null(test$parameter)) as.numeric(test$parameter) else NA,
  p_value = format(test$p.value, scientific = FALSE, digits = 6), estimate = as.numeric(test$estimate), conf_low = test$conf.int[1], conf_high = test$conf.int[2], alternative = test$alternative
)
readr::write_tsv(test_out, "data_clean/Comparison_HetCandidates_ContaminationCorrelation.tsv")

#RColorBrewer::brewer.pal(8, "Dark2")
plot_contam_line <- ggplot(combined_diffvcontam, aes(x = Microbiome_Percentage, y = Diff)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "#1B9E77") +  # Add trendline
  annotate("text", x = max(combined_diffvcontam$Microbiome_Percentage) * 0.7, 
                 y = max(combined_diffvcontam$Diff) * 0.9, 
                 label = equation, parse = TRUE, size = 6) +  # Add formula
  labs(x = "Estimated Microbiome Sequence Presence (%)", 
       y = "Reduction in Candidate Variants\nAfter Cleanup") +
  ggplot_theme_basic
ggsave(plot=plot_contam_line, filename="graphs/Validation_CleanUp_Comp_Contamination_vs_Effect.pdf", device="pdf", width = 30, height = 30, units = "cm", dpi = 640)

####################################################################################################################
###                                              Multiple Plot Graph                                            ####
####################################################################################################################

plot_group_boxplot <- plot_group_boxplot + labs(tag = "A)") + theme(plot.tag = element_text(size = 28, face = "bold"))
plot_contam_line <- plot_contam_line + labs(tag = "B)") + theme(plot.tag = element_text(size = 28, face = "bold"))
plot <- plot_group_boxplot + plot_contam_line  + plot_layout(ncol = 2, axes = "collect", widths = c(2, 3))

ggsave(plot=plot, filename="graphs/Validation_Combined_Graph.pdf", device="pdf", width = 55, height = 30, units = "cm", dpi = 640)

####################################################################################################################
###                                              	 fin			                                            ####
####################################################################################################################

quit()
