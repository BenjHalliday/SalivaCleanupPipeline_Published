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
topdir <- paste0(PROJECT_DIR,"/Training/")
dir <- paste0(topdir,"Training_Analyses")
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
###                                               Meta Data	                                                    ####
####################################################################################################################

fqd <- fread(file = paste0(topdir,"/Meta_Data/Training_Sequenced_Samples_Info.tsv"), sep = '\t')
fqd_main <- fread(file = paste0(META_DIR,"/Comparison_Sequenced_Samples_Cleanup.tsv"), sep = '\t')

####################################################################################################################
###                                               Kraken2 Data	                                                ####
####################################################################################################################

for(number in fqd$Order) {
	line <- fqd %>% filter(Order == number)
	file <- paste0(topdir,"/Kraken2/",line$Individual_ID,"/",line$Individual,"_",line$DNA_Number,"_",line$Library_Number,"_",line$Run_Number,"/",line$Filename1,"_Report.tsv")
	if(file.exists(file)) {
	    if(file.info(file)$size > 0){
	        data <- fread(file, sep = '\t') %>% filter(V6 %in% c("unclassified","root","Eukaryota","Bacteria","Archaea","Viruses","other entries") ) %>% 
	        	select(V2,V6) %>% pivot_wider( names_from = V6, values_from = V2 ) %>% mutate(FastQ = paste0(line$Individual,"_",line$DNA_Number,"_",line$Library_Number,"_",line$Run_Number)) %>%
	        	mutate(Individual_ID = line$Individual)
	        assign(paste0("kraken_",number),data)
	    } else {
	        print(paste0(file," is empty!!!"))
	    } 
	} else {
	    print(paste0(file," is missing!!!"))
	}
}
ls(pattern=paste0("^kraken_*")) -> listOfDataFrames
bind_rows(mget(listOfDataFrames)) -> Combined_Kraken

write.table(Combined_Kraken, file=paste0("data_clean/Kraken_Stats.tsv"), quote=FALSE, sep='\t', row.names=F)

Combined_Kraken_Summary <- Combined_Kraken %>% group_by(Individual_ID) %>% summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>% 
	mutate(total = unclassified + root) %>% mutate(unclassified_percentage = (unclassified / total)*100) %>% mutate(Eukaryota_percentage = (Eukaryota / total)*100) %>% 
	mutate(Bacteria_percentage = (Bacteria / total)*100)

write.table(Combined_Kraken_Summary, file=paste0("data_clean/Kraken_Stats_Summary.tsv"), quote=FALSE, sep='\t', row.names=F)

Combined_Kraken_Summary_CI <- Combined_Kraken_Summary %>% select(unclassified_percentage,Eukaryota_percentage,Bacteria_percentage) %>% summarise(
    across(
      where(is.numeric),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE),
        lower_CI = ~mean(.x, na.rm = TRUE) - qt(0.975, df = length(.x) - 1) * sd(.x, na.rm = TRUE) / sqrt(length(.x)),
        upper_CI = ~mean(.x, na.rm = TRUE) + qt(0.975, df = length(.x) - 1) * sd(.x, na.rm = TRUE) / sqrt(length(.x))
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

write.table(Combined_Kraken_Summary_CI, file=paste0("data_clean/Kraken_Stats_Summary_CI.tsv"), quote=FALSE, sep='\t', row.names=F)

cor_res <- cor.test(~ unclassified_percentage + Bacteria_percentage, data = Combined_Kraken_Summary,method = "pearson")

tibble(
  method = cor_res$method,
  estimate = cor_res$estimate,
  statistic = cor_res$statistic,
  p_value = cor_res$p.value,
  conf_low = cor_res$conf.int[1],
  conf_high = cor_res$conf.int[2],
  n = cor_res$parameter + 2  # df = n - 2
) %>% write.table(file=paste0("data_clean/Kraken_Stats_Summary_CorTest.tsv"), quote=FALSE, sep='\t', row.names=F)

####################################################################################################################
###                                               Megahit Data	                                                ####
####################################################################################################################

for(id in unique(fqd$Individual_ID)) {
	file <- paste0(topdir,"/MEGAHIT/",id,"/final.contigs.fa")
	if(file.exists(file)) {
	    if(file.info(file)$size > 0){
        
        	data <- read_lines(file) %>% keep(~ grepl("^>", .x)) %>% sub("^>", "", .) %>% as.data.frame() %>% rename(header = ".") %>% 
        		separate(header, into = paste0("col", 1:4), sep = " ", fill = "right") %>% mutate(across(everything(), ~ sub(".*=", "", .))) %>%
        		rename(Contig=col1,Flag=col2,Multiplicity=col3,Length=col4) %>% mutate(Individual_ID = id)
        
	        assign(paste0("megahit_",id),data)
	    } else {
	        print(paste0(file," is empty!!!"))
	    } 
	} else {
	    print(paste0(file," is missing!!!"))
	}
}
ls(pattern=paste0("^megahit_*")) -> listOfDataFrames
bind_rows(mget(listOfDataFrames)) -> Combined_Megahit

write.table(Combined_Megahit, file=paste0("data_clean/Megahit_Stats.tsv"), quote=FALSE, sep='\t', row.names=F)

summary_fn <- function(x) {
    x <- as.numeric(x)
    n      <- sum(!is.na(x))
    m      <- mean(x, na.rm = TRUE)
    z      <- median(x, na.rm = TRUE)
    s      <- sd(x, na.rm = TRUE)
    q25    <- quantile(x, 0.25, na.rm = TRUE)
    q75    <- quantile(x, 0.75, na.rm = TRUE)
    tibble(n = n, mean = m, sd = s, median = z, q25 = q25, q75 = q75,
           min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
    }

# per individual
Combined_Megahit_Summary <- Combined_Megahit %>% mutate(Individual_ID = as.character(Individual_ID)) %>% mutate(across(c(Multiplicity, Length), as.numeric)) %>%
    group_by(Individual_ID) %>% summarise(across(c(Multiplicity, Length), summary_fn, .unpack = TRUE), .groups = "drop")

# total summary
Combined_Megahit_Summary_Total <- Combined_Megahit %>% mutate(across(c(Multiplicity, Length), as.numeric)) %>%
    mutate(Individual_ID = "Total") %>% group_by(Individual_ID) %>%
    summarise(across(c(Multiplicity, Length), summary_fn, .unpack = TRUE), .groups = "drop")

# combined
Combined_Megahit_Summary_All <- bind_rows(Combined_Megahit_Summary, Combined_Megahit_Summary_Total)

write.table(Combined_Megahit_Summary_All, file=paste0("data_clean/Megahit_Stats_Summary.tsv"), quote=FALSE, sep='\t', row.names=F)

####################################################################################################################
###                                               BLAST Data	                                                ####
####################################################################################################################

for(id in unique(fqd$Individual_ID)) {
	file <- paste0(topdir,"/BLAST/",id,"/megahit_filtered_contigs_blastn-nt_hits_annotated.tsv")
	if(file.exists(file)) {
	    if(file.info(file)$size > 0){
			data <- fread(file, sep = '\t') %>% mutate(Individual_ID = id)
	        assign(paste0("blast_",id),data)
	    } else {
	        print(paste0(file," is empty!!!"))
	    } 
	} else {
	    print(paste0(file," is missing!!!"))
	}
}
ls(pattern=paste0("^blast_*")) -> listOfDataFrames
bind_rows(mget(listOfDataFrames)) %>% rename(Contig = "Query_Contig") %>%
	left_join(Combined_Megahit,., by = c("Contig","Individual_ID" ),relationship = "many-to-many") %>%
	mutate(Length = as.numeric(Length)) %>% mutate(Multiplicity = as.numeric(Multiplicity)) %>% mutate(Hit = !is.na(Subject)) -> Combined_BLAST

write.table(Combined_BLAST, file=paste0("data_clean/BLAST_Stats.tsv"), quote=FALSE, sep='\t', row.names=F)

res <- wilcox.test(Length ~ Hit, data = Combined_BLAST)
write.table(data.frame(statistic=res$statistic, p_value=res$p.value, method=res$method, alternative=res$alternative), file=paste0("data_clean/BLAST_Stats_Wilcoxon_HitvLength.tsv"), quote=FALSE, sep='\t', row.names=F)
Combined_BLAST %>% group_by(Hit) %>% summarise(mean = mean(Length, na.rm = TRUE),median = median(Length, na.rm = TRUE))

res <- wilcox.test(Multiplicity ~ Hit, data = Combined_BLAST)
write.table(data.frame(statistic=res$statistic, p_value=res$p.value, method=res$method, alternative=res$alternative), file=paste0("data_clean/BLAST_Stats_Wilcoxon_HitvDepth.tsv"), quote=FALSE, sep='\t', row.names=F)
Combined_BLAST %>% group_by(Hit) %>% summarise(mean = mean(Multiplicity, na.rm = TRUE),median = median(Multiplicity, na.rm = TRUE))

# Work around for some poorly formatted Subject_Organism fields 
Combined_BLAST_Taxon <- Combined_BLAST %>% filter(!is.na(Subject_Taxon_ID)) %>% arrange(Individual_ID) %>% group_by(Individual_ID) %>% mutate(Subject_Organism = first(Subject_Organism)) %>% ungroup() %>%
  select(Subject_Taxon_ID, Subject_Organism, Individual_ID) %>% distinct(Individual_ID,Subject_Taxon_ID, .keep_all = TRUE)

write.table(Combined_BLAST_Taxon, file=paste0("data_clean/BLAST_Stats_Taxon.tsv"), quote=FALSE, sep='\t', row.names=F)

Combined_BLAST_Taxon_List <- Combined_BLAST_Taxon %>% filter(!is.na(Subject_Taxon_ID)) %>% pull(Subject_Taxon_ID) %>% unique() 

write.table(Combined_BLAST_Taxon_List, file=paste0("data_clean/BLAST_Stats_Taxon_list.tsv"), quote=FALSE, sep='\t', row.names=F, col.names=F)

for (taxon in Combined_BLAST_Taxon_List) {
  system2("Rscript", args = c(paste0(script_dir,"/BuildTaxonomicLineage.R"), taxon, paste0(topdir,"/NCBI_Taxonomy/")))
}
#need cleanup step

for(taxon in Combined_BLAST_Taxon_List) {
	file <- paste0(topdir,"/NCBI_Taxonomy/Lineages/",taxon,"_full.tsv")
	if(file.exists(file)) {
	    if(file.info(file)$size > 0){
			data <- fread(file, sep = '\t') 
	        assign(paste0("taxon_",taxon),data)
	    } else {
	        print(paste0(file," is empty!!!"))
	    } 
	} else {
	    print(paste0(file," is missing!!!"))
	}
}

Combined_BLAST_Taxon_Domain <- data.frame()
for (taxon in Combined_BLAST_Taxon_List) {
	file <- paste0(topdir,"/NCBI_Taxonomy/Lineages/",taxon,"_full.tsv")
	if(file.exists(file)) {
	    if(length(readLines(file)) > 2){
	    	# Viruses aren't a domain, but making them fit into the df structure under domain
	    	data <- fread(file, sep = '\t') %>% filter(rank %in% c("acellular root","kingdom","domain","phylum","genus","class")) %>% mutate(rank = ifelse(rank=="acellular root","domain",rank)) %>%
	    		select(rank,name) %>% pivot_wider(names_from = rank, values_from = name) 
	    		
			for (col in c("domain","kingdom","phylum","class","genus")) {
  				if (!col %in% names(data)) data[[col]] <- NA
			}
	    	data$taxid <- taxon
	    		
			Combined_BLAST_Taxon_Domain <- rbind(Combined_BLAST_Taxon_Domain, data %>% select(domain, kingdom, phylum, class, genus, taxid) )
	    } else {
	        print(paste0(file," is empty!!!"))
	    } 
	} else {
	    print(paste0(file," is missing!!!"))
	}
}
write.table(Combined_BLAST_Taxon_Domain, file=paste0("data_clean/BLAST_Stats_Taxon_Domain.tsv"), quote=FALSE, sep='\t', row.names=F, col.names=T)

####################################################################################################################
###                                               Joint Data	                                                ####
####################################################################################################################

Combined_Summary <- inner_join(Combined_Kraken_Summary %>% mutate(Individual_ID = as.character(Individual_ID)), Combined_Megahit_Summary, by = "Individual_ID")
Combined_Summary <- Combined_Summary %>% mutate(Estimate_contamination_percentage = 100 - Eukaryota_percentage)

cor_res <- cor.test(~ Multiplicity_n + Eukaryota_percentage, data = Combined_Summary, method = "pearson")

write.table(Combined_Summary, file=paste0("data_clean/Combined_Stats_Summary.tsv"), quote=FALSE, sep='\t', row.names=F)

####################################################################################################################
###                                               Blood Comp	                                                ####
####################################################################################################################

for(number in fqd_main$Order) {
	line <- fqd_main %>% filter(Order == number)
	filetag <- gsub("_R1.*","",line$Filename1)
	file <- paste0(dirname(topdir),"/Alignments/",line$Individual_ID,"_NoFastP_Kraken_trained_taxid_k19/Kraken2/",filetag,"_subsampled_Kraken2Screened_Report.tsv")
	if(file.exists(file)) {
	    if(file.info(file)$size > 0){
	        data <- fread(file, sep = '\t') %>% filter(V6 %in% c("unclassified","root","Eukaryota","Bacteria","Archaea","Viruses","other entries") ) %>% 
	        	select(V2,V6) %>% pivot_wider( names_from = V6, values_from = V2 ) %>% mutate(FastQ = paste0(line$Individual,"_",line$DNA_Number,"_",line$Library_Number,"_",line$Run_Number)) %>%
	        	mutate(Individual_ID = line$Individual)
	        assign(paste0("kraken_matched1_",number),data)
	    } else {
	        print(paste0(file," is empty!!!"))
	    } 
	} else {
	    print(paste0(file," is missing!!!"))
	}
	file <- paste0(dirname(topdir),"/Alignments/",line$Individual_ID,"_NoFastP_Krakentaxid_k19/Kraken2/",filetag,"_subsampled_Kraken2Screened_Report.tsv")
	if(file.exists(file)) {
	    if(file.info(file)$size > 0){
	        data <- fread(file, sep = '\t') %>% filter(V6 %in% c("unclassified","root","Eukaryota","Bacteria","Archaea","Viruses","other entries") ) %>% 
	        	select(V2,V6) %>% pivot_wider( names_from = V6, values_from = V2 ) %>% mutate(FastQ = paste0(line$Individual,"_",line$DNA_Number,"_",line$Library_Number,"_",line$Run_Number)) %>%
	        	mutate(Individual_ID = line$Individual)
	        assign(paste0("kraken_matched2_",number),data)
	    } else {
	        print(paste0(file," is empty!!!"))
	    } 
	} else {
	    print(paste0(file," is missing!!!"))
	}
}
ls(pattern=paste0("^kraken_matched1_*")) -> listOfDataFrames
bind_rows(mget(listOfDataFrames)) %>% mutate(Group="Trained_Database") -> Combined_Kraken_Matched1
ls(pattern=paste0("^kraken_matched2_*")) -> listOfDataFrames
bind_rows(mget(listOfDataFrames)) %>% mutate(Group="Standard_Database") -> Combined_Kraken_Matched2

Combined_Kraken_Matched <- rbind(Combined_Kraken_Matched1,Combined_Kraken_Matched2) %>% mutate(SampleType=ifelse(substr(Individual_ID, nchar(Individual_ID), nchar(Individual_ID)) == "S","Saliva",ifelse(substr(Individual_ID, nchar(Individual_ID), nchar(Individual_ID)) == "B","Blood","ERROR")))
write.table(Combined_Kraken_Matched, file=paste0("data_clean/Kraken_Stats_Matched.tsv"), quote=FALSE, sep='\t', row.names=F)

Combined_Kraken_Matched_Summary <- Combined_Kraken_Matched %>% filter(SampleType == "Blood") %>% group_by(Individual_ID,Group) %>% summarise( across( where(is.numeric), ~ sum(.x, na.rm = TRUE),.names = "{.col}" ),.groups = "drop") %>%
	mutate(total = unclassified + root) %>% mutate(unclassified_percentage = (unclassified / total) * 100, Eukaryota_percentage = (Eukaryota / total) * 100, Bacteria_Percentage = (Bacteria / total) * 100, Viruses_Percentage = (Viruses / total) * 100 )
write.table(Combined_Kraken_Matched_Summary, file=paste0("data_clean/Kraken_Stats_Matched_Summary.tsv"), quote=FALSE, sep='\t', row.names=F)

Combined_Kraken_Matched_Summary_Cohort <- Combined_Kraken_Matched %>% filter(SampleType == "Blood") %>% mutate(total = unclassified + root, unclassified_percentage = (unclassified / total) * 100, Eukaryota_percentage = (Eukaryota / total) * 100, Bacteria_Percentage = (Bacteria / total) * 100, Viruses_Percentage = (Viruses / total) * 100) %>% group_by(Group) %>% 
	summarise(across(where(is.numeric), list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)), .names = "{.col}_{.fn}"), .groups = "drop")
write.table(Combined_Kraken_Matched_Summary_Cohort, file = "data_clean/Kraken_Stats_Matched_Summary_Cohort.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

results <- data.frame()
for (col in c("unclassified_percentage", "Eukaryota_percentage", "Bacteria_Percentage","Viruses_Percentage")) {
	g1 <- Combined_Kraken_Matched_Summary[[col]][Combined_Kraken_Matched_Summary$Group == unique(Combined_Kraken_Matched_Summary$Group)[1]]
	g2 <- Combined_Kraken_Matched_Summary[[col]][Combined_Kraken_Matched_Summary$Group == unique(Combined_Kraken_Matched_Summary$Group)[2]]
  	t <- t.test(g1, g2, paired = TRUE)
	res <- data.frame(Variable = col,
    	Mean_Group1 = mean(Combined_Kraken_Matched_Summary[[col]][Combined_Kraken_Matched_Summary$Group == unique(Combined_Kraken_Matched_Summary$Group)[1]], na.rm = TRUE),
    	Mean_Group2 = mean(Combined_Kraken_Matched_Summary[[col]][Combined_Kraken_Matched_Summary$Group == unique(Combined_Kraken_Matched_Summary$Group)[2]], na.rm = TRUE),
    	t_statistic = unname(t$statistic),
    	p_value = t$p.value,
    	conf_low = t$conf.int[1],
    	conf_high = t$conf.int[2],
    	df = unname(t$parameter)
    )
    results <- rbind(results, res)
}
write.table(results, file=paste0("data_clean/Kraken_Stats_Matched_T-Tests.tsv"), quote=FALSE, sep='\t', row.names=F)

####################################################################################################################
###                                               Depth Stats	                                                ####
####################################################################################################################

for(id in unique(fqd$Individual_ID)) {
	file <- paste0(T_SEQ_DIR,"/hg38/Alignments/Genomic/",id,"/PicardStats/collect_wgs_metrics_",id,".txt")
	if(file.exists(file)) {
	    if(file.info(file)$size > 0){
			data <- fread(file, sep = '\t', skip = 6, nrows = 2) 
	        assign(paste0("picard_",id),data)
	    } else {
	        print(paste0(file," is empty!!!"))
	    } 
	} else {
	    print(paste0(file," is missing!!!"))
	}
}
ls(pattern=paste0("^picard_*")) -> listOfDataFrames
bind_rows(mget(listOfDataFrames)) -> Combined_Picard

mean(Combined_Picard$MEAN_COVERAGE)
sd(Combined_Picard$MEAN_COVERAGE)

####################################################################################################################
###                                               	 Fin		                                                ####
####################################################################################################################

quit()
