####################################################################################################################
###                                               Start_Up                                                      ####
####################################################################################################################

rm(list=ls())

args <- commandArgs()
script_path <- sub("^--file=", "", args[grep("^--file=", args)])

root <- args[6]
tag <- args[7]
db1 <- args[8]
db2 <- args[9]

if (length(script_path) == 0) {
  stop("Unable to determine script path")
}

script_path <- normalizePath(script_path, mustWork = TRUE)
script_dir <- dirname(script_path)
parent_dir <- dirname(script_dir)

source(file.path(parent_dir, "parameters.R"))

###   Directory creation  ### 
dir <- PROJECT_DIR
aligndir <- paste0(PROJECT_DIR,"/Alignments/")
dir <- paste0(PROJECT_DIR,"/Analyses/Kmer_Comparison_Analysis/data_raw/",root,"/",tag)
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

if (!require(foreach)) install.packages("foreach")
library(foreach)

if (!require(BiocManager)) install.packages("BiocManager")
library(BiocManager)

if (!require(doParallel)) BiocManager::install("doParallel")
library(doParallel)

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
###                                               Set up Parallel                                               ####
####################################################################################################################

num_cores <- 2 #detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)
registerDoParallel(cl)

####################################################################################################################
###                                               Data Read In                                                  ####
####################################################################################################################

db1_name <- gsub("_dataset_kmc","",db1)
db2_name <- gsub("_dataset_kmc","",db2)

files <- c(paste0("data_raw/unique_to_",db1_name,"_dataset_kmc"),paste0("data_raw/unique_to_",db2_name,"_dataset_kmc"),paste0(aligndir,db1_name,"/KMC/",db1_name,"_dataset_kmc"),paste0(aligndir,db2_name,"/KMC/",db2_name,"_dataset_kmc"),"data_raw/common_kmers")

for( ID in files ) {
	basename <- basename(ID)
	if (file.exists(paste0("data_clean/",basename,"_k-mer_Counts.csv")) && file.info(paste0("data_clean/",basename,"_k-mer_Counts.csv"))$size > 0) {
		print(paste0(basename,"_k-mer_Counts.csv already exists, continuing..."))
	} else {
		# Set Variables
		file_path <- paste0(ID,".txt")
		file_size <- fread(paste0(ID,"_size.txt")) %>% pull(V1) %>% as.numeric()
		chunk_size <- 1e9
		
		# Set Functions
		read_chunk <- function(file, start_row, num_rows) {
		  fread(file, skip = start_row - 1, select = 2, nrows = num_rows, header = FALSE, col.names = c("count"))
		}
		summarise_chunk <- function(chunk) {
		  chunk %>%
		    group_by(count) %>%
		    summarise(Group_Count = n(), .groups = 'drop') %>%
		    as.data.table()
		}
		combine_summaries <- function(summary1, summary2) {
		  rbindlist(list(summary1, summary2))[, .(Group_Count = sum(Group_Count)), by = count]
		}
		
		# Read in and run functions
		num_chunks <- ceiling( file_size / chunk_size )
		summaries <- foreach(i = 1:num_chunks, .combine = combine_summaries, .packages = c('data.table', 'tidyverse')) %dopar% {
		  start_row <- (i - 1) * chunk_size + 1
		  chunk <- read_chunk(file_path, start_row, chunk_size)
		  summarise_chunk(chunk)
		}
		fwrite(summaries, paste0("data_clean/",basename,"_k-mer_Counts.csv"))
		
		plot <- ggplot(summaries, aes(x = count, y = Group_Count)) +
		  geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.8) +
		  labs(x = "k-mer Count", y = "k-mer Count Frequency (Millions)") +
		  scale_y_continuous(labels = scales::comma_format(scale = 1/1e6)) +
		  scale_x_continuous(expand = c(0, 0),limits = c(0, 150),breaks=seq(0, 140, by = 10)) + 
		  ggplot_theme_basic
		ggsave(plot=plot,filename=paste0("graphs/",basename,"_k-mer_histogram.pdf"),device="pdf", width = 40, height = 30, units = "cm",dpi = 640)
	}
}

stopCluster(cl)

####################################################################################################################
###                                              	 fin			                                            ####
####################################################################################################################

quit()
