####################################################################################################################
###                                               Start_Up                                                      ####
####################################################################################################################

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
TOP_DIR <- PROJECT_DIR
WRK_DIR <- file.path(TOP_DIR, "Analyses/Bracken_Analysis")

TRA_META <- file.path(TOP_DIR, "Training/Meta_Data/Training_Sequenced_Samples_Info.tsv")
TRA_ALG_DIR <- file.path(TOP_DIR, "Training/Kraken2")

VAL_META <- file.path(TOP_DIR, "Validation/Meta_Data/Validation_Sequenced_Samples_Info.tsv")
VAL_ALG_DIR <- file.path(TOP_DIR, "Validation/Alignments/CMM_Alignments")

COM_META <- file.path(TOP_DIR, "MetaData/Comparison_Sequenced_Samples_Cleanup.tsv")
COM_ALG_DIR <- file.path(TOP_DIR, "Alignments")

Dirs <- c("graphs","data_clean","data_raw")

if (!dir.exists(WRK_DIR)) dir.create(WRK_DIR, recursive = TRUE)
setwd(WRK_DIR)
for (subDir in Dirs) if (!dir.exists(file.path(WRK_DIR, subDir))) dir.create(file.path(WRK_DIR, subDir))

message("[INFO] Directory: ", WRK_DIR)

####################################################################################################################
###                                              Packages                                                       ####
####################################################################################################################

if(!require(tidyverse)) install.packages("tidyverse")
library(tidyverse)
if(!require(data.table)) install.packages("data.table")
library(data.table)
if(!require(vegan)) install.packages("vegan")
library(vegan)
if(!require(patchwork)) install.packages("patchwork")
library(patchwork)

####################################################################################################################
###                                           ggplot theme                                                      ####
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
###                                           Set Genus Groups 	                                                ####
####################################################################################################################

# HOMD (oral) and Salter et al. 2014 (reagent)
# Accessed 10/10/26
oral    <- c("Abiotrophia", "Acinetobacter", "Actinomyces", "Aggregatibacter", "Alloprevotella", "Alloscardovia", "Anaeroglobus", "Anaerosphaera", "Arachnia", "Atopobium", "Bacteroides", "Bifidobacterium", "Brooklawnia", "Bulleidia", "Butyrivibrio", "Campylobacter", "Capnocytophaga", "Cardiobacterium", "Catonella", "Centipeda", "Cloacibacterium", "Corynebacterium", "Cryptobacterium", "Desulfobulbus", "Desulfomicrobium", "Desulfovibrio", "Dialister", "Eggerthia", "Eikenella", "Fannyhessea", "Fastidiosipila", "Filifactor", "Fretibacterium", "Fusobacterium", "Gallibacter", "Gemella", "Granulicatella", "Haemophilus", "Hallella", "Helicobacter", "Hornefia", "Hoylesella", "Johnsonella", "Jonquetella", "Kingella", "Lachnoanaerobaculum", "Lactobacillus", "Lancefieldella", "Lautropia", "Leptotrichia", "Ligilactobacillus", "Limosilactobacillus", "Megasphaera", "Metamycoplasma", "Methanobrevibacter", "Mitsuokella", "Mogibacterium", "Moraxella", "Moryella", "Mycoplasmopsis", "Nanogingivalis", "Nanoperiomorbus", "Nanosynbacter", "Nanosyncoccus", "Neisseria", "Odoribacter", "Olsenella", "Oribacterium", "Ottowia", "Parascardovia", "Parvimonas", "Parvisynbacter", "Peptidiphaga", "Peptoanaerobacter", "Peptococcus", "Peptoniphilus", "Peptostreptococcus", "Phocaeicola", "Porphyromonas", "Prevotella", "Propionibacterium", "Pseudoleptotrichia", "Pseudomonas", "Pseudoramibacter", "Pyramidobacter", "Ralstonia", "Riemerella", "Rothia", "Saccharimonas", "Scardovia", "Schaalia", "Segatella", "Selenomonas", "Serratia", "Shuttleworthia", "Simonsiella", "Slackia", "Solobacterium", "Stomatobaculum", "Streptococcus", "Tannerella", "Treponema", "Veillonella")
reagent <- c("Afipia", "Aquabacterium", "Asticcacaulis", "Aurantimonas", "Beijerinckia", "Bosea", "Bradyrhizobium", "Brevundimonas", "Caulobacter", "Craurococcus", "Devosia", "Hoeflea", "Mesorhizobium", "Methylobacterium", "Novosphingobium", "Ochrobactrum", "Paracoccus", "Pedomicrobium", "Phyllobacterium", "Rhizobium", "Roseomonas", "Sphingobium", "Sphingomonas", "Sphingopyxis", "Acidovorax", "Azoarcus", "Azospira", "Burkholderia", "Comamonas", "Cupriavidus", "Curvibacter", "Delftia", "Duganella", "Herbaspirillum", "Janthinobacterium", "Kingella", "Leptothrix", "Limnobacter", "Massilia", "Methylophilus", "Methyloversatilis", "Oxalobacter", "Pelomonas", "Polaromonas", "Ralstonia", "Schlegelella", "Sulfuritalea", "Undibacterium", "Variovorax", "Acinetobacter", "Enhydrobacter", "Enterobacter", "Escherichia", "Nevskia", "Pseudomonas", "Pseudoxanthomonas", "Psychrobacter", "Stenotrophomonas", "Xanthomonas", "Aeromicrobium", "Arthrobacter", "Beutenbergia", "Brevibacterium", "Corynebacterium", "Curtobacterium", "Dietzia", "Geodermatophilus", "Janibacter", "Kocuria", "Microbacterium", "Micrococcus", "Microlunatus", "Patulibacter", "Propionibacterium", "Rhodococcus", "Tsukamurella", "Abiotrophia", "Bacillus", "Brevibacillus", "Brochothrix", "Facklamia", "Paenibacillus", "Streptococcus", "Chryseobacterium", "Dyadobacter", "Flavobacterium", "Hydrotalea", "Niastella", "Olivibacter", "Pedobacter", "Wautersiella", "Deinococcus")
# Note 8 overlapping genera removed from reagent (retained in oral) 
reagent <- setdiff(reagent, oral)

src <- function(x) case_when(x %in% oral ~ "Oral commensal", x %in% reagent ~ "Reagent/environmental", TRUE ~ "Other/unassigned")

####################################################################################################################
###                                           Read Bracken outputs                                              ####
####################################################################################################################

safe <- function(x) tryCatch(print(x), error = function(e) message("Skipped: ", conditionMessage(e)))
cohorts <- tribble(~Cohort, ~Meta, ~Dir, "Training", TRA_META, TRA_ALG_DIR, "Validation", VAL_META, VAL_ALG_DIR, "Comparison", COM_META, COM_ALG_DIR)

read_bracken <- function(Cohort, Dir, ...) {
	f <- list.files(Dir, pattern = "_Report_(Species|Genus)\\.bracken$", recursive = TRUE, full.names = TRUE)
	if (!length(f)) { warning("No Bracken files for ", Cohort, " under ", Dir); return(NULL) }
	rel <- str_remove(f, fixed(paste0(sub("/+$", "", Dir), "/")))
	map2_dfr(f, rel, ~ read_tsv(.x, show_col_types = FALSE) %>% mutate(Sample = str_split_i(.y, "/", 1), Level = if_else(str_detect(.x, "_Species\\.bracken$"), "Species", "Genus"), Cohort = Cohort, File = .x))
}
bracken <- pmap_dfr(cohorts, read_bracken) %>% select(Sample, Cohort, Level, everything()) %>% mutate(Sample = sub("_.*","",Sample)) %>% mutate(Sample_Type = if_else(str_ends(Sample, "B"), "Blood", "Saliva")) %>%
	mutate(DataBaseType = if_else(str_ends(Cohort, "Validation"), "Standard_Trained", "Standard"))

bracken %>% distinct(Cohort, Sample, Level) %>% count(Cohort, Level)
multi <- bracken %>% distinct(Sample, Level, File) %>% count(Sample, Level) %>% filter(n > 1)
if (nrow(multi)) warning(nrow(multi), " sample/level combination(s) have more than one Bracken file; rows would be double counted: ", paste(head(unique(multi$Sample), 5), collapse = ", "))

bracken <- select(bracken, -File)
write_tsv(bracken, "data_raw/bracken_tidy.tsv")

pair_key <- function(x) str_remove(x, "[-_.]?(B|S)$") 
samp <- bracken %>% distinct(Sample, Cohort,Sample_Type) %>% add_count(Cohort, Sample_Type, name = "N_Grp")
sinfo <- distinct(samp, Sample, Cohort, Sample_Type) %>% group_by(Cohort) %>% arrange(Sample, .by_group = TRUE) %>%
	mutate(Sample_Label = case_when( Cohort == "Training"   ~ paste0("Train_", row_number()), Cohort == "Validation" ~ paste0("Val_", row_number()), Cohort == "Comparison" ~ paste0("Comp_", dense_rank(pair_key(Sample)), if_else(Sample_Type == "Blood", "B", "S")))) %>%
	ungroup() %>% mutate(Sample_Label = factor(Sample_Label, levels = str_sort(unique(Sample_Label), numeric = TRUE)))

genus <- bracken %>% filter(taxonomy_lvl == "G", name != "Homo") %>%
	transmute(Sample, Genus = name, Reads = new_est_reads, Frac_Total = fraction_total_reads) %>%
	left_join(sinfo, by = "Sample") %>% group_by(Sample) %>% mutate(Rel = Reads / sum(Reads), Source = src(Genus)) %>% ungroup()

####################################################################################################################
###                                        1. Source Composition                                                ####
####################################################################################################################

n_top <- 3; n_s <- n_distinct(genus$Sample)
md  <- requireNamespace("ggtext", quietly = TRUE)
hdr <- function(x) if (md) paste0("**", x, "**") else x

cohort_lv <- c("Comparison", "Training", "Validation"); type_lv <- c("Blood","Saliva")
cat_cols  <- c("Oral commensal" = "#1B9E77", "Reagent/environmental" = "#D95F02", "Other/unassigned" = "#666666")
cat_label <- c("Oral commensal" = "Oral", "Reagent/environmental" = "Reagent/Environmental", "Other/unassigned" = "Other")

cat_order <- intersect(names(cat_cols), unique(genus$Source))

sp_hdr <- TRUE; sp_grp <- TRUE        
plot_gap <- 0.35; leg_txt <- 28; bar_alpha <- 0.9; stack_reverse <- FALSE

top_g <- genus %>% group_by(Source, Genus) %>% summarise(M = sum(Rel) / n_s, .groups = "drop") %>% group_by(Source) %>% slice_max(M, n = n_top, with_ties = FALSE) %>% ungroup()

key <- imap_dfr(cat_order, function(cs, i) {
	g  <- top_g %>% filter(Source == cs) %>% arrange(desc(M)) %>% pull(Genus)
	sh <- colorRampPalette(c(cat_cols[[cs]], "white"))(length(g) + 3)[seq_len(length(g) + 1)]
	tibble(Key    = c(if (sp_grp && i > 1) paste0("SPG_", cs), paste0("HDR_", cs), if (sp_hdr) paste0("SPH_", cs), g, paste0(cs, "_other")),
	       Label  = c(if (sp_grp && i > 1) "",                 hdr(cat_label[[cs]]), if (sp_hdr) "",              g, "Other Genera"),
	       Colour = c(if (sp_grp && i > 1) "transparent",      "transparent",        if (sp_hdr) "transparent",   sh))
})

comp <- genus %>% left_join(top_g %>% select(Source, Genus) %>% mutate(Top = TRUE), by = c("Source", "Genus")) %>% mutate(Key = factor(if_else(!is.na(Top), Genus, paste0(Source, "_other")), levels = key$Key)) %>%
	group_by(Sample, Sample_Label, Cohort, Sample_Type, Key) %>% summarise(Frac = sum(Rel), .groups = "drop")

human <- bracken %>% filter(taxonomy_lvl == "G", name == "Homo") %>% group_by(Sample) %>% summarise(Human = sum(fraction_total_reads) * 100)
ord <- sinfo %>% left_join(human, by = "Sample") %>% mutate(Human = replace_na(Human, 0), Microbial = 100 - Human) %>% arrange(factor(Cohort, cohort_lv), factor(Sample_Type, type_lv), Human, Sample_Label)
lev <- as.character(ord$Sample_Label)
fix <- function(d) d %>% mutate(Cohort = factor(Cohort, cohort_lv), Sample_Type = factor(Sample_Type, type_lv), Sample_Label = factor(as.character(Sample_Label), levels = lev))
ord %>% select(Sample, Cohort, Sample_Type, Human) %>% write_tsv("data_clean/human_fraction.tsv")   

txt_theme <- theme(
  axis.text.y = element_text(size = 24, margin = margin(r = 10), hjust = 0),
  axis.title.y = element_text(size = 28, margin = margin(r = 20)),
  strip.text = element_text(size = 28, margin = margin(t = 10, b = 24, l = 8, r = 8)),
  strip.background = element_blank(),
  legend.title = element_text(size = 30, hjust = 0),
  legend.text = element_text(hjust = 0),
  legend.key = element_blank(),
  legend.background = element_rect(fill = "white"),
  legend.box.just = "left",
  legend.justification = "left",
  legend.title.align = 0,
  legend.text.align = 0,
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)

p_hum <- fix(ord) %>% select(Sample_Label, Cohort, Sample_Type, Human, Microbial) %>%
	pivot_longer(c(Human, Microbial), names_to = "Read_Type", values_to = "Pct") %>% mutate(Read_Type = factor(Read_Type, c("Microbial", "Human"))) %>%
	ggplot(aes(Sample_Label, Pct, fill = Read_Type, colour = Read_Type)) + geom_col(width = 1, alpha = bar_alpha) +
		facet_grid(~ Cohort + Sample_Type, scales = "free_x", space = "free_x") +
		scale_fill_manual(values = c(Microbial = "#1B9E77", Human = "#7570B3")) + 
		scale_colour_manual(values = c(Microbial = "#1B9E77", Human = "#7570B3"), guide = "none") +
		scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
		labs(x = NULL, y = "Percentage of\nClassified Reads", fill = NULL) + theme_classic() + txt_theme + theme(legend.text = element_text(size = leg_txt))

p_comp <- ggplot(fix(comp), aes(Sample_Label, Frac * 100, fill = Key, colour = Key)) +
	geom_col(width = 1, alpha = bar_alpha, position = position_stack(reverse = stack_reverse)) +
	facet_grid(~ Cohort + Sample_Type, scales = "free_x", space = "free_x") +
	scale_fill_manual(values = setNames(key$Colour, key$Key), limits = key$Key, breaks = key$Key, labels = key$Label) +
	scale_colour_manual(values = setNames(key$Colour, key$Key), guide = "none") +
	scale_x_discrete(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
	labs(x = NULL, y = "Percentage of\nMicrobial Reads", fill = NULL) + theme_classic() + txt_theme +
	theme(strip.text = element_blank(), legend.key.size = unit(0.35, "cm"),
	      legend.text = if (md) ggtext::element_markdown(size = leg_txt) else element_text(size = leg_txt))

p <- p_hum / plot_spacer() / p_comp + plot_layout(heights = c(1, plot_gap, 4))
ggsave("graphs/Source_Composition.pdf", p, width = 32, height = 18)

source_frac <- genus %>% group_by(Sample, Sample_Label, Cohort, Sample_Type, Source) %>% summarise(Frac = sum(Rel), .groups = "drop") %>%
	complete(nesting(Sample, Sample_Label, Cohort, Sample_Type), Source, fill = list(Frac = 0))

p <- ggplot(source_frac, aes(Sample_Label, Frac * 100, fill = Source)) + 
    scale_color_brewer(palette = "Dark2") +
  	scale_fill_brewer(palette = "Dark2") +
	geom_col() + 
	facet_grid(~ Cohort + Sample_Type, scales = "free_x", space = "free_x") +
	labs(x = NULL, y = "% of microbial reads") + 
	theme_classic() + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("graphs/Source_Composition_Simple.pdf", p, width = 16, height = 7)

####################################################################################################################
###                      3. Shared vs sample-specific genera, within and across groups                          ####
####################################################################################################################

prev <- genus %>% group_by(Cohort, Sample_Type, Genus, Source) %>% summarise(N_Samples = n_distinct(Sample), Sum_Rel = sum(Rel), .groups = "drop") %>%
	left_join(distinct(samp, Cohort, Sample_Type, N_Grp), by = c("Cohort", "Sample_Type")) %>%
	mutate(Prevalence = N_Samples / N_Grp, Mean_Rel = Sum_Rel / N_Grp, Class = case_when(N_Samples == 1 ~ "Sample-specific", Prevalence >= 0.9 ~ "Core (>=90%)", Prevalence >= 0.5 ~ "Common (50-90%)", TRUE ~ "Variable"))
write_tsv(prev, "data_clean/genus_prevalence.tsv")
count(prev, Cohort, Sample_Type, Class)

top30 <- prev %>% group_by(Cohort, Sample_Type) %>% slice_max(Mean_Rel, n = 30, with_ties = FALSE) %>%
	summarise(N_Oral_Of_Top30 = sum(Source == "Oral commensal"), Non_Oral = paste(Genus[Source != "Oral commensal"], collapse = ", "), .groups = "drop")
write_tsv(top30, "data_clean/top30_oral_check.tsv"); top30

shared <- prev %>% filter(Sample_Type == "Saliva", Prevalence >= 0.5) %>% group_by(Genus, Source) %>% summarise(N_Cohorts = n_distinct(Cohort), Cohorts = paste(sort(unique(Cohort)), collapse = ","), .groups = "drop")
write_tsv(shared, "data_clean/genus_shared_across_cohorts.tsv"); count(shared, N_Cohorts, Source)
mean_rel <- prev %>% filter(Sample_Type == "Saliva") %>% select(Cohort, Genus, Mean_Rel) %>% pivot_wider(names_from = Cohort, values_from = Mean_Rel, values_fill = 0) %>% column_to_rownames("Genus")
round(cor(mean_rel, method = "spearman"), 2)     # between-cohort agreement in genus abundance (saliva only)

top_g <- prev %>% group_by(Genus) %>% summarise(M = mean(Mean_Rel)) %>% slice_max(M, n = 30, with_ties = FALSE) %>% arrange(desc(M)) %>% pull(Genus)
p <- prev %>% filter(Genus %in% top_g) %>% mutate(Genus = factor(Genus, rev(top_g)), Group = paste(Cohort, Sample_Type)) %>%
	ggplot(aes(Group, Genus, fill = Prevalence)) + geom_tile() + facet_grid(Source ~ ., scales = "free_y", space = "free_y") + scale_fill_viridis_c(limits = c(0, 1), na.value = "#440154") +
	labs(x = NULL, y = NULL) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("graphs/Genus_Prevalence_Top30.pdf", p, width = 10, height = 9)

####################################################################################################################
###                        					4. Alpha diversity         									        ####
####################################################################################################################

alpha <- genus %>% group_by(Sample, Cohort, Sample_Type) %>% summarise(Microbial_Reads = sum(Reads), Genera = n(), Shannon = -sum(Rel * log(Rel)), .groups = "drop")
write_tsv(alpha, "data_clean/alpha_diversity.tsv")
alpha %>% group_by(Cohort, Sample_Type) %>% summarise(N = n(), across(c(Genera, Shannon, Microbial_Reads), median), .groups = "drop")
safe(kruskal.test(Shannon ~ Cohort, data = filter(alpha, Sample_Type == "Saliva")))
safe(wilcox.test(Shannon ~ Sample_Type, data = alpha))
p <- ggplot(alpha, aes(Cohort, Shannon, fill = Sample_Type)) + geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_point(position = position_jitterdodge(jitter.width = 0.15), size = 1.5) + theme_classic()
ggsave("graphs/Alpha_Shannon_ByGroup.pdf", p, width = 8, height = 6)
p <- ggplot(alpha, aes(Microbial_Reads, Genera, colour = Cohort, shape = Sample_Type)) + geom_point(size = 2) + scale_x_log10() + theme_classic()   # does depth, not biology, drive genus counts?
ggsave("graphs/Genera_vs_Depth.pdf", p, width = 8, height = 6)

####################################################################################################################
###                        5. Beta diversity: similarity within and between groups						        ####
####################################################################################################################

mat <- genus %>% select(Sample, Genus, Rel) %>% pivot_wider(names_from = Genus, values_from = Rel, values_fill = 0) %>% column_to_rownames("Sample") %>% as.matrix()
bm <- tibble(Sample = rownames(mat)) %>% left_join(sinfo, by = "Sample")
bc <- vegdist(mat, "bray")
grp <- setNames(paste(bm$Cohort, bm$Sample_Type), bm$Sample)
sim <- as.data.frame(as.table(as.matrix(1 - bc))) %>% mutate(A = as.character(Var1), B = as.character(Var2)) %>% filter(A < B) %>%
	mutate(G1 = pmin(unname(grp[A]), unname(grp[B])), G2 = pmax(unname(grp[A]), unname(grp[B])))
pair_sim <- sim %>% group_by(G1, G2) %>% summarise(N_Pairs = n(), Median_Similarity = round(median(Freq), 3), Min = round(min(Freq), 3), .groups = "drop")
write_tsv(pair_sim, "data_clean/pairwise_similarity_by_group.tsv"); pair_sim    # G1 == G2 rows = within-group

pc <- cmdscale(bc, eig = TRUE, k = 2); ve <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
p <- bind_cols(bm, PC1 = pc$points[, 1], PC2 = pc$points[, 2]) %>% ggplot(aes(PC1, PC2, colour = Cohort, shape = Sample_Type)) + geom_point(size = 3) +
	labs(x = paste0("PCoA1 (", ve[1], "%)"), y = paste0("PCoA2 (", ve[2], "%)")) + theme_classic()
ggsave("graphs/Genus_PCoA_BrayCurtis.pdf", p, width = 8, height = 6)
k <- !is.na(bm$Sample_Type); bck <- as.dist(as.matrix(bc)[k, k]); bmk <- bm[k, ]
safe(adonis2(bck ~ Cohort + Sample_Type, data = bmk, permutations = 9999, by = "margin"))
safe(anova(betadisper(bck, paste(bmk$Cohort, bmk$Sample_Type)), permutations = 9999))

####################################################################################################################
###                    6. Biomass relationship: does the reagent fraction rise as microbial load falls?	        ####
####################################################################################################################
 
bio <- source_frac %>% left_join(genus %>% group_by(Sample) %>% summarise(Microbial_Frac = sum(Frac_Total)), by = "Sample")
bio %>% filter(Sample_Type == "Saliva") %>% { bind_rows(., mutate(., Cohort = "All")) } %>% group_by(Cohort, Source) %>% filter(n() >= 4) %>%
	summarise(N = n(), rho = cor(Microbial_Frac, Frac, method = "spearman"), P = cor.test(Microbial_Frac, Frac, method = "spearman", exact = FALSE)$p.value, .groups = "drop") %>% write_tsv("data_clean/source_vs_biomass.tsv")
p <- ggplot(bio, aes(Microbial_Frac, Frac * 100, colour = Cohort, shape = Sample_Type)) + geom_point() +
	geom_smooth(data = filter(bio, Sample_Type == "Saliva"), aes(group = Cohort), method = "lm", se = FALSE) + scale_x_log10() + facet_wrap(~ Source, scales = "free_y") +
	labs(x = "Microbial fraction of classified reads (log10)", y = "% of microbial reads") + theme_classic()
ggsave("graphs/Source_vs_Biomass.pdf", p, width = 12, height = 7)
p <- bio %>% distinct(Sample, Cohort, Sample_Type, Microbial_Frac) %>% ggplot(aes(Cohort, Microbial_Frac, fill = Sample_Type)) + geom_boxplot(alpha = 0.6) + scale_y_log10() + theme_classic()
ggsave("graphs/MicrobialFraction_ByGroup.pdf", p, width = 8, height = 6)

####################################################################################################################
###                                               Finish                                                       ####
####################################################################################################################

quit()
