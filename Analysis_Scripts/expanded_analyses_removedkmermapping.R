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
TOP_DIR <- PROJECT_DIR
WRK_DIR <- file.path(TOP_DIR, "Analyses/Removed_Read_Analysis")

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
###                                           Parameters and Helpers                                            ####
####################################################################################################################

set.seed(1)
conf     <- list(mapq = 20, aligned = 0.9)    # "high confidence" read
chrs     <- c(paste0("chr", 1:22), "chrX", "chrY")
bin_w    <- 1e4                               # localisation bin size (later parts)
#tracks   <- c(SegDup = "segdup.bed", SimpleRepeat = "simple_repeats.bed", CenSat = "censat.bed", Blacklist = "blacklist.bed", LowMappability = "low_mappability.bed")   # EDIT paths
#sv_bed   <- "sv_breakpoints.bed"              # EDIT
chr_cols <- c("name", "chr", "mate_chr")      # character columns in the feature files; everything else is numeric

# id and seed from .../kmer_<id>_NoFastP_NoKraken/<lost|ctrl>_features_k<seed>.tsv.gz
parse_path <- function(f) tibble( id = str_match(f, "/kmer_([^/]+)_NoFastP_NoKraken/")[, 2], kmer = as.integer(str_match(f, "_k(\\d+)\\.tsv\\.gz$")[, 2]))

add_class <- function(d) d %>% mutate( confident = mapq >= conf$mapq & aligned_frac >= conf$aligned, Class = case_when(mapq0 == 1 ~ "MAPQ0", aligned_frac < 0.5 ~ "Poorly aligned (<50%)", confident ~ "High confidence", TRUE ~ "Intermediate"))

# --- helpers used by the later annotation / enrichment parts ---
read_bed <- function(f) read_tsv(f, col_names = FALSE, col_types = cols_only(X1 = col_character(), X2 = col_integer(), X3 = col_integer()), comment = "#", progress = FALSE) %>% transmute(chr = X1, bs = X2 + 1L, be = X3)

# TRUE for rows of d (needs chr, lo, hi; 1-based inclusive) overlapping any interval in bed
in_track <- function(d, bed) {
	d <- mutate(d, .i = row_number())
	d$.i %in% semi_join(select(d, .i, chr, lo, hi), bed, by = join_by(chr, overlaps(lo, hi, bs, be)))$.i
}

# Mantel-Haenszel across individuals; OR > 1 = cases over-represented in the track
mh_test <- function(li, lo, ci, co) {
	if (length(li) < 2) return(tibble(OR = NA_real_, P = NA_real_))
	r <- tryCatch(mantelhaen.test(array(as.numeric(rbind(li, ci, lo, co)), c(2, 2, length(li)))), error = function(e) NULL)
	if (is.null(r)) tibble(OR = NA_real_, P = NA_real_) else tibble(OR = unname(r$estimate), P = r$p.value)
}

# d: id, kmer, case (lgl), in_tr (lgl)
mh_enrich <- function(d) {
	d %>% group_by(kmer, id) %>% summarise(li = sum(case & in_tr), lo = sum(case & !in_tr), ci = sum(!case & in_tr), co = sum(!case & !in_tr), .groups = "drop") %>% group_by(kmer) %>% group_modify(~ bind_cols(tibble(Pct_Case_In = 100 * sum(.x$li) / sum(.x$li + .x$lo), Pct_Other_In = 100 * sum(.x$ci) / sum(.x$ci + .x$co)), mh_test(.x$li, .x$lo, .x$ci, .x$co))) %>% ungroup()
}

####################################################################################################################
###                                  Read In Lost/Ctrl Read Features                                            ####
####################################################################################################################

lost_files <- Sys.glob(file.path(WRK_DIR, "data_raw", "kmer_*_NoFastP_NoKraken", "lost_features_k*.tsv.gz"))
ctrl_files <- gsub("lost_features", "ctrl_features", lost_files)

stopifnot(length(lost_files) > 0, all(file.exists(ctrl_files)))

all_k <- sort(unique(parse_path(lost_files)$kmer))
bit   <- setNames(bitwShiftL(1L, seq_along(all_k) - 1L), all_k)

# Read one feature file
read_feat <- function(file, grp, n_keep = Inf, exclude = character()) {
	d <- read_tsv(file, col_types = cols(.default = col_character()), na = c("", "NA"), progress = FALSE)
	if (!"chr" %in% names(d)) { message("  SKIPPED (empty / no chr column): ", basename(file)); return(NULL) }
	d <- d %>% filter(chr %in% chrs, !name %in% exclude) %>% mutate(across(-any_of(chr_cols), as.numeric))
	if (nrow(d) > n_keep) d <- slice_sample(d, n = n_keep)
	p <- parse_path(file)
	mutate(d, id = p$id, kmer = p$kmer, grp = grp)
}

read_files <- function(files, grp) {
	out <- vector("list", length(files))
	for (j in seq_along(files)) {
		t0 <- Sys.time()
		d  <- read_feat(files[j], grp)
		message(sprintf("  [%d/%d] %s: %s rows (%.0f s)", j, length(files), basename(files[j]), format(if (is.null(d)) 0L else nrow(d), big.mark = ","), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
		out[j] <- list(d)
	}
	bind_rows(out)
}

# Collapse all reads from one ID into one row (lowest-seed kept)
collapse_id <- function(x) {
	id_i <- x$id[1]
	x <- x %>% distinct(kmer, name, mate, .keep_all = TRUE) %>% arrange(kmer)
	ks <- sort(unique(x$kmer))

	mem <- x %>% mutate(b = unname(bit[as.character(kmer)])) %>% group_by(name, mate) %>% summarise(mask = sum(b), n_seeds = n(), first_k = min(kmer), .groups = "drop")

	# Percent of reads that are lost in subsequent kmers
	nest <- map_dfr(seq_len(length(ks) - 1), function(j) {
		a <- bitwAnd(mem$mask, bit[[as.character(ks[j])]]) > 0
		b <- bitwAnd(mem$mask, bit[[as.character(ks[j + 1])]]) > 0
		tibble(id = id_i, k_from = ks[j], k_to = ks[j + 1], n_from = sum(a), n_new = sum(b & !a), pct_from_still_lost = 100 * sum(a & b) / sum(a))
	})

	n_exp <- length(ks) - match(mem$first_k, ks) + 1L
	nm <- unique(x$name); nm <- sample(nm, min(1e5, length(nm)))
	same <- x %>% filter(name %in% nm) %>% group_by(name, mate) %>% filter(n() > 1) %>% summarise(same = n_distinct(chr, pos, mapq, aligned) == 1, .groups = "drop")
	chk <- tibble(id = id_i, n_reads = nrow(mem), pct_noncontiguous = 100 * mean(mem$n_seeds != n_exp), pct_features_same = 100 * mean(same$same))

	reads <- x %>% distinct(name, mate, .keep_all = TRUE) %>% left_join(select(mem, name, mate, mask, n_seeds), by = c("name", "mate"))
	list(reads = reads, nest = nest, chk = chk)
}

## Lost reads (read and collapsed one id at a time)
if (!file.exists("data_clean/lost_read_dedup.tsv")) {
	lost_by_id <- split(lost_files, parse_path(lost_files)$id)
	res <- list()
	for (id_i in names(lost_by_id)) {
		message(id_i, ": ", length(lost_by_id[[id_i]]), " seed file(s)")
		x <- read_files(lost_by_id[[id_i]], "lost")
		if (nrow(x) == 0) { message("  no reads, skipping"); next }
		res[[id_i]] <- collapse_id(x)
		rm(x); gc(verbose = FALSE)
	}
	reads_lost <- bind_rows(map(res, "reads")) %>% add_class()
	write_tsv(bind_rows(map(res, "nest")), "data_clean/lost_read_overlap_full.tsv")
	write_tsv(bind_rows(map(res, "chk")),  "data_clean/lost_read_overlap_summary.tsv")
	write_tsv(reads_lost, "data_clean/lost_read_dedup.tsv")
	rm(res)
} else {
	reads_lost <- read_tsv("data_clean/lost_read_dedup.tsv", show_col_types = FALSE)
	stopifnot(all(c("kmer", "grp", "Class") %in% names(reads_lost)))    # fails on a cache written by the old code: delete it
}

## Control reads (sampled reads matched to the number lost by sample x kmer)
n_need <- reads_lost %>% count(id, kmer, name = "n_lost")

if (!file.exists("data_clean/ctrl_read_subset.tsv")) {
	res <- vector("list", length(ctrl_files))
	for (j in seq_along(ctrl_files)) {
		p <- parse_path(ctrl_files[j])
		n_i <- n_need$n_lost[n_need$id == p$id & n_need$kmer == p$kmer]
		if (length(n_i) == 0 || n_i == 0) next      # no reads first lost at this seed
		t0 <- Sys.time()
		d <- read_feat(ctrl_files[j], "control", n_keep = n_i, exclude = reads_lost$name[reads_lost$id == p$id])
		message(sprintf("[%d/%d] %s k%d: kept %s of %s needed (%.0f s)", j, length(ctrl_files), p$id, p$kmer, format(if (is.null(d)) 0L else nrow(d), big.mark = ","), format(n_i, big.mark = ","), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
		res[j] <- list(d)
	}
	reads_ctrl <- bind_rows(res) %>% add_class()
	write_tsv(reads_ctrl, "data_clean/ctrl_read_subset.tsv")
	rm(res)
} else {
	reads_ctrl <- read_tsv("data_clean/ctrl_read_subset.tsv", show_col_types = FALSE)
}

reads <- bind_rows(reads_lost, reads_ctrl)

####################################################################################################################
###                                  Basic QC Lost/Ctrl Read Features                                           ####
####################################################################################################################

comp <- reads %>% count(grp, id, kmer, Class) %>% group_by(grp, id, kmer) %>% mutate(Pct = 100 * n / sum(n)) %>% ungroup()
write_tsv(comp, "data_clean/lost_read_class_composition.tsv")

avg  <- function(x) mean(x, na.rm = TRUE)
feat <- reads %>% group_by(grp, kmer) %>% summarise(N = n(), MAPQ0 = avg(mapq0), Aligned_Frac = median(aligned_frac, na.rm = TRUE), Double_Clip = avg(dbl_clip), Split_SA = avg(has_sa), Multi_Hit = avg(n_xa > 0), Longest_Exact = median(longest_exact, na.rm = TRUE), Proper_Pair = avg(proper), Mate_Mapped = avg(mate_mapped), Mate_MAPQ0 = avg(mate_mapq0), Mate_Also_Lost = avg(mate_also_lost), .groups = "drop")
write_tsv(feat, "data_clean/lost_vs_control_features.tsv")

# Check to ensure that lost reads aren't well aligned
reads %>% filter(grp == "lost", !is.na(longest_exact)) %>% group_by(kmer) %>% summarise(Pct_Exact_Run_ge_Seed = 100 * mean(longest_exact >= kmer), .groups = "drop")

p <- comp %>% filter(grp == "lost") %>% ggplot(aes(id, Pct, fill = Class)) + 
	geom_col() +
	facet_grid(~ kmer, scales = "free_x", space = "free_x") + 
	labs(x = NULL, y = "% of lost reads") +
	theme_classic() + 
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("graphs/Lost_Read_Classes.pdf", p, width = 12, height = 6)

####################################################################################################################
###                                 	     Localised or Sporadic	                                            ####
####################################################################################################################

min_ids <- 3    # a hotspot recurs if the same bin is hot in at least this many individuals (same as r in recur())

bins <- reads %>% count(id, kmer, grp, chr, bin = pos %/% bin_w) %>% pivot_wider(names_from = grp, values_from = n, values_fill = 0L) %>% group_by(id, kmer) %>%
	mutate(pr = (control + 0.5) / sum(control + 0.5), expct = sum(lost) * pr, enrich = (lost + 0.5) / (expct + 0.5),
	       pval = pbinom(lost - 1, sum(lost), pr, lower.tail = FALSE), padj = p.adjust(pval, "BH"),
	       hot = padj < 0.01 & enrich >= 3 & lost >= 10,
	       bin_start = bin * bin_w, bin_end = (bin + 1) * bin_w - 1) %>% ungroup()  

# optional: which annotation tracks a region overlaps (tracks that don't exist are skipped)
annotate <- function(d) {
	d$Overlaps <- ""
	for (tr in names(tracks)) {
		if (!file.exists(tracks[[tr]])) next
		hit <- in_track(transmute(d, chr, lo = start, hi = end), read_bed(tracks[[tr]]))
		d$Overlaps[hit] <- paste0(d$Overlaps[hit], if_else(d$Overlaps[hit] == "", "", ","), tr)
	}
	d
}

## Per-individual hotspots: bin level, then adjacent hot bins merged into regions
hot_bins <- bins %>% filter(hot) %>% select(id, kmer, chr, bin, bin_start, bin_end, lost, control, expct, enrich, pval, padj) %>% arrange(id, kmer, factor(chr, chrs), bin)
write_tsv(hot_bins, "data_clean/lost_read_hotspot_bins.tsv")

hot_regions <- hot_bins %>% group_by(id, kmer, chr) %>% mutate(region = cumsum(bin != lag(bin, default = -Inf) + 1)) %>% group_by(id, kmer, chr, region) %>%
	summarise(start = min(bin_start), end = max(bin_end), n_bins = n(), lost = sum(lost), control = sum(control), expected = sum(expct), enrich = (lost + 0.5) / (expected + 0.5), min_padj = min(padj), .groups = "drop") %>%
	select(-region) %>% arrange(id, kmer, factor(chr, chrs), start) %>% annotate()
write_tsv(hot_regions, "data_clean/lost_read_hotspot_regions.tsv")

## Hotspots recurring across individuals (same bin hot in >= min_ids individuals), then merged into regions
recur_bins <- bins %>% filter(hot) %>% group_by(kmer, chr, bin, bin_start, bin_end) %>% summarise(n_ids = n_distinct(id), ids = paste(sort(unique(id)), collapse = ","), lost = sum(lost), control = sum(control), .groups = "drop") %>%
	filter(n_ids >= min_ids) %>% arrange(kmer, factor(chr, chrs), bin)

recur_regions <- recur_bins %>% group_by(kmer, chr) %>% mutate(region = cumsum(bin != lag(bin, default = -Inf) + 1)) %>% group_by(kmer, chr, region) %>%
	summarise(start = min(bin_start), end = max(bin_end), n_bins = n(), max_ids_in_a_bin = max(n_ids), ids = paste(sort(unique(unlist(strsplit(ids, ",")))), collapse = ","), lost = sum(lost), control = sum(control), .groups = "drop") %>%
	mutate(n_ids = lengths(strsplit(ids, ","))) %>% select(-region) %>% arrange(kmer, factor(chr, chrs), start) %>% annotate()
write_tsv(recur_regions, "data_clean/lost_read_hotspot_recurrent_regions.tsv")

# BED for IGV / bedtools (0-based start)
recur_regions %>% transmute(chr, start = pmax(start - 1, 0), end, name = paste0("k", kmer, "_n", n_ids), score = lost) %>% write_tsv("data_clean/lost_read_hotspot_recurrent_regions.bed", col_names = FALSE)

count(hot_regions, kmer, name = "regions"); count(recur_regions, kmer, name = "recurrent_regions")




## ---------- 1D. Enrichment in annotated regions (read level, Mantel-Haenszel across individuals) ----------
enr <- imap_dfr(tracks, function(f, tr) {
	bed <- read_bed(f)
	d   <- mutate(reads, case = grp == "lost", in_tr = in_track(reads, bed))
	bind_rows(mutate(mh_enrich(d), Subset = "all reads"), mutate(mh_enrich(filter(d, confident)), Subset = "high confidence")) %>% mutate(Track = tr, .before = 1)
})
write_tsv(enr, "data_clean/lost_read_annotation_enrichment.tsv")

## ---------- 1E. Distance to SV breakpoints (nearest breakpoint on either side) ----------
bp <- read_bed(sv_bed) %>% transmute(chr, bp_pos = (bs + be) %/% 2L) %>% distinct()
x  <- select(reads, chr, pos)
up <- left_join(x, bp, by = join_by(chr, closest(pos >= bp_pos)))
dn <- left_join(x, bp, by = join_by(chr, closest(pos <= bp_pos)))
stopifnot(nrow(up) == nrow(x), nrow(dn) == nrow(x))
reads$bp_dist <- pmin(abs(up$pos - up$bp_pos), abs(dn$bp_pos - dn$pos), na.rm = TRUE)

reads %>% group_by(grp, type, kmer) %>% summarise(Within_1kb = 100 * mean(bp_dist <= 1e3, na.rm = TRUE), Within_10kb = 100 * mean(bp_dist <= 1e4, na.rm = TRUE), Median_Dist = median(bp_dist, na.rm = TRUE), .groups = "drop") %>%
	write_tsv("data_clean/lost_read_bp_distance.tsv")

####################################################################################################################
###                                     Region-Level SV Evidence and Depth                                      ####
####################################################################################################################

## ---------- 2.0 Load: sv_evidence_*.tsv ----------
sv <- Sys.glob(file.path(base, "kmer_*", "features_*", "sv_evidence_*.tsv")) %>% map_dfr(function(f) read_tsv(f, col_types = cols(.default = col_integer(), id = col_character(), type = col_character(), window = col_character()), progress = FALSE) %>% mutate(kmer = parse_path(f)$kmer))

## ---------- 2A. Change in SV evidence at breakpoint windows, k19 vs k (paired within each seed's files) ----------
sv_cmp <- sv %>% filter(aln == 19 | aln == kmer) %>% mutate(aln = if_else(aln == 19, "k19", "kk")) %>% pivot_wider(names_from = aln, values_from = c(total, q20, mapq0, clipped, split, discordant), names_sep = "_") %>%
	mutate(ev19 = split_k19 + discordant_k19, evk = split_kk + discordant_kk, Set = if_else(str_starts(window, "rand_"), "Matched random windows", "SV breakpoint windows"))   # name control windows with a rand_ prefix

sv_sum <- sv_cmp %>% filter(ev19 >= 5) %>% group_by(Set, type, kmer) %>% summarise(Windows = n(), Median_Evidence_Retained = median(evk / ev19), Pct_Lose_25pct = 100 * mean(evk / ev19 < 0.75), Median_Q20_Depth_Retained = median(q20_kk / q20_k19, na.rm = TRUE), .groups = "drop")
write_tsv(sv_sum, "data_clean/sv_evidence_retention.tsv")

## ---------- 2B. Regional depth dropout (mosdepth 10 kb bins) ----------
depth_files <- Sys.glob(file.path(base, "kmer_*", "features_*", "depth_k*_Q*.regions.bed.gz"))
dp <- map_dfr(depth_files, function(f) {
	m  <- parse_path(f)
	fn <- str_match(basename(f), "^depth_k(\\d+)_Q(\\d+)\\.")
	read_tsv(f, col_names = c("chr", "s", "e", "depth"), col_types = cols(chr = col_character(), s = col_integer(), e = col_integer(), depth = col_double()), progress = FALSE) %>% filter(chr %in% chrs) %>% mutate(id = m$id, type = m$type, kmer = m$kmer, aln = as.integer(fn[2]), Q = as.integer(fn[3]))
})
dw <- dp %>% mutate(src = if_else(aln == 19, "d19", "dk")) %>% select(-aln) %>% pivot_wider(names_from = src, values_from = depth) %>% group_by(id, type, kmer, Q) %>% mutate(ratio = (dk / d19) / median((dk / d19)[d19 >= 10], na.rm = TRUE), dropout = d19 >= 10 & ratio < 0.8) %>% ungroup()   # >20% depth loss vs genome-wide median ratio

dropout_sum <- dw %>% filter(d19 >= 10) %>% group_by(id, type, kmer, Q) %>% summarise(Pct_Dropout_Bins = 100 * mean(dropout), Dropout_Mb = sum(dropout) * 0.01, .groups = "drop") %>% group_by(type, kmer, Q) %>% summarise(N_Samples = n(), across(c(Pct_Dropout_Bins, Dropout_Mb), median), .groups = "drop")
write_tsv(dropout_sum, "data_clean/depth_dropout_summary.tsv")

dropout_runs <- dw %>% filter(dropout) %>% arrange(id, type, kmer, Q, chr, s) %>% group_by(id, type, kmer, Q) %>% mutate(bin = s %/% 1e4, run = cumsum(is.na(lag(bin)) | chr != lag(chr) | bin != lag(bin) + 1)) %>% group_by(id, type, kmer, Q, run) %>% summarise(chr = first(chr), start = min(s), end = max(e), Len_kb = n() * 10, .groups = "drop")
write_tsv(dropout_runs, "data_clean/depth_dropout_regions.tsv")

## ---------- 2C. Do dropout bins sit in annotated regions? (reuses the shared helpers) ----------
dq <- 20                                          # MAPQ threshold of the depth track to use
drop_bins <- dw %>% filter(Q == dq, d19 >= 10) %>% mutate(lo = s + 1L, hi = e)
enr_depth <- imap_dfr(tracks, function(f, tr) {
	bed <- read_bed(f)
	mh_enrich(mutate(drop_bins, case = dropout, in_tr = in_track(drop_bins, bed))) %>% mutate(Track = tr, .before = 1)
})
write_tsv(enr_depth, "data_clean/depth_dropout_annotation_enrichment.tsv")

####################################################################################################################
###                                               Finish                                                        ####
####################################################################################################################

quit()
