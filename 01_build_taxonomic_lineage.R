#!/usr/bin/env Rscript
####################################################################################################################
###                                    BuildTaxonomicLineage.R                                                  ####
####################################################################################################################
#
# Downloads the NCBI taxonomy database and builds the full ancestor lineage
# for a given taxid, writing out a plain list of taxids (root -> ... ->
# target) that can be used to define an "inclusion list" for kraken2filter.sh
# (i.e. the set of taxids that should be RETAINED, e.g. the full human
# lineage under taxid 9606, when filtering out everything else as
# contamination).
#
# Usage:
#   Rscript BuildTaxonomicLineage.R <taxid> <output_directory>
#
# Example:
#   Rscript BuildTaxonomicLineage.R 9606 ./resources/NCBI_Taxonomy
#
# Requirements: R packages `tidyverse` and `data.table` (installed
# automatically on first run if missing).
#
# Output:
#   <output_directory>/Lineages/<taxid>_full.tsv    -- taxid, name, rank for every ancestor
#   <output_directory>/Lineages/<taxid>_taxids.tsv  -- one taxid per line, for use with kraken2filter.sh -t
#
####################################################################################################################
###                                          Parse arguments                                                    ####
####################################################################################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2 || args[1] %in% c("-h", "--help")) {
	cat(
		"Usage: Rscript BuildTaxonomicLineage.R <taxid> <output_directory>\n",
		"  <taxid>             NCBI taxid to build the ancestor lineage for (e.g. 9606 for human)\n",
		"  <output_directory>  Directory to store the downloaded taxonomy dump and output lineage files\n",
		sep = ""
	)
	quit(status = if (length(args) != 2) 1 else 0)
}

taxid <- args[1]
dest <- args[2]

if (!grepl("^[0-9]+$", taxid)) {
	stop("taxid must be numeric, got: '", taxid, "'")
}

####################################################################################################################
###                                          Directory creation                                                 ####
####################################################################################################################

if (!dir.exists(dest)) dir.create(dest, recursive = TRUE)
setwd(dest)
if (!dir.exists("Lineages")) dir.create("Lineages")

####################################################################################################################
###                                              Packages                                                       ####
####################################################################################################################

required_packages <- c("tidyverse", "data.table")
for (pkg in required_packages) {
	if (!requireNamespace(pkg, quietly = TRUE)) {
		message("Installing missing package: ", pkg)
		install.packages(pkg, repos = "https://cloud.r-project.org")
	}
}
suppressPackageStartupMessages({
	library(tidyverse)
	library(data.table)
})

message("Packages loaded.")

####################################################################################################################
###                                          Pull NCBI taxonomy                                                 ####
####################################################################################################################

url <- "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
destfile <- "taxdump.tar.gz"

if (!file.exists(destfile)) {
	download.file(url, destfile, mode = "wb")
	message("Download complete.")
} else {
	message("File already exists, skipping download: ", destfile)
}

if (!file.exists("nodes.dmp") || !file.exists("names.dmp")) {
	untar(destfile)
	message("Extraction complete.")
} else {
	message("nodes.dmp / names.dmp already present, skipping extraction.")
}

####################################################################################################################
###                                          Build lineage                                                      ####
####################################################################################################################

taxids_out <- file.path("Lineages", paste0(taxid, "_taxids.tsv"))
full_out <- file.path("Lineages", paste0(taxid, "_full.tsv"))

if (!file.exists(taxids_out)) {

	nodes <- read_delim("nodes.dmp", delim = "|", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE) %>%
		select(taxid = X1, parent_taxid = X2, rank = X3) %>%
		mutate(across(everything(), trimws)) %>%
		mutate(taxid = as.character(taxid))

	names_tbl <- read_delim("names.dmp", delim = "|", col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE) %>%
		select(taxid = X1, name = X2, name_class = X4) %>%
		filter(name_class == "scientific name") %>%
		select(taxid, name) %>%
		mutate(taxid = as.character(taxid))

	nodes_with_names <- nodes %>% left_join(names_tbl, by = "taxid")

	if (!(taxid %in% nodes_with_names$taxid)) {
		stop("taxid '", taxid, "' was not found in the NCBI taxonomy dump (nodes.dmp). ",
		     "Check the taxid is correct and current.")
	}

	get_lineage <- function(start_taxid) {
		lineage <- tibble(taxid = character(), name = character(), rank = character())
		current_taxid <- start_taxid
		visited <- character()
		while (!is.na(current_taxid) && current_taxid != "1") {
			if (current_taxid %in% visited) {
				warning("Cycle detected in taxonomy graph at taxid ", current_taxid, "; stopping lineage walk.")
				break
			}
			visited <- c(visited, current_taxid)
			current_node <- nodes_with_names %>%
				filter(taxid == current_taxid) %>%
				select(taxid, name, rank, parent_taxid)
			if (nrow(current_node) == 0) break
			lineage <- bind_rows(lineage, select(current_node, taxid, name, rank))
			current_taxid <- current_node$parent_taxid
		}
		lineage <- lineage %>% bind_rows(tibble(taxid = "1", name = "root", rank = "no rank"))
		return(lineage)
	}

	lineage <- get_lineage(taxid)
	print(lineage, n = 100)

	fwrite(lineage, full_out, sep = "\t")
	write_lines(lineage$taxid, taxids_out)

	message("Wrote ", nrow(lineage), " lineage entries to:")
	message("  ", full_out)
	message("  ", taxids_out)

} else {
	message("Lineage already exists, skipping: ", taxids_out)
}

####################################################################################################################
###                                                 fin                                                         ####
####################################################################################################################

quit(status = 0)
