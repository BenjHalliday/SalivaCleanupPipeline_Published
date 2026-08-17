# Kraken2-based host-microbiome contamination filtering

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Scripts for taxonomically classifying and removing non-target (e.g.
microbial) contaminant reads from short-read sequencing data, as used in
[Citation - add on publication].

Repository: https://github.com/BenjHalliday/SalivaCleanupPipeline_Published

## Pipeline overview

1. **`01_build_taxonomic_lineage.R`** — downloads the NCBI taxonomy dump and
   builds the full ancestor lineage for a target taxid (e.g. `9606` for
   human). Produces a taxid list that can be used as an inclusion list in
   step 3, so that filtering retains an entire taxonomic lineage rather than
   only exact matches to the target taxid.
2. **`02_kraken2_classify.sh`** — runs Kraken2 on paired-end (and optionally
   unpaired) fastq reads against a chosen database, producing classified and
   unclassified read sets plus a taxonomic report.
3. **`03_kraken2_filter.sh`** — uses step 2's classified output to identify
   reads not belonging to the retained taxon/lineage, and filters them out
   of the original fastq files, producing decontaminated fastqs.

```
(optional) 01_build_taxonomic_lineage.R --> taxid list
                                                 v         
raw fastqs --> 02_kraken2_classify.sh --> 03_kraken2_filter.sh --> filtered fastqs
```

## Requirements

- [Kraken2](https://github.com/DerrickWood/kraken2) (tested with v2.1.3)
- [seqkit](https://github.com/shenwei356/seqkit) (tested with v2.4.0 and v2.8.2)
- R (≥4.0) with `tidyverse` and `data.table`
- A Kraken2-compatible database (e.g. the [standard database](https://benlangmead.github.io/aws-indexes/k2), or a custom database, as detailed in the manuscript)

Both shell scripts can be run as standalones (`bash script.sh ...`) or can be submitted
via `sbatch` on a Slurm cluster. If `kraken2`/`seqkit` are not already
in `PATH`, set the `KRAKEN2_MODULE` / `SEQKIT_MODULE` environment variables
to have the scripts `module load` them.

## Usage

```bash
# 1. (optional) Build a taxid inclusion list, e.g. the full human lineage
Rscript 01_build_taxonomic_lineage.R 9606 ./resources/NCBI_Taxonomy

# 2. Classify reads
bash 02_kraken2_classify.sh \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -d /path/to/kraken2_db \
  -o results/sample/kraken2 \
  -s sample01 -p 8

# 3. Filter out non-retained reads (defaults to retaining only taxid 9606)
bash 03_kraken2_filter.sh \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  -k results/sample/kraken2 \
  -s sample01 \
  -t ./resources/NCBI_Taxonomy/Lineages/9606_taxids.tsv
```

Filtered, decontaminated reads are written as
`<sample>_kraken_filtered_R1.fastq.gz` / `_R2.fastq.gz` (and
`_unpaired.fastq.gz` if applicable) in the specified output directory.

## Testing on public data

For testing this pipeline against real publicly available human saliva WGS data, we'd suggest:

- **SRA accessions ERX1462737, ERX1462740, SRX2830683, SRX2830684, SRX2830689**
  — saliva-derived WGS datasets used as public test cases in a
   related study (Samson et al., *Sci Rep* 2020,
  [10.1038/s41598-020-76022-4](https://doi.org/10.1038/s41598-020-76022-4)),
  chosen because of known non-human contamination.
  Download via `sra-tools` (`prefetch` then `fasterq-dump`).

## Notes on adapting this pipeline

- `01_build_taxonomic_lineage.R` downloads a large taxonomy dump
  (`taxdump.tar.gz`, several hundred MB) from NCBI on first run — this is
  cached in the given output directory and not re-downloaded on subsequent runs.
- The default filtering behaviour for `03_kraken2_filter.sh` (no `-t` given) retains only reads
  classified as human (taxid `9606`). Change the hardcoded `9606` or pass `-t` if your target organism differs.

## License

MIT — see [LICENSE](LICENSE).

## Citation

If you use these scripts, please cite:
[Citation - add on publication]

