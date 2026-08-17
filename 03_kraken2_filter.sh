#!/usr/bin/env bash
#SBATCH --cpus-per-task 2
#SBATCH --mem 4G
#SBATCH --error slurm/KF_%j.out
#SBATCH --output slurm/KF_%j.out
#SBATCH --time 4:00:00
#
# 03_kraken2_filter.sh
#
# Uses the output of 02_kraken2_classify.sh to remove reads NOT matching a
# given set of "retain" taxids from the original fastq files. By default
# (no -t given) this retains only reads classified as human (NCBI taxid
# 9606) and drops everything else as contamination -- pass a taxid lineage
# file (see 01_build_taxonomic_lineage.R) via -t to retain a different
# taxon (or its full lineage) instead.
#
# Run directly (bash 03_kraken2_filter.sh ...) or via `sbatch` on Slurm.
#
# Requirements: seqkit must be available on PATH, or loadable via
# `module load $SEQKIT_MODULE` if the SEQKIT_MODULE environment variable is
# set (for Lmod/Environment-Modules based clusters).

set -euo pipefail

(echo "$(date) on $(hostname)" 1>&2)
(echo "$0 $*" 1>&2)

usage() {
(echo -e "\
*************************************
* Filters Kraken2-classified reads  *
* to remove non-target contaminants *
*************************************

Usage: $0 -1 R1.fastq.gz -2 R2.fastq.gz -k /path/to/classify_output -o /path/to/output [options]

Required:
  -1 [FILE]   Original Read 1 fastq(.gz) (the same file given to 02_kraken2_classify.sh)
  -2 [FILE]   Original Read 2 fastq(.gz)
  -k [DIR]    Directory containing 02_kraken2_classify.sh's output for this sample

Optional:
  -u [FILE]   Original unpaired reads fastq(.gz), if classified in step 2
  -o [DIR]    Output directory for filtered fastqs (default: same as -k)
  -s [STR]    Sample ID, used only for log readability (default: 'sample')
  -t [FILE]   Taxid lineage file (one taxid per line) to RETAIN during
              filtering -- see 01_build_taxonomic_lineage.R. If omitted,
              defaults to retaining only taxid 9606 (human).

Environment:
  SEQKIT_MODULE   If set, 'module load \$SEQKIT_MODULE' is run before
                  invoking seqkit. If unset, seqkit is expected on PATH.
*********************************" 1>&2)
}

SAMPLE="sample"
R1=""; R2=""; UNPAIRED=""; KRAKENDIR=""; OUTDIR=""; TAXIDFILE=""

while getopts "s:1:2:u:k:o:t:h" OPTION
do
	case $OPTION in
		s) SAMPLE=${OPTARG} ;;
		1) R1=${OPTARG} ;;
		2) R2=${OPTARG} ;;
		u) UNPAIRED=${OPTARG} ;;
		k) KRAKENDIR=${OPTARG} ;;
		o) OUTDIR=${OPTARG} ;;
		t) TAXIDFILE=${OPTARG} ;;
		h) usage; exit 0 ;;
		?) usage; exit 1 ;;
	esac
done

if [[ -z "${R1}" || -z "${R2}" || -z "${KRAKENDIR}" ]]; then
	echo "FAIL: Missing required parameter(s)." 1>&2
	usage
	exit 1
fi
for f in "${R1}" "${R2}"; do
	if [[ ! -e "${f}" ]]; then
		echo "FAIL: Input file ${f} does not exist!" 1>&2
		exit 1
	fi
done
if [[ -n "${UNPAIRED}" && ! -e "${UNPAIRED}" ]]; then
	echo "FAIL: Unpaired input file ${UNPAIRED} does not exist!" 1>&2
	exit 1
fi
if [[ -n "${TAXIDFILE}" && ! -e "${TAXIDFILE}" ]]; then
	echo "FAIL: Taxid lineage file ${TAXIDFILE} does not exist!" 1>&2
	exit 1
fi
: "${OUTDIR:=${KRAKENDIR}}"

READNAME="${SAMPLE}"
CLASSIFIED_R1="${KRAKENDIR}/${READNAME}_Kraken2Screened_Classified_1.fastq"
CLASSIFIED_R2="${KRAKENDIR}/${READNAME}_Kraken2Screened_Classified_2.fastq"
CLASSIFIED_UNPAIRED="${KRAKENDIR}/${READNAME}_Kraken2Screened_Unpaired_Classified.fastq"

for f in "${CLASSIFIED_R1}" "${CLASSIFIED_R2}"; do
	if [[ ! -e "${f}" ]]; then
		echo "FAIL: Expected classify-step output not found: ${f}" 1>&2
		echo "      Has 02_kraken2_classify.sh been run for this sample?" 1>&2
		exit 1
	fi
done

printf "%-22s%s\n" "Sample ID" "${SAMPLE}" 1>&2
printf "%-22s%s\n" "Kraken2 output dir" "${KRAKENDIR}" 1>&2
printf "%-22s%s\n" "Output directory" "${OUTDIR}" 1>&2
[[ -n "${TAXIDFILE}" ]] && printf "%-22s%s\n" "Retain taxid file" "${TAXIDFILE}" || printf "%-22s%s\n" "Retain taxid" "9606 (human, default)"

if [[ -n "${SEQKIT_MODULE:-}" ]]; then
	module purge
	module load "${SEQKIT_MODULE}"
elif command -v module >/dev/null 2>&1 && module avail seqkit/2.8.2 &>/dev/null; then
	# Fall back to the version this pipeline was tested against, if available
	module purge
	module load seqkit/2.8.2
fi
if ! command -v seqkit >/dev/null 2>&1; then
	echo "FAIL: seqkit not found on PATH. Set SEQKIT_MODULE or add seqkit to PATH." 1>&2
	exit 1
fi

mkdir -p "${OUTDIR}"
HEADER="KF"
CONTAM_LIST="${OUTDIR}/${READNAME}_Readnames_Contaminants.list"

if [[ ! -f "${OUTDIR}/kraken2readlist.done" ]]; then
	TMP_LIST="${OUTDIR}/${READNAME}_Readnames_Contaminants_temp.list"
	: > "${TMP_LIST}"

	extract_non_retained() {
		local fastq="$1"
		if [[ -n "${TAXIDFILE}" ]]; then
			seqkit seq "${fastq}" -n | grep -v -f <(sed 's/^/kraken:taxid|/' "${TAXIDFILE}") | cut -d' ' -f1 >> "${TMP_LIST}"
		else
			seqkit seq "${fastq}" -n | grep -v 'kraken:taxid|9606' | cut -d' ' -f1 >> "${TMP_LIST}"
		fi
	}

	echo "${HEADER}: building contaminant read list from ${CLASSIFIED_R1}" 1>&2
	extract_non_retained "${CLASSIFIED_R1}"
	echo "${HEADER}: building contaminant read list from ${CLASSIFIED_R2}" 1>&2
	extract_non_retained "${CLASSIFIED_R2}"
	if [[ -e "${CLASSIFIED_UNPAIRED}" ]]; then
		echo "${HEADER}: building contaminant read list from ${CLASSIFIED_UNPAIRED}" 1>&2
		extract_non_retained "${CLASSIFIED_UNPAIRED}"
	fi

	sort --unique "${TMP_LIST}" -o "${CONTAM_LIST}"
	rm -f "${TMP_LIST}"
	touch "${OUTDIR}/kraken2readlist.done"
fi

if [[ ! -f "${OUTDIR}/kraken2readfilter1.done" ]]; then
	echo "${HEADER}: filtering ${R1}" 1>&2
	seqkit grep -v -f "${CONTAM_LIST}" "${R1}" -o "${OUTDIR}/${READNAME}_kraken_filtered_R1.fastq.gz"
	touch "${OUTDIR}/kraken2readfilter1.done"
fi

if [[ ! -f "${OUTDIR}/kraken2readfilter2.done" ]]; then
	echo "${HEADER}: filtering ${R2}" 1>&2
	seqkit grep -v -f "${CONTAM_LIST}" "${R2}" -o "${OUTDIR}/${READNAME}_kraken_filtered_R2.fastq.gz"
	touch "${OUTDIR}/kraken2readfilter2.done"
fi

if [[ -n "${UNPAIRED}" && -e "${CLASSIFIED_UNPAIRED}" && ! -f "${OUTDIR}/kraken2readfilterunpaired.done" ]]; then
	echo "${HEADER}: filtering ${UNPAIRED}" 1>&2
	seqkit grep -v -f "${CONTAM_LIST}" "${UNPAIRED}" -o "${OUTDIR}/${READNAME}_kraken_filtered_unpaired.fastq.gz"
	touch "${OUTDIR}/kraken2readfilterunpaired.done"
fi

# Remove the large intermediate classified fastqs now that filtering is done
rm -f "${CLASSIFIED_R1}" "${CLASSIFIED_R2}" "${CLASSIFIED_UNPAIRED}"

touch "${OUTDIR}/kraken2filter.done"
echo "${HEADER}: Done. Filtered reads in ${OUTDIR}" 1>&2

exit 0
