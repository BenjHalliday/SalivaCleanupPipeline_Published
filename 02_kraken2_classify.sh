#!/bin/bash
#SBATCH --cpus-per-task 6
#SBATCH --mem 74G
#SBATCH --error slurm/KC_%j.out
#SBATCH --output slurm/KC_%j.out
#SBATCH --time 4:00:00

# 02_kraken2_classify.sh
# Runs Kraken2 taxonomic classification on paired-end (and optionally unpaired) fastq reads. Intended as step 1 of a two-step contamination screening pipeline; step 2 is 03_kraken2_filter.sh.
# This script can be run directly (bash 02_kraken2_classify.sh ...) or submitted via sbatch on a Slurm cluster
# Requirements: kraken2 must be available on PATH, or loadable via `module load $KRAKEN2_MODULE` if the KRAKEN2_MODULE environment variable is set.

set -euo pipefail

(echo "$(date) on $(hostname)" 1>&2)
(echo "$0 $*" 1>&2)

usage() {
(echo -e "\
*************************************
* This script uses Kraken2 to       *
* categorise reads taxonomically    *
*************************************

Usage: $0 -1 R1.fastq.gz -2 R2.fastq.gz -d /path/to/kraken2_db -o /path/to/output [options]

Required:
  -1 [FILE]   Read 1 fastq(.gz)
  -2 [FILE]   Read 2 fastq(.gz)
  -d [DIR]    Full path to the Kraken2 database directory
  -o [DIR]    Output directory (created if it doesn't exist)

Optional:
  -u [FILE]   Unpaired reads fastq(.gz) (e.g. surviving singletons after
              adapter/quality trimming), classified as an additional step
  -s [STR]    Sample ID, used only for log/job-name readability (default: 'sample')
  -p [INT]    Threads (default: 6)

Environment:
  KRAKEN2_MODULE   If set, 'module load \$KRAKEN2_MODULE' is run before
                   invoking kraken2 (for Lmod/Environment-Modules clusters).
                   If unset, kraken2 is expected to already be in PATH.
*********************************" 1>&2)
}

SAMPLE="sample"
THREADS=6
R1=""; R2=""; UNPAIRED=""; KRAKENDB=""; OUTDIR=""

while getopts "s:1:2:u:d:o:p:h" OPTION
do
	case $OPTION in
		s) SAMPLE=${OPTARG} ;;
		1) R1=${OPTARG} ;;
		2) R2=${OPTARG} ;;
		u) UNPAIRED=${OPTARG} ;;
		d) KRAKENDB=${OPTARG} ;;
		o) OUTDIR=${OPTARG} ;;
		p) THREADS=${OPTARG} ;;
		h) usage; exit 0 ;;
		?) usage; exit 1 ;;
	esac
done

if [[ -z "${R1}" || -z "${R2}" || -z "${KRAKENDB}" || -z "${OUTDIR}" ]]; then
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
if [[ ! -e "${KRAKENDB}" ]]; then
	echo "FAIL: Kraken2 database ${KRAKENDB} does not exist!" 1>&2
	exit 1
fi
if [[ -n "${UNPAIRED}" && ! -e "${UNPAIRED}" ]]; then
	echo "FAIL: Unpaired input file ${UNPAIRED} does not exist!" 1>&2
	exit 1
fi

printf "%-22s%s\n" "Sample ID" "${SAMPLE}" 1>&2
printf "%-22s%s\n" "Read 1" "${R1}" 1>&2
printf "%-22s%s\n" "Read 2" "${R2}" 1>&2
[[ -n "${UNPAIRED}" ]] && printf "%-22s%s\n" "Unpaired reads" "${UNPAIRED}" 1>&2
printf "%-22s%s\n" "Kraken2 database" "${KRAKENDB}" 1>&2
printf "%-22s%s\n" "Output directory" "${OUTDIR}" 1>&2

if [[ -n "${KRAKEN2_MODULE:-}" ]]; then
	module purge
	module load "${KRAKEN2_MODULE}"
fi
if ! command -v kraken2 >/dev/null 2>&1; then
	echo "FAIL: kraken2 not found in PATH. Set KRAKEN2_MODULE or add kraken2 to PATH." 1>&2
	exit 1
fi

mkdir -p "${OUTDIR}"
READNAME="${SAMPLE}"
HEADER="KC"

if [[ ! -f "${OUTDIR}/kraken2paired_classified.done" ]]; then
	echo "${HEADER}: ${KRAKENDB} ${R1} ${R2}" 1>&2
	kraken2 --db "${KRAKENDB}" --threads "${THREADS}" \
		--output "${OUTDIR}/${READNAME}_Output.tsv" \
		--classified-out "${OUTDIR}/${READNAME}_Kraken2Screened_Classified#.fastq" \
		--unclassified-out "${OUTDIR}/${READNAME}_Kraken2Screened_Unclassified#.fastq" \
		--report "${OUTDIR}/${READNAME}_Kraken2Screened_Report.tsv" \
		--paired --gzip-compressed --use-names "${R1}" "${R2}"
	touch "${OUTDIR}/kraken2paired_classified.done"
fi

if [[ -n "${UNPAIRED}" && ! -f "${OUTDIR}/kraken2unpaired_classified.done" ]]; then
	echo "${HEADER}: ${KRAKENDB} ${UNPAIRED}" 1>&2
	kraken2 --db "${KRAKENDB}" --threads "${THREADS}" \
		--output "${OUTDIR}/${READNAME}_Unpaired_Output.tsv" \
		--classified-out "${OUTDIR}/${READNAME}_Kraken2Screened_Unpaired_Classified.fastq" \
		--unclassified-out "${OUTDIR}/${READNAME}_Kraken2Screened_Unpaired_Unclassified.fastq" \
		--report "${OUTDIR}/${READNAME}_Unpaired_Kraken2Screened_Report.tsv" \
		--gzip-compressed --use-names "${UNPAIRED}"
	touch "${OUTDIR}/kraken2unpaired_classified.done"
fi

touch "${OUTDIR}/kraken2classify.done"
echo "${HEADER}: Done. Outputs in ${OUTDIR}" 1>&2

exit 0
