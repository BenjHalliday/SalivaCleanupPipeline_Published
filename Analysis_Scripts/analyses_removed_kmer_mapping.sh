#!/bin/bash 
#SBATCH -c 2    # cores requested
#SBATCH --mem=10G  # memory in Mb
#SBATCH -e slurm_kmer_alignmnet_mapping_outfile_%A_%a.txt
#SBATCH -o slurm_kmer_alignmnet_mapping_outfile_%A_%a.txt
#SBATCH -t 24:00:00  # time requested in hour:minute:second

### Extracting reads that are newly unmapped when the aligner seed length is increased from k19 to k${kmer}
### Newly unmapped = unmapped at k${kmer}, not unmapped at k19. Only the unmapped mate of the pair is extracted.

set -eo pipefail
export LC_ALL=C    # same collation for sort and comm

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load samtools

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

id=${1:-}
kmer=${2:-}

## Config
if [ -n "${SLURM_JOB_ID:-}" ]
then
	script_path=$(scontrol show job "${SLURM_JOB_ID}" | awk -F= '/Command=/{print $2}')
else
	script_path=$(realpath "$0")
fi
script_dir=$(dirname "${script_path}")
parent_dir=$(dirname "${script_dir}")

source ${parent_dir}/parameters.config 

## Inputs
wkdir=${PROJECT_DIR}
andir=${wkdir}/Analyses/Removed_Read_Analysis/kmer_${id}_NoFastP_NoKraken

bam=${wkdir}/Alignments/${id}B_NoFastP_NoKraken_k19/${id}B.cram
bamnew=${wkdir}/Alignments/${id}B_NoFastP_NoKraken_k${kmer}/${id}B.cram

fqd=${wkdir}/MetaData/Comparison_Sequenced_Samples_Cleanup.tsv
ref=${GENRES_DIR}/Resources/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

export TMP_DIR="/weka/users/${USER}/${SLURM_JOB_ID:-$$}"
export MAMBA_ROOT_PREFIX="${TMP_DIR}"
export XDG_CACHE_HOME="${TMP_DIR}"
mkdir -p "${TMP_DIR}"

trap 'rm -rf "$TMP_DIR"' EXIT

if [[ -z "$id" ]]
then 
    echo "No id specified, exiting"
    exit 1
else
	if ! awk -F'\t' -v id="$id" '$1 ~ "^" id {found=1} END{exit !found}' "$fqd"
	then
    	echo "'${id}' not found in $fqd"
    	echo "Exiting..."
    	exit 1
    else
      	echo "'${id}' found in $fqd"  	
	fi
fi

if [[ -z "$kmer" ]]
then 
    echo "No kmers file specified, exiting"
    exit 1
else
	if [[ "$kmer" =~ ^[0-9]+$ ]] && (( $kmer > 19 && $kmer <= 27 ))
	then
    	echo "kmer is between 20 and 27"
	else
    	echo "kmer is NOT between 20 and 27 (k19 is the baseline), exiting..."
    	exit 1
	fi
fi

for f in "${bam}" "${bamnew}"
do
	if [[ ! -s ${f} ]]
	then
		echo "Missing input: ${f}"
		exit 1
	fi
done

if [ -n "${SLURM_JOB_ID:-}" ]
then
	scontrol update jobid=${SLURM_JOBID} jobname="kmer_alignment_mapping_${id}_${kmer}"
	echo "This job is running on node: $SLURMD_NODENAME"
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=${SLURM_CPUS_PER_TASK} 
else
    threads=4
fi

###########################################################################
#########    				ReadName + Mate Approach		     ##########
###########################################################################

mkdir -p ${andir}

unmapped_mates=${andir}/unmapped_mates
lost=${unmapped_mates}_k${kmer}_lost.txt

# One line per unmapped record: read name <tab> mate (1 or 2, from flags 0x40/0x80)
list_unmapped() {
	samtools view -@ ${threads} -T ${ref} -f 4 $1 | awk '{print $1 "\t" (int($2/64)%2 ? 1 : 2)}' | sort -u -T ${MYTMPDIR}
}

if [[ ! -s ${unmapped_mates}_k19.txt ]]
then
	list_unmapped ${bam} > ${unmapped_mates}_k19.txt
fi

if [[ ! -s ${unmapped_mates}_k${kmer}.txt ]]
then
	list_unmapped ${bamnew} > ${unmapped_mates}_k${kmer}.txt
fi

# Newly unmapped mates: unmapped at k${kmer} but not at k19
if [[ ! -e ${lost} ]]
then
	comm -13 ${unmapped_mates}_k19.txt ${unmapped_mates}_k${kmer}.txt > ${lost}
fi
echo "Newly unmapped mates (k19 -> k${kmer}): $(wc -l < ${lost})"

cut -f1 ${lost} | sort -u > ${lost%.txt}.names

# Pull ONLY the newly unmapped mate: -N narrows the scan by read name, awk keeps name + mate matches, primary records only
extract_mate() {   # usage: extract_mate <cram> <out.bam> [extra samtools view flags]
	local in=$1 out=$2; shift 2
	samtools view -@ ${threads} -T ${ref} -h -F 0x900 "$@" -N ${lost%.txt}.names ${in} | awk -F'\t' 'FNR==NR {L[$1 SUBSEP $2]=1; next}; /^@/ {print; next}; {m = (int($2/64)%2 ? 1 : 2)}; (($1 SUBSEP m) in L)' ${lost} - | samtools view -@ ${threads} -b -o ${out} -
}

if [[ -s ${lost} ]]
then
	# these mates as they were aligned at k19
	if [[ ! -s ${lost%.txt}.bam ]]
	then
		extract_mate ${bam} ${lost%.txt}.bam
		samtools index ${lost%.txt}.bam
	fi
	# the same mates as unmapped records at k${kmer} (use these for FASTQ / Kraken2)
	if [[ ! -s ${lost%.txt}.unmapped.bam ]]
	then
		extract_mate ${bamnew} ${lost%.txt}.unmapped.bam -f 4
	fi
else
	echo "No newly unmapped mates for ${id} at k${kmer}, nothing to extract."
fi

exit 0

###########################################################################
#########    				Mapp Kmers to BAM			         ##########
###########################################################################

if [[ ! -s ${andir}/removed_from_kmer_${id}_NoFastP_NoKraken_k19vk${kmer}.kmc_pre ]]
then
	kmc_tools -t${threads} simple "${basekmer}" "${compkmer}" kmers_subtract "${andir}/removed_from_kmer_${id}_NoFastP_NoKraken_k19vk${kmer}"
fi
if [[ ! -s ${andir}/removed_from_kmer_${id}_NoFastP_NoKraken_k19vk${kmer}.txt ]]
then
	kmc_tools -t${threads} transform "${andir}/removed_from_kmer_${id}_NoFastP_NoKraken_k19vk${kmer}" dump "${andir}/removed_from_kmer_${id}_NoFastP_NoKraken_k19vk${kmer}.txt"
fi

cut -f1 "${andir}/removed_from_kmer_${id}_NoFastP_NoKraken_k19vk${kmer}.txt" > "${andir}/removed_from_kmer_${id}_NoFastP_NoKraken_k19vk${kmer}_clean.txt"

if [[ ! -s ${dname}/${bname}_${bamtag} ]]
then
	samtools view -b -T ${ref} "$bam" | katcher --input - \
        --kmers "${andir}/removed_from_kmer_${id}_NoFastP_NoKraken_k19vk${kmer}_clean.txt" \
        --output "${andir}/removed_from_kmer_${id}_NoFastP_NoKraken_k19vk${kmer}_clean.bam" \
        --threads "${threads}"
fi

if [[ ! -s ${dname}/${bname}_${sorttag}.bai ]]
then
	cmd="samtools sort --write-index --threads ${threads} -O BAM -o ${dname}/${bname}_${sorttag}##idx##${dname}/${bname}_${sorttag}.bai ${dname}/${bname}_${bamtag}"
	echo -e "Running samtools sort for $bam for kmers ${kmers}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
fi

rm -f ${dname}/${bname}_temp.txt
rm -f ${dname}/${bname}_${bamtag}

###########################################################################
#########    					Fin						         ##########
###########################################################################

exit 0
