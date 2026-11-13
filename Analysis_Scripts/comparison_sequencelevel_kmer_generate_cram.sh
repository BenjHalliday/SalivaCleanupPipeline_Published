#!/bin/bash 
#SBATCH -c 5    # cores requested
#SBATCH --mem=45G  # memory in Mb
#SBATCH -e slurm_kmer_generate_cram_outfile_%A_%a.txt
#SBATCH -o slurm_kmer_generate_cram_outfile_%A_%a.txt
#SBATCH -t 24:00:00  # time requested in hour:minute:second

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load samtools
module load kmc/3.2.4

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

if [ ! -z ${1+x} ]
then 
	file=$1
fi
if [[ ! -n "$file" ]]
then 
    echo "No CRAM specified, exiting"
    exit 1
else
	if [[ ! -s "$file" && "$file" == *.cram  ]]
	then
		echo "${file} doesn't exist or is not a CRAM file, exiting"
		exit 1
	fi
fi

if [ ! -z ${2+x} ]
then 
	type=$2
fi

andir=$( dirname ${file} )
id=$( basename ${andir} )

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

export TMPDIR=${HOME}/mamba_tmp/$SLURM_JOB_ID
export APPTAINER_TMPDIR=$TMPDIR
export APPTAINER_CACHEDIR=$TMPDIR
export MAMBA_ROOT_PREFIX=$TMPDIR
export XDG_CACHE_HOME=$TMPDIR
mkdir -p $TMPDIR
chmod 770 $TMPDIR

if [ -n "${SLURM_MEM_PER_NODE:-}" ]
then
    mem=$(echo "scale=0; $SLURM_MEM_PER_NODE / 1024" | bc)
else
    mem=64
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=$((SLURM_CPUS_PER_TASK / 2))
else
    threads=8
fi

## Inputs
wkdir=${PROJECT_DIR}
fqd=${wkdir}/MetaData/Comparison_Sequenced_Samples_Cleanup.tsv
if [[ "$type" == "decoy"  ]]
then
	ref=${GENRES_DIR}/Resources/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_plus_MT169739_analysis_set.fna
	echo "${file} is aligned using a decoy, using: ${ref}"
else
	ref=${GENRES_DIR}/Resources/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
	echo "${file} is aligned using the standard reference, using: ${ref}"
fi

minlen=3
maxlen=150
kmerlen=31

###########################################################################
#########    				Make FastQ					         ##########
###########################################################################

mkdir -p ${andir}/KMC

if [[ ! -s ${andir}/KMC/${id}_dataset_kmc_size.txt ]]
then
	rm -f -r ${andir}/KMC/*
fi

if [[ ! -s ${andir}/KMC/${id}.fastq ]] && [[ ! -s ${andir}/KMC/${id}_dataset_kmc_size.txt ]]
then
	#cmd="samtools fastq -F 4 -@ ${threads} --reference ${ref} -o ${andir}/KMC/${id}.fastq ${file}"
	cmd="samtools view --reference ${ref} -@ ${threads} -h ${file} $(printf "chr%s " {1..22} X Y) | samtools fastq -@ ${threads} -F 4 -o ${andir}/KMC/${id}.fastq"
	echo -e "Running SAMTools to produce a FASTQ for ${id}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
fi

###########################################################################
#########    				Run KMC						         ##########
###########################################################################

if [[ ! -s ${andir}/KMC/${id}_dataset_kmc.kmc_pre ]]
then
	mkdir -p ${andir}/KMC/${id}_kmc_temp
	cmd="kmc -t${threads} -m${mem} -k${kmerlen} -ci${minlen} -cs${maxlen} -fq ${andir}/KMC/${id}.fastq ${andir}/KMC/${id}_dataset_kmc ${andir}/KMC/${id}_kmc_temp"
	echo -e "Running KMC creating 31-mer database for ${id}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
else
	echo "KMC dataset already generated for ${id}"
fi

if [[ ! -s ${andir}/KMC/${id}_dataset_kmc.txt ]]
then
	cmd="kmc_tools -t${threads} transform "${andir}/KMC/${id}_dataset_kmc" dump "${andir}/KMC/${id}_dataset_kmc.txt""
	echo -e "Running KMC dumping 31-mer database for ${id}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
else
	echo "KMC dump already generated for ${id}"
fi

if [[ ! -s ${andir}/KMC/${id}_dataset_kmc_size.txt ]]
then
	awk 'END { print NR }' "${andir}/KMC/${id}_dataset_kmc.txt" > "${andir}/KMC/${id}_dataset_kmc_size.txt"
else
	echo "KMC dump already linecounted for ${id}"
fi

rm -f -r ${andir}/KMC/${id}_kmc_temp/
rm -f ${andir}/KMC/${id}.fastq

rm -r $TMPDIR
sleep $((1 + $RANDOM % 5))

exit 0
