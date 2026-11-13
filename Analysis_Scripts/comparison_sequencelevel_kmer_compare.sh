#!/bin/bash 
#SBATCH -c 5    # cores requested
#SBATCH --mem=30G  # memory in Mb
#SBATCH -e slurm_kmer_compare_outfile_%A_%a.txt
#SBATCH -o slurm_kmer_compare_outfile_%A_%a.txt
#SBATCH -t 24:00:00  # time requested in hour:minute:second

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load r/4.4.3-k6zf4xz
module load kmc/3.2.4

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

if [ ! -z ${1+x} ]
then 
	root=$1
fi
if [[ ! -n "$root" ]]
then 
    echo "No root specified, exiting"
    exit 1
fi

if [ ! -z ${2+x} ]
then 
	kmc1=$2
fi
if [[ ! -n "$kmc1" ]]
then 
    echo "No KMC DB 1 specified, exiting"
    exit 1
fi

if [ ! -z ${3+x} ]
then 
	kmc2=$3
fi
if [[ ! -n "$kmc2" ]]
then 
    echo "No KMC DB 2 specified, exiting"
    exit 1
fi

if [ ! -z ${4+x} ]
then 
	tag=$4
fi
if [[ ! -n "$tag" ]]
then 
    echo "No tag specified, exiting"
    exit 1
fi

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
    mem=140
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=${SLURM_CPUS_PER_TASK}
else
    threads=16
fi

## Inputs
wkdir=${PROJECT_DIR}
andir=${wkdir}/Analyses/Kmer_Comparison_Analysis/data_raw/${root}/kmer_${tag}/data_raw
fqd=${wkdir}/MetaData/Comparison_Sequenced_Samples_Cleanup.tsv

###########################################################################
#########    				Run KMC						         ##########
###########################################################################

mkdir -p ${andir}

kmc1_dir=$(echo $kmc1 | sed 's/\.kmc.*//g')
kmc2_dir=$(echo $kmc2 | sed 's/\.kmc.*//g')

kmc1_name=$(basename $kmc1_dir)
kmc2_name=$(basename $kmc2_dir)

if [[ ! -s ${andir}/common_kmers.kmc_pre ]]
then
	kmc_tools -t${threads} simple "${kmc1_dir}" "${kmc2_dir}" intersect "${andir}/common_kmers"
fi

if [[ ! -s ${andir}/unique_to_${kmc1_name}.kmc_pre ]]
then
	kmc_tools -t${threads} simple "${kmc1_dir}" "${kmc2_dir}" kmers_subtract "${andir}/unique_to_${kmc1_name}"
fi

if [[ ! -s ${andir}/unique_to_${kmc2_name}.kmc_pre ]]
then
	kmc_tools -t${threads} simple "${kmc2_dir}" "${kmc1_dir}" kmers_subtract "${andir}/unique_to_${kmc2_name}"
fi

if [[ ! -s ${andir}/common_kmers.txt ]]
then
	kmc_tools -t${threads} transform "${andir}/common_kmers" dump "${andir}/common_kmers.txt"
fi

if [[ ! -s "${andir}/unique_to_${kmc1_name}.txt" ]]
then
	kmc_tools -t${threads} transform "${andir}/unique_to_${kmc1_name}" dump "${andir}/unique_to_${kmc1_name}.txt"
fi

if [[ ! -s "${andir}/unique_to_${kmc2_name}.txt" ]]
then
	kmc_tools -t${threads} transform "${andir}/unique_to_${kmc2_name}" dump "${andir}/unique_to_${kmc2_name}.txt"
fi

if [[ ! -s "${andir}/common_kmers_size.txt" ]]
then
	awk 'END { print NR }' "${andir}/common_kmers.txt" > "${andir}/common_kmers_size.txt"
fi

if [[ ! -s "${andir}/unique_to_${kmc1_name}_size.txt" ]]
then
	awk 'END { print NR }' "${andir}/unique_to_${kmc1_name}.txt" > "${andir}/unique_to_${kmc1_name}_size.txt"
fi

if [[ ! -s "${andir}/unique_to_${kmc2_name}_size.txt" ]]
then
	awk 'END { print NR }' "${andir}/unique_to_${kmc2_name}.txt" > "${andir}/unique_to_${kmc2_name}_size.txt"
fi

cmd="Rscript ${script_dir}/comparison_sequencelevel_kmer_compare_plot.R ${root} kmer_${tag} ${kmc1_name} ${kmc2_name}"
echo $cmd
eval $cmd

rm -r $TMPDIR
sleep $((1 + $RANDOM % 5))

exit 0
