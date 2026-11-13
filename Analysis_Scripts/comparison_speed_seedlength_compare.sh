#!/bin/bash 
#SBATCH -c 2    # cores requested
#SBATCH --mem=4G  # memory in Mb
#SBATCH -e slurm_speed_compare_outfile_%A_%a.txt
#SBATCH -o slurm_speed_compare_outfile_%A_%a.txt
#SBATCH -t 24:00:00  # time requested in hour:minute:second

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load r/4.4.3-k6zf4xz

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

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
    mem=4
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=${SLURM_CPUS_PER_TASK}
else
    threads=2
fi

## Inputs
wkdir=${PROJECT_DIR}
spdir=${wkdir}/Analyses/Compute_Comparison

fqd=${wkdir}/MetaData/Comparison_Sequenced_Samples_Cleanup.tsv

###########################################################################
#########    				Run untar					         ##########
###########################################################################

mkdir -p ${spdir}
mkdir -p ${spdir}/data_raw

ls ${wkdir}/Alignments/V*/V*_*_k*_SACCT.txt | while read -r stat
do
	id=$( basename ${stat} | sed 's/_SACCT\.txt//g' )
	awk '{$1=$1; gsub(/[ ]+/, "\t"); print}' ${stat} | sed 2d | \
		awk -F'\t' '!($2 ~ /batch/ || $2 ~ /extern/) { print }' OFS='\t' | \
		awk -F'\t' '!($1 ~ /\./ || $1 ~ /\+/) { print }' OFS='\t' > ${spdir}/data_raw/${id}_sacct_output.tsv
done

cmd="Rscript ${script_dir}/analyses_speed_seedlength.R"
echo $cmd
eval $cmd

exit 0
