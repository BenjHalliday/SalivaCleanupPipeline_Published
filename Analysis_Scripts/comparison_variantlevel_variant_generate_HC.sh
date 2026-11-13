#!/bin/bash 
#SBATCH -c 1    # cores requested
#SBATCH --mem=1G  # memory in Mb
#SBATCH -e slurm_variant_generate_HC_outfile_%A_%a.txt
#SBATCH -o slurm_variant_generate_HC_outfile_%A_%a.txt
#SBATCH -t 24:00:00  # time requested in hour:minute:second

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load bcftools 

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

# Ensure all variables are provided
if [[ -z ${1+x} || -z ${2+x} || -z ${3+x} ]]; then
  echo "Usage: $0 <Sample_ID> <Group> <SeedLength>"
fi

# Check if id is provided and not empty
if [ ! -z ${1+x} ]; then
	id=$1
else
	echo "Error: ID must be provided"
	exit 1
fi

# Check if type is one of the allowed values (Decoy, Kraken, KrakenTrained, SeedLength, Minimap2)
if [ ! -z ${2+x} ]; then
	type=$2
	if [[ ! "$type" =~ ^(Decoy|Kraken|KrakenTrained|SeedLength|Minimap2)$ ]]; then
		echo "Error: Type must be one of 'Decoy', 'Kraken', 'KrakenTrained', 'SeedLength', or 'Minimap2'"
		exit 1
	fi
else
 	echo "Error: Type must be provided (one of 'Decoy', 'Kraken', 'KrakenTrained', 'SeedLength', or 'Minimap2')"
  	exit 1
fi

if [[ "$type" == "Decoy" ]]; then
	name="Decoy"
elif [[ "$type" == "SeedLength" ]]; then
	name="NoKraken"
elif [[ "$type" == "Kraken" ]]; then
	name="Krakentaxid"
elif [[ "$type" == "KrakenTrained" ]]; then
	name="Kraken_trained_taxid"
fi

# Check if seedlength is a valid number
if [ ! -z ${3+x} ]; then
	seedlength=$3
	if ! [[ "$seedlength" =~ ^[0-9]+$ ]]; then
		echo "Error: Seedlength must be a valid number"
		exit 1
	fi
else
	echo "Error: Seedlength must be provided"
	exit 1
fi

# Print out the valid values for confirmation
echo "ID: $id"
echo "Type: $type"
echo "SeedLength: $seedlength"

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
    mem=8
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=$((SLURM_CPUS_PER_TASK / 2))
else
    threads=2
fi

bloodid=$(echo $id | sed 's/S/B/g')

## Inputs
dir=${PROJECT_DIR}
alndir=${dir}/Alignments/${id}_NoFastP_${name}_k${seedlength}
balndir=${dir}/Alignments/${bloodid}_NoFastP_${name}_k${seedlength}
wkdir=${dir}/Analyses/Variant_Comparison_Analysis

regions=${dir}/Analyses/SoftClip_Comparison_Analysis/data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Masked_Saliva.tsv #${dir}/Analyses/SoftClip_Comparison_Analysis/data_clean/Combined_DoubleSoftClipped_Repeat_Filtered_Mappability_Combined_Masked_Saliva.tsv
fqd=${wkdir}/MetaData/Comparison_Sequenced_Samples_Cleanup.tsv

annscript=${GENRES_DIR}/Pipelines/VCFAnnotation/Annotate_my_VCF.sh
anndefault=${GENRES_DIR}/Pipelines/VCFAnnotation/defaults_hg38.txt

freq=4 #5
buffer=20 #5

if [ -d "$alndir" ]; then
    echo "Directory $alndir exists."
else
    echo "Directory $alndir does not exist, exiting..."
    exit 1
fi

###########################################################################
#########    				Pull Variants				         ##########
###########################################################################

mkdir -p ${wkdir}/data_raw


if [[ ! -s ${wkdir}/data_raw/${id}_NoFastP_NoKraken_k19_SoftClip_FreqTable.tsv ]]
then
	awk -v id="$id" -v type="SeedLength" -v seedlength="19" -v freq="$freq" '
	BEGIN {
	    FS = OFS = "\t";  # Set field separator to tab
	}
	NR == 1 {
    	# Find the column indices for Sample_ID, SeedLength, Group, and Frequency
    	for (i = 1; i <= NF; i++) {
        	if ($i == "Sample_ID") col1 = i;
        	if ($i == "SeedLength") col2 = i;
        	if ($i == "Group") col3 = i;
        	if ($i == "Frequency") col4 = i;
    	}
    	print $0;  # Print the header
	}
	NR > 1 {
    	# Filter rows based on input values
    	if ($col1 == id && $col3 == type && $col2 == seedlength && $col4 >= freq) {
        	print $0;
    	}
	}' "$regions" > "${wkdir}/data_raw/${id}_NoFastP_NoKraken_k19_SoftClip_FreqTable.tsv"
fi

if [[ ! -s ${wkdir}/data_raw/${id}_NoFastP_NoKraken_k19_SoftClip_FreqTable.bed ]]
then
	awk -v buffer="$buffer" 'BEGIN { FS = OFS = "\t" } NR > 1 { print $2,$3-buffer,$4+buffer}' ${wkdir}/data_raw/${id}_NoFastP_NoKraken_k19_SoftClip_FreqTable.tsv > ${wkdir}/data_raw/${id}_NoFastP_NoKraken_k19_SoftClip_FreqTable.bed
fi

##############
##	Saliva	##
##############

if [[ ! -s ${wkdir}/data_raw/${id}_NoFastP_${name}_k${seedlength}_SoftClip.vcf.gz ]]
then
	bcftools view --output ${wkdir}/data_raw/${id}_NoFastP_${name}_k${seedlength}_SoftClip.vcf.gz --output-type z --regions-file ${wkdir}/data_raw/${id}_NoFastP_NoKraken_k19_SoftClip_FreqTable.bed ${alndir}/${id}_HC.vcf.gz
fi

if [[ ! -s ${wkdir}/data_raw/${id}_NoFastP_${name}_k${seedlength}_SoftClip.vcf.gz.tbi ]]
then
	bcftools index -t ${wkdir}/data_raw/${id}_NoFastP_${name}_k${seedlength}_SoftClip.vcf.gz
fi

${script_dir}/vcf_human_readable.sh ${wkdir}/data_raw/${id}_NoFastP_${name}_k${seedlength}_SoftClip.vcf.gz

##############
##	Blood	##
##############

if [[ ! -s ${wkdir}/data_raw/${bloodid}_NoFastP_${name}_k${seedlength}_SoftClip.vcf.gz ]]
then
	bcftools view --output ${wkdir}/data_raw/${bloodid}_NoFastP_${name}_k${seedlength}_SoftClip.vcf.gz --output-type z --regions-file ${wkdir}/data_raw/${id}_NoFastP_NoKraken_k19_SoftClip_FreqTable.bed ${balndir}/${bloodid}_HC.vcf.gz
fi

if [[ ! -s ${wkdir}/data_raw/${bloodid}_NoFastP_${name}_k${seedlength}_SoftClip.vcf.gz.tbi ]]
then
	bcftools index -t ${wkdir}/data_raw/${bloodid}_NoFastP_${name}_k${seedlength}_SoftClip.vcf.gz
fi

${script_dir}/vcf_human_readable.sh ${wkdir}/data_raw/${bloodid}_NoFastP_${name}_k${seedlength}_SoftClip.vcf.gz

###########################################################################
#########    						Fin					         ##########
###########################################################################

exit 0
