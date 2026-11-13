#!/bin/bash 
#SBATCH -c 8    # cores requested
#SBATCH --mem=42G  # memory in Mb
#SBATCH -e slurm_validation_variant_annotation_outfile_%A_%a.txt
#SBATCH -o slurm_validation_variant_annotation_outfile_%A_%a.txt
#SBATCH -t 64:00:00  # time requested in hour:minute:second

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load bcftools

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

# Ensure all variables are provided
if [[ -z ${1+x} || -z ${2+x}  ]]; then
  echo "Usage: $0 <Group> <SeedLength>"
fi

# Check if type is one of the allowed values (Decoy, Kraken, KrakenTrained, SeedLength, Minimap2)
if [ ! -z ${1+x} ]; then
	type=$1
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
if [ ! -z ${2+x} ]; then
	seedlength=$2
	if ! [[ "$seedlength" =~ ^[0-9]+$ ]]; then
		echo "Error: Seedlength must be a valid number"
		exit 1
	fi
else
	echo "Error: Seedlength must be provided"
	exit 1
fi

# Print out the valid values for confirmation
echo "Type: $type"
echo "SeedLength: $seedlength"

if [ -n "${SLURM_JOB_ID:-}" ]
then
	scontrol update jobid=${SLURM_JOBID} jobname="validation_variant_annotation_"${type}"_"${seedlength}""
	echo "This job is running on node: $SLURMD_NODENAME"
fi

if [ -n "${SLURM_MEM_PER_NODE:-}" ]
then
    mem=$(echo "scale=0; $SLURM_MEM_PER_NODE / 1024" | bc)
else
    mem=42
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=$((SLURM_CPUS_PER_TASK / 2))
else
    threads=8
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

## Inputs
dir=${PROJECT_DIR}
alndir=${dir}/Validation/Alignments/CMM_Alignments
vardir=${dir}/Validation/VariantCalls

fqd=${dir}/Validation/Meta_Data/Validation_Sequenced_Samples_Info.tsv

annscript=${GENRES_DIR}/Pipelines/VCFAnnotation/Annotate_my_VCF.sh
anndefault=${GENRES_DIR}/Pipelines/Variant_Annotation/defaults_hg38.txt

ref=${GENRES_DIR}/Resources/hg38/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
regions=${GENRES_DIR}/Resources/hg38/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_nodecoyandnounlocalised.fna.bed

if [ -d "$alndir" ]; then
    echo "Directory $alndir exists."
else
    echo "Directory $alndir does not exist, exiting..."
    exit 1
fi

###########################################################################
#########    				Merge Variants				         ##########
###########################################################################

mkdir -p ${vardir}

rm -f "${vardir}/NoFastP_${name}_k${seedlength}_Files.txt"; touch "${vardir}/NoFastP_${name}_k${seedlength}_Files.txt"
cat $fqd | sed 1d | cut -f1 | sort -u | while read -r id
do
	if [[ -f "${alndir}/${id}_NoFastP_${name}_k${seedlength}/${id}_HC.vcf.gz" ]]
	then
        echo "${alndir}/${id}_NoFastP_${name}_k${seedlength}/${id}_HC.vcf.gz" >> "${vardir}/NoFastP_${name}_k${seedlength}_Files.txt"
    else
    	echo "Could not find "${alndir}/${id}_NoFastP_${name}_k${seedlength}/${id}_HC.vcf.gz""
    fi
done

if [[ ! -s ${vardir}/NoFastP_${name}_k${seedlength}_Merged.vcf.gz.csi ]]
then
	bcftools merge --threads ${threads} --merge none --file-list "${vardir}/NoFastP_${name}_k${seedlength}_Files.txt" -O z -o ${vardir}/NoFastP_${name}_k${seedlength}_Merged.vcf.gz
	bcftools index ${vardir}/NoFastP_${name}_k${seedlength}_Merged.vcf.gz
fi

if [[ ! -s ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm.vcf.gz.csi ]]
then
	bcftools view ${vardir}/NoFastP_${name}_k${seedlength}_Merged.vcf.gz --regions-file $regions | bcftools norm -m - --multi-overlaps . --fasta-ref ${ref} --threads ${threads} -Oz -o ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm.vcf.gz
	bcftools index ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm.vcf.gz
fi

if [[ ! -s ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz && ! -s ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_Annotated_1/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz ]]
then
	${annscript} -i ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm.vcf.gz -d ${anndefault}
fi

max_retries=20
wait_time=600

attempt=0
while [[ ! -s ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz && ! -s ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_Annotated_1/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz && $attempt -lt $max_retries ]]; do
    echo "Attempt $((attempt+1))/$max_retries: Waiting for $target..."
    sleep "$wait_time"
    ((attempt++))
done

if [[ -f "${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_Annotated_1/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz" ]]; then
    echo "${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_Annotated_1/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz found! Continuing..."
    mv ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_Annotated_1/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz
    mv ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_Annotated_1/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz.tbi ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz.tbi
	rm -r ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_Annotated_1
else
	if [[ ! -s ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf.gz ]]
	then
    	echo "ERROR: $target not found after $max_retries attempts. Exiting."
    	exit 1
	fi
fi

###########################################################################
#########    						Fin					         ##########
###########################################################################

exit 0
