#!/bin/bash 
#SBATCH -c 8    # cores requested
#SBATCH --mem=42G  # memory in Mb
#SBATCH -e slurm_validation_variant_filtering_outfile_%A_%a.txt
#SBATCH -o slurm_validation_variant_filtering_outfile_%A_%a.txt
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
	scontrol update jobid=${SLURM_JOBID} jobname="validation_variant_filtering_"${type}"_"${seedlength}""
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
fltdir=${dir}/Validation/VariantFiltering/data_raw

fqd=${dir}/Validation/Meta_Data/Validation_Sequenced_Samples_Info.tsv

annscript=${GENRES_DIR}/Pipelines/VCFAnnotation/Annotate_my_VCF.sh
anndefault=${GENRES_DIR}/Pipelines/Variant_Annotation/defaults_hg38.txt

ref=${GENRES_DIR}/Resources/hg38/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
regions=${GENRES_DIR}/Resources/hg38/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set_nodecoyandnounlocalised.fna.bed

af=0.001

if [ -d "$alndir" ]; then
    echo "Directory $alndir exists."
else
    echo "Directory $alndir does not exist, exiting..."
    exit 1
fi

###########################################################################
#########    				Filter Variants				         ##########
###########################################################################

mkdir -p ${fltdir}

input=${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf
zcat ${input}.gz | head -n 10000 | grep '^#CHROM' | cut -f10- | tr '\t' '\n' > ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_ids.txt

cat ${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_ids.txt | while read -r id
do
	echo "Loop started for ${id}"
	
    input=${vardir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann.vcf
    output=${fltdir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_${id}.vcf
	idind=$(bcftools query -l ${input}.gz | nl -v 0 | grep ${id} | cut -f1 | sed 's/^[ \t]*//;s/[ \t]*$//')
    if [[ ! -s ${output}.gz.tbi ]]; then
        echo "${id} Extracting"
        rm -r -f ${output}*
        bcftools view --regions-file=${regions} -i 'FORMAT/DP['${idind}']>=8 && GQ['${idind}']>=10 && (GT['${idind}']=="alt")' -Oz -o ${output}.gz ${input}.gz
    else
    	echo "${id} extraction already complete, continuing..."
    fi
    if [[ $? != 0 ]]; then
        rm -r -f ${output}*
        echo "Failed to split out variants for "${id}""
        exit 1  
    else
    	if [[ ! -s ${output}.gz.tbi ]]; then
        	bcftools index --tbi ${output}.gz
        fi
    fi

    ### Rarity Filter    ###
    input=${fltdir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_${id}.vcf
    output=${fltdir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_${id}_rare.vcf
    if [[ ! -s ${output}.gz.tbi ]]; then
        echo "${id} Rarity Filter ${af}"
        rm -r -f -f ${output}*
        expression='( gnomADv4_1_AF_grpmax_joint[*] <= '${af}' || gnomADv4_1_AF_grpmax_joint[*] = "." )'
    	command="bcftools view -i '${expression}' ${input}.gz -Oz -o ${output}.gz"
        eval "${command}"
    fi
    if [[ $? != 0 ]]; then
        rm -r -f -f ${output}*
        echo "Failed to filter rare variants for "${id}""
        exit 1  
    else
    	if [[ ! -s ${output}.gz.tbi ]]; then
        	bcftools index --tbi ${output}.gz
        fi
    fi

    ### Coding Filter    ###
    input=${fltdir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_${id}_rare.vcf
    output=${fltdir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_${id}_rare_coding.vcf   
    if [[ ! -s ${output}.gz.tbi ]]; then
        echo "${id} Coding Filter"
        rm -r -f -f ${output}*
	    bcftools filter -i 'INFO/CSQ~"HIGH" || INFO/CSQ~"MODERATE"' -O z -o ${output}.gz ${input}.gz
    fi
    if [[ $? != 0 ]]; then
        rm -r -f -f ${output}*
        echo "Failed to filter coding variants for "${id}""
        exit 1  
    else
    	sleep 5
        vcf_human_readable.sh ${output}.gz    
		wait
    fi

		### Het Filter    ###
		input=${fltdir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_${id}_rare_coding.vcf   
		output=${fltdir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_${id}_rare_coding_het.vcf    
		idind=$(bcftools query -l ${input}.gz | nl -v 0 | grep ${id} | cut -f1 | sed 's/^[ \t]*//;s/[ \t]*$//')
		if [[ ! -s ${output}.gz.tbi ]]; then
		    echo "${id} Het Filter"
		    rm -r -f -f ${output}*
	        bcftools filter -i 'GT['${idind}']="RA"' -O z -o ${output}.gz ${input}.gz
		    sleep 1
		fi
	    if [[ $? != 0 ]]; then
	        rm -r -f -f ${output}*
	        echo "Failed to filter coding variants for "${id}""
	        exit 1  
	    else
			bcftools index --tbi ${output}.gz
	    	sleep 5
	        vcf_human_readable.sh ${output}.gz    
			wait
	    fi	
	    
		### Hom Filter    ###
		input=${fltdir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_${id}_rare_coding.vcf   
		output=${fltdir}/NoFastP_${name}_k${seedlength}_Merged_Norm_ann_${id}_rare_coding_hom.vcf    
		idind=$(bcftools query -l ${input}.gz | nl -v 0 | grep ${id} | cut -f1 | sed 's/^[ \t]*//;s/[ \t]*$//')
		if [[ ! -s ${output}.gz.tbi ]]; then
		    echo "${id} Hom Filter"
		    rm -r -f -f ${output}*
	        bcftools filter -i 'GT['${idind}']="AA"' -O z -o ${output}.gz ${input}.gz
		    sleep 1
		fi
	    if [[ $? != 0 ]]; then
	        rm -r -f -f ${output}*
	        echo "Failed to filter coding variants for "${id}""
	        exit 1  
	    else
			bcftools index --tbi ${output}.gz
	    	sleep 5
	        vcf_human_readable.sh ${output}.gz    
			wait
	    fi	
done

###########################################################################
#########    						Fin					         ##########
###########################################################################

exit 0
