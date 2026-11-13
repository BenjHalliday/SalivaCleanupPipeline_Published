#!/bin/bash 
#SBATCH -c 6     # cores requested
#SBATCH --mem=24G  # memory in Mb
#SBATCH --error     slurm/Kraken2/Kraken2_IL_MEGAHIT_out_%A_%a.txt
#SBATCH --output    slurm/Kraken2/Kraken2_IL_MEGAHIT_out_%A_%a.txt
#SBATCH -t 24:00:00  # time requested in hour:minute:second
 
set -e
set -u
set -o pipefail

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load megahit/1.2.9

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

if [ ! -z ${1+x} ]
then 
	id=$1
else
	IDARRAYTRUE=($(echo $IDARRAYFALSE | tr '-' " " ))
	id=${IDARRAYTRUE[$(( ${SLURM_ARRAY_TASK_ID} - 1 ))]} 
	scontrol update jobid=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} jobname="Kraken2_MEGAHIT_"${id}""
fi

if [[ ! -n "$id" ]]
then 
    echo "No ID specified, exiting"
    exit 1
fi

if [ -n "${SLURM_MEM_PER_NODE:-}" ]
then
    mem=$(echo "scale=0; $SLURM_MEM_PER_NODE / 1024" | bc)
else
    mem=24
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=$SLURM_CPUS_PER_TASK
else
    threads=6
fi

indir=${krkdir}/Kraken2/${id}/
outdir=${krkdir}/MEGAHIT/${id}
temp=${SCRATCH}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

###########################################################################
#########    						Code				         ##########
###########################################################################

mkdir -p ${krkdir}/MEGAHIT/

if [[ -f ${indir}/Missing.txt ]]
then
	rm ${indir}/Missing.txt
fi

awk 'BEGIN{FS="\t";OFS="\t"}$10 == "Genomic" && $5 == "Illumina" { print $0 }' ${fqd} | \
awk 'BEGIN{FS="\t";OFS="\t"}($14~/\.*[0-9]$/) { print $0 }' | awk 'BEGIN{FS="\t";OFS="\t"}$1 == "'${id}'" { print $0 }' | \
cut -f6 | sort -u | while read -r read
do
	Iid=$(awk -F"\t" '$6 == "'${read}'" { print $1 }' ${fqd})
	Did=$(awk -F"\t" '$6 == "'${read}'" { print $2 }' ${fqd})
	lib=$(awk -F"\t" '$6 == "'${read}'" { print $3 }' ${fqd})
	run=$(awk -F"\t" '$6 == "'${read}'" { print $4 }' ${fqd})
	R1=$(awk -F"\t" '$6 == "'${read}'" { print $6 }' ${fqd})
	R2=$(awk -F"\t" '$6 == "'${read}'" { print $7 }' ${fqd})

	readname=$(echo $R1 | sed 's/_R1\.fastq\.gz/_Kraken2/g')
	prefix=${Iid}_${Did}_${lib}_${run}

	Uread1=${indir}/${prefix}/${readname}_Unclassified_1.fastq.gz
	Uread2=${indir}/${prefix}/${readname}_Unclassified_2.fastq.gz
	
	if [[ -z "${fastq1:-}" ]]
    then    
        echo
        if [[ -f ${Uread1} && ${Uread2} ]]
        then
            fastq1="${Uread1}"
            fastq2="${Uread2}"
        else
            fastq1="${Uread1}"
            fastq2="${Uread2}"
    		echo "${fastq1}" > ${indir}/Missing.txt
            echo "${Uread1} does not exist, exiting..."
            exit 1
        fi
    else
        if [[ -f ${Uread1} && ${Uread2} ]]
        then
            fastq1=${fastq1}" ${Uread1}"
            fastq2=${fastq2}" ${Uread2}"
        else
            fastq1="${Uread1}"
            fastq2="${Uread2}"
    		echo "${fastq1}" > ${indir}/Missing.txt
            echo "${Uread1} does not exist, exiting..."
            exit 1
        fi
    fi
    echo "${fastq1}" > ${indir}/Inputs1.txt
    echo "${fastq2}" > ${indir}/Inputs2.txt
done

if [[ -s ${indir}/Missing.txt ]]
then
	echo "Missing FASTQ file detected, exiting..."
	rm ${indir}/Inputs1.txt
	rm ${indir}/Inputs2.txt
	exit 2
else
	fastq1=`cat ${indir}/Inputs1.txt`
	fastq2=`cat ${indir}/Inputs2.txt`
fi

mkdir -p ${temp}

if [[ ! -s ${indir}/${id}_Kraken2_Unclassified_1.fastq.gz ]]
then
	cat ${fastq1} > ${indir}/${id}_Kraken2_Unclassified_1.fastq.gz
else
	echo "Individual level combination already complete for R1, finishing..."
fi

if [[ ! -s ${indir}/${id}_Kraken2_Unclassified_2.fastq.gz ]]
then
	cat ${fastq2} > ${indir}/${id}_Kraken2_Unclassified_2.fastq.gz
else
	echo "Individual level combination already complete for R2, finishing..."
fi

### Remove NNNN sequences

if [[ ! -f ${indir}/${id}_Kraken2_Unclassified_2_Trimmed.fastq.gz ]] 
then
	module purge
	module load cutadapt/2.10
	cutadapt -j $SLURM_CPUS_PER_TASK --trim-n -o ${indir}/${id}_Kraken2_Unclassified_1_Trimmed.fastq.gz \
		-p ${indir}/${id}_Kraken2_Unclassified_2_Trimmed.fastq.gz ${indir}/${id}_Kraken2_Unclassified_1.fastq.gz ${indir}/${id}_Kraken2_Unclassified_2.fastq.gz
	if [[ $? != 0 ]]; then
      	echo "Trimming output fastq files failed, exiting..."
      	rm ${indir}/${id}_Kraken2_Unclassified_*_Trimmed.fastq.gz
      	exit 7  
    fi
else
	echo "CutAdapt trimming complete for ${indir}/${id}_Kraken2_Unclassified_[1-2].fastq.gz, finishing..."
	sleep 1
fi

if [[ ! -s ${outdir}/final.contigs.fa ]]
then
	module purge
	module load megahit/1.2.9
	if [[ -d ${outdir}/ ]]
	then
		rm -r ${outdir}
	fi
	megahit -o ${outdir} --tmp-dir ${temp} --k-list 59,79,99,119,141 --no-mercy \
		--bubble-level 0 --prune-level 3 --prune-depth 5 --cleaning-rounds 3 --no-local --kmin-1pass \
		--continue -t $SLURM_CPUS_PER_TASK  --mem-flag 0 --min-contig-len 500 --min-count 30 \
		-1 ${indir}/${id}_Kraken2_Unclassified_1_Trimmed.fastq.gz -2 ${indir}/${id}_Kraken2_Unclassified_2_Trimmed.fastq.gz
else
	echo "Individual level MEGAHIT already complete, finishing..."
fi

rm -r ${temp}
rm ${indir}/Inputs1.txt
rm ${indir}/Inputs2.txt

exit 0
