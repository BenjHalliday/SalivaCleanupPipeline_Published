#!/bin/bash 
#SBATCH -c 8     # cores requested
#SBATCH --mem=84G  # memory in Mb
#SBATCH --error     slurm/Kraken2/Kraken2_RL_out_%A_%a.txt
#SBATCH --output    slurm/Kraken2/Kraken2_RL_out_%A_%a.txt
#SBATCH -t 36:00:00  # time requested in hour:minute:second

set -e
set -u
set -o pipefail

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load kraken2
module load pigz

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

if [ ! -z ${1+x} ]
then 
	read=$1
else
	RDARRAYTRUE=($(echo $RDARRAYFALSE | tr '-' " " ))
	read=${RDARRAYTRUE[$(( ${SLURM_ARRAY_TASK_ID} - 1 ))]} 
	scontrol update jobid=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} jobname="Kraken2_"${read}""
fi

if [[ ! -n "$read" ]]
then 
    echo "No Read specified, exiting"
    exit 1
fi

if [ -n "${SLURM_MEM_PER_NODE:-}" ]
then
    mem=$(echo "scale=0; $SLURM_MEM_PER_NODE / 1024" | bc)
else
    mem=84
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=$SLURM_CPUS_PER_TASK
else
    threads=8
fi

Iid=$(awk 'BEGIN {FS = "\t";OFS = "\t" } ($14~/\.*[0-9]$/) { print $0 }' $fqd | awk -F"\t" '$6 == "'${read}'" { print $1 }' )
Did=$(awk 'BEGIN {FS = "\t";OFS = "\t" } ($14~/\.*[0-9]$/) { print $0 }' $fqd | awk -F"\t" '$6 == "'${read}'" { print $2 }' )
lib=$(awk 'BEGIN {FS = "\t";OFS = "\t" } ($14~/\.*[0-9]$/) { print $0 }' $fqd | awk -F"\t" '$6 == "'${read}'" { print $3 }' )
run=$(awk 'BEGIN {FS = "\t";OFS = "\t" } ($14~/\.*[0-9]$/) { print $0 }' $fqd | awk -F"\t" '$6 == "'${read}'" { print $4 }' )
R1=$(awk 'BEGIN {FS = "\t";OFS = "\t" } ($14~/\.*[0-9]$/) { print $0 }' $fqd | awk -F"\t" '$6 == "'${read}'" { print $6 }' )
R2=$(awk 'BEGIN {FS = "\t";OFS = "\t" } ($14~/\.*[0-9]$/) { print $0 }' $fqd | awk -F"\t" '$6 == "'${read}'" { print $7 }' )
readdir=$(awk 'BEGIN {FS = "\t";OFS = "\t" } ($14~/\.*[0-9]$/) { print $0 }' $fqd | awk -F"\t" '$6 == "'${read}'" { print $8 }' )

outdir=${krkdir}/Kraken2/${Iid}/
temp=${SCRATCH}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
prefix=${Iid}_${Did}_${lib}_${run}
readname=$(echo $R1 | sed 's/_R1\.fastq\.gz/_Kraken2/g')

###########################################################################
#########    						Code				         ##########
###########################################################################


if [[ ! -f ${outdir}/${prefix}/${readname}_Report.tsv ]] 
then
	mkdir -p ${outdir}
	mkdir -p ${temp}/${prefix}
	
	if [[ ! -f ${temp}/${R1} ]]
	then
	    echo "Copying "${R1}" to scratch..."
		MAX_TRIES=50
		attempts=0
		until rsync -aP ${readdir}/$R1 ${temp}/
		do
			sleep 2
			if [ $((++attempts)) -gt $MAX_TRIES ]
			then
				(echo "FAILED to sync "$R1"!" 1>&2)
				exit 1
			fi
		done
	    echo "Copying "${R1}" complete..."
	fi
	if [[ ! -f ${temp}/${R2} ]]
	then
	    echo "Copying "${R2}" to scratch..."
		MAX_TRIES=50
		attempts=0
		until rsync -aP ${readdir}/$R2 ${temp}/
		do
			sleep 2
			if [ $((++attempts)) -gt $MAX_TRIES ]
			then
				(echo "FAILED to sync "$R2"!" 1>&2)
				exit 1
			fi
		done
	    echo "Copying "${R2}" complete..."
	fi
	
    if [[ -f ${temp}/${R1} && ${temp}/${R2} ]] 
    then
        echo "Processing "${R1}"..."
        kraken2 --db ${krakendb}/${db} --memory-mapping --threads 8 --output ${temp}/${prefix}/${readname}_Output.tsv --unclassified-out ${temp}/${prefix}/${readname}_Unclassified#.fastq --report ${temp}/${prefix}/${readname}_Report.tsv --paired --gzip-compressed --use-names ${temp}/$R1 ${temp}/$R2
		if [[ $? != 0 ]]; then
        	echo "Kraken failed, exiting..."
        	exit 2  
    	fi

		echo "Kraken2 complete, compressing..."
		pigz -p 6 ${temp}/${prefix}/${readname}_*.fastq
		if [[ $? != 0 ]]; then
        	echo "Compressing fasta files failed, exiting..."
        	exit 3  
    	fi
		
		echo "Compression complete, copying to "${outdir}"..."
		MAX_TRIES=5
		attempts=0
		until rsync -aPr ${temp}/${prefix} ${outdir}/
		do
			sleep 2
			if [ $((++attempts)) -gt $MAX_TRIES ]
			then
				(echo "FAILED to sync "${temp}/${prefix}"!" 1>&2)
				exit 1
			fi
		done
		if [[ $? != 0 ]]; then
        	echo "Moving files to output directory failed, exiting..."
        	exit 4  
    	fi
    	
		echo "Kraken2 pipeline successfully completed for "$R1"!"
		touch ${temp}/.done
    else 
        echo "FastQ "${R1}" is missing, exiting..."
        exit 5     
    fi
    
	if [[ -f ${temp}/.done ]]
	then
		rm -r ${temp}/
	else
		echo "Failure detected, exiting..."
		exit 6
	fi
else
    echo "Kraken2 Screen complete for ${R1}, continuing..."
    sleep 1
fi

exit 0
