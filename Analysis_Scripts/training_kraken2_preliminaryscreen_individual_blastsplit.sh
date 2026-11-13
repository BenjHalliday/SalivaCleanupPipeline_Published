#!/bin/bash
#SBATCH -c 1     # cores requested
#SBATCH --mem=1G  # memory in Mb
#SBATCH --error     slurm/Kraken2/Kraken2_IL_BLASTsplit_out_%A_%a.txt
#SBATCH --output    slurm/Kraken2/Kraken2_IL_BLASTsplit_out_%A_%a.txt
#SBATCH -t 128:00:00  # time requested in hour:minute:second

set -e
set -u
set -o pipefail

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load BLAST+/2.13.0-gompi-2021a

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

if [ ! -z ${1+x} ]
then
	id=$1
else
	IDARRAYTRUE=($(echo $IDARRAYFALSE | tr '-' " " ))
	id=${IDARRAYTRUE[$(( ${SLURM_ARRAY_TASK_ID} - 1 ))]}
	scontrol update jobid=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} jobname="Kraken2_BLAST_"${id}""
fi

if [[ ! -n "$id" ]]
then
    echo "No ID specified, exiting"
    exit 1
fi

if [ -n "${SLURM_JOB_ID:-}" ]
then
	script_path=$(scontrol show job "${SLURM_JOB_ID}" | awk -F= '/Command=/{print $2}')
else
	script_path=$(realpath "$0")
fi
script_dir=$(dirname "${script_path}")
parent_dir=$(dirname "${script_dir}")

if [ -n "${SLURM_MEM_PER_NODE:-}" ]
then
    mem=$(echo "scale=0; $SLURM_MEM_PER_NODE / 1024" | bc)
else
    mem=1
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=$SLURM_CPUS_PER_TASK
else
    threads=1
fi

assdir=${krkdir}/MEGAHIT/${id}
blastdir=${krkdir}/BLAST/${id}
chunkdir=${blastdir}/chunks
chunkthreads=20
chunksize=1000 #500 entries of two lines

###########################################################################
#########    						Code				         ##########
###########################################################################

mkdir -p ${chunkdir}

if [[ -f ${assdir}/done ]]
then
    if [[ ! -f ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.done ]]
    then
    	cat ${assdir}/final.contigs.fa | sed -e 's/multi=//g' | sed -e 's/len=//g' | sed -e 's/^>//g' | awk 'BEGIN {FS=" ";OFS="\t"} $3 >= '${depthmin}' && $4 >= '${lengthmin}' { print $1,$3,$4 }'  > ${blastdir}/megahit_filtered_contignames.tsv
    	grep -w --no-group-separator -A1 ${assdir}/final.contigs.fa -f <(cut -f1  ${blastdir}/megahit_filtered_contignames.tsv) > ${blastdir}/megahit_filtered_contigs.fa

    	if [[ ! -s ${blastdir}/megahit_filtered_contigs.fa ]]
    	then
    	    echo "No MEGAHIT contigs were left after filtering for ${id}, finishing..."
    	    touch ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt
    		touch ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.done
    	    exit 0
    	fi

    	# Split the FASTA file into chunks
    	split -l 1000 -d --additional-suffix=.fasta ${blastdir}/megahit_filtered_contigs.fa ${chunkdir}/chunk_

    	# Submit BLAST jobs for each chunk
    	chunk_files=(${chunkdir}/chunk_*.fasta)
    	chunk_count=${#chunk_files[@]}
    	job_ids=()

    	for ((i=0; i<chunk_count; i++))
    	do
    		chunk=${chunk_files[$i]}
    		chunk_name=$(basename $chunk .fasta)
    		output=${chunkdir}/${chunk_name}_blastn-nt_hits.txt

    		# Submit a separate SLURM job for each chunk and store the job ID
    		job_id=$(sbatch --export=chunk=${chunk},output=${output},blastdb=${blastdb},percentmatch=${percentmatch},querycover=${querycover},chunkthreads=${chunkthreads} --parsable ${script_dir}/training_kraken2_preliminaryscreen_individual_blastsplit_chunk.sh)
    		job_ids+=($job_id)
    	done

    	# Wait for all chunk jobs to complete
    	for job_id in "${job_ids[@]}"
    	do
    		scontrol wait job=$job_id
    		job_state=$(scontrol show job $job_id | grep 'JobState=' | awk -F= '{print $2}')
    		if [[ $job_state != "COMPLETED" ]]
    		then
    			echo "Error: Job $job_id failed."
    			exit 1
    		fi
    	done

    	# Concatenate all chunk outputs into the final output file
    	head -n1 ${chunkdir}/chunk_00_blastn-nt_hits.txt > ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt
    	tail -n +2 -q ${chunkdir}/chunk_*_blastn-nt_hits.txt >> ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt
    	echo "BLASTn analysis complete for ${id}, finishing..."
    	touch ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.done
    	rm -r -f ${chunkdir}
    else
    	echo "BLASTn analysis already complete for ${id}, finishing..."
    fi
else
    echo "MEGAHIT failed for ${id}, exiting..."
    exit 2
fi

exit 0