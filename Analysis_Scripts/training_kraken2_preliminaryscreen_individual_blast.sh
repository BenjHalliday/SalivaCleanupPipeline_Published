#!/bin/bash 
#SBATCH -c 12     # cores requested
#SBATCH --mem=200G  # memory in Mb
#SBATCH --error     slurm/Kraken2/Kraken2_IL_BLAST_out_%A_%a.txt
#SBATCH --output    slurm/Kraken2/Kraken2_IL_BLAST_out_%A_%a.txt
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

if [ -n "${SLURM_MEM_PER_NODE:-}" ]
then
    mem=$(echo "scale=0; $SLURM_MEM_PER_NODE / 1024" | bc)
else
    mem=230
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=$SLURM_CPUS_PER_TASK
else
    threads=24
fi

assdir=${krkdir}/MEGAHIT/${id}
blastdir=${krkdir}/BLAST/${id}

###########################################################################
#########    						Code				         ##########
###########################################################################

mkdir -p  ${blastdir}

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
		
    	blastn -task megablast -num_threads $SLURM_CPUS_PER_TASK -db ${blastdb}/nt -query ${blastdir}/megahit_filtered_contigs.fa -out ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt -perc_identity ${percentmatch} -qcov_hsp_perc ${querycover} -max_target_seqs 1 -outfmt "6 sqseqid    sqgi       sqacc      sqaccver   sqlen      ssseqid    ssallseqid ssgi       ssallgi    ssacc      ssaccver   ssallacc   sslen      sqstart    sqend      ssstart    ssend      sqseq      ssseq      sevalue    sbitscore  sscore     slength    spident    snident    smismatch  spositive  sgapopen   sgaps      sppos      sframes    sqframe    ssframe    sbtop      sstaxids   ssscinames sscomnames ssblastnamessskingdomssstitle    ssalltitlesssstrand  sqcovs    sqcovhsp " 
    	echo "BLASTn query complete, continuing..."
    	sed -i '1i Query_Contig\tSubject\tIdentity\tMax_Identity\tMismatches\tGaps\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_Value\tScore' ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt 
    	echo "BLASTn analysis complete for ${id}, finishing..."
    	touch ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.done
    else
    	echo "BLASTn analysis already complete for ${id}, finishing..."    	
    fi
else
    echo "MEGAHIT failed for ${id}, exiting..."
    exit 2
fi

exit 0
