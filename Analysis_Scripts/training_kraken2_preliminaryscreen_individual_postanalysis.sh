#!/bin/bash 
#SBATCH -c 1     # cores requested
#SBATCH --mem=2G  # memory in Mb
#SBATCH --error     slurm/Kraken2/Kraken2_IL_Post_out_%A_%a.txt
#SBATCH --output    slurm/Kraken2/Kraken2_IL_Post_out_%A_%a.txt
#SBATCH -t 48:00:00  # time requested in hour:minute:second

set -e
set -u
set -o pipefail

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge

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

blastdir=${krkdir}/BLAST/${id}
outdir=${krkdir}/DatabaseCreation

###########################################################################
#########    						Code				         ##########
###########################################################################

test_apps="esearch"
for cmd in $test_apps; do
    type &>/dev/null $cmd && echo "edirect installed continuing.." ||
        { echo >&2 "edirect not installed, downloading..."; sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"; export PATH=${PATH}:${HOME}/edirect; }
done

if [[ ! -f ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_annotated.done ]]
then

	if [[ ! -s /${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt ]]
	then
		echo "No MEGAHIT contigs were left after filtering for ${id}, finishing..."
    	touch ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_annotated.tsv
    	touch ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_annotated.done
    	touch ${blastdir}/megahit_filtered_contigs_blastn-nt_nohits.fa
    	exit 0
	fi
		### Non-hits ### 
	grep --no-group-separator -A1 ${blastdir}/megahit_filtered_contigs.fa -w -f <( grep -v -w -f <( cat ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt | cut -f1 | sort -u ) ${blastdir}/megahit_filtered_contignames.tsv | cut -f1) > ${blastdir}/megahit_filtered_contigs_blastn-nt_nohits.fa
	
		###  Retrieve Hit Fasta and format for Kraken DataBase addition ### 
	echo "Retrieving and annoting FASTA for ${id}, processing..."
	cat ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt | sed 1d | cut -f2 | sort -u 
	
	grep -f <(cat ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt | sed 1d | cut -f2 | sort -u) ${outdir}/Combined_filtered_contigs_blastn-nt_hits_acc-annotations.tsv > ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_acc-annotations.tsv
	
		### Annotate Hits ### 
	echo "BLAST annotation for ${id}, processing..."
	cat ${blastdir}/megahit_filtered_contigs.fa | grep '^>' | cut -d' ' -f 1,3,4 | sed 's/^>//g' | sed 's/multi=//g' |  sed 's/len=//g' | tr ' ' '\t' > ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_contig-annotations.tsv
	paste --delimiters='\t' <( head -n1 ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt ) <(echo -e "Query_Depth\tQuery_Length\tSubject_Taxon_ID\tSubject_Organism\tSubject_Fragment_Name\tQuery_Coverage") > ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_ann.tsv
	
	join -t $'\t' -1 2 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,2.2,2.3,2.4 \
	 <( join -t $'\t' -1 1 -2 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.2,2.3 <( sed 1d ${blastdir}/megahit_filtered_contigs_blastn-nt_hits.txt | sort -k1,1 ) <( sort -k1,1 ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_contig-annotations.tsv ) | sort -k 2,2 ) \
	 <( sort -k1,1 ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_acc-annotations.tsv ) | awk 'BEGIN {FS=" ";OFS="\t"} {print $0, ($8-$7+1) / $14 }' | tr -d $'\r' >> ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_ann.tsv
	
	awk 'BEGIN {FS = "\t";OFS = "\t" } { print $2,$15,$16,$17, $1,$13,$14,$18,$3,$9,$10,$7,$8, $4,$5,$6,$11,$12 }' ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_ann.tsv > ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_annotated.tsv
	
	rm ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_ann.tsv
	rm ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_acc-annotations.tsv
	rm ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_contig-annotations.tsv

	touch ${blastdir}/megahit_filtered_contigs_blastn-nt_hits_annotated.done
else 
	echo "BLAST annotation already complete for ${id}, finishing..."
fi

exit 0

