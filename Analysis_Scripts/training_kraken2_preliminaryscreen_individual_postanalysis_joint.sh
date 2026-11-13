#!/bin/bash 
#SBATCH -c 1     # cores requested
#SBATCH --mem=2G  # memory in Mb
#SBATCH --error     slurm/Kraken2/Kraken2_IL_Post_Joint_out_%A_%a.txt
#SBATCH --output    slurm/Kraken2/Kraken2_IL_Post_Joint_out_%A_%a.txt
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

outdir=${krkdir}/DatabaseCreation

if [ -n "${SLURM_JOB_ID:-}" ]
then
	script_path=$(scontrol show job "${SLURM_JOB_ID}" | awk -F= '/Command=/{print $2}')
	scontrol update jobid=$SLURM_JOBID jobname="Kraken2_IndividualLevel_PostAnalysis_Joint"
	IDARRAYTRUE=($(echo $IDARRAYFALSE | tr '-' " " ))
else
	echo "Run using slurm launch script"
	exit 1
fi
script_dir=$(dirname "${script_path}")
parent_dir=$(dirname "${script_dir}")


###   Setup    ###
test_apps="esearch"
for cmd in $test_apps; do
    type &>/dev/null $cmd && echo "edirect installed continuing.." ||
        { echo >&2 "edirect not installed, downloading..."; sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"; export PATH=${PATH}:${HOME}/edirect; }
done

###########################################################################
#########    						Code				         ##########
###########################################################################

mkdir -p ${outdir}

rm -f ${outdir}/Combined_filtered_contigs_blastn-nt_hits.txt 2>&1; touch ${outdir}/Combined_filtered_contigs_blastn-nt_hits.txt 
for (( i=1; i<=$(( ${#IDARRAYTRUE[*]} -1 )); i++ ))
do 
	id="${IDARRAYTRUE[$i]//i}"
    cat ${krkdir}/BLAST/${id}/megahit_filtered_contigs_blastn-nt_hits.txt | cut -f2 | sed 1d | while read line
    do
        # Append the ID as a new column
        echo -e "${line}\t${id}" >> ${outdir}/Combined_filtered_contigs_blastn-nt_hits.txt
    done
done

sort -u -o ${outdir}/Combined_filtered_contigs_blastn-nt_hits.txt ${outdir}/Combined_filtered_contigs_blastn-nt_hits.txt 

mkdir -p ${outdir}/Combined_filtered_contigs_blastn-nt_hits
rm -f ${outdir}/Combined_filtered_contigs_blastn-nt_hits_acc-annotations.tsv 2>&1; touch ${outdir}/Combined_filtered_contigs_blastn-nt_hits_acc-annotations.tsv
cat ${outdir}/Combined_filtered_contigs_blastn-nt_hits.txt | cut -f1 | sort -u | while read -r acc
do
	echo $acc
	taxinfo=$( esearch -db nuccore -query $acc < /dev/null | efetch -format xml | xtract -pattern Org-ref -element Object-id_id,Org-ref_taxname | uniq )
	taxid=$( echo $taxinfo | cut -f1 -d" "  )
	taxname=$( echo $taxinfo | cut -f2- -d" "  )
	if [[ ! -s ${outdir}/Combined_filtered_contigs_blastn-nt_hits/${acc}.fasta ]]
	then
		esearch -db nuccore -query $acc < /dev/null | efetch -format fasta | awk 'NR==1 {gsub("'$acc'", "'$acc'|kraken:taxid|'$taxid'", $0); quit};1' > ${outdir}/Combined_filtered_contigs_blastn-nt_hits/${acc}.fasta
	fi
	name=$( head -n1 ${outdir}/Combined_filtered_contigs_blastn-nt_hits/${acc}.fasta | sed 's/ /'$'\t''/' | cut -f2 )
	echo -e "$acc\t$taxid\t$taxname\t$name" >> ${outdir}/Combined_filtered_contigs_blastn-nt_hits_acc-annotations.tsv
	sleep 1
done

exit 0

