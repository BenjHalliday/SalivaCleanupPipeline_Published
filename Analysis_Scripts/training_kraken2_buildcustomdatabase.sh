#!/bin/bash 

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
stddb=${GENRES_DIR}/Resources/Kraken2_Databases/Standard_20240605/
rtdir=${PROJECT_DIR}
dbcreate=${rtdir}/Training/DatabaseCreation
customdb=${rtdir}/Training/Databases
dbname=Standard_PlusTraining_PlusCHM13
blacklist=${customdb}/blacklist.txt # blacklist for problematic accession numbers

###########################################################################
#########    				Taxonomy Tree Function 		         ##########
###########################################################################

# Written by Frédéric Mahé https://www.biostars.org/p/13452/#google_vignette
# Modified by Ben Halliday

TAXDIR="${dbcreate}/NCBI_Taxonomy/"
NAMES=${TAXDIR}/"names.dmp"
NODES=${TAXDIR}/"nodes.dmp"

if [[ ! -f ${NAMES} ]]
then
	mkdir -p ${TAXDIR}
	if [[ ! -f ${TAXDIR}/taxdump.tar.gz ]]
	then
		wget -P ${TAXDIR} https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
	fi
	tar -xvf ${TAXDIR}/taxdump.tar.gz --directory ${TAXDIR}/
fi

# Obtain the name corresponding to a taxid or the taxid of the parent taxa
get_name_or_taxid() {
    grep --max-count=1 "^${1}"$'\t' "${2}" | cut --fields="${3}"
}

# Loop until you reach the root of the taxonomy (i.e. taxid = 1)
create_taxtree() {
	TAXONOMY_NAME=""
	TAXONOMY_ID=""
	TAXID="${1}"
	while [[ "${TAXID}" -gt 1 ]] ; do
		    	# Obtain the scientific name corresponding to a taxid
    	NAME=$(get_name_or_taxid "${TAXID}" "${NAMES}" "3")
    	# Obtain the parent taxa taxid
    	PARENT=$(get_name_or_taxid "${TAXID}" "${NODES}" "3")
    	# Build the taxonomy path
    	TAXONOMY_NAME="${NAME};${TAXONOMY_NAME}"
    	TAXONOMY_ID="${TAXID};${TAXONOMY_ID}"
    	TAXID="${PARENT}"
	done
	echo -e "${TAXONOMY_ID}"
}

###########################################################################
#########    				Add Saliva Contaminants		         ##########
###########################################################################

mkdir -p slurm/Kraken2_BuildCustomDataBase/
mkdir -p $customdb/$dbname

if [[ ! -f ${blacklist} ]]
then
	touch ${blacklist}
fi

### Add custom fragments, removing any fragments with taxid that are non-human Eukaryota (9606=homo sapiens; 2759=Eukaryota)
cat ${dbcreate}/Combined_filtered_contigs_blastn-nt_hits_acc-annotations.tsv | cut -f1,2 | sort -u -k1,1 | grep -v -f $blacklist | while read acc taxid
do
	TAXTREE=$(create_taxtree "${taxid}")
	if echo "${TAXTREE}" | grep -qw "9606"
	then
		echo "Human, adding..."
		head -n1 ${dbcreate}/Combined_filtered_contigs_blastn-nt_hits/${acc}.fasta
		kraken2-build --add-to-library ${dbcreate}/Combined_filtered_contigs_blastn-nt_hits/${acc}.fasta --db $customdb/$dbname
	else
		if echo "${TAXTREE}" | grep -qw "2759"
		then
			echo "Non-Human Eukaryota, ignoring..."
		else
			echo "Non-Eukaryota, adding"
			head -n1 ${dbcreate}/Combined_filtered_contigs_blastn-nt_hits/${acc}.fasta
			kraken2-build --add-to-library ${dbcreate}/Combined_filtered_contigs_blastn-nt_hits/${acc}.fasta --db $customdb/$dbname
		fi
	fi	
done

###########################################################################
#########    						Add CHM13 			         ##########
###########################################################################

if [[ ! -f ${dbcreate}/CHM13_v2.0_Kraken2Formatted.fna ]]
then

	curl "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz" -o "${dbcreate}/CHM13_v2.0.fna.gz"
	pigz -d "${dbcreate}/CHM13_v2.0.fna.gz"

	awk -v taxid=9606 '/^>/ {sub(/(>NC_[0-9]+\.[0-9]+)/, "&|kraken:taxid|" taxid); print; next} 1' ${dbcreate}/CHM13_v2.0.fna > ${dbcreate}/CHM13_v2.0_Kraken2Formatted.fna
	rm -f ${dbcreate}/CHM13_v2.0.fna 
fi

echo "Adding CHM13..."
kraken2-build --add-to-library ${dbcreate}/CHM13_v2.0_Kraken2Formatted.fna --db $customdb/$dbname

###########################################################################
#########    						Download Tax		         ##########
###########################################################################

### rsync connectivity issue, might take a couple of tries to get running
if [[ ! -f $customdb/$dbname/taxonomy/taxdump.untarflag ]]
then
	MAX_TRIES=20
	attempts=0
	until kraken2-build --download-taxonomy --db $customdb/$dbname
	do
		sleep 5
		if [[ $((++attempts)) -gt $MAX_TRIES ]]
		then
			echo $attempts
			exit 1
		fi
	done
fi

###########################################################################
##	    						Copy over Standard	        			 ##
###########################################################################

# Copying over the same database that was used for the initial screen (not through the standard download due to version conflicts)

cp -r ${stddb}/library/* $customdb/$dbname/library

### Once custom fragments are added and taxonomy is download, build database
srun -c24 --mem 240G kraken2-build --build --threads 24 --max-db-size 100000000000 --db $customdb/$dbname

### Clean up database construction 
#kraken2-build --clean --db $customdb/$dbname

exit 0
