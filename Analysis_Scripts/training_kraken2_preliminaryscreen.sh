#!/bin/bash 

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
rtdir=${PROJECT_DIR}
krkdir=${rtdir}/Training
fqd=${rtdir}/Training/Meta_Data/Training_Sequenced_Samples_Info.tsv

krakendb=${GENRES_DIR}/Resources/Kraken2_Databases/
blastdb=${GENRES_DIR}/Resources/BLAST_Databases/
db=Standard_20240605

percentmatch=90
querycover=80
depthmin=30
lengthmin=500

###########################################################################
#########    						Code				         ##########
###########################################################################

mkdir -p slurm/Kraken2/

RDARRAY=()
IDARRAY=()

while read -r read readdir id
do
    if [[ -f ${readdir}/${read} ]]
    then
        #echo "${read} does exist"
    	RDARRAY+=(${read})
		IDARRAY+=(${id})
    else
        echo "${read} does not exist"
    fi
done < <(awk 'BEGIN {FS = "\t";OFS = "\t" } $10 == "Genomic" && $5 == "Illumina" { print $0 }' ${fqd} | awk 'BEGIN {FS = "\t";OFS = "\t" } ($14~/\.*[0-9]$/) { print $0 }' | cut -f6,8,14 | sort -u)

RDARRAY=($(echo ${RDARRAY[@]} | tr ' ' "\n" | sort -u | tr '\n' " " ))
RDARRAYNUM=${#RDARRAY[@]}
RDARRAYFALSE=($(echo ${RDARRAY[@]} | tr ' ' "-" ))

IDARRAY=($(echo ${IDARRAY[@]} | tr ' ' "\n" | sort -u | tr '\n' " " ))
IDARRAYNUM=${#IDARRAY[@]}
IDARRAYFALSE=($(echo ${IDARRAY[@]} | tr ' ' "-" ))

###  Submit Call Jobs  ###

cmd="sbatch -J ReadLevel_Kraken2 \
	--export=ALL,fqd=${fqd},krakendb=${krakendb},blastdb=${blastdb},db=${db},krkdir=${krkdir},RDARRAYFALSE=${RDARRAYFALSE},RDARRAYNUM=${RDARRAYNUM} \
	--array=1-${RDARRAYNUM}%16 ${script_dir}/training_kraken2_preliminaryscreen_read.sh"
#echo $cmd
jobcall=$(eval $cmd | awk '{print $4}')
echo "Kraken2 job submitted: ${jobcall}"

cmd="sbatch -J IndividualLevel_MEGAHIT --dependency=afterok:${jobcall} \
	--export=ALL,fqd=${fqd},krakendb=${krakendb},blastdb=${blastdb},db=${db},krkdir=${krkdir},IDARRAYFALSE=${IDARRAYFALSE},IDARRAYNUM=${IDARRAYNUM} \
	--array=1-${IDARRAYNUM}%12 ${script_dir}/training_kraken2_preliminaryscreen_individual_megahit.sh"
#echo $cmd
jobhitIL=$(eval $cmd | awk '{print $4}')
echo "MEGAHIT individual level job submitted: ${jobhitIL}"

cmd="sbatch -J IndividualLevel_BLAST --dependency=aftercorr:${jobhitIL} \
	--export=ALL,fqd=${fqd},krakendb=${krakendb},blastdb=${blastdb},db=${db},krkdir=${krkdir},IDARRAYFALSE=${IDARRAYFALSE},IDARRAYNUM=${IDARRAYNUM},percentmatch=${percentmatch},querycover=${querycover},depthmin=${depthmin},lengthmin=${lengthmin} \
	--array=1-${IDARRAYNUM}%2 ${script_dir}/training_kraken2_preliminaryscreen_individual_blastsplit.sh"
#echo $cmd
jobhitILB=$(eval $cmd | awk '{print $4}')
echo "BLAST individual level job submitted: ${jobhitILB}"

cmd="sbatch -J Joint_PostAnalysis --dependency=afterok:${jobhitILB} \
	--export=ALL,fqd=${fqd},krakendb=${krakendb},blastdb=${blastdb},db=${db},krkdir=${krkdir},IDARRAYFALSE=${IDARRAYFALSE},IDARRAYNUM=${IDARRAYNUM},percentmatch=${percentmatch},querycover=${querycover},depthmin=${depthmin},lengthmin=${lengthmin} \
	${script_dir}/training_kraken2_preliminaryscreen_individual_postanalysis_joint.sh"
#echo $cmd
jobhitJoP=$(eval $cmd | awk '{print $4}')
echo "BLAST Post Analysis individual level job submitted: ${jobhitJoP}"

cmd="sbatch -J Individual_PostAnalysis --dependency=afterok:${jobhitJoP} \
	--export=ALL,fqd=${fqd},krakendb=${krakendb},blastdb=${blastdb},db=${db},krkdir=${krkdir},IDARRAYFALSE=${IDARRAYFALSE},IDARRAYNUM=${IDARRAYNUM},percentmatch=${percentmatch},querycover=${querycover},depthmin=${depthmin},lengthmin=${lengthmin} \
	--array=1-${IDARRAYNUM}%4 ${script_dir}/training_kraken2_preliminaryscreen_individual_postanalysis.sh"
#echo $cmd
jobhitILP=$(eval $cmd | awk '{print $4}')
echo "BLAST Post Analysis individual level job submitted: ${jobhitILP}"


exit 0
