#!/bin/bash 
#SBATCH -c 2    # cores requested
#SBATCH --mem=4G  # memory in Mb
#SBATCH -e slurm_fastq_subsampling_outfile_%A_%a.txt
#SBATCH -o slurm_fastq_subsampling_outfile_%A_%a.txt
#SBATCH -t 24:00:00  # time requested in hour:minute:second

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load fastq-screen_multi
module load seqtk/1.3

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

if [ ! -z ${1+x} ]
then 
	id=$1
fi
if [[ ! -n "$id" ]]
then 
    echo "No id specified, exiting"
    exit 1
fi

id1=${id}B
id2=${id}S

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
wkdir=${PROJECT_DIR}
illdir=${wkdir}/IlluminaReads

fqd=${META_DIR}/Comparison_Sequenced_Samples_Cleanup.tsv
fastqscreenconfig=${META_DIR}/FastQ_Screen.config

###########################################################################
#########    				Run Prep					         ##########
###########################################################################

mkdir -p ${illdir}/RawReads
mkdir -p ${illdir}/FastQStats
mkdir -p ${illdir}/SubsampledReads

rm -f ${illdir}/FastQStats/${id1}/${id1}_stats.txt
rm -f ${illdir}/FastQStats/${id2}/${id2}_stats.txt

# Generate Read Stats
while read -r id R1 R2 DIR
do
	mkdir -p ${illdir}/RawReads/${id}
 	mkdir -p ${illdir}/FastQStats/${id}  
 	mkdir -p ${illdir}/SubsampledReads/${id}  
 	
 	if [[ ! -s "${illdir}/RawReads/${id}/${R1}" ]]
	then
    	cp "${DIR}/${R1}" "${illdir}/RawReads/${id}"
    fi
	if [[ ! -s "${illdir}/RawReads/${id}/${R2}" ]]
	then
    	cp "${DIR}/${R2}" "${illdir}/RawReads/${id}"
    fi
 
    fq1=$(echo ${R1} | sed 's/\.fastq\.gz//g')
    fq2=$(echo ${R2} | sed 's/\.fastq\.gz//g')
    
	if [[ ! -s "${illdir}/FastQStats/${id}/${fq1}_screen.html" ]]
	then
    	fastq_screen --aligner BWA --outdir ${illdir}/FastQStats/${id} --conf ${fastqscreenconfig} --force --subset 10000 "${illdir}/RawReads/${id}/${R1}"    	
    fi
  	if [[ ! -s "${illdir}/FastQStats/${id}/${fq2}_screen.html" ]]
	then
    	fastq_screen --aligner BWA --outdir ${illdir}/FastQStats/${id} --conf ${fastqscreenconfig} --force --subset 10000 "${illdir}/RawReads/${id}/${R2}"    	
    fi
    
    if [[ ! -s "${illdir}/FastQStats/${id}/${fq1}_stats.txt" ]]
	then
   	    awk '{s++} END {printf "%.0f\n", s/4}' "${illdir}/RawReads/${id}/${R1}" > ${illdir}/FastQStats/${id}/${fq1}_stats.txt
    fi
    if [[ ! -s "${illdir}/FastQStats/${id}/${fq2}_stats.txt" ]]
	then
   	    awk '{s++} END {printf "%.0f\n", s/4}' "${illdir}/RawReads/${id}/${R2}" > ${illdir}/FastQStats/${id}/${fq2}_stats.txt
    fi
    
    human_content=$( awk 'BEGIN {FS = "\t";OFS = "\t" } $1 == "Human" { print 1 - ( $3/$2 ) }' "${illdir}/FastQStats/${id}/${fq1}_screen.txt" )
    read_number=$( cat "${illdir}/FastQStats/${id}/${fq1}_stats.txt" ) 
    human_number=$( echo "${human_content}*${read_number}" | bc )
    echo -e ""${illdir}/RawReads/${id}/${R1}"\t${read_number}\t${human_content}\t${human_number}" >> ${illdir}/FastQStats/${id}/${id}_stats.txt
    
done < <(awk 'BEGIN {FS = "\t";OFS = "\t" } ( $1 == "'${id1}'" || $1 == "'${id2}'" ) && $5 == "Illumina" { print $0 }' ${fqd} | cut -f1,6-8)

id1_total_sum=$( cut -f4 ${illdir}/FastQStats/${id1}/${id1}_stats.txt | awk '{s+=$1} END {printf "%.0f\n", s}' )
id2_total_sum=$( cut -f4 ${illdir}/FastQStats/${id2}/${id2}_stats.txt | awk '{s+=$1} END {printf "%.0f\n", s}' )

if [ "$id1_total_sum" -gt "$id2_total_sum" ]
then
	echo "${id1} human read number ($id1_total_sum) is greater than ${id2} human read number ($id2_total_sum)"
	scaling=$(echo "$id2_total_sum/$id1_total_sum" | bc -l)
	scal_id=${id1}
	echo "Subsampling ${id1} reads by ${scaling}"
elif [ "$id1_total_sum" -lt "$id2_total_sum" ]
then
	echo "${id1} human read number ($id1_total_sum) is less than ${id2} human read number ($id2_total_sum)"
	scaling=$(echo "$id1_total_sum/$id2_total_sum" | bc -l)
	scal_id=${id2}
	echo "Subsampling ${id2} reads by ${scaling}"
else
	echo "${id1} human read number ($id1_total_sum) is equal to ${id2} human read number ($id2_total_sum), or something has gone wrong. That sucks."
fi

while read -r id R1 R2 DIR
do
    fq1=$(echo ${R1} | sed 's/_R1\.fastq\.gz//g')
    fq2=$(echo ${R2} | sed 's/_R2\.fastq\.gz//g')

	if [[ "${scal_id}" == "${id}" ]]
	then
		if [[ ! -f ${illdir}/SubsampledReads/${id}/${fq1}_subsampled_R1.fastq.gz.done ]]
		then
			echo "${id} read file ${R1} is being subsampled by ${scaling}, subsampling..."
			sbatch -c1 --mem 1G --wrap="seqtk sample -s123 ${illdir}/RawReads/${id}/${R1} ${scaling} | gzip > ${illdir}/SubsampledReads/${id}/${fq1}_subsampled_R1.fastq.gz"
			touch ${illdir}/SubsampledReads/${id}/${fq1}_subsampled_R1.fastq.gz.done
		fi
		if [[ ! -f ${illdir}/SubsampledReads/${id}/${fq2}_subsampled_R2.fastq.gz.done ]]
		then
			echo "${id} read file ${R2} is being subsampled by ${scaling}, subsampling..."
			sbatch -c1 --mem 1G --wrap="seqtk sample -s123 ${illdir}/RawReads/${id}/${R2} ${scaling} | gzip > ${illdir}/SubsampledReads/${id}/${fq2}_subsampled_R2.fastq.gz"
			touch ${illdir}/SubsampledReads/${id}/${fq2}_subsampled_R2.fastq.gz.done
		fi
	else
		if [[ ! -f ${illdir}/SubsampledReads/${id}/${fq1}_subsampled_R1.fastq.gz.done ]]
		then
			echo "${id} read file ${R1} does not need scaling, copying.."	
			cp ${illdir}/RawReads/${id}/${R1} ${illdir}/SubsampledReads/${id}/${fq1}_subsampled_R1.fastq.gz
			touch ${illdir}/SubsampledReads/${id}/${fq1}_subsampled_R1.fastq.gz.done
		fi
		if [[ ! -f ${illdir}/SubsampledReads/${id}/${fq2}_subsampled_R2.fastq.gz.done ]]
		then
			echo "${id} read file ${R2} does not need scaling, copying.."	
			cp ${illdir}/RawReads/${id}/${R2} ${illdir}/SubsampledReads/${id}/${fq2}_subsampled_R2.fastq.gz
			touch ${illdir}/SubsampledReads/${id}/${fq2}_subsampled_R2.fastq.gz.done
		fi
	fi   
done < <(awk 'BEGIN {FS = "\t";OFS = "\t" } ( $1 == "'${id1}'" || $1 == "'${id2}'" ) && $5 == "Illumina" { print $0 }' ${fqd} | cut -f1,6-8)

exit 0
