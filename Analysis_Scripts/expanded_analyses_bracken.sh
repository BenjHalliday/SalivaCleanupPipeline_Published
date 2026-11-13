#!/bin/bash 
#SBATCH -c 1    # cores requested
#SBATCH --mem=2G  # memory in Mb
#SBATCH -e slurm_saliva_expanded_analyses_bracken_%A.txt
#SBATCH -o slurm_saliva_expanded_analyses_bracken_%A.txt
#SBATCH -t 2:00:00  # time requested in hour:minute:second

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load bracken
module load krakentools

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
wkdir=${PROJECT_DIR}
db=${wkdir}/Training/Databases/Standard_PlusTraining_PlusCHM13/

train_meta=${wkdir}/Training/Meta_Data/Training_Sequenced_Samples_Info.tsv
train_algdir=${wkdir}/Training/Kraken2

val_meta=${wkdir}/Validation/Meta_Data/Validation_Sequenced_Samples_Info.tsv
val_algdir=${wkdir}/Validation/Alignments/CMM_Alignments

comp_meta=${wkdir}/MetaData/Comparison_Sequenced_Samples_Cleanup.tsv
comp_algdir=${wkdir}/Alignments

if [ -n "${SLURM_JOB_ID:-}" ]
then
	scontrol update jobid=${SLURM_JOBID} jobname="expanded_analyses_bracken"
	echo "This job is running on node: $SLURMD_NODENAME"
fi

if [ -n "${SLURM_MEM_PER_NODE:-}" ]
then
    mem=$(echo "scale=0; $SLURM_MEM_PER_NODE / 1024" | bc)
else
    mem=2
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=${SLURM_CPUS_PER_TASK} 
else
    threads=1
fi

###########################################################################
#########    				ReadName Approach			         ##########
###########################################################################

if [[ ! -s ${db}/database150mers.kmer_distrib ]]
then
	echo "build bracken database using: bracken-build -d ${db} -t ${threads} -k 35 -l 150"
	exit 1
fi

cat $train_meta | sed 1d | cut -f1 | sort -u | while read -r id
do
	if [[ ! -s ${train_algdir}/${id}/${id}_Kraken2_Report.tsv ]]
	then
		combine_kreports.py --only-combined -r ${train_algdir}/${id}/*/*_Report.tsv -o ${train_algdir}/${id}/${id}_Kraken2_Report.tsv
	fi
	
	if [[ ! -s ${train_algdir}/${id}/${id}_Kraken2_Report_Species.bracken ]]
	then
		bracken -d ${db} -i ${train_algdir}/${id}/${id}_Kraken2_Report.tsv \
			-o ${train_algdir}/${id}/${id}_Kraken2_Report_Species.bracken \
			-w ${train_algdir}/${id}/${id}_Kraken2_Report_Species.tsv -r 150
	fi
	
	if [[ ! -s ${train_algdir}/${id}/${id}_Kraken2_Report_Species_NonHuman.bracken ]]
	then
		filter_bracken.out.py -i ${train_algdir}/${id}/${id}_Kraken2_Report_Species.bracken \
			-o ${train_algdir}/${id}/${id}_Kraken2_Report_Species_NonHuman.bracken --exclude 9606
	fi
	
	if [[ ! -s ${train_algdir}/${id}/${id}_Kraken2_Report_Genus.bracken ]]
	then
		bracken -d ${db} -i ${train_algdir}/${id}/${id}_Kraken2_Report.tsv \
			-o ${train_algdir}/${id}/${id}_Kraken2_Report_Genus.bracken \
			-w ${train_algdir}/${id}/${id}_Kraken2_Report_Genus.tsv -r 150 -l G
	fi

	if [[ ! -s ${train_algdir}/${id}/${id}_Kraken2_Report_Genus_NonHuman.bracken ]]
	then
		filter_bracken.out.py -i ${train_algdir}/${id}/${id}_Kraken2_Report_Genus.bracken \
			-o ${train_algdir}/${id}/${id}_Kraken2_Report_Genus_NonHuman.bracken --exclude 9605
	fi
done

cat $val_meta | sed 1d | cut -f1 | sort -u | while read -r id
do
	if [[ ! -s ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report ]]
	then
		combine_kreports.py --only-combined -r ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/*_P001_Kraken2Screened_Report.tsv -o ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report.tsv
	fi
	
	if [[ ! -s ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Species.bracken ]]
	then
		bracken -d ${db} -i ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report.tsv \
			-o ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Species.bracken \
			-w ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Species.tsv -r 150
	fi
	
	if [[ ! -s ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Species_NonHuman.bracken ]]
	then
		filter_bracken.out.py -i ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Species.bracken \
			-o ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Species_NonHuman.bracken --exclude 9606
	fi
	
	if [[ ! -s ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Genus.bracken ]]
	then
		bracken -d ${db} -i ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report.tsv \
			-o ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Genus.bracken \
			-w ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Genus.tsv -r 150 -l G
	fi

	if [[ ! -s ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Genus_NonHuman.bracken ]]
	then
		filter_bracken.out.py -i ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Genus.bracken \
			-o ${val_algdir}/${id}_NoFastP_Kraken_trained_taxid_k24/Kraken2/${id}_NoFastP_Kraken_trained_taxid_k24_Kraken2Screened_Report_Genus_NonHuman.bracken --exclude 9605
	fi
done

cat $comp_meta | sed 1d | cut -f1 | sort -u | while read -r id
do
	if [[ ! -s ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report.tsv ]]
	then
		combine_kreports.py --only-combined -r ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/*_subsampled_Kraken2Screened_Report.tsv -o ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report.tsv
	fi
	
	if [[ ! -s ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Species.bracken ]]
	then
		bracken -d ${db} -i ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report.tsv \
			-o ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Species.bracken \
			-w ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Species.tsv -r 150
	fi
	
	if [[ ! -s ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Species_NonHuman.bracken ]]
	then
		filter_bracken.out.py -i ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Species.bracken \
			-o ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Species_NonHuman.bracken --exclude 9606
	fi
	
	if [[ ! -s ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Genus.bracken ]]
	then
		bracken -d ${db} -i ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report.tsv \
			-o ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Genus.bracken \
			-w ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Genus.tsv -r 150 -l G
	fi

	if [[ ! -s ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Genus_NonHuman.bracken ]]
	then
		filter_bracken.out.py -i ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Genus.bracken \
			-o ${comp_algdir}/${id}_NoFastP_Krakentaxid_k19/Kraken2/${id}_NoFastP_Krakentaxid_k19_Kraken2_Report_Genus_NonHuman.bracken --exclude 9605
	fi
done

###########################################################################
#########    					 Finish 				         ##########
###########################################################################

exit 0
