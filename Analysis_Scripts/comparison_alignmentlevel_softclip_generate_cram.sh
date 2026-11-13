#!/bin/bash 
#SBATCH -c 8    # cores requested
#SBATCH --mem=40G  # memory in Mb
#SBATCH -e slurm_softclip_generate_cram_outfile_%A_%a.txt
#SBATCH -o slurm_softclip_generate_cram_outfile_%A_%a.txt
#SBATCH -t 24:00:00  # time requested in hour:minute:second

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
#module load gridss
module load jvarkit
module load samtools
module load bedtools

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

if [ ! -z ${1+x} ]
then 
	file=$1
fi
if [[ ! -n "$file" ]]
then 
    echo "No CRAM specified, exiting"
    exit 1
else
	if [[ ! -s "$file" && "$file" == *.cram  ]]
	then
		echo "${file} doesn't exist or is not a CRAM file, exiting"
		exit 1
	fi
fi

if [ ! -z ${2+x} ]
then 
	type=$2
fi

andir=$( dirname ${file} )
id=$( basename ${andir} )

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

export TMPDIR=${HOME}/mamba_tmp/$SLURM_JOB_ID
export APPTAINER_TMPDIR=$TMPDIR
export APPTAINER_CACHEDIR=$TMPDIR
export MAMBA_ROOT_PREFIX=$TMPDIR
export XDG_CACHE_HOME=$TMPDIR
mkdir -p $TMPDIR
chmod 770 $TMPDIR

if [ -n "${SLURM_MEM_PER_NODE:-}" ]
then
    mem=$(echo "scale=0; $SLURM_MEM_PER_NODE / 1024" | bc)
else
    mem=40
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=${SLURM_CPUS_PER_TASK}
else
    threads=8
fi

## Inputs
wkdir=${PROJECT_DIR}
fqd=${wkdir}/MetaData/Comparison_Sequenced_Samples_Cleanup.tsv
resdir=${GENRES_DIR}/Resources/hg38
repeats=${resdir}/Repeats/rmsk.bed.gz
cytobands=${resdir}/Cytobands/cytoBand.bed.gz
umap=${resdir}/Mappability/k24.Umap.MultiTrackMappability.bed.gz

if [[ "$type" == "decoy"  ]]
then
	ref=${GENRES_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_plus_MT169739_analysis_set.fna
	regions=${GENRES_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_plus_MT169739_analysis_set.fna.bed
	echo "${file} is aligned using a decoy, using: ${ref}"
else
	ref=${GENRES_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
	regions=${GENRES_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bed
	echo "${file} is aligned using the standard reference, using: ${ref}"
fi

refname=$( echo $ref | sed 's/\.[^.]*$//' )

###########################################################################
#########    			Extract Double SoftClipped		         ##########
###########################################################################

mkdir -p ${andir}/CollectCigarMetrics/

if [[ ! -f ${refname}.genome ]]
then
	cut -f1,2 ${ref}.fai > ${refname}.genome
fi

if [[ ! -f ${andir}/CollectCigarMetrics/RepeatMasker.bed ]]
then
	zcat ${repeats} | awk '$1 !~ /_/' > ${andir}/CollectCigarMetrics/RepeatMasker.bed
fi

if [[ ! -s ${andir}/CollectCigarMetrics/RepeatMasker+Buffer150.bed ]]
then
	bedtools slop -i ${andir}/CollectCigarMetrics/RepeatMasker.bed -g ${refname}.genome -b 150 > ${andir}/CollectCigarMetrics/RepeatMasker+Buffer150.bed
fi

if [[ ! -s ${andir}/CollectCigarMetrics/TargetCytobands.bed ]]
then
	zcat ${cytobands} | awk -F'\t' '$5 ~ /gneg|gpos/' | awk '$1 !~ /_/' | awk '$1 !~ /chrY/' > ${andir}/CollectCigarMetrics/TargetCytobands.bed
fi

if [[ ! -s ${andir}/CollectCigarMetrics/Mappability_k24.bed ]]
then
	zcat ${umap} | bedtools merge -i - > ${andir}/CollectCigarMetrics/MappableRegions_k24.bed
fi

if [[ ! -s ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped.bam.csi ]]
then
	cmd="samjdk --reference ${ref} --expression 'if(record.getReadUnmappedFlag()) return false;final Cigar c=record.getCigar();if(c==null || c.numCigarElements()<2) return false; return (c.getFirstCigarElement().getOperator().isClipping() && c.getLastCigarElement().getOperator().isClipping()) ;' ${file} | samtools view --write-index -b -h -L ${regions} -O BAM -o ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped.bam"
	echo -e "Extracting double softclipped reads for ${file}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
else
	echo "Double softclipped reads already generated for ${file}"
fi

if [[ ! -s ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeats.bam.csi ]]
then
	cmd="bedtools intersect -a ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped.bam -b ${andir}/CollectCigarMetrics/RepeatMasker+Buffer150.bed -wa | samtools view -h -@ ${threads} --write-index -O BAM -o ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeats.bam"
	echo -e "Extracting repeat overlapping double softclipped reads for ${file}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
else
	echo "Double softclipped repeat overlapping reads already generated for ${file}"
fi

if [[ ! -s ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered.bam.csi ]]
then
	cmd="bedtools subtract -a ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped.bam -b ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeats.bam | samtools view -h -@ ${threads} --write-index -O BAM -o ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered.bam"
	echo -e "Extracting non repeat overlapping double softclipped reads for ${file}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
else
	echo "Double softclipped non repeat overlapping reads already generated for ${file}"
fi

if [[ ! -s ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_TargetCytobands.bam.csi ]]
then
	cmd="samtools view -h -@ ${threads} --write-index -O BAM -L ${andir}/CollectCigarMetrics/TargetCytobands.bed -o ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_TargetCytobands.bam ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered.bam"
	echo -e "Extracting target cytoband softclipped reads for ${file}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
else
	echo "Double softclipped non repeat overlapping reads already generated for ${file}"
fi

if [[ ! -s ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_stats.tsv ]]
then
	cmd="samtools stats -@ ${threads} ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped.bam > ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_stats.tsv"
	echo -e "Extracting stats for raw softclipped reads for ${file}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
else
	echo "Stats for raw softclipped target reads already generated for ${file}"
fi

if [[ ! -s ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_TargetCytobands_stats.tsv ]]
then
	cmd="samtools stats -@ ${threads} ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_TargetCytobands.bam > ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_TargetCytobands_stats.tsv"
	echo -e "Extracting stats for cytoband softclipped reads for ${file}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
else
	echo "Stats for cytoband softclipped target reads already generated for ${file}"
fi

if [[ ! -s ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_TargetCytobands_FreqTable.tsv ]]
then
	#cmd="samtools view ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_TargetCytobands.bam | cut -f 3,4 | uniq -c | sort -k1,1 -n > ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_TargetCytobands_FreqTable.tsv"
	#echo -e "   $cmd"
	#eval "$cmd" || exit $?
	
	samtools view ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_TargetCytobands.bam | cut -f 3,4,6 | awk '{ $3 = gensub(/^[0-9]+S|[0-9]+S$/, "", "g", $3); print $0 }' | awk '{ col = $3; sum = 0; while (match(col, /[0-9]+/)) { sum += substr(col, RSTART, RLENGTH); col = substr(col, RSTART + RLENGTH);} print $1,$2,$2+sum,sum }' | sort -k1,1V -k2,2n -k3,3n | uniq -c > ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_TargetCytobands_FreqTable.tsv
	
	echo -e "Extracting Freq table for cytoband softclipped reads for ${file}..."
else
	echo "Freq table for cytoband softclipped target reads already generated for ${file}"
fi

if [[ ! -s ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_Mappability.bam.csi ]]
then
	cmd="bedtools intersect -f 1.0 -abam ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered.bam -b ${andir}/CollectCigarMetrics/MappableRegions_k24.bed > ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_Mappability.bam"
	echo -e "Extracting mappable softclipped reads for ${file}..."
	echo -e "   $cmd"
	eval "$cmd" || exit $?
else
	echo "Double softclipped non repeat overlapping reads already generated for ${file}"
fi

if [[ ! -s ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_Mappability_FreqTable.tsv ]]
then
	samtools view ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_Mappability.bam | cut -f 3,4,6 | awk '{ $3 = gensub(/^[0-9]+S|[0-9]+S$/, "", "g", $3); print $0 }' | awk '{ col = $3; sum = 0; while (match(col, /[0-9]+/)) { sum += substr(col, RSTART, RLENGTH); col = substr(col, RSTART + RLENGTH);} print $1,$2,$2+sum,sum }' | sort -k1,1V -k2,2n -k3,3n | uniq -c > ${andir}/CollectCigarMetrics/${id}_DoubleSoftClipped_Repeat_Filtered_Mappability_FreqTable.tsv	
	echo -e "Extracting Freq table for cytoband softclipped reads for ${file}..."
else
	echo "Freq table for cytoband softclipped target reads already generated for ${file}"
fi

exit 0
