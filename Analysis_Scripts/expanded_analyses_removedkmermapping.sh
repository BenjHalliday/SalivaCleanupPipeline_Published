#!/bin/bash 
#SBATCH -c 2    # cores requested
#SBATCH --mem=40G  # memory in Mb
#SBATCH -e slurm_kmer_alignmnet_mapping_outfile_%A_%a.txt
#SBATCH -o slurm_kmer_alignmnet_mapping_outfile_%A_%a.txt
#SBATCH -t 24:00:00  # time requested in hour:minute:second

### Extracting reads that are newly unmapped when the aligner seed length is increased from k19 to k${kmer}
### Newly unmapped = unmapped at k${kmer}, not unmapped at k19. Only the unmapped mate of the pair is extracted.

set -eo pipefail
export LC_ALL=C    # same collation for sort and comm

###########################################################################
#########    					Modules					         ##########
###########################################################################

module purge
module load samtools

###########################################################################
#########    				Inputs and Setup			         ##########
###########################################################################

id=${1:-}
kmer=${2:-}
type=${3:-B}

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
andir=${wkdir}/Analyses/Removed_Read_Analysis/data_raw/kmer_${id}_NoFastP_NoKraken

bam=${wkdir}/Alignments/${id}${type}_NoFastP_NoKraken_k19/${id}${type}.cram
bamnew=${wkdir}/Alignments/${id}${type}_NoFastP_NoKraken_k${kmer}/${id}${type}.cram

fqd=${wkdir}/MetaData/Comparison_Sequenced_Samples_Cleanup.tsv
ref=${GENRES_DIR}/Resources/hg38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

unmapped_mates=${andir}/unmapped_mates
lost=${unmapped_mates}_k${kmer}_lost.txt

export TMP_DIR="/weka/users/${USER}/${SLURM_JOB_ID:-$$}"
export MAMBA_ROOT_PREFIX="${TMP_DIR}"
export XDG_CACHE_HOME="${TMP_DIR}"
mkdir -p "${TMP_DIR}"

trap 'rm -rf "$TMP_DIR"' EXIT

if [[ -z "$id" ]]
then 
    echo "No id specified, exiting"
    exit 1
else
	if ! awk -F'\t' -v id="$id" '$1 ~ "^" id {found=1} END{exit !found}' "$fqd"
	then
    	echo "'${id}' not found in $fqd"
    	echo "Exiting..."
    	exit 1
    else
      	echo "'${id}' found in $fqd"  	
	fi
fi

if [[ -z "$kmer" ]]
then 
    echo "No kmers file specified, exiting"
    exit 1
else
	if [[ "$kmer" =~ ^[0-9]+$ ]] && (( $kmer > 19 && $kmer <= 27 ))
	then
    	echo "kmer is between 20 and 27"
	else
    	echo "kmer is NOT between 20 and 27 (k19 is the baseline), exiting..."
    	exit 1
	fi
fi

for f in "${bam}" "${bamnew}"
do
	if [[ ! -s ${f} ]]
	then
		echo "Missing input: ${f}"
		exit 1
	fi
done

if [ -n "${SLURM_JOB_ID:-}" ]
then
	scontrol update jobid=${SLURM_JOBID} jobname="kmer_alignment_mapping_${id}_${kmer}"
	echo "This job is running on node: $SLURMD_NODENAME"
fi

if [ -n "${SLURM_CPUS_PER_TASK:-}" ]
then
    threads=${SLURM_CPUS_PER_TASK} 
else
    threads=4
fi

###########################################################################
#########    				ReadName + Mate Approach		     ##########
###########################################################################

mkdir -p ${andir}

# One line per unmapped record: read name <tab> mate (1 or 2, from flags 0x40/0x80)
list_unmapped() {
	samtools view -@ ${threads} -T ${ref} -f 4 $1 | awk '{print $1 "\t" (int($2/64)%2 ? 1 : 2)}' | sort -u -T ${MYTMPDIR}
}

if [[ ! -s ${unmapped_mates}_k19.txt ]]
then
	list_unmapped ${bam} > ${unmapped_mates}_k19.txt
fi

if [[ ! -s ${unmapped_mates}_k${kmer}.txt ]]
then
	list_unmapped ${bamnew} > ${unmapped_mates}_k${kmer}.txt
fi

if [[ ! -e ${lost} ]]
then
	comm -13 ${unmapped_mates}_k19.txt ${unmapped_mates}_k${kmer}.txt > ${lost}
fi
echo "Newly unmapped mates (k19 -> k${kmer}): $(wc -l < ${lost})"

if [[ ! -s ${unmapped_mates}_k${kmer}.txt ]]
then
	cut -f1 ${lost} | sort -u > ${lost%.txt}.names
fi

extract_mate() {   # usage: extract_mate <cram> <out.bam> [extra samtools view flags]
	local in=$1 out=$2; shift 2
	samtools view -@ ${threads} -T ${ref} -h -F 0x900 "$@" -N ${lost%.txt}.names ${in} | awk -F'\t' 'FNR==NR {L[$1 SUBSEP $2]=1; next}; /^@/ {print; next}; {m = (int($2/64)%2 ? 1 : 2)}; (($1 SUBSEP m) in L)' ${lost} - | samtools view -@ ${threads} -b -o ${out} -
}

if [[ -s ${lost} ]]
then
	# Lost reads as represented in k19
	if [[ ! -s ${lost%.txt}.bam ]]
	then
		extract_mate ${bam} ${lost%.txt}.bam
		samtools index ${lost%.txt}.bam
	fi
	# Lost reads as represented in k${kmer}
	if [[ ! -s ${lost%.txt}.unmapped.bam ]]
	then
		extract_mate ${bamnew} ${lost%.txt}.unmapped.bam -f 4
	fi
else
	echo "No newly unmapped mates for ${id} at k${kmer}, nothing to extract."
fi

###########################################################################
#########    					Finish					         ##########
###########################################################################

ctrl_frac=005	# sampled control read %
minclip=1       # min soft-clipped to call a read double-soft-clipped
seed=42

#samtools idxstats $bam | awk '{s+=$3} END{print s}' > $andir/n_mapped_k19.txt   # denominator for per-million rates

cat > $TMP_DIR/feat.awk <<'EOF'
function has(f,b){return int(f/b)%2}
function na(x){return x==""?"NA":x}
function cig(c,   n,op,first){aligned=0;qlen=0;hard=0;sl=0;sr=0;first=1
  while(match(c,/^[0-9]+[MIDNSHP=X]/)){n=substr(c,1,RLENGTH-1)+0;op=substr(c,RLENGTH,1);c=substr(c,RLENGTH+1)
    if(op=="M"||op=="="||op=="X"){aligned+=n;qlen+=n}else if(op=="I")qlen+=n
    else if(op=="S"){qlen+=n;if(first)sl=n;else sr=n}else if(op=="H")hard+=n
    if(op!="H")first=0}}
BEGIN{FS=OFS="\t"
  if(lostfile!=""){while((getline l<lostfile)>0){split(l,a,"\t");L[a[1] SUBSEP a[2]]=1};close(lostfile)}
  if(skipfile!=""){while((getline l<skipfile)>0)S[l]=1;close(skipfile)}
  print "name","mate","chr","pos","mapq","mapq0","readlen","aligned","aligned_frac","clip_frac","dbl_clip","nm","longest_exact","has_sa","n_xa","proper","tlen","mate_mapped","mate_chr","mate_pos","mate_mapq","mate_mapq0","mate_same_chr","mate_also_lost"}
/^@/{next}
{f=$2; if(has(f,256)||has(f,2048))next
 m=has(f,64)?1:2; n=$1
 mf[n,m]=f;mq[n,m]=$5;mc[n,m]=$3;mp[n,m]=$4
 if(has(f,4))next
 if(lostfile!=""){if(!((n SUBSEP m) in L))next} else if(n in S)next
 cig($6); rl=qlen+hard
 nm="";md="";sa=0;xa=0
 for(i=12;i<=NF;i++){t=substr($i,1,5)
  if(t=="NM:i:")nm=substr($i,6)+0; else if(t=="MD:Z:")md=substr($i,6); else if(t=="SA:Z:")sa=1
  else if(t=="XA:Z:"){x=substr($i,6);xa=gsub(/;/,";",x)}}
 le="";if(md!=""){le=0;while(match(md,/[0-9]+/)){v=substr(md,RSTART,RLENGTH)+0;if(v>le)le=v;md=substr(md,RSTART+RLENGTH)}}
 row[n,m]=n OFS m OFS $3 OFS $4 OFS $5 OFS ($5==0) OFS rl OFS aligned OFS (rl?aligned/rl:0) OFS (rl?(sl+sr)/rl:0) OFS (sl>=minclip&&sr>=minclip) OFS na(nm) OFS na(le) OFS sa OFS xa OFS has(f,2) OFS $9
 rch[n,m]=$3}
END{for(k in row){split(k,a,SUBSEP);n=a[1];m=a[2];o=(m==1?2:1)
  if((n,o) in mf){mm=!has(mf[n,o],4);mch=mc[n,o];mpo=mp[n,o];mmq=(mm?mq[n,o]:"NA")}else{mm=0;mch="NA";mpo="NA";mmq="NA"}
  print row[k],mm,mch,mpo,mmq,(mmq=="NA"?"NA":(mmq+0==0)),(mm&&mch==rch[k]),((n SUBSEP o) in L)}}
EOF

if [[ ! -s $andir/lost_features_k${kmer}.tsv.gz ]]
then
	# lost
	samtools view -@ $threads -T $ref -F 0x900 -N ${lost%.txt}.names $bam | awk -v lostfile=$lost -v minclip=$minclip -f $TMP_DIR/feat.awk | gzip > $andir/lost_features_k${kmer}.tsv.gz
fi

if [[ ! -s $andir/ctrl_features_k${kmer}.tsv.gz ]]
then
	# control (excluding lost)
	samtools view -@ $threads -T $ref -F 0x900 -s ${seed}.${ctrl_frac} $bam | awk -v skipfile=${lost%.txt}.names -v minclip=$minclip -f $TMP_DIR/feat.awk | gzip > $andir/ctrl_features_k${kmer}.tsv.gz
fi

###########################################################################
#########    					Finish					         ##########
###########################################################################

exit 0
