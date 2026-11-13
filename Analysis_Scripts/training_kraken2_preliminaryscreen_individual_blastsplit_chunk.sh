#!/bin/bash
#SBATCH -c 20     # cores requested
#SBATCH --mem=220G  # memory in Mb
#SBATCH --error     slurm/Kraken2/Kraken2_IL_BLASTsplit_chunk_%A_%a.txt
#SBATCH --output    slurm/Kraken2/Kraken2_IL_BLASTsplit_chunk_%A_%a.txt
#SBATCH -t 128:00:00  # time requested in hour:minute:second

set -e
set -u
set -o pipefail

module purge
module load BLAST+/2.13.0-gompi-2021a

# Check if the output file already exists
if [[ -s $output ]]
then
    echo "Output for chunk $chunk already exists, skipping..."
    exit 0
fi

echo ${chunk}

# Run BLAST command
blastn -task megablast -num_threads $chunkthreads -db $blastdb/nt -query $chunk -out $output -perc_identity $percentmatch -qcov_hsp_perc $querycover -max_target_seqs 1 -outfmt "6 sqseqid    sqgi       sqacc      sqaccver   sqlen      ssseqid    ssallseqid ssgi       ssallgi    ssacc      ssaccver   ssallacc   sslen      sqstart    sqend      ssstart    ssend      sqseq      ssseq      sevalue    sbitscore  sscore     slength    spident    snident    smismatch  spositive  sgapopen   sgaps      sppos      sframes    sqframe    ssframe    sbtop      sstaxids   ssscinames sscomnames ssblastnamessskingdomssstitle    ssalltitlesssstrand  sqcovs    sqcovhsp "

# Check if BLAST command was successful
if [[ $? -ne 0 ]]
then
    echo "Error: BLAST command failed for chunk $chunk"
    exit 1
fi

# Append the header to the first chunk's output
if [[ ! -f ${output}.done ]]
then
    sed -i '1i Query_Contig\tSubject\tIdentity\tMax_Identity\tMismatches\tGaps\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_Value\tScore' $output
    touch ${output}.done
fi

exit 0