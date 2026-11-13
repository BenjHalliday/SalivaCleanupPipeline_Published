#!/bin/bash 

set -e
set -u
set -o pipefail

###  Modules  ###
module purge
module load bcftools
module load htslib

###  Inputs   ###
vcf=$1

vcf_human_readable () {
    input=$1
    temp=$(echo ${input} | sed "s/\.gz//g")
    out=$(echo ${input} | sed "s/\.vcf.gz/_out\.tsv/g")
    bgzip -d $input -c > ${temp}
    info=; for var in `less ${temp} | grep "##INFO" | cut -d\, -f1 | sed "s/##INFO=<ID=//g" | grep -v "-"`; do if [ -z "$info" ]; then info="%INFO/${var}"; else info="${info}\t%INFO/${var}"; fi; done;
    format=; for var in `less ${temp} | grep "^##FORMAT" | cut -d\, -f1 | sed "s/##FORMAT=<ID=//g"`; do if [ -z "$format" ]; then format="%${var}"; else format="${format}\t%${var}"; fi; done; 
    bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t'$info'\t['$format'\t]\n' ${temp} > ${out}
	wait < <(jobs -p)
    rm ${temp}
}


if [[ ! -f "$vcf" ]]
then
    echo "Input VCF file ${vcf} doesn't exist, exiting..."
else
    vcf_human_readable ${vcf}
    wait
fi

exit 0
