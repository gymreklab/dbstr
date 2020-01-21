#!/bin/bash

for chrom in $(seq 1 22)
do
    VCFFILE=/storage/mgymrek/ssc-imputation/filtered_vcfs/hipstr.chr${chrom}.allfilters.vcf.gz
    cmd="./parse_ssc_afreqs.py ${VCFFILE} > /storage/mgymrek/webstr_files/afreqs/ssc_hipstr_afreqs_${chrom}.csv"
    echo $cmd
done | xargs -I% -n1 -P10 sh -c "%"
