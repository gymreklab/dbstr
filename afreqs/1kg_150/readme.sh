#!/bin/bash

for pop in eur afr eas
do
    echo ${pop}
    for chrom in $(seq 1 22)
    do
	VCFFILE=/storage/s1saini/manuscript_strsnp/fig3/hipstr.1kg/${pop}/chr${chrom}/hipstr.chr${chrom}.${pop}.vcf.gz
	./parse_1kg_afreqs.py ${VCFFILE}
    done >  /storage/mgymrek/webstr_files/afreqs/1kg_${pop}_afreqs.csv
done
