Chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

for chr in $Chroms
do
(bcftools query  -f "%CHROM|%POS|%ID|%REF|%FILTER| %ALT\n" "/storage/mgymrek/ssc-imputation/filtered_vcfs/hipstr.chr"$chr".allfilters.vcf.gz" > "/storage/resources/sscchr"$chr"alt.txt";
Rscript --vanilla "extractalt.R" "/storage/resources/sscchr"$chr"alt.txt" "/storage/resources/Altschr"$chr".csv") > "mout"$chr".out"  &

done


