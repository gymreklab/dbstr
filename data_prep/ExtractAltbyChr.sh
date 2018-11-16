Chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

SrcLoc="/storage/s1saini/manuscript_strsnp/fig3/hipstr.1kg/afr/chr"

for chr in $Chroms
do
(bcftools query  -f "%CHROM|%POS|%ID|%REF|%FILTER| %ALT\n" $SrcLoc$chr"/hipstr.chr"$chr".afr.vcf.gz" > "/storage/resources/sscchr"$chr"alt.txt";
Rscript --vanilla "extractalt.R" "/storage/resources/sscchr"$chr"alt.txt" "/storage/resources/Altschr"$chr".csv") > "mout"$chr".out"  &

done


