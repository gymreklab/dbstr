

SrcLoc="/storage/mgymrek/ssc-imputation/filtered_vcfs/"
SrcLoc="/storage/s1saini/manuscript_strsnp/fig3/hipstr.1kg/afr/chr"
DestHetro="/storage/resources/het"
DestGTData="/storage/resources/GTdata"
Chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

for Chr in $Chroms
do
nohup python3 vcf_meltgtver2.py "$SrcLoc""$Chr"/hipstr.chr"$Chr".afr.vcf.gz "$DestHetro$Chr".csv "$DestGTData$Chr".csv > melt"$Chr".out &
done

