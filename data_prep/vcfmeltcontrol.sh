

SrcLoc="/storage/mgymrek/ssc-imputation/filtered_vcfs/"
DestHetro="/storage/resources/het"
DestGTData="/storage/resources/GTdata"
Chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

for Chr in $Chroms
do
nohup python3 vcf_meltgtver2.py "$SrcLoc"hipstr.chr"$Chr".allfilters.vcf.gz "$DestHetro$Chr".csv "$DestGTData$Chr".csv > melt"$Chr".out &
done

