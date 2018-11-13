Chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"

for chr in $Chroms
do
cp script_of_commands2.txt "script_of_commands"$chr".txt"
echo ".import ../../GTdata"$chr".csv tempGT" >> "script_of_commands"$chr".txt"
echo ".import ../../het"$chr".csv" vcfhomozyg >> "script_of_commands"$chr".txt"
echo ".import Altschr"$chr".csv vcfAlt" >> "script_of_commands"$chr".txt" 
#echo ".import Vepdet"$chr".csv vepdet" >> "script_of_commands"$chr".txt"
echo ".read file2.sql" >> "script_of_commands"$chr".txt"
echo ".exit" >> "script_of_commands"$chr".txt"

nohup sqlite3 "dbSTR"$chr".db" ".read script_of_commands"$chr".txt" &

done
