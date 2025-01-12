
rm summary.txt
for name in *_f.pdb
do
python caculate_interaction.py  $name  >> summary.txt

done
