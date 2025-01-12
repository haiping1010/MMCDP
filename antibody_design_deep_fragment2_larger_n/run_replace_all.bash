
for name  in  replace_fragement_identify_point_??
do

cd  $name

rm  -rf dict_*.npy
rm -rf *.pkl
python   read_mutat.py

cd  ../

done
