
cd ../
for name in replace_fragement_identify_point_??
do
base=${name:18 }

echo $base

cd   $name

python   read_mutat.py

cd ../

cp  -r  $name/possible_mutations.pkl    all_generate_single_m/$base'_possible_mutations.pkl'


done
