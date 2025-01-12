
for name in   new_antibody_*L_???.pdb
do

	base=${name%.pdb}

nohup pdbfixer  $base'.pdb'  --add-atoms=all --output=$base'_f.pdb' &
sleep 4s



done
