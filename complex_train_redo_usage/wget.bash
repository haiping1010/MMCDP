#cat    final_pdbid_rest.txt  | while read  line

cat   final_tail30000.list  | while read line

do


#base=${line:1:4}
base=${line:0:4}

wget 'https://files.rcsb.org/download/'$base'.pdb'


done

