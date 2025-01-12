#cat    final_pdbid_rest.txt  | while read  line

cat   final_6789.list  | while read line

do


#base=${line:1:4}
base=${line:0:4}

wget 'https://files.rcsb.org/download/'$base'.pdb' 


done

