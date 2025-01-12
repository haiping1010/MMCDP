cat    summary_sort_sel_cutoff6.txt   | while read line
do

IFS=' ' read -r -a array <<< $line
##wget 'http://zinc15.docking.org/substances/'${array[0]}'.sdf'

cp -r    ../static_eval_cutoff10/${array[0]}  .

done

cat    summary_sort_sel_cutoff10.txt   | while read line
do

IFS=' ' read -r -a array <<< $line
##wget 'http://zinc15.docking.org/substances/'${array[0]}'.sdf'

cp -r    ../static_eval_cutoff10/${array[0]}  .

done

