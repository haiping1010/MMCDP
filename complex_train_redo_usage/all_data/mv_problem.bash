mkdir problem_fold
cat    need_remove.txt  | while read line
do

IFS=' ' read -r -a array <<< $line
##wget 'http://zinc15.docking.org/substances/'${array[0]}'.sdf'

mv  positive/*${array[0]}*  problem_fold
mv  negative/*${array[0]}*  problem_fold


done
