
file=$1

echo $file
makeblastdb -in  $file   -dbtype prot

