



for name in  L{1..2}_query.fasta
do
second_char=${name:1:1}

base=${name%.fasta}
basex=${name%_query.fasta}

id=${name:1:1}

echo $id
blastp -query $name -db  'all_cdr'$id'_aa_light_n_uniq.fasta'    -evalue 1000  -out  $base'_out.txt'

grep  '^>id_\|^Query \|^Sbjct'  $base'_out.txt'  > 'replace_fragement_identify_point_'$basex'/output_file_f.txt'


done

