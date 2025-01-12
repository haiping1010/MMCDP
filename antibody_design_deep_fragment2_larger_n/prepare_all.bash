
for name in all_cdr?_aa_heavy_n.txt  all_cdr?_aa_light_n.txt

do

	base=${name%.txt}
	sort $name | uniq | grep -v "\"\"" > $base'_uniq.txt'
        python  convert_fasta.py  $base'_uniq.txt'

        makeblastdb -in  $base'_uniq.fasta'   -dbtype prot





done
