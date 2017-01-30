script_dir="$(dirname "$(readlink -f "$(which "$0")")")"

fastacmd -d refseq_genomic -p F -D 2 > ./refseq_gis.txt
fastacmd -d refseq_genomic -i ./refseq_gis.txt -T > ./refseq_gis.tax

$script_dir/format_gi_tax_list.py ./refseq_gis.tax > ./refseq_tax_by_gi.txt

$script_dir/build_tree_refseq.sh ./refseq_tax_by_gi.txt
du -h -a -b . | grep "\.fasta" | sort -n -r > ./files.list

#add some from nt
fastacmd -d nt -p F -D 2 > ./nt_gis.txt
$script_dir/remove_from_db.py nt_gis.txt refseq_gis.txt > nt_gis.clean.txt 2>nt_gis.removed.txt
fastacmd -d nt -i ./nt_gis.clean.txt -T > ./nt_gis.tax
$script_dir/format_gi_tax_list.py ./nt_gis.tax > nt_tax_by_gi.txt

cat ./files.list | grep "/Viruses/" > ./files.list.viruses
$script_dir/filter_tax_by_gi.py nt_tax_by_gi.txt files.list.viruses > nt_tax_by_gi.txt.viruses 
$script_dir/build_tree_nt.sh ./nt_tax_by_gi.txt.viruses

du -h -a -b . | grep "\.fasta" | sort -n -r > ./files.list

$script_dir/build_tax_id_parents.py ./files.list > ./tax.parents

