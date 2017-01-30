script_dir="$(dirname "$(readlink -f "$(which "$0")")")"

fastacmd -d nt -p F -D 2 > ./nt_gis.txt
fastacmd -d nt -i ./nt_gis.txt -T > ./nt_gis.tax
$script_dir/format_gi_tax_list.py ./nt_gis.tax > ./nt_tax_by_gi.txt
$script_dir/build_tree_nt.sh ./nt_tax_by_gi.txt

du -h -a -b . | grep "\.fasta" | sort -n -r > ./files.list

$script_dir/build_tax_id_parents.py ./files.list > ./tax.parents

