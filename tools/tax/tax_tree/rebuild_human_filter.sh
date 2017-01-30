script_dir="$(dirname "$(readlink -f "$(which "$0")")")"

cat ./files.list | grep "/Homo sapiens/" > ./files.list.human
cat ./files.list | grep "/Mus musculus/" > ./files.list.mouse
cat ./files.list | grep -v "/Eukaryota/" > ./files.list.non_eukaryota

$script_dir/search/build_index ./files.list.human ./tax.parents 1000 32 > ./human_filter 2>./human_filter.log
$script_dir/search/check_index ./files.list.mouse ./tax.parents ./human_filter > ./human_filter.mouse_checked 2> ./human_filter.mouse_checked.log
$script_dir/search/check_index ./files.list.non_eukaryota ./tax.parents ./human_filter.mouse_checked > ./human_filter.mouse_checked.non_euk_checked 2> ./human_filter.mouse_checked.non_euk_checked.log
$script_dir/filter_db ./human_filter.mouse_checked.non_euk_checked -only_tax 9606 > ./human_filter.mouse_checked.non_euk_checked.clean 2> ./human_filter.mouse_checked.non_euk_checked.clean.log
$script_dir/db_fasta_to_bin ./human_filter.mouse_checked.non_euk_checked.clean ./human_filter.db

