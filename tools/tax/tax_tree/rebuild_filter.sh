script_dir="$(dirname "$(readlink -f "$(which "$0")")")"

$script_dir/search/build_index ./files.list ./tax.parents 1000 32 > ./tree_filter 2>./tree_filter.log
$script_dir/search/check_index ./files.list ./tax.parents ./tree_filter > ./tree_filter.checked 2> ./tree_filter.checked.log
$script_dir/filter_db ./tree_filter.checked > ./tree_filter.checked.clean 2> ./tree_filter.checked.clean.log
$script_dir/db_fasta_to_bin ./tree_filter.checked.clean ./tree_filter.checked.clean.dbs
$script_dir/sort_dbs ./tree_filter.checked.clean.dbs ./tree_filter.dbss

