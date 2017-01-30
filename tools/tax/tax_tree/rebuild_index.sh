script_dir="$(dirname "$(readlink -f "$(which "$0")")")"

$script_dir/search/build_index ./files.list ./tax.parents 1 32 > ./tree_index 2>./tree_index.log
$script_dir/search/check_index ./files.list ./tax.parents ./tree_index > ./tree_index.checked 2> ./tree_index.checked.log
$script_dir/filter_db ./tree_index.checked > ./tree_index.checked.clean 2> ./tree_index.checked.clean.log
$script_dir/db_fasta_to_bin ./tree_index.checked.clean ./tree_index.dbs

