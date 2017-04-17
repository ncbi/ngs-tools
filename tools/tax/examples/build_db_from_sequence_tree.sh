set -e
bin_dir=../bin
name=example
du -h -a -b . | grep "sequence_tree.*\.fasta$" | sort -n -r > ./files.list.$name
$bin_dir/build_index ./files.list.$name ./example_data/tax.parents 1000 32 > ./$name.kmers 2>./$name.kmers.log
$bin_dir/check_index ./files.list.$name ./example_data/tax.parents ./$name.kmers > ./$name.kmers.checked 2> ./$name.kmers.checked.log
$bin_dir/filter_db ./$name.kmers.checked > ./$name.kmers.checked.clean 2> ./$name.kmers.checked.clean.log
$bin_dir/db_fasta_to_bin ./$name.kmers.checked.clean ./$name.dbs
$bin_dir/sort_dbs ./$name.dbs ./$name.dbss
