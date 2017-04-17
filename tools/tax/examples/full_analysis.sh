set -e
bin_dir=../bin
acc=SRR1553418

#todo: fetch dbs and dbss from public site
$bin_dir/aligns_to -dbs /panfs/pan1.be-md.ncbi.nlm.nih.gov/tax_analysis/tax_tree_13_feb/tree_index.dbs $acc > $acc.fast_hits 2>$acc.fast_hits.log
$bin_dir/get_tax_list_from_hits.py $acc.fast_hits > $acc.tax_list
$bin_dir/aligns_to -dbss /panfs/pan1.be-md.ncbi.nlm.nih.gov/tax_analysis/tax_tree_13_feb/tree_filter.dbss -tax_list $acc.tax_list $acc > $acc.full_hits 2>$acc.full_hits.log
./hits_to_report.sh $acc.full_hits > $acc.report
