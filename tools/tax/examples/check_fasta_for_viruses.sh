set -e
bin_dir=../bin/
$bin_dir/aligns_to -dbs ./example_data/viruses.dbs -print_counts ./example_data/SRR1553418.fasta > SRR1553418.fasta.hits
