set -e
bin_dir=../bin
$bin_dir/aligns_to -dbs ./example_data/viruses.dbs ./example_data/SRR1553418.fasta > SRR1553418.fasta.hits
./hits_to_report.sh ./SRR1553418.fasta.hits
