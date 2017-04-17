set -e
bin_dir=../bin
$bin_dir/aligns_to -dbs ./example.dbs -print_counts ./example_data/SRR4841604.fasta > SRR4841604.fasta.hits 2>SRR4841604.fasta.hits.log
echo "full database analysis"
./hits_to_report.sh ./SRR4841604.fasta.hits
echo "checking only for particular viruses from database"
$bin_dir/aligns_to -dbss ./example.dbss -tax_list ./example_data/some_particular_viruses.tax_list -print_counts ./example_data/SRR4841604.fasta > SRR4841604.fasta.particular_viruses.hits 2>SRR4841604.fasta.hits.log
./hits_to_report.sh ./SRR4841604.fasta.particular_viruses.hits
