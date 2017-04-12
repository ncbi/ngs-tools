set -e
bin_dir=../bin/
$bin_dir/aligns_to -dbs ./example_data/viruses.dbs -print_counts ./example_data/SRR1553418.fasta > SRR1553418.fasta.hits
sort ./SRR1553418.fasta.hits > ./SRR1553418.fasta.hits.sorted
if [ ! -f ./taxdump.tar.gz ]; then
	wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
fi
python ../bin/tax_analysis_parser.py -t ./taxdump.tar.gz -c ./gettax.cache.sqlite ./SRR1553418.fasta.hits.sorted > ./SRR1553418.fasta.hits.sorted.report.xml
python ../bin/tax_analysis_report.py ./SRR1553418.fasta.hits.sorted.report.xml