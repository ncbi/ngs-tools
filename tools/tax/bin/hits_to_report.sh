set -e

script_dir="$(dirname "$(readlink -f "$(which "$0")")")"

hits=$1
sort $hits > $hits.sorted
if [ ! -f ./taxdump.tar.gz ]; then
	wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
fi
python2 $script_dir/tax_analysis_parser.py -t ./taxdump.tar.gz -c ./gettax.cache.sqlite $hits.sorted > $hits.sorted.report.xml
python2 $script_dir/tax_analysis_report.py $hits.sorted.report.xml