set -e

script_dir="$(dirname "$(readlink -f "$(which "$0")")")"

hits=$1
sort $hits > $hits.sorted
if [ ! -f ./gettax.sqlite ]; then
	wget https://sra-download.ncbi.nlm.nih.gov/traces/sra_references/tax_analysis/gettax.sqlite
fi

_tmpdir=$(mktemp -d)
python3 -m venv $_tmpdir/venv
source $_tmpdir/venv/bin/activate
pip install -r $script_dir/requirements.txt
python3 $script_dir/tax_analysis_parser.py  -c ./gettax.sqlite $hits.sorted > $hits.sorted.report.xml
python3 $script_dir/tax_analysis_report.py $hits.sorted.report.xml
deactivate
