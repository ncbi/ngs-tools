#/bin/bash
set -e # exit on error
make clean # force rebuild on linker depends, just in case
make -j32
make check -j32
git push
cp -v bin/aligns_to_dbss.py bin/aligns_to bin/gettax.py bin/tax_analysis_parser.py /panfs/traces01.be-md.ncbi.nlm.nih.gov/trace_software/tax_analysis/
echo copied
