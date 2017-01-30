#/bin/bash
set -e # exit on error
make clean # force rebuild on linker depends, just in case
make -j32
make check -j32
git push
cp aligns_to_dbss.py aligns_to /panfs/traces01.be-md.ncbi.nlm.nih.gov/trace_software/tax_analysis/
echo copied
