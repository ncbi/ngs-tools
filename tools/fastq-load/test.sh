export PYTHONPATH=`cd ../../shared/python; pwd`

R1=R1.10.fastq
R2=R2.10.fastq

# single file test
python3 fastq-load.py --output=foo ${R1} >/dev/null

# two file test
python3 fastq-load.py --output=foo --read1PairFiles=${R1} --read2PairFiles=${R2} ${R1} ${R2} >/dev/null
