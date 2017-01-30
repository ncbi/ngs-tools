#!/bin/bash
k=1
while [ $k -le 30 ]; do
	../../kmer_distribution.py ../../eng.txt $k >./eng"$k".txt
	let k=k+1
done

#./start.sh SRR1564836 &
