#!/bin/bash
k=1
#1553610
while [ $k -le 28 ]; do
	../kmer_distribution.py ../rus.txt $k 100 >./rus"$k".txt
	let k=k+1
done

#./start.sh SRR1564836 &
