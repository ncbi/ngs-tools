#!/opt/python-2.7/bin/python
#!/usr/bin/python
import os
import sys
import perf

def weight(_seq):
	seq, cov = _seq
	return len(seq)*cov

def get_seqs(data):
	lines = data.split('>')[1:]
	lines = [l.split() for l in lines]

	return [(l[1], int(l[0])) for l in lines]

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <log file>" # [-align]"
		return

	seqs = get_seqs(open(sys.argv[1]).read())
	print len(seqs)
		
	seqs.sort(key = lambda seq: weight(seq), reverse = True)
	for seq, cov in seqs:
		print ">cov", cov, "weight", weight((seq, cov))
		print seq

main()
