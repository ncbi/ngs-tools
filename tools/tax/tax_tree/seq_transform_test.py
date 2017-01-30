#!/usr/bin/python
import seq_transform

def equal(a, b):
	if a!=b:
		print str(a) + " != " + str(b)
	else:
		print "ok"

def main():
	equal(seq_transform.reverse_complement("ACTG"), "CAGT")
	equal(seq_transform.reverse_complement("AAATATTCATATAAAAGGAATCTCGGCCCTCT"), "AGAGGGCCGAGATTCCTTTTATATGAATATTT")
	equal(seq_transform.reverse_complement("TAAAAAAAAACTGGGG"), "CCCCAGTTTTTTTTTA")
	equal(seq_transform.reverse_complement("AAAAAAAAAAAAAAAAAAAAAAAAAAATCGTC"), "GACGATTTTTTTTTTTTTTTTTTTTTTTTTTT")
	equal(seq_transform.reverse_complement("TAGAAAATAAATATTACAATGTTAAATTAAAT"), "ATTTAATTTAACATTGTAATATTTATTTTCTA")
	
                                                                                 
main()