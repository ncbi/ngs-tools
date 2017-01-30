#!/usr/bin/python
import subprocess
import os
import stat
import glob
import sys
import shell

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <acc list file> [-skip_reads]"
		return

	script_dir = os.path.dirname(os.path.realpath(__file__))
	f = open(sys.argv[1])

	skip_reads = len(sys.argv) >= 3 and sys.argv[2] == "-skip_reads"

	for acc in f:
		acc = acc.rstrip()
		cmdline = script_dir + "/start_farm.py " + acc
		if skip_reads:
			cmdline += " ./" + acc + ".matches"

		print shell.execute_local(cmdline)
		
main()
