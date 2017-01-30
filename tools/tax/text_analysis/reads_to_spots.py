#!/usr/bin/python
import sys

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <reads file>"
		return

	f = open(sys.argv[1])

	prev_spot = 0
	for line in f:
		if line[0] == '>':
			spot = int(line.split()[0].split('.')[1])
		else:
			spot = int(line.split('.')[0])

		if spot != prev_spot:
			print spot
		prev_spot = spot
	
main()
