#!/opt/python-2.7/bin/python
import sys
#import argparse

def median_of(v):
	v = list(v)
	v.sort()
	return v[len(v)/2]

def main():
	if len(sys.argv) < 2:
		print >> sys.stderr, "need <connectivity file>"
		return

	f = open(sys.argv[1])
	level = None		
	for line in f:
		line = line.split()
		line = [int(x) for x in line]
		if level == None:
			level = median_of(line)
			print >> sys.stderr, "level is ", level

		break_points = []
		for pos in xrange(len(line)):
			if line[pos] < 70 and line[pos] < level / 2 and pos < len(line) - level: # todo: probably -level*2/3
				break_points.append((pos, line[pos]))

		print len(break_points),
		for pos, conn in break_points:
			print pos, conn,

		print ""
			
if __name__ == "__main__":
	main()
