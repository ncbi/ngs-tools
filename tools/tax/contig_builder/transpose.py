#!/opt/python-2.7/bin/python
#!/usr/bin/python
import os
import sys

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print "need <coverage file>" # [-align]"
		return

	v = eval(open(sys.argv[1]).read())
	for x in v:
		print x
#		if isinstance(x, int):
#			print 100
#		else:
#			cov, sum_cov = x
#			print 100*cov/sum_cov

main()