#!/opt/python-2.7/bin/python
#!/usr/bin/python
import os
import sys
import glob

SEPARATOR = '\t'

def get(field, text, default_value = ""):
	text = text.split(field)
	if len(text) < 2:
		return default_value

	return text[1].splitlines()[0].lstrip()

def add(report_file, field):
	report_file.write(field)
	report_file.write(SEPARATOR)

report = []

def add_report(text):
	if get("total time (ms)", text, None) == None:
		return False

	x = get("accession:", text), int(get("run spots:", text, "0"))/1000, int(get("total time (ms)", text, "0"))/60000, get("mem usage before build_contigs (G)", text), get("reported contamination contigs count", text), get("reported contamination contigs %", text), get("reported contigs count", text), get("reported sum %", text), get("16 real size:", text), get("16 optimized size:", text)  
	# get("reported contamination contigs count", text), 
	report.append(x)
	return True

#	add(report_file, get("accession:", text))
#	add(report_file, str(int(get("run spots:", text, "0"))/1000))
#	add(report_file, str(int(get("total time (ms)", text, "0"))/60000))
#	add(report_file, get("mem usage before build_contigs (G)", text))
#	add(report_file, get("reported contamination contigs count", text))
#	add(report_file, get("reported contigs count", text))
#	add(report_file, get("reported sum %", text))
#	report_file.write('\n')

def two_digits(x):
	return "%.2f" % float(x)

def main():
	if __name__ != "__main__":
		return

	file_mask = "*.log"
	if len(sys.argv) > 1:
		file_mask = sys.argv[1]

#	report_file = open("./report.txt", "w")
	print "accession" + SEPARATOR + "spots (K)" + SEPARATOR + "time (min)" + SEPARATOR + "mem_usage (G)" + SEPARATOR + "contamination contigs" + SEPARATOR + "%contamination" + SEPARATOR + "good contigs" + SEPARATOR + "sum % reported"

	bad_logs = 0
	reads_too_short = 0

	for filename in glob.glob(file_mask):
#		print filename
		text = open(filename).read()
		if "Reads are too short" in text:
			reads_too_short += 1
			continue

		if not add_report(text):
			bad_logs += 1

	report.sort(key = lambda x: x[1])
	for r in report:
		print r[0], SEPARATOR, r[1], SEPARATOR, r[2], SEPARATOR, r[3], SEPARATOR, r[4], SEPARATOR, two_digits(r[5]), SEPARATOR, r[6], SEPARATOR, two_digits(r[7]), SEPARATOR, r[8], SEPARATOR, r[9]

	print "too short reads", SEPARATOR, reads_too_short
	print "bad logs", SEPARATOR, bad_logs

main()
