#!/opt/python-2.7/bin/python
import sys
import argparse
import subprocess
import xml.etree.ElementTree as ET

def fetch_tax_analysis_xml(acc):
	p = subprocess.Popen(
		["sqsh-ms -S SRA_PRIMARY -D SRA_Main -U sra_sa -P sra_sa_pw -b -h -w 1000000000"],
		stdout = subprocess.PIPE, 
		stdin = subprocess.PIPE, 
		shell = True)

	input = "select result from SRA_Track..LoadResults(nolock) lr where name = '" + acc + "' and job_type = 'tax_analysis' and ok = 1 and result <> ''"
	input += "\ngo\nexit\n"
	out, err = p.communicate(input)
	if err:
		raise Exception("fetching sql data error: " + err)

	return out

def parse_tax_ids(xml):
	if not xml:
		return []
	root = ET.fromstring(xml)

	return [ ( int(taxon.attrib["tax_id"]), int(taxon.attrib["self_count"]) ) for taxon in root.iter("taxon") ]

def main():
	parser = argparse.ArgumentParser(description='fetch tax list of analyzed run')
	parser.add_argument('accession', help='accession')
	parser.add_argument('--min-count', metavar='N', type=int, help='cut items with less than given amount of hits', default = 0)

	args = parser.parse_args()
	xml = fetch_tax_analysis_xml(args.accession)
	ids = parse_tax_ids(xml)
	for tax, count in ids:
		if count >= args.min_count:
			print tax#, count

if __name__ == "__main__":
	main()
