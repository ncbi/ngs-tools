#!/usr/bin/python
import os
import shell
import sys

class Taxonomy:
	def __init__(self, tax_id, scientific_name):
		self.tax_id = tax_id
#		self.common_name = common_name
		self.scientific_name = scientific_name

	def equal(self, tax):
		return self.tax_id == tax.tax_id

UNKNOWN_TAXONOMY = Taxonomy(0, "UNKNOWN")

def get_taxonomy_of_gi(gi, database_name):
	cmdline = "/netopt/ncbi_tools64/bin/fastacmd -d " + database_name + " -T -s " + str(gi)
	out = ""
	try:
		out = shell.execute_local(cmdline)
	except Exception:
		return Taxonomy(0, "UNKNOWN")

#	print out
	return Taxonomy(parse_tax_id(out), parse_scientific_name(out))

def get_taxonomy_of_tax_id(tax_id):
	out = ""
	try:
		out = shell.execute_local("/netopt/ncbi_tools64/bin/gettax " + str(tax_id))
	except Exception:
		return Taxonomy(tax_id, "UNKNOWN")

#	print out
	return Taxonomy(tax_id, parse_scientific_name_tax_id(out))

def parse_scientific_name_tax_id(out):
	return out.split("scientific name: ")[1].split('\n')[0]

def parse_scientific_name(out):
	return out.split("Scientific name: ")[1].split('\n')[0]

def parse_tax_id(out):
	out = out.split("NCBI taxonomy id: ")[1]
	out = out.split('\n')
	return int(out[0])

def parse_common_name(out):
	return out.split("Common name: ")[1].split('\n')[0]

def get_tax_lineage(tax_id):
	if tax_id == 2:
		return ['Bacteria'] # todo: fix
	if tax_id == 2759:
		return ['Eukaryota']
	if tax_id == 10239:
		return ['Viruses']
	if tax_id == 2157:
		return ['Archaea']

	out = shell.execute_local("/netopt/ncbi_tools64/bin/gettax " + str(tax_id))
	return parse_lineages(out) + [parse_scientific_name_tax_id(out)]

def parse_lineages(out):
	out = out.split("lineage: ")
	if len(out) == 1:
		return []
	out = out[1].split('id_gc:')[0]	
	out = out.split(';')
	return [o.strip() for o in out]

def get_tax_id_of_species(tax):
	out = shell.execute_local('/netopt/ncbi_tools64/bin/gettax "' + tax + '"', False)
	return parse_lineage_tax_id(out)

def get_parent_tax_id_of_species(tax):
	out = shell.execute_local('/netopt/ncbi_tools64/bin/gettax "' + tax + '"', False)
	return parse_lineage_parent_tax_id(out)
	
def parse_lineage_tax_id(out_):
	out = out_.split("tax id: ")
	if len(out) < 2:
#		print >> sys.stderr, "bad tax out : ", out_
		return None
	out = out[1]
	out = out.split('\n')
	return int(out[0])
	
def parse_lineage_parent_tax_id(out_):
	out = out_.split("parent tax id: ")
	if len(out) < 2:
#		print >> sys.stderr, "bad tax out : ", out_
		return None
	out = out[1]
	out = out.split('\n')
	return int(out[0])

def main():
	if len(sys.argv) < 2:
		print >> sys.stderr, "need <filename with taxonomies>"
		return

	f = open(sys.argv[1])
	for line in f:
		if not line:
			break

		tax = int(line.split()[0])
		print line.rstrip() + '\t' + str(get_tax_lineage(tax))
	
if __name__ == "__main__":
	main()
