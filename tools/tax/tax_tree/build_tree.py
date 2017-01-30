#!/opt/python-2.7/bin/python
#!/usr/bin/python
import os
import sys
import taxonomy

def load_tax_by_gi(filename):
	f = open(filename)
	tax_by_gi = dict()
	for line in f:
		line = line.split()
		gi = int(line[0])
		tax_id = int(line[1])
		tax_by_gi[gi] = tax_id

	return tax_by_gi
	
def make_folders(filename):
	directory = os.path.dirname(filename)
	if not os.path.exists(directory):
		os.makedirs(directory)	

def get_tree_node_filename(tax_id, lineage):
	name = "./tree/"
	for l in lineage:
		name += l + '/'

	return name + str(tax_id) + ".fasta"

lineages = dict()
def get_tax_lineage(tax_id):
	if not tax_id in lineages:
		lineages[tax_id] = taxonomy.get_tax_lineage(tax_id)

	return lineages[tax_id]

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 2:
		print >> sys.stderr, "need <tax by gi file>"
		print >> sys.stderr, "taxes input from stdin"
		return

	tax_by_gi = load_tax_by_gi(sys.argv[1])
#	detected_gis = set()
#	detected_taxes = set()

	f = None
	gi = None
	f_filename = None
	written = 0
	gis = set()
	for line in sys.stdin:
		if line[0] == '>':
			gi = int(line.split("gi|")[1].split('|')[0])
			if not gi in tax_by_gi:
				print >> sys.stderr, "tax for gi", gi, "not found"
				f = None
				f_filename = None
				continue

			gis.add(gi)
			tax_id = tax_by_gi[gi]
			lineage = get_tax_lineage(tax_id)
#			print >> sys.stderr, gi, tax_id, lineage
#			detected_gis.add(gi)
#			detected_taxes.add(tax_id)
			filename = get_tree_node_filename(tax_id, lineage)
#			print >> sys.stderr, filename
			print "gis:", len(gis), str(written/1000000) + "M written", filename
			if f_filename != filename:
				if f != None:
					f.close() # not necessary ?
					f = None

				make_folders(filename)
				f = open(filename, "a")
				f_filename = filename

			f.write(line)
		else:
			if f != None:			
				f.write(line)

		if f != None:
			written += len(line)
		
#	print >> sys.stderr, "gis: ", len(detected_gis)
#	print >> sys.stderr, "taxonomies: ", len(detected_taxes)


main()
