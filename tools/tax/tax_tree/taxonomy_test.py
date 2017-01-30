#!/usr/bin/python
import taxonomy
import os

def equal(a, b):
	if a!=b:
		print str(a) + " != " + str(b)
	else:
		print "ok"

def main():
	equal(taxonomy.get_tax_lineage(186536), ['Viruses', 'ssRNA viruses', 'ssRNA negative-strand viruses', 'Mononegavirales', 'Filoviridae', 'Ebolavirus'])
	equal(taxonomy.get_tax_id_of_species('Ebolavirus'), 186536)
	equal(taxonomy.get_parent_tax_id_of_species('Ebolavirus'), 11266)
	equal(taxonomy.get_taxonomy_of_gi(410427591, "nr").tax_id, 458326)
	equal(taxonomy.get_taxonomy_of_tax_id(428406).scientific_name, "Ralstonia pickettii 12D")


main()