#!/usr/bin/python
import tree_index_merge_taxes
import taxonomy
import os

def equal(a, b):
	if a!=b:
		print str(a) + " != " + str(b)
	else:
		print "ok"

def main():
	equal(tree_index_merge_taxes.consensus_tax([363253, 1234378]), 29546)

	return
	line_30195 =  ['Eukaryota', 'Metazoa', 'Ecdysozoa', 'Arthropoda', 'Hexapoda', 'Insecta', 'Pterygota', 'Neoptera', 'Endopterygota', 'Hymenoptera', 'Apocrita', 'Aculeata', 'Apoidea', 'Apidae', 'Bombus', 'Bombus', 'Bombus terrestris']
	line_132113 = ['Eukaryota', 'Metazoa', 'Ecdysozoa', 'Arthropoda', 'Hexapoda', 'Insecta', 'Pterygota', 'Neoptera', 'Endopterygota', 'Hymenoptera', 'Apocrita', 'Aculeata', 'Apoidea', 'Apidae', 'Bombus', 'Pyrobombus', 'Bombus impatiens']
	line_7458   = ['Eukaryota', 'Metazoa', 'Ecdysozoa', 'Arthropoda', 'Hexapoda', 'Insecta', 'Pterygota', 'Neoptera', 'Endopterygota', 'Hymenoptera', 'Apocrita', 'Aculeata', 'Apoidea', 'Apidae']


	equal(taxonomy.get_tax_lineage(30195), line_30195)
	equal(taxonomy.get_tax_lineage(132113), line_132113)
	equal(taxonomy.get_tax_lineage(7458), line_7458)
	equal(tree_index_merge_taxes.get_common_lineage([line_30195, line_132113]), ['Eukaryota', 'Metazoa', 'Ecdysozoa', 'Arthropoda', 'Hexapoda', 'Insecta', 'Pterygota', 'Neoptera', 'Endopterygota', 'Hymenoptera', 'Apocrita', 'Aculeata', 'Apoidea', 'Apidae', 'Bombus'])
	equal(tree_index_merge_taxes.get_tax_id_of_lineage(['Eukaryota', 'Metazoa', 'Ecdysozoa', 'Arthropoda', 'Hexapoda', 'Insecta', 'Pterygota', 'Neoptera', 'Endopterygota', 'Hymenoptera', 'Apocrita', 'Aculeata', 'Apoidea', 'Apidae', 'Bombus'], line_30195), 7458)

	equal(tree_index_merge_taxes.consensus_tax([30195,	132113]), 7458)
	
	equal(tree_index_merge_taxes.consensus_tax([633,273123,502801,	748672,	1286084,	1286089]), 1649845)


	equal(tree_index_merge_taxes.get_tax_id_of_lineage(['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacteriales', 'Enterobacteriaceae', 'Salmonella', 'Salmonella enterica subsp. enterica serovar Bareilly str. 0'], ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacteriales', 'Enterobacteriaceae', 'Salmonella', 'Salmonella enterica subsp. enterica serovar Bareilly str. 0']), 590)
	equal(tree_index_merge_taxes.get_tax_id_of_lineage(['Bacteria', 'Firmicutes', 'Bacilli', 'Bacillales', 'Bacillaceae', 'Bacillus'], ['Bacteria', 'Firmicutes', 'Bacilli', 'Bacillales', 'Bacillaceae', 'Bacillus', 'Bacillus atrophaeus 1013-1']), 1452)

	line_9531 = ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini', 'Catarrhini', 'Cercopithecidae', 'Cercopithecinae', 'Cercocebus', 'Cercocebus atys']
	line_9545 = ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini', 'Catarrhini', 'Cercopithecidae', 'Cercopithecinae', 'Macaca', 'Macaca nemestrina'] 
	line_9837 = ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Laurasiatheria', 'Cetartiodactyla', 'Tylopoda', 'Camelidae', 'Camelus', 'Camelus bactrianus']
	line_9233 = ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Archelosauria', 'Archosauria', 'Dinosauria', 'Saurischia', 'Theropoda', 'Coelurosauria', 'Aves', 'Neognathae', 'Sphenisciformes', 'Spheniscidae', 'Aptenodytes', 'Aptenodytes forsteri']

	equal(tree_index_merge_taxes.get_common_lineage([line_9531, line_9545]), ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Euarchontoglires', 'Primates', 'Haplorrhini', 'Catarrhini', 'Cercopithecidae', 'Cercopithecinae'])
	equal(tree_index_merge_taxes.get_common_lineage([line_9531, line_9545, line_9233]), ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi'])
	equal(tree_index_merge_taxes.get_common_lineage([line_9531, line_9545, line_9837, line_9233]), ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi'])

	equal(tree_index_merge_taxes.consensus_tax([9233, 9531, 9545, 9837]), 117571)

	db = tree_index_merge_taxes.load_db("./test/testdb")
	equal(len(db), 3)
	equal(db["AAAAAAAAAAAAAAAAAAAAAAAAAAATCGTC"], [9233, 9531, 9837, 9545])
	equal(db["AAAAAAAAAAAAAAAAAAAAAAAAGTTGAAAC"], [9545])
	equal(db["ATTTAATTTAACATTGTAATATTTATTTTCTA"], [9233, 85066])

	equal(tree_index_merge_taxes.min_transform_of("TAAAAAAAAACTGGGG"), "CCCCAGTTTTTTTTTA")
	equal(tree_index_merge_taxes.min_transform_of("CCCCAGTTTTTTTTTA"), "CCCCAGTTTTTTTTTA")
	
main()