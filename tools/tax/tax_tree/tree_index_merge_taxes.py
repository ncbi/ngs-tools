#!/usr/bin/python
import sys
import seq_transform
import taxonomy
import filter_db

def any_of(v):
	for x in v:
		return x

	return None

def min_transform_of(kmer):
	rev_compl = seq_transform.reverse_complement(kmer)
	if rev_compl < kmer:
		return rev_compl

	return kmer

KMER_LEN = 32

def parse_db_line(line):
	line = line.split()
	if len(line[0]) != KMER_LEN:
		print >> sys.stderr, "bad kmer", line[0]

	kmer = min_transform_of(line[0])
	taxes = [int(tax) for tax in line[1:]]

	return kmer, taxes

def load_db(filename):
	f = open(filename)
	db = dict()
	line_number = 0
	for line in f:
		kmer, taxes = parse_db_line(line)

		if not kmer in db:
			db[kmer] = taxes
		else:
			db[kmer] = list( set(list(db[kmer]) + taxes) )

#		if not kmer in db:
#			db[kmer] = set()

#		db[kmer].update(set(taxes))

		line_number += 1
		if line_number % 1024*128 == 0:
			print >> sys.stderr, line_number

	return db

def level_up(lineage):
	return lineage[:-1]

tax_id_by_lineage = dict()
lineage_by_tax_id = dict()

f_common_lineages_cache = None
LINEAGES_CACHE = "lineages.tax.cache"

def get_tax_id_of_lineage(lineage, bigger_lineage):
	if not lineage:
		print >> sys.stderr, "no lineage!"
		return None

	lineage_hashable = tuple(lineage)

	if lineage_hashable in tax_id_by_lineage:
		return tax_id_by_lineage[lineage_hashable]

	print >> sys.stderr, "getting tax of", lineage, bigger_lineage

	tax_id = taxonomy.get_tax_id_of_species(lineage[-1])
	if tax_id == None:
		print >> sys.stderr, "cannot get tax_id of", lineage, bigger_lineage
		if len(bigger_lineage) > len(lineage):
			tax_id = taxonomy.get_parent_tax_id_of_species(bigger_lineage[len(lineage)])
		else:
			print >> sys.stderr, "level up", lineage, bigger_lineage
			tax_id = get_tax_id_of_lineage(level_up(lineage), lineage)
	
	if tax_id == None:
		print >> sys.stderr, "completely cannot get tax_id of", lineage, bigger_lineage
		tax_id = get_tax_id_of_lineage(level_up(lineage), lineage)
#		raise "fail"

	tax_id_by_lineage[lineage_hashable] = tax_id	

	global f_common_lineages_cache
	if f_common_lineages_cache == None:
		f_common_lineages_cache = open(LINEAGES_CACHE, 'a')

	if tax_id != None: # todo: think ?
		f_common_lineages_cache.write(str(tax_id) + '\t' + '/'.join(lineage) + '\n')
	else:
		print >> sys.stderr, "tax_id is none!", lineage, bigger_lineage

	return tax_id

def consensus_tax(taxes):
	if len(taxes) >= 100:
		return None

	if len(taxes) == 0:
		raise "not taxes?" 

	if len(taxes) == 1:
		return any_of(taxes)

	lineages = [get_tax_lineage(tax) for tax in taxes]
	common_lineage = get_common_lineage(lineages)
	if not common_lineage:
		return None

	return get_tax_id_of_lineage(common_lineage, lineages[0])

def get_common_lineage(lineages):
	if not lineages:
		return None

	main_lineage = lineages[0]

	checked_common = []
	index = 0
	while index < len(main_lineage):
		common = checked_common + [ main_lineage[index] ]
		if not has_common_lineage(lineages, common):
			return checked_common
		
		checked_common = common
		index += 1

	return checked_common	

def has_common_lineage(lineages, common):
	if not lineages or not common:
		return False

	for l in lineages:
		if len(l) < len(common):
			return False

		if l[:len(common)]  != common:
			return False

	return True

#lineages = dict()
def get_tax_lineage(tax_id):
	if not tax_id in lineage_by_tax_id:
		lineage_by_tax_id[tax_id] = taxonomy.get_tax_lineage(tax_id)
		print >> sys.stderr, "unknown tax_id", tax_id

	return lineage_by_tax_id[tax_id]
	
def low_complexity(kmer):
	if len(kmer) != 32:
		print >> sys.stderr, "not implemented yet for kmer", kmer, len(kmer)
		raise "not implemented yet"
	return filter_db.predicted(kmer) >= 13 # todo: tune const. and it works just for kmer == 32

def load_tax_lineages(filename):
	lineage_by_tax_id = dict()
	tax_id_by_lineage = dict()
	f = open(filename)
	for line in f:
		line = line.split('/')[2:]
		tax_id = int(line[-1].split('.')[0])
		lineage = line[:-1]

		lineage_by_tax_id[tax_id] = lineage
		tax_id_by_lineage[tuple(lineage)] = tax_id

	return lineage_by_tax_id, tax_id_by_lineage

def load_tax_cache(filename, tax_id_by_lineage):
	f = open(filename)
	for line in f:
		line = line.split('\t')
		tax_id = int(line[0])
		lineage = line[1].rstrip().split('/')
#		print tax_id
#		print lineage
		tax_id_by_lineage[tuple(lineage)] = tax_id

line_number = 0
def process(kmer, taxes):
	global line_number

	if low_complexity(kmer):
		print >> sys.stderr, kmer, "low complexity"
		return
		
	tax_id = consensus_tax(taxes)
	if tax_id != None:
		print kmer, tax_id
	else:
		print >> sys.stderr, kmer, "no consensus"

	line_number += 1
	if line_number % 1024 == 0:
		print >> sys.stderr, "lines processed", line_number

def main():
	if __name__ != "__main__":
		return

	if len(sys.argv) < 3:
		print "need <tree index file> <files.list> [-line_by_line]"
		return

	line_by_line = False
	if len(sys.argv) == 4 and sys.argv[3] == "-line_by_line":
		line_by_line = True
		print >> sys.stderr, "processing line by line"

	global lineage_by_tax_id
	global tax_id_by_lineage
	lineage_by_tax_id, tax_id_by_lineage = load_tax_lineages(sys.argv[2])
	load_tax_cache(LINEAGES_CACHE, tax_id_by_lineage)

	if not line_by_line:
		db = load_db(sys.argv[1])
		print >> sys.stderr, "db len is", len(db)
		for kmer in db:
			process(kmer, db[kmer])

	else:
		f = open(sys.argv[1])
		for line in f:
			kmer, taxes = parse_db_line(line)
			process(kmer, taxes)

main()
