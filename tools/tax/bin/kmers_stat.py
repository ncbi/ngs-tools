#!/opt/python-2.7/bin/python
import sys
from collections import Counter
import sqlite3

def print_stat(tax_id, tax_name, total_kmer_for_tax_id, kmers):
    weight = get_weight(kmers)
    kmers = [(kmer, kmers[kmer]) for kmer in kmers]
    kmers.sort(key = lambda(kmer, count) : - count)
    print tax_id,
    if tax_name != None:
        print tax_name,

    kmers_used = len(kmers)
    print "\t", weight, "kmers", kmers_used, "unique kmers",
    if total_kmer_for_tax_id != None:
        print "of", total_kmer_for_tax_id, "(", int( 0.5 + 100.0 * kmers_used / total_kmer_for_tax_id ), "% )",
    print ""

    for kmer, count in kmers:
        print kmer, '\t', count    
    print ""

def load_tax_names(filename):
    tax_names = dict()
    conn = sqlite3.connect(filename)
    cur = conn.cursor()
    cur.execute('select tax_id, scientific_name from taxons')
    for tax_id, scientific_name in cur:
        tax_names[int(tax_id)] = scientific_name

    return tax_names    

def load_annotation(filename):
    f = open(filename)
    ann = dict()
    for line in f:
        line = line.split()
        ann[int(line[0])] = int(line[1])

    return ann         

def get_weight(kmers):
    return sum([kmers[kmer] for kmer in kmers])

def main():
    if len(sys.argv) < 2:
        print >> sys.stderr, "need <.hits file> [dbss .annotation] [gettax.sqlite for naming]"
        return

    annotation = dict()
    if len(sys.argv) >= 3:
        annotation = load_annotation(sys.argv[2])

    tax_names = dict()        
    if len(sys.argv) >= 4:
        tax_names = load_tax_names(sys.argv[3])

    kmers = dict()

    f = open(sys.argv[1])
    for line in f:
        line = line.rstrip()
        if not line:
            break
        line = line.split()
        kmer = line[0]
        tax_id = int(line[1])

        if not tax_id in kmers:
            kmers[tax_id] = Counter()

        kmers[tax_id][kmer] += 1;        

    kmers = [ (tax_id, kmers[tax_id], get_weight(kmers[tax_id])) for tax_id in kmers ]
    kmers.sort(key = lambda(tax_id, kmers, weight) : - weight)

    for tax_id, kmer_set, weight in kmers:
        total_kmer_for_tax_id = None
        if tax_id in annotation:
            total_kmer_for_tax_id = annotation[tax_id]

        tax_id_name = None
        if tax_id in tax_names:
            tax_id_name = tax_names[tax_id]

        print_stat(tax_id, tax_id_name, total_kmer_for_tax_id, kmer_set)

if __name__ == "__main__":
    main()
