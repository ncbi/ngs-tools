import sys

tax_ids = set()

f = open(sys.argv[1])
for line in f:
    line = line.split('\t')[1].split()
    for hit in line:
        tax_id = int(hit.split('x')[0])
        tax_ids.add(tax_id)

for tax_id in tax_ids:
    print tax_id