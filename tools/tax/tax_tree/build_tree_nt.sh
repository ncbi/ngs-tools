#fastacmd -d refseq_genomic -p F -D 1 | ~shutovo/blast/tax_tree/build_tree.py $1 2>build_tree.err
blastdbcmd -db nt -entry all | ~shutovo/blast/tax_tree/build_tree.py $1 2>build_tree_nt.err