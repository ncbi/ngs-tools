~/blast/tax_tree/db_fasta_to_bin ./test/data/test_db ./test/tmp/test_db.dbs
./sort_dbs ./test/tmp/test_db.dbs ./test/tmp/test_db.dbss
diff ./test/tmp/test_db.dbss.annotation ./test/data/test_db.dbss.annotation

~/blast/tax_tree/db_fasta_to_bin ./test/data/test_db2 ./test/tmp/test_db2.dbs
./sort_dbs ./test/tmp/test_db2.dbs ./test/tmp/test_db2.dbss
diff ./test/tmp/test_db2.dbss.annotation ./test/data/test_db2.dbss.annotation

#~/blast/tax_tree/db_fasta_to_bin ./test/data/test_db2 ./test/tmp/test_db2.dbs

