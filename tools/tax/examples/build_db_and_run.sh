set -e
bin_dir=../bin
# To build the example database we will be using example files from the ./sequence_tree/ folder
# Each file has to have <tax id>.fasta name. The folder structure is irrelevant
du -h -a -b . | grep "sequence_tree.*\.fasta$" | sort -n -r > ./files.list.example

# The taxonomy structure is sets of pairs (tax id, parent tax id) in the ./example_data/tax.parents text file
# We usually generate it from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

KMER_LEN=32
DENSE_WINDOW=4 # 1 kmer of 4 for dense db (just for example)
SPARSE_WINDOW=128 # 1 kmer of 128 for sparse db (just for example)

# build_index_of_each_file will create <filename>.dense.db for every file we want to index
# in the real life we split files.list into ~100 parts and run multiple build_index_of_each_file processes in the cloud
# Please not that build_index_of_each_file does not overwrite existing .db files so you have to manually delete them to rebuild everything from scratch

$bin_dir/build_index_of_each_file ./files.list.example $DENSE_WINDOW $KMER_LEN .dense.db >./build_index_of_each_file.dense.log
$bin_dir/build_index_of_each_file ./files.list.example $SPARSE_WINDOW $KMER_LEN .sparse.db >./build_index_of_each_file.sparse.log

# then we merge all results into a single file
# We use a basic hashtable with predefined size to fit into 128G ram for the largest production database
MAX_KMER_DICTIONARY_SIZE=5000000 # This number should be roughly as max kmers expected * 2. 
$bin_dir/merge_db ./files.list.example .dense.db ./example.dense.db $MAX_KMER_DICTIONARY_SIZE >./merge.dense.log
$bin_dir/merge_db ./files.list.example .sparse.db ./example.sparse.db $MAX_KMER_DICTIONARY_SIZE >./merge.sparse.log

# Then we need to go through the data one more time and identify tax id for every kmer and store in the .db.tax_ids file
# This process can be also run in multiple instances and then the tool called bin/merge_kingdoms can be used to combine results into a single file
$bin_dir/identify_tax_ids ./files.list.example ./example_data/tax.parents ./example.dense.db ./example.dense.db.tax_ids 2>identify_tax_ids.dense.log
$bin_dir/identify_tax_ids ./files.list.example ./example_data/tax.parents ./example.sparse.db ./example.sparse.db.tax_ids 2>identify_tax_ids.sparse.log

# Combine .db and .db.tax_id file into a single .dbs file
$bin_dir/db_tax_id_to_dbs ./example.dense.db ./example.dense.db.tax_ids ./example.dense.dbs
$bin_dir/db_tax_id_to_dbs ./example.sparse.db ./example.sparse.db.tax_ids ./example.sparse.dbs

# Sort dense db by tax id for 2 step processing
$bin_dir/sort_dbs ./example.dense.dbs ./example.dense.dbss

# We can try to analyze some short read fasta file, like one from SRR4841604
# This is an example of straightforward 1 step processing using full dense database:
$bin_dir/aligns_to -dbs ./example.dense.dbs ./example_data/SRR4841604.fasta > ./SRR4841604.fasta.hits

# This is also going to work: (useful when unpacking data on the fly)
# cat ./example_data/SRR4841604.fasta | $bin_dir/aligns_to -dbs ./example.dense.dbs stdin > ./SRR4841604.fasta.hits

# Create a xml report and display it in some human readable format:
$bin_dir/hits_to_report.sh ./SRR4841604.fasta.hits

# Or we can do 2 step processing, using sparse database first to identify tax id diversity
# then use only parts of the dense database to do the rest. 
$bin_dir/aligns_to -dbs ./example.sparse.dbs ./example_data/SRR4841604.fasta > ./SRR4841604.fasta.1ststep.hits
python $bin_dir/hits_to_tax_list.py ./SRR4841604.fasta.1ststep.hits > ./tax_list
$bin_dir/aligns_to -dbss ./example.dense.dbss -tax_list ./tax_list ./example_data/SRR4841604.fasta > ./SRR4841604.fasta.2ndstep.hits

# Create a xml report and display it in some human readable format:
$bin_dir/hits_to_report.sh ./SRR4841604.fasta.2ndstep.hits

