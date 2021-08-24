ngs-tools::Sequence Taxonomic Analysis Tool (STAT)
===

### Description 
* ### STAT is a collection of tools for building `k-mer` databases, and querying against them.

### Quick start guide
* ### Run quickbuild.sh to build only the files for minimum functionality.
* ### Go to ./examples folder and run build_db_and_run.sh shell script. 
* ### Read build_db_and_run.sh: it is richly annotated with comments about each step performed.
* ### If you need NGS library support, or want all the tools you can run ./build.sh from any folder. 

### For Reproducing results found in `STAT: A fast, scalable, MinHash-based k-mer tool to assess Sequence Read Archive next generation sequence submission`
* ### Run build.sh 
* ### Genereate results for each of the files (`accuracy_1.fasta, accuracy_2.fasta`)
    * #### NOTE: This test will require approximately 120GB available memory
    * #### `{path}/ngs-tools/tools/tax/bin/aligns_to -dbss 20200518_tree_filter.dbss -tax_list TaxID_file -out accuracy_1.hits accuracy_1.fasta`
* ### To see the accuracy counts
    * #### Go to `{path}/ngs-tools/tools/tax/bin/`
    * #### Create python virtual environment, e.g.
        * ####`python3 -m ./venv`
    * #### Activate python virtual environment
        * #### `source v./venv/bin/activate`
    * #### Add requirements to the environment `pip install -r requirements.txt`
        * #### `pip install -r requirements.txt`
    * #### For each result file execute stat_accuracy.py
        * #### `./stat_accuracy accuracy_1.hits excluded_taxids`
    * #### Output to stdout is 
        * #### First three columns for each read: its excluded taxid, read name,  resolved taxid
            * #####
                    ...
                    10804 taxid_10804.62 4098
                    10804 taxid_10804.324 4098
                    10804 taxid_10804.589 9528
                    10804 taxid_10804.699 9528
                    10804 taxid_10804.805 239935
                    11801 taxid_11801.10 10090
                    11801 taxid_11801.17 91061
                    ...
        * #### Last ten lines summarizes each category total for the file
            * ##### 
                    ...
                    Virus:species: 2032
                    Virus:genus: 3347
                    Virus:ancestor: 700
                    Virus:false positive: 233
                    Virus:false negative: 274
                    Bacteria:species: 11437
                    Bacteria:genus: 11156
                    Bacteria:ancestor: 14680
                    Bacteria:false positive: 0
                    Bacteria:false negative: 0
