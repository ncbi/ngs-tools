/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

SparkGenes is a tool that counts how many alignments of a given SRA-accession
intersect with each feature in a GTF-file.

The tool needs at least:
(1) a GTF-file
(2) a SRA-accession

The tool is compiled into the file "genes.jar".

The external dependencies of the tool are:
(a) spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar
(b) libncbi-vdb.so
(c) libngs-sdk.so


This is an example shell-script to execute this task:

count.sh:
---------------------------------------------------------------------------------------------
#!/bin/bash

SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
TOOLJAR="./genes.jar"
TOOL="SparkGenes"
GTF_SRC="-f human.gtf"
SRA_ACC="-a SRR1237963"
VERBOSE="-p"

java -cp $SPARKJAR:$TOOLJAR $TOOL $GTF_SRC $SRA_ACC $VERBOSE
---------------------------------------------------------------------------------------------

To execute this script you have to be in the directory where the file "genes.jar" is located.
The 2 libraries "libncbi-vdb.so" and "libngs-sdk.so" have to be in this directory too.
And you need a correct installation of the Java-Runtime and Apache-Spark.

The tool will print it's progress on the terminal. After it is done, a file SRR1237963.counts.txt
will be created in the current directory.

This command makes use of some default-values.
(1) The run-mode is per default local/multi-threaded, without Spark
(2) It uses 8 threads
(3) It filters the GTF-file for "exons".
(4) It assembles a feature by the Id "FeatureID".
(5) It does not count unaligned reads
(6) It uses the count-mode "UNION"


Useful helper-functions:


prepare.sh:
---------------------------------------------------------------------------------------------
#!/bin/bash

SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
TOOLJAR="./genes.jar"
TOOL="SparkGenes"
GTF_SRC="-f human.gtf"
OUTFILE="-o filtered.gtf"
MODE="-r PRE_GTF"
VERBOSE="-p"

java -cp $SPARKJAR:$TOOLJAR $TOOL $GTF_SRC $OUTFILE $MODE $VERBOSE
---------------------------------------------------------------------------------------------

This script processes the original GTF-file and creates a filtered gtf-file.
If you use this file (filtered.gtf) as input instead of the original GTF-file all further
scripts will execute faster. The filtered file will be much smaller.



compare_refs.sh:
---------------------------------------------------------------------------------------------
#!/bin/bash

SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
TOOLJAR="./genes.jar"
TOOL="SparkGenes"
GTF_SRC="-f filtered.gtf"
SRA_ACC="-a SRR1237963"
MODE="-r REFS"
VERBOSE="-p"

java -cp $SPARKJAR:$TOOLJAR $TOOL $GTF_SRC $SRA_ACC $MODE $VERBOSE
---------------------------------------------------------------------------------------------

This script compares the names of the references used in the filtered or original gtf-file and the
SRA-accession. It gives you an overview about which, if any, reference-names match between
GTF-file and Accession - and which references do not match. This script executes fast and helps
you avoiding the counting of a SRA-accession versus a GTF-file that do not match.
In the given example the GTF-file is human but the SRA-accession is from a mouse.
Even if the species do match, there can be problems with the names of the references. There is no
standard how to name a chromosome. "chr1", "c1", "1", "chromosome1" "I" etc. are examples for that.
The chromosome-names inside a SRA-accession are the ones the submitter has chosen.

For our example of SRR1237963 (mouse), you can for instance download a GTF-file from
"ftp://ftp.ensembl.org/pub/release-79/gtf/mus_musculus". If you run the compare_refs.sh script of
this gtf-file vs the SRA-accession you will notice that none of the references match, because the
creators of this GTF-file choose to uses "1" for the first chromosome and the submitter of the
SRA used "chr1". You cannot changed the SRA-accession. You can try to fix it by editing the GTF-file,
it is just a text-file after all (about 700 MByte !). But this tool offers a alternative: the
translation option. Let us create a simple text-file:

from_1_to_chr1.txt:
---------------------------------------------------------------------------------------------
1=chr1
2=chr2
3=chr3
...
19=chr19
X=chrX
Y=chrY
MT=chrM
---------------------------------------------------------------------------------------------

You are writing this file to translate from GTF into Accession: on the left the reference-
name used in the GTF-file, on the right the one used in the Accession.

Now you can repeat the comparison with the translation-file inserted like this:

compare_refs_trans.sh:
---------------------------------------------------------------------------------------------
#!/bin/bash

SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
TOOLJAR="./genes.jar"
TOOL="SparkGenes"
GTF_SRC="-f mouse.gtf"
SRA_ACC="-a SRR1237963"
TRANSLATION="-x from_1_to_chr1.txt"
MODE="-r REFS"
VERBOSE="-p"

java -cp $SPARKJAR:$TOOLJAR $TOOL $GTF_SRC $SRA_ACC $TRANSLATION $MODE $VERBOSE
---------------------------------------------------------------------------------------------


Now the main chromosomes match up, and it actually makes sense to perform the counting!
But before you do that you can speed up further scripts by repeating the prepare-script
with the translation inserted and against the mouse-GTF.

prepare_trans.sh:
---------------------------------------------------------------------------------------------
#!/bin/bash

SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
TOOLJAR="./genes.jar"
TOOL="SparkGenes"
GTF_SRC="-f mouse.gtf"
OUTFILE="-o filtered.gtf"
TRANSLATION="-x from_1_to_chr1.txt"
MODE="-r REFS"
VERBOSE="-p"

java -cp $SPARKJAR:$TOOLJAR $TOOL $GTF_SRC $OUTFILE $TRANSLATION $MODE $VERBOSE
---------------------------------------------------------------------------------------------


Now we can start to do the counting! The tool offers 3 different 'run-modes'.

(1)	LOCAL_S ... local, sequentially
	does not use Spark, only 1 thread, slow, but low memory footprint
	
(2)	LOCAL_P ... local, parallel
	does not use Spark, uses multiple threads, fast, but needs lots of memory
	
(3) SPARKED ... local or remote, parallel
	does use Spark, can run locally for testing, or remote on nodes for scaling

Notice that in the following examples there is no mention of the translation-file because we
are using the prepared file "filtered.gtf", it is already "translated".

run_mode_local_s.sh
---------------------------------------------------------------------------------------------
#!/bin/bash

SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
TOOLJAR="./genes.jar"
TOOL="SparkGenes"
GTF_SRC="-f filtered.gtf"
SRA_ACC="-a SRR1237963"
OUTFILE="-o local_s_counters.txt"
MODE="-r LOCAL_S"
VERBOSE="-p"

java -cp $SPARKJAR:$TOOLJAR $TOOL $GTF_SRC $SRA_ACC $OUTFILE $MODE $VERBOSE
---------------------------------------------------------------------------------------------


run_mode_local_p.sh
---------------------------------------------------------------------------------------------
#!/bin/bash

SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
TOOLJAR="./genes.jar"
TOOL="SparkGenes"
GTF_SRC="-f filtered.gtf"
SRA_ACC="-a SRR1237963"
OUTFILE="-o local_p_counters.txt"
MODE="-r LOCAL_P"
THREADS="-s 8"
VERBOSE="-p"

java -cp $SPARKJAR:$TOOLJAR $TOOL $GTF_SRC $SRA_ACC $OUTFILE $MODE $THREADS $VERBOSE
---------------------------------------------------------------------------------------------


run_mode_sparked_local.sh
---------------------------------------------------------------------------------------------
#!/bin/bash

SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
TOOLJAR="./genes.jar"
TOOL="SparkGenes"
GTF_SRC="-f filtered.gtf"
SRA_ACC="-a SRR1237963"
OUTFILE="-o local_x_counters.txt"
MODE="-r SPARKED -d local[16]"
SLICES="-s 16"
VERBOSE="-p"

java -cp $SPARKJAR:$TOOLJAR $TOOL $GTF_SRC $SRA_ACC $OUTFILE $MODE $SLICES $VERBOSE
---------------------------------------------------------------------------------------------

This script lets the tool run in Spark-mode but locally in 16 simulated nodes. That will need
some memory. Yes that will be slower than the 16 parallel threads, but that is not the point.
This mode allows you to verify that everything works fine for spark. Especially that you do
not see any warnings in the exhaustive Spark-Log's on the screen. If Spark complains that
the data transmitted to the nodes is too big, you have to increase the number of slices.
( -s 16 ---> -s 32 ) These warnings do not stop Spark in local mode, but in remote-mode
this can lead to canceling the jobs.


Before we can execute the tool remotely via Spark we have to make sure to have access to
a cluster. There are many different way's to do that, read the Spark-documentation for details.
It also depends on your IT infrastructure how access is granted. The given example uses
yarn to access a hadoop cluster. After starting the cluster a job-id is obtained...

run_mode_sparked_remote.sh
---------------------------------------------------------------------------------------------
#!/bin/bash

TOOLJAR="./genes.jar"
TOOL="SparkGenes"
GTF_SRC="-f filtered.gtf"
SRA_ACC="-a SRR1237963"
OUTFILE="-o local_r_counters.txt"
LIBS="libncbi-vdb.so,libngs-sdk.so"
MODE="-r SPARKED"
SLICES="-s 16"
VERBOSE="-p"

ARGS="$GTF_SRC $SRA_ACC $OUTFILE $MODE $SLICES $VERBOSE"

JOBID="XXXXXX"

export YARN_CONF_DIR=/home/user/hadoop/conf.$JOBID
export HADOOP_CONF_DIR=$YARN_CONF_DIR
execute "spark-submit --class $TOOL --master yarn-client --files $LIBS $TOOLJAR $ARGS"
---------------------------------------------------------------------------------------------

For remote execution we do not call java directly. We have to call a script "spark-submit".
This script is told what class inside the jar to call, and what master to connect to. We also
have to tell the script to transmit the two ".so"-files to the nodes. The jar-file is followed
by the actual arguments. The more alignments the SRA-accession has the more slices you should
ask for. The given Accession has about 25 Million alignments, 16 slices work fine for that.
Increase the number of slices or increase the amount of memory for the worker-nodes if the
accession is bigger. How do you know how many alignments an accession has?

	"vdb-dump SRR1237963 -T PRIMARY_ALIGNMENT -r"

The way a job has to be submitted can change in detail as Spark is under development.

The last 4 scripts have produced an outputfile each:

local_s_counters.txt, local_p_counters.txt, local_x_counters.txt and local_r_counters.txt

All 4 files should have identical content.



More Options:


(A)	the overlap-count-mode:

	UNION 		-c UNION
	STRICT		-c STRICT
	NONEMPTY	-c NONEMPTY

	if the parameter is omitted, UNION is used as default-mode.
	
	For each reference-position of an alginment build a set of features at this position.
	Merge all of these "position-sets" into a new set.
	
	UNION 		... merge by using union
	STRICT		... merge by using intersection
	NONEMPTY	...	merge by using intersection of the non-empty position-sets
	
	If this new set contains only one feature, the alginment counts for this feature.
	If this new set contains more than one feature, the alignment counts as ambiguous.
	If this new set is empty, the alignment counts as no-feature.
	
	UNION is the most useful overlap-mode.
	
	

(B) ignore orientation

	-g
	
	If an alignment is in forward orientation and a feature has reverse orientation
	( or vice versa ), the alignment is not considered for this feature because of the
	mismatch of orientation. That is the default behavior. But you can turn that off and
	let the feature be considered for the alignment.
	
	

(C)	feature_id

	-i transcript_id
	
	When analyzing the GTF-file, this is the string used as a key to merge multiple lines
	of the GTF-file	into one feature. The default value is "gene_id", an optinal value
	could be "transcript-id".
	

(D) minimal mapping quality

	-m 20
	
	Alignments with a mapping-quality lower than this value are rejected. The default-value
	is zero, that means no alignments are rejected.
	
	

(E) feature-type

	-t gene
	
	When analyzing the GTF-file, this is the string used as a filter for the 3rd column.
	GTF-lines without it are filtered out. The default-value is "exon", alternative
	values could be "gene" or "transcript".
	
	

(F)	unaligned

	-u
	
	Without this option unaligned or paritally unaligned reads are ignored. If you want to
	add the count of these reads too, you have to supply this option. Of course they will
	not count against any feature - they have no alignment information after all, the will
	only increase the no-feature count.
	
