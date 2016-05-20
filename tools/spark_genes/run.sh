#!/bin/bash

execute()
{
    echo "------------------------------------------------------"
    echo $1
    eval $1
    echo "."
}

SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
MAINJAR="./genes.jar"

MOUSE_GTF="$HOME/spark/spark_genes/mouse.gtf"

GTF_M="-f $MOUSE_GTF"
GTFPRE="-f preprocessed.txt"
ACC="-a SRR1273905"
UNION="-c UNION"
SLICES="-s 16"
VERBOSE="-p"
LOCAL_P="-r LOCAL_P"
SPARKED="-r SPARKED"
PRE_GTF="-r PRE_GTF"
COMPARE_REFS="-r REFS"
SPARK_LOCAL="-d local[16]"
TIMED="-e"
MINMAPQ="-m 20"
OUTPUT="-o features.txt"
TRANSLATE="-x from_1_to_chr1.txt"

ARGS="$GTF_M $ACC $UNION $SLICES $VERBOSE $TIMED $TRANSLATE"

run_local_without_spark()
{
	execute "java -cp $MAINJAR:$SPARKJAR SparkGenes $ARGS $LOCAL_P"
}

run_local_with_spark()
{
	execute "java -cp $MAINJAR:$SPARKJAR SparkGenes $ARGS $SPARKED $SPARK_LOCAL"
}

run_on_cluster()
{
	#submit the jar to spark via spark-submit for execution on the SGE-cluster
	# !!! before that can be done: cd into ~/hadoop, execute start.sh
	# !!! wait for the cluster to be ready via "tail -f ./hadoop.o<jobid>"
	# !!! after job has been run : qdel <jobid>
	JOBID="2148495"

	export YARN_CONF_DIR="/$HOME/hadoop/conf.$JOBID"
	export HADOOP_CONF_DIR=$YARN_CONF_DIR
	execute "spark-submit --class SparkGenes --master yarn-client $MAINJAR $ARGS $SPARKED"
}


run_local_without_spark
#run_local_with_spark
#run_on_cluster
