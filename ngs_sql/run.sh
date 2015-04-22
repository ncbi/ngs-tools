#!/bin/bash


SPARKJAR="/usr/local/spark/1.2.0/lib/spark-assembly-1.2.0-hadoop2.5.0-cdh5.3.0.jar"
MAINJAR="./ngs_sql.jar"
MAINCLASS="NGS_SQL"
ARGS="sql.txt"

run_local()
{
	java -cp $MAINJAR:$SPARKJAR $MAINCLASS $ARGS
}

submit_remote()
{
	#submit the jar to spark via spark-submit for execution on the SGE-cluster
	# !!! before that can be done: cd into ~/hadoop, execute start.sh
	# !!! wait for the cluster to be ready via "tail -f ./hadoop.o<jobid>"
	# !!! after job has been run : qdel <jobid>
	JOBID="3561173"

	export YARN_CONF_DIR=/home/raetzw/hadoop/conf.3561173
	export HADOOP_CONF_DIR=$YARN_CONF_DIR
	spark-submit --class $MAINCLASS --master yarn-client $MAINJAR $ARGS
}

run_local
