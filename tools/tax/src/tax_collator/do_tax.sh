#!/bin/bash
DIR=$1
ACC=$1
MODE=$2
COMPACT=''
if [ "$3" == 'compact' ]; then
    COMPACT='-compact'
    
fi

[ -d "$DIR" ] || mkdir $DIR
./prefetch -L err --max-size 100000000 -O $DIR $ACC > $DIR/prefetch.log 2>&1
f=''

FD="../fasterq-dump --seq-defline >\$si --stdout --skip-technical"
times=''
if [ "$MODE" == "collate" ]; then

    times="$DIR/times_new"
    [ -e $times ] && rm $times    

    log="$DIR/log_new"
    [ -e $log ] && rm $log

    name="${ACC}_new.fasta"
    f="$DIR/$name"
    if [ ! -s "$f" ]; then
        cd $DIR
        /usr/bin/time -f 'fasterq   : %e\tmem: %M' -o ../$times --append $FD --fasta-unsorted $ACC > $name 2>> ../$log
        cd ..
        [ ! -s $f ] && echo "Failed to generate fasta" && exit 1
    fi
    
elif [ "$MODE" == "vc" ]; then

    times="$DIR/times_msh"
    [ -e $times ] && rm $times    

    log="$DIR/log_msh"
    [ -e $log ] && rm $log
    
    name="${ACC}_msh.fasta"
    f="$DIR/$name"
    if [ ! -s "$f" ]; then
        cd $DIR
        /usr/bin/time -f 'fasterq   : %e\tmem: %M' -o ../$times --append $FD --fasta-unsorted $ACC > $name 2>> ../$log
        cd ..
        [ ! -s $f ] && echo "Failed to generate fasta" && exit 1
    fi

else
    times="$DIR/times_old"
    [ -e $times ] && rm $times    

    log="$DIR/log_old"
    [ -e $log ] && rm $log
    
    name="${ACC}_old.fasta"
    f="$DIR/$name"
    if [ ! -s "$f" ]; then
        cd $DIR
        /usr/bin/time -f 'fasterq   : %e\tmem: %M' -o ../$times --append $FD --fasta $ACC > $name  2>> ../$log
        cd ..
        [ ! -s $f ] && echo "Failed to generate fasta" && exit 1
    fi
fi

hits="${f}.hits"
if [ ! -s "$hits" ]; then
    /usr/bin/time -f 'first hits: %e\tmem: %M' -o $times --append ./aligns_to -dbs ./tree_index.dbs -hide_counts -num_threads 16 $f > $hits 2>> $log
    [ ! -s $hits ] && echo "Failed to generate hits" && exit 1
fi    

tax_list="${f}.tax_list"
if [ ! -s "$tax_list" ]; then
    /usr/bin/time -f 'hits_tax  : %e\tmem: %M' -o $times --append python2 ./hits_to_tax_list.py $hits > $tax_list 2>> $log
    [ ! -s $tax_list ] && echo "Failed to generate tax_list" && exit 1
fi    

TA="./tax_analysis_parser.py -c ./gettax.sqlite"
if [ "$COMPACT" != '' ]; then 
    TA="$TA --compact"
fi

if [ "$MODE" == 'collate' ]; then
    /usr/bin/time -f 'aligns_to : %e\tmem: %M' -o $times --append ./aligns_to -dbss ./tree_filter.dbss -tax_list $tax_list $f -hide_counts -collate $COMPACT -num_threads 16 > "${f}.out" 2>> $log
    /usr/bin/time -f 'xml_build : %e\tmem: %M' -o $times --append python2 $TA --collated "${f}.out" > ${f}.xml 2>> $log
elif [ "$MODE" == "vc" ]; then
    /usr/bin/time -f 'aligns_to : %e\tmem: %M' -o $times --append ./aligns_to -dbss ./tree_filter.dbss -tax_list $tax_list $f -collate -vectorize $COMPACT -num_threads 16 > "${f}.out" 2>> $log
    export LD_PRELOAD=libmesh.so
    #export MALLOCSTATS=1 
    /usr/bin/time -f 'collator  : %e\tmem: %M' -o $times --append ./tax_collator $f $COMPACT -num_threads 16 > "${f}.out" 2>> $log   
    /usr/bin/time -f 'xml_build : %e\tmem: %M' -o $times --append python2 $TA --collated "${f}.out" > ${f}.xml 2>> $log

else
    /usr/bin/time -f 'aligns_to : %e\tmem: %M' -o $times --append ./aligns_to -dbss ./tree_filter.dbss -tax_list $tax_list $f -hide_counts $COMPACT -num_threads 16 > "${f}.out" 2>> $log
    /usr/bin/time -f 'sort      : %e\tmem: %M' -o $times --append sort "${f}.out" 2>> $log | /usr/bin/time -f 'xml_build : %e\tmem: %M' -o $times --append python2 $TA > ${f}.xml 2>> $log
fi

ALL_TIME=`(cut -f 1 $times | cut -d ':' -f 1 --complement)`
RUNTIME=0
for i in $ALL_TIME; do RUNTIME=`echo "$RUNTIME + $i"|bc`; done

ALL_MEM=`(cut -f 2 $times | cut -d ':' -f 1 --complement)`
RUN_MEM=0
for i in $ALL_MEM; do 
    RUN_MEM=$(( $i > $RUN_MEM ? $i : $RUN_MEM ))
done

TOTAL_MEM=0
for i in $ALL_MEM; do 
    TOTAL_MEM=`echo "$TOTAL_MEM + $i * 1024"|bc`
done


RUN_MEM=`echo "$RUN_MEM * 1024"|bc`; 
FMT_RUN_MEM=`(numfmt --to=si $RUN_MEM)`
FMT_TOTAL_MEM=`(numfmt --to=si $TOTAL_MEM)`
PREFIX=''
if [ "$MODE" == "collate" ]; then
    PREFIX='New'
elif [ "$MODE" == "vc" ]; then
    PREFIX='Msh'
else
    PREFIX='Old'
fi    
echo "$PREFIX Total time(s) : $RUNTIME" 
echo "$PREFIX Max memory    : $FMT_RUN_MEM" 
echo "$PREFIX Total memory  : $FMT_TOTAL_MEM" 
echo "-------------------------------------"

echo "$PREFIX Total time(s) : $RUNTIME" >> $DIR/totals
echo "$PREFIX Max memory    : $FMT_RUN_MEM" >> $DIR/totals
echo "$PREFIX Total memory  : $FMT_TOTAL_MEM" >> $DIR/totals
echo "-------------------------------------" >> $DIR/totals

