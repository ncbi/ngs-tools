#!/bin/bash

./do_tax.sh $1 no-collate $2
if [ $? -ne 0 ]; then 
    cat $1/log_old
    exit 1
fi

./do_tax.sh $1 collate $2
if [ $? -ne 0 ]; then
    cat $1/log_new
    exit 1
fi    

#./getf.sh $1 vc
#if [ $? -ne 0 ]; then
#    cat $1/log_msh
#    exit 1
#fi    

#cat $1/totals

DIFF=$(./xmldiffs "$1/$1_new.fasta.xml" "$1/$1_old.fasta.xml")
if [ "$DIFF" == "" ] 
then
    echo "Checked New Old!"
else    
    echo $DIFF
fi

#DIFF=$(./xmldiffs "$1/$1_msh.fasta.xml" "$1/$1_old.fasta.xml")
#if [ "$DIFF" == "" ] 
#then
#    echo "Checked Msh Old!"
#else    
#    echo $DIFF
#fi
