#!/bin/bash
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#echo "$0 $*"

#
# A generic script to run a command line and diff the output against previously saved
#
# $1 - the executable
# $2 - test case ID 
# $3 - work directory (expected results under expected/, actual results and temporaries created under actual/)
#       expected/$2.stdout is expected to exists and will be diffed against the sdtout of the run
#       if expected/$2.stderr exists, it is diffed against the sdterr of the run
# $4 - expected return code
# $5, $6, ... - command line options for the executable
#
# return codes:
# 0 - passed
# 1 - could not create actual/$2/
# 2 - unexpected return code from the executable
# 3 - stdouts differ
# 4 - stderrs differ

EXE=$1
CASEID=$2
WORKDIR=$3
RC=$4
shift 4
ARGS=$*

TEMPDIR=$WORKDIR/actual
EXPECTED_STDOUT=$WORKDIR/expected/$CASEID.stdout
EXPECTED_STDERR=$WORKDIR/expected/$CASEID.stderr
ACTUAL_STDOUT=$TEMPDIR/$CASEID.stdout
ACTUAL_STDERR=$TEMPDIR/$CASEID.stderr

echo -n "running test case $CASEID... "

mkdir -p $TEMPDIR
if [ "$?" != "0" ] ; then
    echo "cannot create "
    exit 1
fi
rm -rf $TEMPDIR/*

CMD="$EXE $ARGS 1>$ACTUAL_STDOUT 2>$ACTUAL_STDERR"
#echo $CMD
eval $CMD
rc="$?"
if [ "$rc" != "$RC" ] ; then
    echo "$EXE returned $rc, expected $RC"
    echo "command executed:"
    echo $CMD
    cat $ACTUAL_STDOUT
    cat $TEMPDIR/$CASEID.stderr
    exit 2
fi

diff $EXPECTED_STDOUT $ACTUAL_STDOUT >$TEMPDIR/$CASEID.stdout.diff
rc="$?"
if [ "$rc" != "0" ] ; then
    cat $TEMPDIR/$CASEID.stdout.diff
    echo "command executed:"
    echo $CMD
    exit 4
fi    

if [ -f $EXPECTED_STDERR ] 
    then
    # clean up stderr:
    #
    # remove timestamps
    sed -i -e 's/^....-..-..T..:..:.. //g' $ACTUAL_STDERR
    # remove pathnames
    sed -i -e 's=/.*/==g' $ACTUAL_STDERR
    # remove source locations
    sed -i -e 's=: .*:[0-9]*:[^ ]*:=:=g' $ACTUAL_STDERR
    # remove version number if present
    sed -i -e 's=$(basename $EXE)\(\.[0-9]*\)*=$(basename $EXE)=g' $ACTUAL_STDERR
    #

    diff $EXPECTED_STDERR $ACTUAL_STDERR >$TEMPDIR/$CASEID.stderr.diff
    rc="$?"
    if [ "$rc" != "0" ] ; then
        cat $TEMPDIR/$CASEID.stderr.diff
        echo "command executed:"
        echo $CMD
        exit 4
    fi    
fi

rm -rf $TEMPDIR

echo "passed"

exit 0
