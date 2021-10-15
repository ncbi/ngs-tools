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

# $1 - path to the tools (loaders and vdb-dump)
# $2 - path to the schema (e.g. ncbi-vdb/interfaces)
# $3 - work directory for temporary files (removed if successful)
# $4 - command line for loader-1
# $5 - output archive of loader-1
# $6 - command line for loader-2
# $7 - output archive of loader-2
#
# return codes:
# 0 - passed
# 1 - coud not create temp dir
# 2 - unexpected return code from loader-1
# 3 - vdb-dump failed on the output of loader-1
# 4 - unexpected return code from loader-2
# 5 - vdb-dump failed on the output of loader-2
# 6 - outputs differ

SRA_BINDIR=$1
VDB_INCDIR=$2
TEMPDIR=$3
CMDLINE1=$4
OUTPUT1=$5
CMDLINE2=$6
OUTPUT2=$7

DUMP="$SRA_BINDIR/vdb-dump"
LOAD1="$SRA_BINDIR/$CMDLINE1"
LOAD2="$SRA_BINDIR/$CMDLINE2"

if [ "$(uname)" == "Darwin" ]; then
    DIFF="diff"
else
    DIFF="diff -Z"
fi

mkdir -p $TEMPDIR
rm -rf $TEMPDIR/*
if [ "$?" != "0" ] ; then
    exit 1
fi

echo "/vdb/schema/paths=\"${VDB_INCDIR}\"" >$TEMPDIR/tmp.kfg
export VDB_CONFIG=$TEMPDIR
#${SRA_BINDIR}/vdb-config --no-user-settings  /vdb/schema/paths

#CMD="$LOAD $CMDLINE -o $TEMPDIR/obj --no-user-settings 1>$TEMPDIR/load.stdout 2>$TEMPDIR/load.stderr"
CMD="$LOAD1 1>$TEMPDIR/load1.stdout 2>$TEMPDIR/load1.stderr"
echo $CMD
eval $CMD
rc="$?"
if [ "$rc" != "0" ] ; then
    echo "$LOAD1 returned $rc"
    echo "command executed:"
    echo $CMD
    cat $TEMPDIR/load1.stderr
    exit 2
fi

CMD="$DUMP $OUTPUT1 >$TEMPDIR/dump1.stdout 2>$TEMPDIR/dump1.stderr"
eval $CMD
rc="$?"
if [ "$rc" != "0" ] ; then
    echo "command executed:"
    echo $CMD
    cat $TEMPDIR/dump1.stderr
    exit 3
fi

CMD="$LOAD2 1>$TEMPDIR/load2.stdout 2>$TEMPDIR/load2.stderr"
echo $CMD
eval $CMD
rc="$?"
if [ "$rc" != "0" ] ; then
    echo "$LOAD2 returned $rc"
    echo "command executed:"
    echo $CMD
    cat $TEMPDIR/load2.stderr
    exit 4
fi

CMD="$DUMP $OUTPUT2 1>$TEMPDIR/dump2.stdout 2>$TEMPDIR/dump2.stderr"
eval $CMD
rc="$?"
if [ "$rc" != "0" ] ; then
    echo "command executed:"
    echo $CMD
    cat $TEMPDIR/dump2.stderr
    exit 5
fi

$DIFF $TEMPDIR/dump1.stdout $TEMPDIR/dump2.stdout >$TEMPDIR/diff
rc="$?"
if [ "$rc" != "0" ] ; then
    echo "outputs differ:"
    cat $TEMPDIR/diff
    exit 5
fi

rm -rf $TEMPDIR/* $OUTPUT1 $OUTPUT2

exit 0
