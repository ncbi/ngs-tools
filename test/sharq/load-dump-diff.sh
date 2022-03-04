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
# Load an archive, vdb-dump and diff the output vs expected
#

# $1 - path to the tools (loader, and vdb-dump)
# $2 - work directory for temporary files (removed if successful)
# $3 - command line for loader
# $4 - output archive
# $5 - expected output from vdb-dump
#
# return codes:
# 0 - passed
# 1 - coud not create temp dir
# 2 - unexpected return code from loader
# 3 - vdb-dump failed on the output of loader
# 6 - output differs from expected

SRA_BINDIR=$1
TEMPDIR=$2
CMDLOAD=$3
OUTPUT=$4

DUMP="$SRA_BINDIR/vdb-dump -f json"
LOAD="$SRA_BINDIR/$CMDLOAD"

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

# Load
CMD="$LOAD 1>$TEMPDIR/load.stdout 2>$TEMPDIR/load.stderr"
eval $CMD
rc="$?"
if [ "$rc" != "0" ] ; then
    echo "$LOAD returned $rc"
    echo "command executed:"
    echo $CMD
    cat $TEMPDIR/load.stderr
    exit 2
fi

# Dump
CMD="$DUMP $OUTPUT >$TEMPDIR/dump.stdout 2>$TEMPDIR/dump.stderr"
eval $CMD
rc="$?"
if [ "$rc" != "0" ] ; then
    echo "command executed:"
    echo $CMD
    cat $TEMPDIR/dump.stderr
    exit 3
fi

# Diff
CMD="$DIFF expected/$OUTPUT.dump.stdout $TEMPDIR/dump.stdout >$TEMPDIR/diff"
eval $CMD
rc="$?"
if [ "$rc" != "0" ] ; then
    echo $CMD
    echo "outputs differ:"
    cat $TEMPDIR/diff
    exit 6
fi

rm -rf $TEMPDIR/* $OUTPUT $OUTPUT2

exit 0
