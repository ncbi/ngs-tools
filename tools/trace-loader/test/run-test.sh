#!/bin/bash

## LYRICS: This script tests trace-loader program:
##
## 1) It will remove all content of directory "WD/tmp"
## 2) It will pass content of directory "WD/trace-data" as
##    input for trace-loader program and produce SRA db
## 3) It will dump content of SRA db to text file and 
##    compare it with content of file
##       "WD/expected/vdb-dump-output.txt"
##
## WD - is a working directory, which is passed as a first
## paramter to script

usage ()
{
    MSS="$@"
    if [ -n "$MSS" ]
    then
        echo >&2
        echo ERROR:$MSS>&2
    fi

    cat <<EOF >&2

This script will run trace-loader and general-loader to produce
SRA database. It will dump database content and compare it with
expected.

Syntax: $( basename $0 ) WDIR TRLD GNLD VDDP SCPT RULS

Where:
        WDIR - working directory
        TRLD - path to trace-loader programm to test
        GNLD - path to general-loader programm
        VDDP - path to vdb-dump program
        SCPT - path to directory which contains schemas
        RULS - path to directory which validation rules

EOF

    exit 1
}

###
##  Exclusives and promotions
#
err_msg ()
{
    MSS="$@"
    if [ -z "$MSS" ]
    then
        echo ERROR: unknown>&2
    else
        echo ERROR: $MSS>&2
    fi
}

warn_msg ()
{
    MSS="$@"
    if [ -z "$MSS" ]
    then
        echo WARNING: unknown>&2
    else
        echo WARNING: $MSS>&2
    fi
}

err_exit ()
{
    MSS="$@"
    if [ -n "$MSS" ]
    then
        echo "ERROR: $MSS">&2
    fi

    echo "ERROR: Exiting">&2

    exit 1
}

run_cmd ()
{
    RCM="$@"

    if [ -n "$RCM" ]
    then
        echo "## $RCM"
        eval $RCM
        if [ $? -ne 0 ]
        then
            err_exit command failed
        fi
    fi
}

###
##  Parameters and environment
#
if [ $# -ne 6 ]
then
    usage Invalid arguments
fi

canonical_path ()
{
    U=$(readlink -e $2)
    if [ $? -ne 0 ]
    then
        err_exit can not stat \'$2\'
    fi

    eval $1=$U
}

canonical_path WORK_DIR $1
canonical_path TRACE_LD $2
canonical_path GENERAL_LD $3
canonical_path VDB_DUMP $4
canonical_path SCHEMA_PATH $5
canonical_path RULES_PATH $6

###
##  Working dir
#
canonical_path TRACE_DIR $WORK_DIR/trace-data
canonical_path TMP_DIR $WORK_DIR/tmp
canonical_path EXPECTED_DIR $WORK_DIR/expected
canonical_path EXPECTED_FILE $EXPECTED_DIR/vdb-dump-output.txt
canonical_path TRACEINFO_FILE $TRACE_DIR/TRACEINFO.tbl
canonical_path SCHEMA_FILE $SCHEMA_PATH/ncbi/trace.vschema

run_cmd rm -rf $TMP_DIR/*

###
##  Running trace-loader
#
GW_PATH=${TMP_DIR}/gw.out
GL_PATH=${TMP_DIR}/gl.out
LOC_LOG="${TMP_DIR}/log.$(basename $TRACE_LD).$( date +%Y_%m_%d-%H:%M:%S ).xml"

CMD="$TRACE_LD --trace-info $TRACEINFO_FILE"
CMD="$CMD --output $GW_PATH"
CMD="$CMD --schema $SCHEMA_FILE"
CMD="$CMD -z $LOC_LOG"
CMD="$CMD --val-rules $RULES_PATH/rule-main.owp"
CMD="$CMD --log-level info"

run_cmd $CMD

###
##  Running general-loader
#
VDB_CFG=$TMP_DIR/vdb_config.kfg
cat <<EOF >$VDB_CFG
/vdb/schema/paths = "$SCHEMA_PATH"
/sra/quality_type = "raw_scores"
EOF

export NCBI_SETTINGS=/
export VDB_CONFIG=$VDB_CFG

run_cmd "$GENERAL_LD -T $GL_PATH <$GW_PATH"

###
##  Running vdb-dump
#
DUMP_FILE=$TMP_DIR/vdb-dump-output.txt
run_cmd "$VDB_DUMP $GL_PATH >$DUMP_FILE"

###
##  Running diff
#
run_cmd diff $DUMP_FILE $EXPECTED_FILE

###
##  Clearing mess
#
run_cmd rm -rf $TMP_DIR/*

echo TEST PASSED
