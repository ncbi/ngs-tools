#!/bin/sh
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

HELP='\n
 CMake configuration script\n
\n
 $1 - source tree root directory\n
 $2 - output directory (will be created if necessary)\n
\n
 Options:\n
 -h|--help                 print usage\n 
 -b|--build Debug|Release  (default Debug) build type\n 
 -a|--args <string>        additional arguments to pass to CMake\n
\n
 return codes:\n
 0 - passed\n
 1 - bad arguments\n
 2 - failed to create output directory\n
 3 - CMake failed\n
 4 - make failed\n
'

if [ "$1" = "-h" ]
then
    echo $HELP
    exit 0
fi
if [ "$1" = "--help" ]
then
    echo $HELP
    exit 0
fi

if [ ! -d "$1" ]
then
    echo "source directory does not exist: $1"
    echo $HELP
    exit 1
fi
case $1 in
  /*) ROOT=$1;;
  *) ROOT=$PWD/$1;;
esac

if [ "$#" -lt "2" ]
then
    echo "output directory not specified"
    echo $HELP
    exit 1
fi
case $2 in
  /*) OUTDIR=$2;;
  *) OUTDIR=$PWD/$2;;
esac

shift 2

BUILD=Debug
RC=0
ARGS=""

while [ "$#" -gt "1" ]
do
    key="$1"
    case $key in
        -h|--help)
            echo $HELP
            exit 0
            ;;
        -b|--build)
            if [ "$2" = "Debug" ]
            then
                BUILD=Debug
            elif [ "$2" = "Release" ]
            then
                BUILD=Release
            else
                echo "for --build, specify Debug or Release"
                exit 1
            fi
            ;;
        -a|--args)
            ARGS="$2"
            shift 
            ;;
        *)
            echo "unknown option " $key
            exit 1
        ;;
    esac
    shift # past argument or value
done

mkdir -p $OUTDIR || ( echo "cannot create $OUTDIR" && exit 2 )
cd $OUTDIR

mkdir -p ${BUILD}
cd ${BUILD}

cmake -DCMAKE_BUILD_TYPE=${BUILD} ${ROOT}  || ( echo "CMake failed" && exit 3 )

make  || ( echo "make failed" && exit 4 )