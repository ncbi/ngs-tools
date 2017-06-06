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

default all std: cmake

#-------------------------------------------------------------------------------
# environment
#
TOP ?= $(CURDIR)
include $(TOP)/build/Makefile.env

test runtests: ctest

.PHONY: default std all test runtests

#-------------------------------------------------------------------------------
# std
#
clean: stdclean
	@ -rm -rf $(ILIBDIR) $(BINDIR)

.PHONY: clean

#-------------------------------------------------------------------------------
# install
#
install: cinstall

uninstall:

.PHONY: install uninstall

#-------------------------------------------------------------------------------
# configuration help
#
help configure:
	@ echo "Before initial build, run './configure --build-prefix=<out>' from"
	@ echo "the project root to set the output directory of your builds."
	@ echo "Run ./configure -h for full description."
	@ echo
	@ echo "Targets:"
	@ echo "all, std        : full build"
	@ echo "clean           : remove build results"
	@ echo "test, runtests  : build and run tests"
	@ echo "                  to control which tests are executed, use 'make test CTEST=<any part of a test's name>"
	@ echo "                  e.g. 'make test CTEST=Slow' will run all tests with 'Slow' in their name"
	@ echo "install         : build, install to $(INST_BINDIR)"
	@ echo

.PHONY: help configure
