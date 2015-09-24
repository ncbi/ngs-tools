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

default: std

TOP = $(abspath ../..)
MODULE = sra-search/test

# the location of Makefile.env is likely to change, for now take it from the parent project's directory
include ../Makefile.env

INCDIRS += -I$(SRCDIR)/.. -I$(NGS_INCDIR)

std: 
	@ $(MAKE) --no-print-directory test-whitebox
	@ $(MAKE) --no-print-directory test-cmdline

clean: 
	rm -f $(TEST_SRA_SEARCH_OBJ) $(TEST_BINDIR)/test-sra-search

.PHONY: clean 

#-------------------------------------------------------------------------------
# test-sra-search
#
TEST_SRA_SEARCH_SRC = \
	test-sra-search \

TEST_SRA_SEARCH_OBJ = \
	$(addprefix $(OBJDIR)/,$(addsuffix .$(OBJX),$(TEST_SRA_SEARCH_SRC)))

TEST_SRA_SEARCH_LIBS = \
    -L$(VDB_LIBDIR)/../ilib     \
    -lktst              \
    -L$(NGS_LIBDIR)     \
    -lncbi-ngs-c++      \
    -L$(VDB_LIBDIR)     \
	-lncbi-vdb-static   \
	-lngs-c++           \

$(TEST_BINDIR)/test-sra-search: $(TEST_SRA_SEARCH_OBJ)
	$(CP) -o $@ $^ $(TEST_SRA_SEARCH_LIBS) $(LDFLAGS)

test-whitebox: $(TEST_BINDIR)/test-sra-search    
	@ echo
	@ echo "Unit-testing $(BINDIR)/sra-search"
	@ $^
	@ echo

vg-whitebox: $(TEST_BINDIR)/test-sra-search    
	valgrind --ncbi --suppressions=$(SRCDIR)/valgrind.suppress $^

#TODO: create a generic runcmdline.sh to check exit code and diff output with expected    
test-cmdline: 
	@ $(MAKE) --no-print-directory -C ..
	@ echo
	@ echo "Testing $(BINDIR)/sra-search"
	@ # help
	@ ./runtestcase.sh $(BINDIR)/sra-search 1.0-h    . 0 -h
	@ ./runtestcase.sh $(BINDIR)/sra-search 1.1-help . 0 --help 
	@ # bad args
	@ ./runtestcase.sh $(BINDIR)/sra-search 2.0-no-args                         . 1 
	@ ./runtestcase.sh $(BINDIR)/sra-search 2.1-no-run                          . 1 AGCT
	@ ./runtestcase.sh $(BINDIR)/sra-search 2.2-missing-algorithm               . 1 AGCT SRR000001 -a
	@ ./runtestcase.sh $(BINDIR)/sra-search 2.3-bad-algorithm                   . 1 --algorithm bad AGCT SRR000001
	@ # single run
	@ ./runtestcase.sh $(BINDIR)/sra-search 3.0-single-run-default              . 0 AGCTAGCTAGCT SRR000001
	@ ./runtestcase.sh $(BINDIR)/sra-search 3.1-single-run-fgrep-dumb           . 0 AGCTAGCTAGCT -a FgrepDumb SRR000001
	@ ./runtestcase.sh $(BINDIR)/sra-search 3.2-single-run-fgrep-boyer-moore    . 0 AGCTAGCTAGCT --algorithm FgrepBoyerMoore SRR000001
	@ # 
	@ echo

.PHONY: test-wb test-cmdline
