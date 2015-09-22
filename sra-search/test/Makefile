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
	-static             \
    -L$(VDB_LIBDIR)/../ilib     \
    -lktst              \
    -L$(NGS_LIBDIR)     \
    -lncbi-ngs-c++      \
    -L$(VDB_LIBDIR)     \
	-lncbi-vdb-static   \
	-lngs-c++           \

$(TEST_BINDIR)/test-sra-search: $(TEST_SRA_SEARCH_OBJ)
	$(CP) -o $@ $^ $(TEST_SRA_SEARCH_LIBS)

test-whitebox: $(TEST_BINDIR)/test-sra-search    
	$^

#TODO: create a generic runcmdline.sh to check exit code and diff output with expected    
test-cmdline: $(BINDIR)/sra-search
	$^ -h 

.PHONY: test-wb test-cmdline
