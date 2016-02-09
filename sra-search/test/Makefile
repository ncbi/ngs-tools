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

runtests std: 
	@ $(MAKE) --no-print-directory test-whitebox
	@ $(MAKE) --no-print-directory test-cmdline
	
slowtests: runtests
	@ $(MAKE) --no-print-directory slowtest-cmdline

clean: 
	rm -f $(TEST_SRA_SEARCH_OBJ) $(TEST_BINDIR)/test-sra-search

.PHONY: clean 

#-------------------------------------------------------------------------------
# test-sra-search
#
TEST_SRA_SEARCH_SRC = \
	test-sra-search \

TEST_SRA_SEARCH_OBJ = \
	$(addprefix $(OBJDIR)/,$(addsuffix .$(OBJX),$(TEST_SRA_SEARCH_SRC) ../vdb-search ../searchblock ../ngs-vdb ../NGS_VDB)) 

TEST_SRA_SEARCH_LIBS = \
    -L$(VDB_LIBDIR)/../ilib     \
    -lktst              \
    -L$(NGS_LIBDIR)     \
    -lncbi-ngs-c++      \
    -L$(VDB_LIBDIR)     \
	-lncbi-vdb-static   \
	-lngs-c++           \

$(TEST_BINDIR)/test-sra-search: $(TEST_SRA_SEARCH_OBJ)
	$(MAKE) -C $(SRCDIR)/..
	$(CP) -o $@ $^ $(TEST_SRA_SEARCH_LIBS) $(LDFLAGS)

test-whitebox: $(TEST_BINDIR)/test-sra-search    
	@ echo
	@ echo "Unit-testing $(BINDIR)/sra-search"
	@ $^
	@ echo

vg-whitebox: $(TEST_BINDIR)/test-sra-search    
	$(VALGRIND) --suppressions=$(SRCDIR)/valgrind.suppress --leak-check=full $^

TIME = `which time` -f"%E"

test-cmdline: 
	@ $(MAKE) --no-print-directory -C ..
	@ echo
	@ echo "Testing $(BINDIR)/sra-search"
	@ # help
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 1.0-h    --args "-h"
	@ $(TIME) -f"%E" ./runtestcase.sh $(BINDIR)/sra-search 1.1-help --args "--help" 
	@ # bad args
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 2.0-no-args          	--rc 1 
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 2.1-no-run              --rc 1 --args AGCT
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 2.2-missing-algorithm   --rc 1 --args "AGCT SRR000001 -a"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 2.3-bad-algorithm       --rc 1 --args "--algorithm bad AGCT SRR000001"
	@ # single run
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.05-single-run-default-WGS       --args "AGCTAGCTA ALAI01"
	@ # query expressions
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 5.1-run-not-nucstrstr-with-expression 	--args "AAAAAAACCCCCCCAAAAAAACCCCCCC\|\|AGCTAGCTAGCT --algorithm FgrepStandard --expression  SRR000001" --rc 1 
	@ # imperfect match
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.0.1-imperfect-match-no-arg           --args "AGCTAGCTAGCT SRR000001 -S" --rc 1 
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.0.2-imperfect-match-bad-arg          --args "AGCTAGCTAGCT SRR000001 -S qq" --rc 1 
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.0.3-imperfect-match-default          --args "AGCTAGCTAGCT -S 90 SRR000001" --rc 1 
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.1-imperfect-match-fgrep-standard     --args "AGCTAGCTAGCT --score 90 -a FgrepStandard SRR000001" --rc 1 
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.2-imperfect-match-fgrep-boyer-moore  --args "AGCTAGCTAGCT --score 90 --algorithm FgrepBoyerMoore    SRR000001" --rc 1 
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.3-imperfect-match-fgrep-aho          --args "AGCTAGCTAGCT --score 90 --algorithm FgrepAho           SRR000001" --rc 1 
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.8-imperfect-match-nucstrstr          --args "AGCTAGCTAGCT --score 90 --algorithm NucStrstr          SRR000001" --rc 1 
	@ # multi-threading
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 7.0.1-threads-no-arg   --rc 1 --args "AGCTAGCTAGCT SRR600094 SRR600095 -T"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 7.0.2-threads-bad-arg  --rc 1 --args "AGCTAGCTAGCT SRR600094 SRR600095 --threads -1"
	@ echo

slowtest-cmdline: 
	@ $(MAKE) --no-print-directory -C ..
	@ echo
	@ echo "Slow-Testing $(BINDIR)/sra-search"
	@ # single run
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.0-single-run-default            --args "AGCTAGCTAGCT SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.1-single-run-fgrep-standard     --args "AGCTAGCTAGCT -a FgrepStandard SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.2-single-run-fgrep-boyer-moore  --args "AGCTAGCTAGCT --algorithm FgrepBoyerMoore    SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.3-single-run-fgrep-aho          --args "AGCTAGCTAGCT --algorithm FgrepAho           SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.4-single-run-agrep-dp           --args "AGCTAGCTAGCT --algorithm AgrepDP            SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.5-single-run-agrep-wumanber     --args "AGCTAGCTAGCT --algorithm AgrepWuManber      SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.6-single-run-agrep-myers        --args "AGCTAGCTAGCT --algorithm AgrepMyers         SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.7-single-run-agrep-myers-unltd  --args "AGCTAGCTAGCT --algorithm AgrepMyersUnltd    SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.8-single-run-nucstrstr          --args "AGCTAGCTAGCT --algorithm NucStrstr          SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.9-single-run-smith-waterman     --args "AGCTAGCTAGCT --algorithm SmithWaterman      SRR000001"
	@ # multiple runs
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 4.0-multiple-runs --args "TACTGGGCGTAAAGCGTGCGCAGGCGGT --algorithm NucStrstr SRR600096 SRR600097"
	@ # query expressions
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 5.0-run-nucstrstr-with-expression 	    --args "AAAAAAACCCCCCCAAAAAAACCCCCCC\|\|AGCTAGCTAGCT --algorithm NucStrstr --expression  SRR000001"
	@ # imperfect match
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.4-imperfect-match-agrep-dp           --args "AGCTAGCTAGCT --score 90 --algorithm AgrepDP            SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.5-imperfect-match-agrep-wumanber     --args "AGCTAGCTAGCT --score 90 --algorithm AgrepWuManber      SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.6-imperfect-match-agrep-myers        --args "AGCTAGCTAGCT --score 90 --algorithm AgrepMyers         SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.7-imperfect-match-agrep-myers-unltd  --args "AGCTAGCTAGCT --score 90 --algorithm AgrepMyersUnltd    SRR000001"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 6.9-imperfect-match-smith-waterman     --args "AGCTAGCTAGCT --score 90 --algorithm SmithWaterman      SRR000001"
	@ # multi-threading
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 7.1-threads-unsorted --sort --args "ACGTAGGGTCC --threads 4 SRR600094 SRR600095 SRR600096 SRR600099"
	@ # nothing found 
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 8.0-nothing-found --rc 1 --args "TTTTTTTTAAAAAAACCCC SRR000001"
	@ # 
	@ echo


onetest: 
	@ $(MAKE) --no-print-directory -C ..
	#@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 3.05-single-run-default-WGS       --args "AGCTAGCTA ALAI01 --blobs"
	@ $(TIME) ./runtestcase.sh $(BINDIR)/sra-search 8.0-nothing-found --rc 1 --args "TTTTTTTTAAAAAAACCCC SRR000001"
    
.PHONY: test-wb test-cmdline onetest
