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

TOP ?= $(abspath ../..)
MODULE = sra-search

include $(TOP)/build/Makefile.config

INCDIRS += -I$(SRCDIR)/.. -I$(NGS_INCDIR) -I$(VDB_INCDIR) -I$(VDB_INCDIR)/os/$(OSINC) -I$(VDB_INCDIR)/cc/gcc/$(ARCH)

# disable default build rules
.SUFFIXES:

# build dir creation
BDIR := .build
OUTDIR := bin
$(shell mkdir -p $(BDIR) $(BDIR)/src $(BDIR)/src/tests >/dev/null)

# compilation settings
export CPATH := $(NGS_PATH)/include:src
#export LIBRARY_PATH := $(NGS_PATH)/lib
#export LD_LIBRARY_PATH := $(NGS_PATH)/lib
CXXFLAGS += -g -std=c++11 -O3 -fopenmp -I$(INCDIRS)
LDFLAGS += -fopenmp -L$(NGS_LIBDIR) -L$(VDB_LIBDIR)
LDLIBS += -lncbi-ngs-c++-static -lngs-c++-static -lncbi-vdb-static -ldl

# compilation rule
SOURCES := $(wildcard src/*.cpp src/tests/*.cpp)
OBJECTS = $(patsubst %.cpp, $(BDIR)/%.o, $(SOURCES))
DEPFLAGS = -MT $@ -MMD -MP -MF $(BDIR)/$*.d.tmp
$(OBJECTS): $(BDIR)/%.o : %.cpp $(BDIR)/%.d
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@
	@mv -f $(BDIR)/$*.d.tmp $(BDIR)/$*.d

# dependency handling
$(BDIR)/%.d: ;
.PRECIOUS: $(BDIR)/%.d
-include $(patsubst %.cpp, $(BDIR)/%.d, $(SOURCES))

# link rules
PROGRAMS := \
	aligns_to \
	build_index \
	check_index \
	db_fasta_to_bin\
	filter_db \
	sort_dbs
PROGRAMS := $(PROGRAMS:%=$(BDIR)/src/%)
TESTS := $(wildcard src/tests/*.cpp)
TESTS := $(TESTS:%.cpp=$(BDIR)/%)

BINARIES := $(PROGRAMS) $(TESTS)
MAIN_OBJECTS = $(BINARIES:%=%.o)
SHARED_OBJECTS := $(filter-out $(MAIN_OBJECTS), $(OBJECTS))

$(BINARIES): %: %.o $(SHARED_OBJECTS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

INSTALLED_BINARIES := $(notdir $(PROGRAMS))
INSTALLED_BINARIES := $(INSTALLED_BINARIES:%=$(OUTDIR)/%)

$(INSTALLED_BINARIES): $(OUTDIR)/%: $(BDIR)/src/%
	cp $^ $@

# phony targets
.PHONY: all clean test

all: $(INSTALLED_BINARIES)

check: $(TESTS)
	@cd src; \
	for t in $(TESTS); do \
		echo "Running $$t:"; \
		../$$t || exit 1; \
		echo ; \
	done

.DEFAULT_GOAL := all

clean:
	rm -rf $(INSTALLED_BINARIES) $(BDIR)
