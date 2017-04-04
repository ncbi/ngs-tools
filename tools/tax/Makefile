# disable default build rules
.SUFFIXES:

# build dir creation
BDIR := .build
OUTDIR := bin
$(shell mkdir -p $(OUTDIR) $(BDIR) $(BDIR)/src $(BDIR)/src/tests >/dev/null)

# compilation settings
NGS_PATH := /panfs/traces01.be-md.ncbi.nlm.nih.gov/trace_software/ngs/ngs-versions/1.2.5/linux/release/x86_64
export CPATH := $(NGS_PATH)/include:src
export LIBRARY_PATH := $(NGS_PATH)/lib
export LD_LIBRARY_PATH := $(NGS_PATH)/lib
CXXFLAGS += -g -std=c++11 -O3 -fopenmp 
LDFLAGS += -fopenmp
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
