# disable default build rules
.SUFFIXES:

# build dir creation
BDIR := .build
$(shell mkdir -p $(BDIR) $(BDIR)/tests >/dev/null)

# compilation settings
NGS_PATH := /panfs/traces01.be-md.ncbi.nlm.nih.gov/trace_software/ngs/ngs-versions/1.2.5/linux/release/x86_64
export CPATH := $(NGS_PATH)/include:.
export LIBRARY_PATH := $(NGS_PATH)/lib
export LD_LIBRARY_PATH := $(NGS_PATH)/lib
CXXFLAGS += -g -std=c++11 -O3 -fopenmp 
LDFLAGS += -fopenmp
LDLIBS += -lncbi-ngs-c++-static -lngs-c++-static -lncbi-vdb-static -ldl

# compilation rule
SOURCES := $(wildcard *.cpp tests/*.cpp)
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
	aligns_to_build_index_db \
	contig_builder \
	kmers_stat
TESTS := $(wildcard tests/*.cpp)
TESTS := $(TESTS:%.cpp=%) # strip extension

BINARIES := $(PROGRAMS) $(TESTS)
MAIN_OBJECTS = $(BINARIES:%=$(BDIR)/%.o)
SHARED_OBJECTS := $(filter-out $(MAIN_OBJECTS), $(OBJECTS))

$(BINARIES): %: $(BDIR)/%.o $(SHARED_OBJECTS)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

# phony targets
.PHONY: all clean test

all: $(BINARIES)

check: $(TESTS)
	@for t in $(TESTS); do \
		echo "Running $$t:"; \
		$$t || exit 1; \
		echo ;\
	done

.DEFAULT_GOAL := all

clean:
	rm -rf $(BINARIES) $(BDIR)
