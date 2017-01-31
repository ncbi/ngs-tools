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

TOP ?= $(abspath ../..)

MODULE = tools/skesa

include $(TOP)/build/Makefile.env

ALL_TOOLS = \
	skesa

# path to pre-built boost
#builddir = $(BUILDDIR)
#include $(builddir)/Makefile.mk


CC = g++ -std=gnu++11  
CFLAGS = -Wall -Wno-format-y2k  -pthread -fPIC -O3 -finline-functions -fstrict-aliasing -fomit-frame-pointer -msse4.2 \
         $(INCDIRS) $(BOOST_INCLUDE)
LIBS = -L $(NGSLIBS) -lncbi-ngs-c++ -lngs-c++ -lncbi-vdb-static \
       $(BOOST_LIBPATH) $(BOOST_PROGRAM_OPTIONS_STATIC_LIBS) $(BOOST_IOSTREAMS_STATIC_LIBS) $(BOOST_REGEX_STATIC_LIBS) -lboost_timer$(BOOST_TAG)-static -lboost_chrono$(BOOST_TAG)-static  $(BOOST_SYSTEM_STATIC_LIBS) \
       -lrt -ldl -lm  -lpthread -lz
# -ltcmalloc

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

all: skesa wgmlst dbgtester guidedassembler

skesa.o: readsgetter.hpp counter.hpp graphdigger.hpp assembler.hpp KmerInit.hpp DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp
skesa: skesa.o 
	$(CC) -o $@ $^ $(LIBS)

wgmlst.o: KmerInit.hpp DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp
glb_align.o: glb_align.hpp
wgmlst: wgmlst.o glb_align.o
	$(CC) -o $@ $^ $(LIBS)

dbgtester.o: graphdigger.hpp KmerInit.hpp DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp
dbgtester: dbgtester.o
	$(CC) -o $@ $< $(LIBS)

guidedassembler.o: guidedpath.hpp readsgetter.hpp counter.hpp graphdigger.hpp KmerInit.hpp DBGraph.hpp Integer.hpp LargeInt.hpp LargeInt1.hpp LargeInt2.hpp Model.hpp config.hpp
guidedassembler: guidedassembler.o glb_align.o
	$(CC) -o $@ $^ $(LIBS)
