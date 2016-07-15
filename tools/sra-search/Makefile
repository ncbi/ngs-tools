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
MODULE = sra-search

include $(TOP)/build/Makefile.config

INCDIRS += -I$(SRCDIR)/.. -I$(NGS_INCDIR) -I$(VDB_INCDIR) -I$(VDB_INCDIR)/os/$(OSINC) -I$(VDB_INCDIR)/cc/gcc/$(ARCH)

std: $(OBJDIR) $(BINDIR)
	$(MAKE) $(BINDIR)/sra-search

clean:
	rm -f $(SRA_SEARCH_OBJ) $(BINDIR)/sra-search*
	$(MAKE) -C test clean

runtests: test

slowtests: slow

.PHONY: clean

#-------------------------------------------------------------------------------
# vers-includes
#
$(TARGDIR)/vers-includes: $(addsuffix .vers.h,sra-search)

.PHONY: $(TARGDIR)/vers-includes

#-------------------------------------------------------------------------------
# sra-search
#
SRA_SEARCH_SRC = \
    searchblock \
	vdb-search \
	main \

SRA_SEARCH_OBJ = \
	$(addprefix $(OBJDIR)/,$(addsuffix .$(OBJX),$(SRA_SEARCH_SRC)))

SRA_SEARCH_LIBS = \
    -L$(ILIBDIR)        \
    -lngs-vdb           \
    -L$(NGS_LIBDIR)     \
    -L$(VDB_LIBDIR)     \
    -lncbi-ngs-c++      \
	-lncbi-vdb-static   \
	-lngs-c++           \
    -lpthread           \
    -ldl

$(BINDIR)/sra-search: $(OBJDIR) $(BINDIR) $(SRA_SEARCH_OBJ)
	$(LP) -o $@$(VERSION_EXEX) $(SRA_SEARCH_OBJ) $(SRA_SEARCH_LIBS)
	ln -sf $(@F)$(VERSION_EXEX) $@$(MAJVERS_EXEX)
	ln -sf $(@F)$(MAJVERS_EXEX) $@

#-------------------------------------------------------------------------------
# testing
#
test slow: std
	$(MAKE) TOP=$(TOP) MODPATH=$(MODPATH)/test -C test $@

