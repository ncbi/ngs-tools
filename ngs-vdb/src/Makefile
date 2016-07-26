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
MODULE = ngs-vdb

include $(TOP)/build/Makefile.config

INCDIRS += -I$(SRCDIR)/.. -I$(NGS_INCDIR) -I$(VDB_INCDIR) -I$(VDB_INCDIR)/os/$(OSINC) -I$(VDB_INCDIR)/cc/gcc/$(ARCH)

std: $(OBJDIR) $(ILIBDIR)
	$(MAKE) $(ILIBDIR)/ngs-vdb

clean:
	rm -f $(ILIBDIR)/$(LPFX)ngs-vdb.$(LIBX)
	#$(MAKE) -C test clean

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
NGS_VDB_SRC = \
    FragmentBlob \
    FragmentBlobIterator \
    NGS-VDB \
	VdbReadCollection \

NGS_VDB_OBJ = \
	$(addprefix $(OBJDIR)/,$(addsuffix .$(OBJX),$(NGS_VDB_SRC)))

$(ILIBDIR)/ngs-vdb: $(NGS_VDB_OBJ)
	$(AR) $(ILIBDIR)/$(LPFX)ngs-vdb.$(LIBX) $(NGS_VDB_OBJ)

.PHONY: $(ILIBDIR)/ngs-vdb

#-------------------------------------------------------------------------------
# testing
#
test slow: std
	$(MAKE) TOP=$(TOP) MODPATH=$(MODPATH)/test -C test $@

