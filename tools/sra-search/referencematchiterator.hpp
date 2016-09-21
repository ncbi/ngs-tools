/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

#ifndef _hpp_reference_match_iterator_
#define _hpp_reference_match_iterator_

#include <vector>
#include <set>

#include <ngs/ReadCollection.hpp>
#include "fragmentmatchiterator.hpp"

class ReferenceMatchIterator : public MatchIterator
{
public:
    typedef std :: vector < std :: string > ReferenceNames;
    typedef std :: set < std :: string > ReportedFragments;

public:
    ReferenceMatchIterator ( SearchBlock :: Factory&    p_factory,
                             const std :: string&       p_accession,
                             const ReferenceNames&      p_references = ReferenceNames(),
                             bool                       p_blobSearch = false );

    virtual ~ReferenceMatchIterator ();

    virtual SearchBuffer* NextBuffer ();

private:
    ngs :: ReadCollection   m_run;

    ReferenceNames                      m_references;
    ReferenceNames :: const_iterator    m_refIt;

    FragmentMatchIterator   m_unalignedReadIt;
    bool                    m_readsDone;

    ReportedFragments m_reported; // used to eliminate double reports
};

#endif
