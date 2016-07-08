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

#include "referencematchiterator.hpp"

#include <ngs/ncbi/NGS.hpp>

#include "searchbuffer.hpp"

//////////////////// ReferenceSearchBuffer

class ReferenceSearchBuffer : public SearchBuffer
{
public:
    ReferenceSearchBuffer ( SearchBlock* p_sb, const std::string& p_accession, const ngs :: Reference& p_reference )
    :   SearchBuffer ( p_sb, p_accession ),
        m_reference ( p_reference ),
        m_done ( false )
    {
    }

    virtual bool NextMatch ( std::string& p_fragmentId )
    {
        return false;
    }

    virtual std::string BufferId () const
    {
        return m_reference . getCanonicalName ();
    }

private:
    ngs :: Reference m_reference;
    bool m_done;
};

//////////////////// ReferenceMatchIterator

ReferenceMatchIterator :: ReferenceMatchIterator ( SearchBlock :: Factory& p_factory, const std::string& p_accession )
:   MatchIterator ( p_factory, p_accession ),
    m_iter ( ncbi :: NGS :: openReadCollection ( p_accession ) . getReferences () )
{
}

ReferenceMatchIterator :: ~ReferenceMatchIterator ()
{
}

SearchBuffer*
ReferenceMatchIterator :: NextBuffer ()
{   // returns a buffer associated with the next reference in the iterator
    if ( ! m_iter . nextReference () )
    {
        return 0;
    }
    return new ReferenceSearchBuffer ( m_factory.MakeSearchBlock(), m_accession, m_iter );
}

