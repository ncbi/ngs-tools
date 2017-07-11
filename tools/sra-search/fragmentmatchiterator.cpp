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

#include "fragmentmatchiterator.hpp"

#include <ngs/ncbi/NGS.hpp>

#include "searchbuffer.hpp"

using namespace std;
using namespace ngs;

//////////////////// FragmentSearchBuffer

class FragmentSearchBuffer : public SearchBuffer
{
public:
    FragmentSearchBuffer ( SearchBlock* p_sb, const std::string& p_accession, const Fragment& p_fragment )
    :   SearchBuffer ( p_sb, p_accession ),
        m_fragment ( p_fragment ),
        m_done ( false )
    {
    }

    virtual SearchBuffer :: Match * NextMatch ()
    {
        if ( ! m_done )
        {
            StringRef bases = m_fragment . getFragmentBases();
            if ( m_searchBlock -> FirstMatch ( bases . data (), bases . size () ) )
            {
                m_done = true;
                return new Match ( m_accession, m_fragment . getFragmentId () . toString (), string ( bases . data (), bases . size () ) );
            }
        }
        return 0;
    }

    virtual std::string BufferId () const
    {
        return m_fragment . getFragmentId () . toString ();
    }

private:
    const Fragment& m_fragment;
    bool m_done;
};

//////////////////// FragmentMatchIterator

FragmentMatchIterator :: FragmentMatchIterator ( SearchBlock :: Factory& p_factory, const std :: string& p_accession, bool p_unalignedOnly )
:   MatchIterator ( p_factory, p_accession ),
    m_coll ( ncbi :: NGS :: openReadCollection ( p_accession ) ),
    m_readIt ( m_coll . getReads ( p_unalignedOnly ?
                                    ( ngs :: Read :: ReadCategory ) ( ngs :: Read :: unaligned | ngs :: Read :: partiallyAligned  ) :
                                    ngs :: Read :: all )
            ),
    m_unalignedOnly ( p_unalignedOnly )
{
    m_readIt . nextRead ();
}

FragmentMatchIterator :: ~FragmentMatchIterator ()
{
}

SearchBuffer*
FragmentMatchIterator :: NextBuffer ()
{
    do
    {
        while ( m_readIt . nextFragment () )
        {
            if ( ! m_unalignedOnly || ! m_readIt . isAligned () )
            {
                // report one match per fragment
                return new FragmentSearchBuffer ( m_factory.MakeSearchBlock(), m_accession, m_readIt );
            }
        }
    }
    while ( m_readIt . nextRead () );
    return 0;
}

