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

#include <iostream>

#include <ngs/ncbi/NGS.hpp>
#include <ngs/Reference.hpp>

#include "searchbuffer.hpp"

using namespace std;
using namespace ngs;

//////////////////// ReferenceSearchBuffer

static
string
ReverseComplementDNA ( const string& p_source)
{
    string ret;
    ret.reserve ( p_source . size () );
    for ( string :: const_reverse_iterator i = p_source . rbegin (); i != p_source . rend (); ++i )
    {
        char ch;
        switch ( *i )
        {
            case 'A' : ch = 'T'; break;
            case 'C' : ch = 'G'; break;
            case 'G' : ch = 'C'; break;
            case 'T' : ch = 'A'; break;
            case 'N' : ch = 'N'; break;
            default:
                assert(false);
        }
        ret += ch;
    }
    return ret;
}

class ReferenceSearchBuffer : public SearchBuffer
{
public:
    ReferenceSearchBuffer ( SearchBlock* p_sb, const string& p_accession, const ReadCollection& p_run, const Reference& p_reference )
    :   SearchBuffer ( p_sb, p_accession ),
        m_run ( p_run ),
        m_reference ( p_reference ),
        m_bases ( m_reference.getReferenceBases ( 0 ) ),
        m_offset ( 0 ),
        m_refSearch ( NewReferenceSearchBlock ( m_searchBlock -> GetQuery () ) ),
        m_refSearchReverse ( 0 ),
        m_alIt ( 0 ),
        m_fragIt ( 0 )
    {
    }

    virtual ~ReferenceSearchBuffer ()
    {
        delete m_refSearch;
        delete m_refSearchReverse;
        delete m_alIt;
        delete m_fragIt;
    }

    virtual bool NextMatch ( string& p_fragmentId )
    {
        while ( true )  // foreach match on the reference
        {
            if ( m_alIt == 0 ) // start searching at m_offset
            {
                uint64_t hitStart;
                uint64_t hitEnd;
                if ( m_refSearch != 0 )
                {
                    if ( ! m_refSearch -> FirstMatch ( m_bases . data () + m_offset, m_bases . size () - m_offset, hitStart, hitEnd ) )
                    {   // no more matches on this reference; switch to reverse search
                        delete m_refSearch;
                        m_refSearch = 0;
                        string reverseQuery = ReverseComplementDNA ( m_searchBlock -> GetQuery () ); //TODO: this calculation takes place for each reference with the same result; push up the call stack
                        //cout << "Swithing to reverse, query=" << reverseQuery << endl;
                        m_refSearchReverse = NewReferenceSearchBlock ( reverseQuery );
                        m_offset = 0;
                        continue;
                    }
                }
                else
                {
                    if ( ! m_refSearchReverse -> FirstMatch ( m_bases . data () + m_offset, m_bases . size () - m_offset, hitStart, hitEnd ) )
                    {   // no more reverse matches on this reference; the end.
                        delete m_refSearchReverse;
                        m_refSearchReverse = 0;
                        break;
                    }
                }
                // cout << "Match on " << BufferId () << " at " << ( m_offset + hitStart ) << "-" << ( m_offset + hitEnd ) << ( m_refSearch == 0 ? " (reverse)" : "" ) << endl;
                m_alIt = new AlignmentIterator ( m_reference . getAlignmentSlice ( ( int64_t ) ( m_offset + hitStart ), hitEnd - hitStart ) );
                m_offset += hitEnd; //TODO: this may be too far
            }

            while ( true ) // for each alignment on the slice
            {
                if ( m_fragIt == 0 )
                {   // start searching alignment on the matched slice
                    if ( ! m_alIt -> nextAlignment() )
                    {
                        delete m_alIt;
                        m_alIt = 0;
                        break;
                    }
                    m_fragIt = new FragmentIterator ( m_run . getRead ( m_alIt -> getReadId () . toString() ) );   //TODO: there may be a shortcut to get to the fragment's bases
                }

                if ( m_fragIt != 0 )
                {
                    while ( m_fragIt -> nextFragment () ) // foreach fragment
                    {
                        StringRef fragBases = m_fragIt -> getFragmentBases ();
                        // cout << "Searching " << m_fragIt -> getFragmentId () . toString () << "'" << fragBases << "'" << endl;
                        if ( m_searchBlock -> FirstMatch ( fragBases . data (), fragBases . size () ) ) // this search is with the original score threshold
                        {
                            p_fragmentId = m_fragIt -> getFragmentId () . toString ();
                            return true;
                        }
                    }
                    // no (more) matches on this read
                    delete m_fragIt;
                    m_fragIt = 0;
                }
            }
        }
        return false;
    }

    virtual string BufferId () const
    {
        return m_reference . getCanonicalName ();
    }

private:
    SearchBlock* NewReferenceSearchBlock ( const string& p_query ) const
    {
        const unsigned int ReferenceMatchTolerancePct = 90; // search on reference has to be looser than in reads
        const unsigned int ThresholdPct = m_searchBlock -> GetScoreThreshold() * ReferenceMatchTolerancePct / 100;
        return new AgrepSearch ( p_query,  AgrepSearch :: AgrepWuManber, ThresholdPct );
    }

private:
    ReadCollection  m_run;
    Reference       m_reference;
    String          m_bases;
    uint64_t        m_offset;

    // loose search on the reference
    SearchBlock*    m_refSearch;
    SearchBlock*    m_refSearchReverse;

    AlignmentIterator* m_alIt;
    FragmentIterator*  m_fragIt;
};

//////////////////// ReferenceMatchIterator

ReferenceMatchIterator :: ReferenceMatchIterator ( SearchBlock :: Factory& p_factory, const string& p_accession, const std :: vector < std :: string >& p_references )
:   MatchIterator ( p_factory, p_accession ),
    m_run ( ncbi :: NGS :: openReadCollection ( p_accession ) ),
    m_references ( p_references ),
    m_unalignedReadIt ( p_factory, p_accession, ( ngs :: Read :: ReadCategory )  ( ngs :: Read :: unaligned | ngs :: Read :: partiallyAligned  ) ),
    m_readsDone ( false )
{
    if ( m_references . empty () )
    {   // search all references if none specified
        ReferenceIterator refIter = m_run . getReferences ();
        while ( refIter . nextReference () )
        {
            m_references . push_back ( refIter . getCanonicalName () );
        }
    }
    m_refIt = m_references . begin ();
}

ReferenceMatchIterator :: ~ReferenceMatchIterator ()
{
}

SearchBuffer*
ReferenceMatchIterator :: NextBuffer ()
{   // returns a buffer associated with the next reference in the iterator, or with the next unaligned read
    if ( m_refIt != m_references . end () )
    {
        //cout << "Searching on " << *m_refIt << endl;
        ReferenceSearchBuffer* ret = new ReferenceSearchBuffer ( m_factory . MakeSearchBlock(), m_accession, m_run, m_run . getReference ( *m_refIt ) );
        ++ m_refIt;
        return ret;
    }

    if ( ! m_readsDone )
    {
        SearchBuffer* ret = m_unalignedReadIt . NextBuffer ();
        if ( ret != 0 )
        {
            return ret;
        }
        m_readsDone = true;
    }

    return 0;
}

