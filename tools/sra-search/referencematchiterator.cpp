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

#include <ngs-vdb/inc/VdbReference.hpp>

#include "searchbuffer.hpp"

using namespace std;
using namespace ngs;
using namespace ncbi::ngs::vdb;

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

class ReferenceSearchBase : public SearchBuffer
{
public:
    ReferenceSearchBase( SearchBlock* p_sb,
                         const ReadCollection& p_run,
                         const Reference& p_reference,
                         ReferenceMatchIterator :: ReportedFragments& p_reported)
    :   SearchBuffer ( p_sb, p_run . getName() ),
        m_run ( p_run ),
        m_reference ( p_reference ),
        m_refSearch ( 0 ),
        m_refSearchReverse ( 0 ),
        m_alIt ( 0 ),
        m_fragIt ( 0 ),
        m_reported ( p_reported ),
        m_reverse ( false )
    {
        const unsigned int ReferenceMatchTolerancePct = 90; // search on reference has to be looser than in reads
        const unsigned int ThresholdPct = m_searchBlock -> GetScoreThreshold() * ReferenceMatchTolerancePct / 100;
        m_refSearch         = new AgrepSearch ( m_searchBlock -> GetQuery (),                           AgrepSearch :: AgrepWuManber, ThresholdPct );
        m_refSearchReverse  = new AgrepSearch ( ReverseComplementDNA ( m_searchBlock -> GetQuery () ),  AgrepSearch :: AgrepWuManber, ThresholdPct );
    }

    virtual ~ReferenceSearchBase()
    {
        delete m_refSearch;
        delete m_refSearchReverse;
        delete m_alIt;
        delete m_fragIt;
    }

    virtual string BufferId () const
    {
        return m_reference . getCanonicalName ();
    }

    Fragment GetCurrentFragment () { assert ( m_fragIt != 0 ); return * m_fragIt; }

protected:
    bool NextFragment ()
    {
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
                        string fragId = m_fragIt -> getFragmentId () . toString ();
                        if ( m_reported . find ( fragId ) == m_reported.end () )
                        {
                            // cout << "Found " << m_fragIt -> getFragmentId () . toString () << endl;
                            m_reported . insert ( fragId );
                            return true;
                        }
                    }
                }
                // no (more) matches on this read
                delete m_fragIt;
                m_fragIt = 0;
            }
        }
        return false;
    }

protected:
    ReadCollection  m_run;
    Reference       m_reference;

    // loose search on the reference
    SearchBlock*    m_refSearch;
    SearchBlock*    m_refSearchReverse;

    AlignmentIterator* m_alIt;
    FragmentIterator*  m_fragIt;

    ReferenceMatchIterator :: ReportedFragments& m_reported; // all fragments reported for the parent ReferenceMatchIterator, to eliminate double reports

    bool m_reverse;
};

//////////////////// ReferenceSearchBuffer

class ReferenceSearchBuffer : public ReferenceSearchBase
{
public:
    ReferenceSearchBuffer ( SearchBlock* p_sb,
                            const ReadCollection& p_run,
                            const Reference& p_reference,
                            ReferenceMatchIterator :: ReportedFragments& p_reported )
    :   ReferenceSearchBase ( p_sb, p_run, p_reference, p_reported ),
        m_start ( 0 ),
        m_bases ( m_reference.getReferenceBases ( 0 ) ),
        m_offset ( 0 )
    {
        if ( m_reference . getIsCircular () )
        {   // append the start of the reference to be used in search for wraparound matches
            m_bases += m_bases . substr ( 0, m_searchBlock -> GetQuery () . size() * 2 );
        }
    }

    ReferenceSearchBuffer ( SearchBlock* p_sb,
                            const ReadCollection& p_run,
                            const Reference& p_reference, // a slice
                            uint64_t p_start,
                            uint64_t p_end,
                            ReferenceMatchIterator :: ReportedFragments& p_reported )
    :   ReferenceSearchBase ( p_sb, p_run, p_reference, p_reported ),
        m_start ( p_start ),
        m_bases ( m_reference.getReferenceBases ( p_start, p_end - p_start ) ),
        m_offset ( 0 )
    {
        //TODO: handle wraparound slices of circular references
    }

    virtual SearchBuffer :: Match * NextMatch ()
    {
        while ( true )  // for each match on the reference
        {
            if ( m_alIt == 0 ) // start searching at m_offset
            {
                //cout << "searching at " << m_offset << endl;
                uint64_t hitStart;
                uint64_t hitEnd;
                if ( ! m_reverse )
                {
                    if ( ! m_refSearch -> FirstMatch ( m_bases . data () + m_offset, m_bases . size () - m_offset, hitStart, hitEnd ) )
                    {   // no more matches on this reference; switch to reverse search
                        m_reverse = true;
                        m_offset = 0;
                        continue;
                    }
                }
                else
                {
                    if ( ! m_refSearchReverse -> FirstMatch ( m_bases . data () + m_offset, m_bases . size () - m_offset, hitStart, hitEnd ) )
                    {   // no more reverse matches on this reference; the end.
                        return false;
                    }
                }
                // cout << "Match on " << BufferId () << " at " << ( m_start + m_offset + hitStart ) << "-" << ( m_start + m_offset + hitEnd ) << ( m_refSearch == 0 ? " (reverse)" : "" ) << endl;
                m_alIt = new AlignmentIterator ( m_reference . getAlignmentSlice ( ( int64_t ) ( m_start + m_offset + hitStart ), hitEnd - hitStart ) );
                m_offset += hitEnd; //TODO: this may be too far
            }

            if ( NextFragment () )
            {
                return new Match ( m_accession, GetCurrentFragment () );
            }
        }
    }

private:
    uint64_t        m_start;
    String          m_bases;
    uint64_t        m_offset;
};

//////////////////// ReferenceBlobSearchBuffer

class ReferenceBlobSearchBuffer : public ReferenceSearchBase
{
public:
    ReferenceBlobSearchBuffer ( SearchBlock* p_sb,
                                const ReadCollection& p_run,
                                const Reference& p_reference,
                                ReferenceMatchIterator :: ReportedFragments& p_reported )
    :   ReferenceSearchBase ( p_sb, p_run, p_reference, p_reported ),
        m_startInBlob ( 0 ),
        m_end ( (uint64_t)-1 ),
        m_blobIter ( VdbReference ( m_reference ) . getBlobs() ),
        m_curBlob ( 0 ),
        m_nextBlob ( 0 ),
        m_offsetInReference ( 0 ),
        m_offsetInBlob ( 0 ),
        m_lastBlob ( false )
    {
        SetupFirstBlob();
    }

    ReferenceBlobSearchBuffer ( SearchBlock* p_sb,
                                const ReadCollection& p_run,
                                const Reference& p_reference, // a slice
                                uint64_t p_start,
                                uint64_t p_end,
                                ReferenceMatchIterator :: ReportedFragments& p_reported )
    :   ReferenceSearchBase ( p_sb, p_run, p_reference, p_reported ),
        m_startInBlob ( 0 ),
        m_end ( p_end ),
        m_blobIter ( VdbReference ( m_reference ) . getBlobs ( p_start, p_end ) ),
        m_curBlob ( 0 ),
        m_nextBlob ( 0 ),
        m_offsetInReference ( p_start ),
        m_offsetInBlob ( 0 ),
        m_lastBlob ( false )
    {
        if ( SetupFirstBlob() && p_start != 0 )
        { // recalculate starting point of the search
            uint64_t inReference;
            uint32_t repeatCount;
            uint64_t increment;
            m_curBlob . ResolveOffset ( 0, inReference, repeatCount, increment );
            assert ( p_start > inReference );
            m_startInBlob = p_start - inReference;
            m_offsetInBlob = m_startInBlob;
// cout << "inReference=" << inReference << " m_startInBlob=" << m_startInBlob << endl;
        }
    }

    virtual SearchBuffer :: Match * NextMatch ()
    {
        if ( m_bases . size () == 0 )
        {
            return 0;
        }

        while ( true ) // for each blob
        {
            m_unpackedBlobSize = m_curBlob . UnpackedSize ();
            m_reverse = false;
            m_offsetInBlob = 0;

            while ( true )  // for each match in the blob
            {
                if ( m_alIt == 0 ) // start searching at m_offsetInBlob
                {
                    uint64_t hitStart;
                    uint64_t hitEnd;
                    if ( ! m_reverse )
                    {
                        int64_t first;
                        uint64_t count;
                        m_curBlob . GetRowRange ( first, count );
                        // cout << "searching at " << m_offsetInReference + m_offsetInBlob
                        //     << " blob size=" << m_bases.size() << " unpacked=" << m_unpackedBlobSize
                        //     << " rows=" <<  first << "-" << ( first + count - 1) << endl;

                        if ( ! m_refSearch -> FirstMatch ( m_bases . data () + m_offsetInBlob, m_bases . size () - m_offsetInBlob, hitStart, hitEnd ) )
                        {   // no more matches in this blob; switch to reverse search
                            m_reverse = true;
                            m_offsetInBlob = m_startInBlob;
                            continue;
                        }
                    }
                    else if ( ! m_refSearchReverse -> FirstMatch ( m_bases . data () + m_offsetInBlob, m_bases . size () - m_offsetInBlob, hitStart, hitEnd ) )
                    {   // no more reverse matches on this reference; the end for this blob
                        break;
                    }

                    // cout << "Match on " << BufferId () << " '" << string ( m_bases . data () + m_offsetInBlob + hitStart, hitEnd - hitStart ) << "'" <<
                    //         " at " << ( m_offsetInReference + m_offsetInBlob + hitStart ) << "-" << ( m_offsetInReference + m_offsetInBlob + hitEnd ) << ( m_reverse ? " (reverse)" : "" ) <<
                    //         " in blob " << ( m_offsetInBlob + hitStart ) <<
                                    // endl;
                    uint64_t inReference;
                    uint32_t repeatCount;
                    uint64_t increment;
                    m_curBlob . ResolveOffset ( m_offsetInBlob + hitStart, inReference, repeatCount, increment );
                    // cout << "Resolved to  " << inReference << " repeat=" << repeatCount << " inc=" << increment << endl;
                    m_offsetInBlob += hitEnd; //TODO: this may be too far
                    if ( m_end != (uint64_t)-1 && inReference >= m_end )
                    {
                        continue;
                    }
                    m_alIt = new AlignmentIterator ( m_reference . getAlignmentSlice ( ( int64_t ) inReference, hitEnd - hitStart ) );
                }

                if ( NextFragment () )
                {
                    return new Match ( m_accession, GetCurrentFragment () );
                }
            }

            m_offsetInReference += m_unpackedBlobSize;

            if ( m_lastBlob )
            {
                return 0;
            }

            m_curBlob = m_nextBlob;

            m_bases . reserve ( m_curBlob . Size() + m_searchBlock -> GetQuery () . size() * 2 );
            m_bases = String ( m_curBlob . Data(), m_curBlob . Size() );
            if ( m_blobIter . hasMore () )
            {   // append querySize bases from the beginning of the next blob, to catch matches across the two blobs' boundary
                m_nextBlob = m_blobIter. nextBlob ();
                m_bases += String ( m_nextBlob . Data(), m_searchBlock -> GetQuery () . size() * 2 );
            }
            else
            {
                m_lastBlob = true;
                m_bases += m_circularStart;
            }
        }
    }

private:
    bool SetupFirstBlob()
    {
        if ( ! m_blobIter . hasMore () )
        {
            return false;
        }

        m_curBlob = m_blobIter. nextBlob();

        if ( m_reference . getIsCircular () )
        {   // save the start of the reference to be used in search for wraparound matches
            m_circularStart = m_reference . getReferenceBases ( 0, m_searchBlock -> GetQuery () . size() * 2 );
        }

        m_bases . reserve ( m_curBlob . Size() + m_searchBlock -> GetQuery () . size() * 2 );
        m_bases = String ( m_curBlob . Data(), m_curBlob . Size() );
        if ( m_blobIter . hasMore () )
        {   // append querySize bases from the beginning of the next blob, to catch matches across the two blobs' boundary
            m_nextBlob = m_blobIter. nextBlob ();
            m_bases += String ( m_nextBlob . Data(), m_searchBlock -> GetQuery () . size() * 2 );
        }
        else
        {
            m_lastBlob = true;
            m_bases += m_circularStart;
        }
        return true;
    }

private:
    uint64_t                m_startInBlob;
    uint64_t                m_end;

    ReferenceBlobIterator   m_blobIter;
    ReferenceBlob           m_curBlob;
    ReferenceBlob           m_nextBlob;

    String                  m_bases;
    String                  m_circularStart;
    uint64_t                m_offsetInReference;
    uint64_t                m_offsetInBlob;
    uint64_t                m_unpackedBlobSize;
    bool                    m_lastBlob;
};

//////////////////// ReferenceMatchIterator

ReferenceMatchIterator :: ReferenceMatchIterator ( SearchBlock :: Factory&  p_factory,
                                                   const string&            p_accession,
                                                   const ReferenceSpecs&    p_references,
                                                   bool                     p_blobSearch )
:   MatchIterator ( p_factory, p_accession ),
    m_run ( ncbi :: NGS :: openReadCollection ( p_accession ) ),
    m_blobSearch ( p_blobSearch ),
    m_references ( p_references ),
    m_unalignedReadIt ( p_factory, p_accession, ( ngs :: Read :: ReadCategory )  ( ngs :: Read :: unaligned | ngs :: Read :: partiallyAligned  ) ),
    m_unalignedDone ( false )
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
    while ( m_refIt != m_references . end () )
    {
        // cout << "Searching on " << m_refIt -> m_name << endl;
        try
        {
            SearchBuffer* ret;
            if ( m_blobSearch )
            {
                if (  m_refIt -> m_full )
                {
                    ret = new ReferenceBlobSearchBuffer ( m_factory . MakeSearchBlock(),
                                                          m_run,
                                                          m_run . getReference ( m_refIt -> m_name ),
                                                          m_reported );
                }
                else
                {
                    ret = new ReferenceBlobSearchBuffer ( m_factory . MakeSearchBlock(),
                                                          m_run,
                                                          m_run . getReference ( m_refIt -> m_name ),
                                                          m_refIt -> m_start,
                                                          m_refIt -> m_end,
                                                          m_reported );
                }
            }
            else
            {
                if (  m_refIt -> m_full )
                {
                    ret = new ReferenceSearchBuffer ( m_factory . MakeSearchBlock(),
                                                      m_run,
                                                      m_run . getReference ( m_refIt -> m_name ),
                                                      m_reported );
                }
                else
                {
                    ret = new ReferenceSearchBuffer ( m_factory . MakeSearchBlock(),
                                                      m_run,
                                                      m_run . getReference ( m_refIt -> m_name ),
                                                      m_refIt -> m_start,
                                                      m_refIt -> m_end,
                                                      m_reported );
                }
            }
            ++ m_refIt;
            return ret;
        }
        catch ( ngs :: ErrorMsg ex )
        {
            const string NotFoundMsg = "Reference not found";
            if ( string ( ex . what (), NotFoundMsg . size () ) == NotFoundMsg )
            {
                ++ m_refIt;
                continue;
            }
            throw;
        }
    }

    if ( m_references . empty () && ! m_unalignedDone ) // unaligned fragments are not searched if reference spec is not empty
    {
        SearchBuffer* ret = m_unalignedReadIt . NextBuffer ();
        if ( ret != 0 )
        {
            return ret;
        }
        m_unalignedDone = true;
    }

    return 0;
}

