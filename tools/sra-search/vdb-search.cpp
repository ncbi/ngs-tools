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

#include "vdb-search.hpp"

#include <sstream>

#include <ngs/ncbi/NGS.hpp>
#include <atomic32.h>

#include <klib/time.h>

#include <kproc/thread.h>
#include <kproc/lock.h>

#include <ngs-vdb/inc/NGS-VDB.hpp>
#include <ngs-vdb/inc/FragmentBlobIterator.hpp>

#include "searchblock.hpp"

using namespace std;
using namespace ngs;
using namespace ncbi::ngs::vdb;

bool VdbSearch :: logResults = false;

//////////////////// VdbSearch :: OutputQueue

class VdbSearch :: OutputQueue
{   // thread safe output queue; one consumer, multiple producers
    // counts producers
    // if there are active producers, Pop will wait for new items to appear or the last producer to go away
public:
    OutputQueue ( unsigned int p_producers ) throw ( ErrorMsg )
    {
        atomic32_set ( & m_producers, p_producers );
        rc_t rc = KLockMake ( & m_outputQueueLock );
        if ( rc != 0 )
        {
            throw ( ErrorMsg ( "KLockMake failed" ) );
        }
    }
    ~OutputQueue()
    {
        queue < Item > empty;
        swap( m_queue, empty );

        KLockRelease ( m_outputQueueLock );
    }

    void ProducerDone () // called by the producers
    {
        assert ( atomic32_read ( & m_producers ) > 0 );
        atomic32_dec ( & m_producers );
    }

    void Push (const string& p_accession, const string& p_fragmentId) // called by the producers
    {
        KLockAcquire ( m_outputQueueLock );
        m_queue . push ( Item ( p_accession , p_fragmentId ) );
        KLockUnlock ( m_outputQueueLock );
    }

    // called by the consumer; will block until items become available or the last producer goes away
    bool Pop (string& p_accession, string& p_fragmentId)
    {
        while ( m_queue . size () == 0 && atomic32_read ( & m_producers ) > 0 )
        {
            KSleepMs(1);
        }

        assert ( m_queue . size () > 0 || atomic32_read ( & m_producers ) == 0 );

        if ( m_queue . size () > 0 )
        {
            KLockAcquire ( m_outputQueueLock );
            p_accession = m_queue . front () . first;
            p_fragmentId = m_queue . front () . second;
            m_queue.pop();
            KLockUnlock ( m_outputQueueLock );
            return true;
        }

        assert ( atomic32_read ( & m_producers ) == 0 );
        return false;
    }

private:
    typedef pair < string, string > Item;
    queue < Item > m_queue;

    KLock* m_outputQueueLock;

    atomic32_t m_producers;
};

////////////////////  VdbSearch :: SearchThreadBlock

struct VdbSearch :: SearchThreadBlock
{
    KLock* m_searchQueueLock;
    VdbSearch :: SearchQueue& m_search;
    VdbSearch :: OutputQueue& m_output;

    SearchThreadBlock ( SearchQueue& p_search, OutputQueue& p_output )
    : m_search ( p_search ), m_output ( p_output )
    {
        rc_t rc = KLockMake ( & m_searchQueueLock );
        if ( rc != 0 )
        {
            throw ( ErrorMsg ( "KLockMake failed" ) );
        }
    }
    ~SearchThreadBlock ()
    {
        KLockRelease ( m_searchQueueLock );
    }
};

//////////////////// VdbSearch :: SearchBuffer

VdbSearch :: SearchBuffer :: SearchBuffer ( SearchBlock* p_sb, const std::string& p_accession )
:   m_searchBlock ( p_sb),
    m_accession ( p_accession )
{
}

VdbSearch :: SearchBuffer :: ~SearchBuffer ()
{
    delete m_searchBlock;
}

class VdbSearch :: FragmentSearchBuffer : public VdbSearch :: SearchBuffer
{
public:
    FragmentSearchBuffer ( SearchBlock* p_sb, const std::string& p_accession, const Fragment& p_fragment )
    :   SearchBuffer ( p_sb, p_accession ),
        m_fragment ( p_fragment ),
        m_done ( false )
    {
    }

    virtual bool NextMatch ( std::string& p_fragmentId )
    {
        if ( ! m_done )
        {
            StringRef bases = m_fragment.getFragmentBases();
            if ( m_searchBlock -> FirstMatch ( bases . data (), bases . size () ) )
            {
                p_fragmentId = m_fragment . getFragmentId () . toString ();
                m_done = true;
                return true;
            }
        }
        return false;
    }

    virtual std::string BufferId () const
    {
        return m_fragment . getFragmentId () . toString ();
    }

private:
    const Fragment& m_fragment;
    bool m_done;
};

class VdbSearch :: BlobSearchBuffer : public VdbSearch :: SearchBuffer
{
public:
    BlobSearchBuffer ( SearchBlock* p_sb, const std::string& p_accession, KLock* p_lock, const FragmentBlob& p_blob )
    :   SearchBuffer ( p_sb, p_accession ),
        m_dbLock ( p_lock ),
        m_blob ( p_blob ),
        m_startInBlob ( 0 )
    {
    }

    virtual bool NextMatch ( std::string& p_fragmentId )
    {
        string id = BufferId();
        if ( VdbSearch :: logResults )
        {
            cout << id << ": m_startInBlob=" << m_startInBlob << " size=" << ( m_blob . Size () - m_startInBlob ) << endl;
        }

        uint64_t hitStart;
        uint64_t hitEnd;
        while ( m_searchBlock -> FirstMatch ( m_blob . Data () + m_startInBlob, m_blob . Size () - m_startInBlob, hitStart, hitEnd  ) )
        {
            if ( VdbSearch :: logResults )
            {
                cout << "hitStart=" << hitStart << " hitEnd=" << hitEnd << endl;
            }

            // convert to offsets from the start of the blob
            hitStart += m_startInBlob;
            hitEnd += m_startInBlob;

            uint64_t startInBlob;
            uint64_t lengthInBases;
            uint64_t fragEnd;
            bool biological;

            KLockAcquire ( m_dbLock );
            m_blob . GetFragmentInfo ( hitStart, p_fragmentId, startInBlob, lengthInBases, biological );
            KLockUnlock ( m_dbLock );

            fragEnd = startInBlob + lengthInBases;
            if ( VdbSearch :: logResults )
            {
                cout << "fragId=" << p_fragmentId << " hitStart=" << hitStart << " hitEnd=" << hitEnd << " fragEnd=" << fragEnd << " biological=" << ( biological ? "true" : "false" ) << endl;
            }
            if ( biological )
            {
                if ( hitEnd < fragEnd ||                                                                    // inside a fragment: report and move to the next fragment; or
                    m_searchBlock -> FirstMatch ( m_blob . Data () + hitStart, fragEnd - hitStart  ) )    // result crosses fragment boundary: retry within the fragment
                {
                    if ( VdbSearch :: logResults )
                    {
                        cout << "secondary startInBlob=" << hitStart << endl;
                    }
                    m_startInBlob = fragEnd; // search will resume with the next fragment
                    if ( VdbSearch :: logResults )
                    {
                        cout << "updated startInBlob=" << m_startInBlob << endl;
                    }
                    return true;
                }
                // false hit
            }
            // move on to the next fragment
            m_startInBlob = fragEnd;
        }
        m_startInBlob = 0;
        return false;
    }

    virtual std::string BufferId () const
    {   // identify by row Id range
        int64_t first;
        uint64_t count;
        m_blob . GetRowRange ( first, count );
        ostringstream ret;
        ret << first << "-" << ( first + count - 1 );
        return ret.str();
    }

private:
    KLock*          m_dbLock;
    FragmentBlob    m_blob;
    uint64_t        m_startInBlob;
};

//////////////////// VdbSearch :: MatchIterator

VdbSearch :: MatchIterator :: MatchIterator ( SearchBlockFactory& p_factory, const std::string& p_accession )
:   m_factory ( p_factory ),
    m_accession ( p_accession )
{
}

VdbSearch :: MatchIterator :: ~MatchIterator ()
{
}

// Searches all reads fragment by fragment
class VdbSearch :: FragmentMatchIterator : public VdbSearch :: MatchIterator
{
public:
    FragmentMatchIterator ( SearchBlockFactory& p_factory, const std::string& p_accession )
    :   MatchIterator ( p_factory, p_accession ),
        m_coll ( ncbi :: NGS :: openReadCollection ( p_accession ) ),
        m_readIt ( m_coll . getReads ( Read :: all ) )
    {
        m_readIt . nextRead ();
    }
    virtual ~FragmentMatchIterator ()
    {
    }

    virtual SearchBuffer* NextBuffer ()
    {
        if ( ! m_readIt . nextFragment () )
        {   // end of read, switch to the next
            bool haveFragment = false;
            while ( m_readIt . nextRead () )
            {
                if ( m_readIt . nextFragment () )
                {
                    haveFragment = true;
                    break;
                }
            }
            if ( ! haveFragment )
            {
                return 0;
            }
        }
        // report one match per fragment
        return new FragmentSearchBuffer ( m_factory.MakeSearchBlock(), m_accession, m_readIt );
    }

private:
    ngs::ReadCollection m_coll;
    ngs::ReadIterator   m_readIt;
};

// Searches blob by blob
class VdbSearch :: BlobMatchIterator : public VdbSearch :: MatchIterator
{
public:
    BlobMatchIterator ( SearchBlockFactory& p_factory, const std::string& p_accession )
    :   MatchIterator ( p_factory, p_accession ),
        m_coll ( NGS_VDB :: openVdbReadCollection ( p_accession ) ),
        m_blobIt ( m_coll . getFragmentBlobs() )
    {
        rc_t rc = KLockMake ( & m_accessionLock );
        if ( rc != 0 )
        {
            throw ( ErrorMsg ( "KLockMake failed" ) );
        }
    }
    virtual ~BlobMatchIterator ()
    {
        KLockRelease ( m_accessionLock );
    }

    virtual SearchBuffer* NextBuffer ()
    {
        if ( m_blobIt . hasMore () )
        {
            return new BlobSearchBuffer ( m_factory.MakeSearchBlock(), m_accession, m_accessionLock, m_blobIt . nextBlob () );
        }
        return 0;
    }

private:
    VdbReadCollection       m_coll;
    KLock*                  m_accessionLock;
    FragmentBlobIterator    m_blobIt;
};

//////////////////// VdbSearch :: Settings

static
const
struct {
    const char* name;
    VdbSearch :: Algorithm value;
} Algorithms[] = {
#define ALG(n) { #n, VdbSearch :: n }
    { "FgrepStandard", VdbSearch :: FgrepDumb },
    ALG ( FgrepBoyerMoore ),
    ALG ( FgrepAho ),
    ALG ( AgrepDP ),
    ALG ( AgrepWuManber ),
    ALG ( AgrepMyers ),
    ALG ( AgrepMyersUnltd ),
    ALG ( NucStrstr ),
    ALG ( SmithWaterman ),
#undef ALG
};

VdbSearch :: Settings :: Settings ()
:   m_algorithm ( VdbSearch :: FgrepDumb ),
    m_isExpression ( false ),
    m_minScorePct ( 100 ),
    m_threads ( 2 ),
    m_useBlobSearch ( true ),
    m_referenceDriven ( false )
{
}

bool
VdbSearch :: Settings :: SetAlgorithm ( const std :: string& p_algStr )
{
    for ( size_t i = 0 ; i < sizeof ( Algorithms ) / sizeof ( Algorithms [ 0 ] ); ++i )
    {
        if ( string ( Algorithms [ i ] . name ) == p_algStr )
        {
            m_algorithm = Algorithms [ i ] . value;
            return true;
        }
    }
    return false;
}

//////////////////// VdbSearch

static
void
CheckArguments ( const VdbSearch :: Settings& p_settings ) throw ( invalid_argument )
{
    if ( p_settings . m_isExpression && p_settings . m_algorithm != VdbSearch :: NucStrstr )
    {
        throw invalid_argument ( "query expressions are only supported for NucStrstr" );
    }
    if ( p_settings . m_minScorePct != 100 )
    {
        switch ( p_settings . m_algorithm )
        {
            case VdbSearch :: FgrepDumb:
            case VdbSearch :: FgrepBoyerMoore:
            case VdbSearch :: FgrepAho:
            case VdbSearch :: NucStrstr:
                throw invalid_argument ( "this algorithm only supports 100% match" );
            default:
                break;
        }
    }
    if ( p_settings . m_referenceDriven )
    {
        throw invalid_argument ( "reference mode is not implemented" );
    }
}

VdbSearch :: VdbSearch ( const Settings& p_settings )
    throw ( invalid_argument )
:   m_settings ( p_settings ),
    m_sbFactory ( m_settings ),
    m_buf ( 0 ),
    m_output ( 0 ),
    m_searchBlock ( 0 )
{
    if ( m_settings . m_useBlobSearch && m_settings . m_algorithm == VdbSearch :: SmithWaterman )
    {
        m_settings . m_useBlobSearch = false; // SW takes too long on big buffers
    }

    CheckArguments ( m_settings );

    for ( vector<string>::const_iterator i = m_settings . m_accessions . begin(); i != m_settings . m_accessions . end(); ++i )
    {
        if ( m_settings . m_useBlobSearch )
        {
            m_searches . push ( new BlobMatchIterator ( m_sbFactory, *i ) );
        }
        else
        {
            m_searches . push ( new FragmentMatchIterator ( m_sbFactory, *i ) );
        }
    }
}

VdbSearch :: ~VdbSearch ()
{
    while ( ! m_searches . empty () )
    {
        delete m_searches . front ();
        m_searches . pop ();
    }

    // make sure all threads are gone before we release shared objects
    for ( ThreadPool :: iterator i = m_threadPool . begin (); i != m_threadPool. end (); ++i )
    {
        KThreadCancel ( *i );
        KThreadWait  ( *i, 0 );
        KThreadRelease ( *i );
    }
    delete m_buf;
    delete m_searchBlock;
    delete m_output;
}

VdbSearch :: SupportedAlgorithms
VdbSearch :: GetSupportedAlgorithms ()
{
    vector < string > ret;
    for ( size_t i = 0 ; i < sizeof ( Algorithms ) / sizeof ( Algorithms [ 0 ] ); ++i )
    {
        ret . push_back ( Algorithms [ i ] . name );
    }
    return ret;
}

rc_t CC VdbSearch :: SearchAccPerThread ( const KThread *self, void *data )
{
    assert ( data );
    SearchThreadBlock& sb = * reinterpret_cast < SearchThreadBlock* > ( data );
    assert ( sb . m_searchQueueLock );

    while ( true )
    {
        KLockAcquire ( sb . m_searchQueueLock );
        if ( sb . m_search . empty () )
        {
            KLockUnlock ( sb . m_searchQueueLock );
            break;
        }

        MatchIterator* it = sb . m_search . front ();
        sb . m_search . pop ();

        KLockUnlock ( sb . m_searchQueueLock );

        string fragmentId;
        SearchBuffer* buf = it -> NextBuffer();
        while ( buf != 0 )
        {
            while ( buf -> NextMatch ( fragmentId ) )
            {
                sb . m_output . Push ( buf -> AccessionName (), fragmentId );
            }
            delete buf;
            buf = it -> NextBuffer();
        }
        delete it;
    }

    sb . m_output . ProducerDone();
    return 0;
}

rc_t CC VdbSearch :: SearchBlobPerThread ( const KThread *self, void *data )
{
    assert ( data );
    SearchThreadBlock& sb = * reinterpret_cast < SearchThreadBlock* > ( data );
    assert ( sb . m_searchQueueLock );
    if ( VdbSearch :: logResults && ! sb . m_search . empty () )
    {
        KLockAcquire ( sb . m_searchQueueLock );
        cout << "Thread " << (void*)self << " sb=" << (void*)sb . m_search . front () << endl;
        KLockUnlock ( sb . m_searchQueueLock );
    }

    while ( true )
    {
        KLockAcquire ( sb . m_searchQueueLock );
        if ( VdbSearch :: logResults )
        {
            cout << "Thread " << (void*)self << " locked" << endl;
        }

        SearchBuffer* buf = 0;
        if ( ! sb . m_search . empty () )
        {
            if ( VdbSearch :: logResults )
            {
                cout << "Thread " << (void*)self << " NextBuffer()" << endl;
            }
            buf = sb . m_search . front () -> NextBuffer();
            while ( buf == 0 )
            {   // no more blobs, discard the accession
                delete sb . m_search . front ();
                sb . m_search . pop ();
                if ( sb . m_search . empty () )
                {
                    break;
                }
                if ( VdbSearch :: logResults )
                {
                    cout << "Thread " << (void*)self << " new accession, NextBuffer()" << endl;
                }
                buf = sb . m_search . front () -> NextBuffer();
            }
        }
        KLockUnlock ( sb . m_searchQueueLock );
        if ( buf == 0 )
        {
            if ( VdbSearch :: logResults )
            {
                cout << "Thread " << (void*)self << " done" << endl;
            }
            break;
        }

        string id;
        id = buf->BufferId();
        if ( VdbSearch :: logResults )
        {
            cout << "Thread " << (void*)self << " buf=" << (void*)buf << " bufId=" << id << " unlocked " << id << endl;
        }

        string fragmentId;
        while ( buf -> NextMatch ( fragmentId ) )
        {
            sb . m_output . Push ( buf -> AccessionName (), fragmentId );
        }

        if ( VdbSearch :: logResults )
        {
            cout << "Thread " << (void*)self << " buf=" << (void*)buf << " bufId=" << id << " deleted" << endl;
        }
        delete buf;
    }

    sb . m_output . ProducerDone();
    return 0;
}

bool
VdbSearch :: NextMatch ( string& p_accession, string& p_fragmentId ) throw ( ErrorMsg )
{
    if ( m_settings . m_threads > 0 )
    {
        if ( m_output == 0 ) // first call to NextMatch() - set up
        {
            size_t threadNum = m_settings . m_threads;

            if ( ! m_settings . m_useBlobSearch && threadNum > m_searches . size () )
            {   // in thread-per-accession mode, no need for more threads than there are accessions
                threadNum = m_searches . size ();
            }

            m_output = new OutputQueue ( threadNum );
            m_searchBlock = new SearchThreadBlock ( m_searches, *m_output );
            for ( unsigned  int i = 0 ; i != threadNum; ++i )
            {
                KThread* t;
                rc_t rc = KThreadMake ( &t, m_settings . m_useBlobSearch ? SearchBlobPerThread : SearchAccPerThread, m_searchBlock );
                if ( rc != 0 )
                {
                    throw ( ErrorMsg ( "KThreadMake failed" ) );
                }
                m_threadPool . push_back ( t );
            }
        }

        // block until a result appears on the output queue or all searches are marked as completed
        bool ret = m_output -> Pop ( p_accession, p_fragmentId );
        return ret;
    }
    else
    {   // no threads, return one result at a time
        while ( ! m_searches . empty () )
        {
            if (m_buf == 0)
            {
                m_buf = m_searches.front()->NextBuffer();
            }

            while (m_buf != 0)
            {
                if (m_buf->NextMatch(p_fragmentId))
                {
                    p_accession = m_buf->AccessionName();
                    return true;
                }
                delete m_buf;
                m_buf = m_searches.front()->NextBuffer();
            }

            delete m_searches . front ();
            m_searches . pop ();
        }

        return false;
    }
}

//////////////////// SearchBlockFactory

VdbSearch :: SearchBlockFactory :: SearchBlockFactory ( const VdbSearch :: Settings& p_settings )
:   m_settings ( p_settings )
{
}

VdbSearch :: SearchBlock*
VdbSearch :: SearchBlockFactory :: MakeSearchBlock () const
{
    switch ( m_settings . m_algorithm )
    {
        case VdbSearch :: FgrepDumb:
        case VdbSearch :: FgrepBoyerMoore:
        case VdbSearch :: FgrepAho:
            return new FgrepSearch ( m_settings . m_query, m_settings . m_algorithm );

        case VdbSearch :: AgrepDP:
        case VdbSearch :: AgrepWuManber:
        case VdbSearch :: AgrepMyers:
        case VdbSearch :: AgrepMyersUnltd:
            return new AgrepSearch ( m_settings . m_query, m_settings . m_algorithm, m_settings . m_minScorePct );

        case VdbSearch :: NucStrstr:
            return new NucStrstrSearch ( m_settings . m_query, m_settings . m_isExpression, m_settings . m_useBlobSearch );

        case VdbSearch :: SmithWaterman:
            return new SmithWatermanSearch ( m_settings . m_query, m_settings . m_minScorePct );

        default:
            throw ( ErrorMsg ( "SearchBlockFactory: unsupported algorithm" ) );
    }
}
