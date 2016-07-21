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

#include <iostream>

#include <atomic32.h>

#include <klib/time.h>

#include <kproc/thread.h>
#include <kproc/lock.h>

#include <ngs/ErrorMsg.hpp>

#include "searchbuffer.hpp"
#include "blobmatchiterator.hpp"
#include "fragmentmatchiterator.hpp"
#include "referencematchiterator.hpp"

using namespace std;
using namespace ngs;

bool VdbSearch :: logResults = false;

//////////////////// VdbSearch :: OutputQueue

class VdbSearch :: OutputQueue
{   // thread safe output queue; one consumer, multiple producers
    // counts producers
    // if there are active producers, Pop will wait for new items to appear or the last producer to go away
public:
    OutputQueue ( unsigned int p_producers )
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
CheckArguments ( const VdbSearch :: Settings& p_settings )
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
        if ( m_settings . m_referenceDriven )
        {
            m_searches . push ( new ReferenceMatchIterator ( m_sbFactory, *i, m_settings . m_references ) );
        }
        else if ( m_settings . m_useBlobSearch )
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
VdbSearch :: NextMatch ( string& p_accession, string& p_fragmentId )
{
    if ( m_settings . m_referenceDriven )
    {   // reference-driven mode, always single threaded for now
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

    // sequence-driven mode
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

//////////////////// VdbSearch :: SearchBlockFactory

VdbSearch :: SearchBlockFactory :: SearchBlockFactory ( const VdbSearch :: Settings& p_settings )
:   m_settings ( p_settings )
{
}

SearchBlock*
VdbSearch :: SearchBlockFactory :: MakeSearchBlock () const
{
    switch ( m_settings . m_algorithm )
    {
        case VdbSearch :: FgrepDumb:
            return new FgrepSearch ( m_settings . m_query, FgrepSearch :: FgrepDumb );
        case VdbSearch :: FgrepBoyerMoore:
            return new FgrepSearch ( m_settings . m_query, FgrepSearch :: FgrepBoyerMoore );
        case VdbSearch :: FgrepAho:
            return new FgrepSearch ( m_settings . m_query, FgrepSearch :: FgrepAho );

        case VdbSearch :: AgrepDP:
            return new AgrepSearch ( m_settings . m_query, AgrepSearch :: AgrepDP, m_settings . m_minScorePct );
        case VdbSearch :: AgrepWuManber:
            return new AgrepSearch ( m_settings . m_query, AgrepSearch :: AgrepWuManber, m_settings . m_minScorePct );
        case VdbSearch :: AgrepMyers:
            return new AgrepSearch ( m_settings . m_query, AgrepSearch :: AgrepMyers, m_settings . m_minScorePct );
        case VdbSearch :: AgrepMyersUnltd:
            return new AgrepSearch ( m_settings . m_query, AgrepSearch :: AgrepMyersUnltd, m_settings . m_minScorePct );

        case VdbSearch :: NucStrstr:
            return new NucStrstrSearch ( m_settings . m_query, m_settings . m_isExpression, m_settings . m_useBlobSearch );

        case VdbSearch :: SmithWaterman:
            return new SmithWatermanSearch ( m_settings . m_query, m_settings . m_minScorePct );

        default:
            throw ( ErrorMsg ( "SearchBlockFactory: unsupported algorithm" ) );
    }
}
