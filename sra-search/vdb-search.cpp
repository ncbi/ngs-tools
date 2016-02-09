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

#include <ngs/ncbi/NGS.hpp>
#include <atomic32.h>

#include <klib/time.h>

#include <kproc/thread.h>
#include <kproc/lock.h>

#include "searchblock.hpp"
#include "ngs-vdb.hpp"

using namespace std;
using namespace ngs;

//////////////////// VdbSearch :: OutputQueue

class VdbSearch :: OutputQueue 
{   // thread safe output queue; one consumer, multiple producers
    // counts producers
    // if there are active producers, Pop will wait for new items to appear or the last producer to go away
public:
    OutputQueue ( unsigned int p_producers ) throw ( ErrorMsg )
    {
        atomic32_set ( & m_producers, p_producers );
        rc_t rc = KLockMake ( & m_lock );
        if ( rc != 0 )
        {
            throw ( ErrorMsg ( "KLockMake failed" ) );
        }  
    }
    ~OutputQueue()
    {
        queue < Item > empty;   
        swap( m_queue, empty );
        
        KLockRelease ( m_lock );
    }
    
    void ProducerDone () // called by the producers
    {
        assert ( atomic32_read ( & m_producers ) > 0 );
        atomic32_dec ( & m_producers );
    }
    
    void Push (const string& p_accession, const string& p_fragmentId) // called by the producers
    {
        KLockAcquire ( m_lock );
        m_queue . push ( Item ( p_accession , p_fragmentId ) );
        KLockUnlock ( m_lock );
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
            KLockAcquire ( m_lock );
            p_accession = m_queue . front () . first;     
            p_fragmentId = m_queue . front () . second;
            m_queue.pop();
            KLockUnlock ( m_lock );
            return true;
        }
        
        assert ( atomic32_read ( & m_producers ) == 0 );
        return false;
    }
    
private:    
    typedef pair < string, string > Item;
    queue < Item > m_queue;
    
    KLock* m_lock;
    
    atomic32_t m_producers; 
};

////////////////////  VdbSearch :: SearchThreadBlock

struct VdbSearch :: SearchThreadBlock
{
    KLock* m_lock;
    VdbSearch :: SearchQueue& m_search;        
    VdbSearch :: OutputQueue& m_output;
    
    SearchThreadBlock ( SearchQueue& p_search, OutputQueue& p_output )
    : m_search ( p_search ), m_output ( p_output )
    {
        rc_t rc = KLockMake ( & m_lock ); 
        if ( rc != 0 )
        {
            throw ( ErrorMsg ( "KLockMake failed" ) );
        }
    }
    ~SearchThreadBlock ()
    {
        KLockRelease ( m_lock );
    }
};

//////////////////// VdbSearch :: MatchIterator

VdbSearch :: MatchIterator :: MatchIterator ( SearchBlock* p_block, const std::string& p_accession )
:   m_searchBlock ( p_block ),
    m_accession ( p_accession )
{
}

VdbSearch :: MatchIterator :: ~MatchIterator ()
{
    delete m_searchBlock;
}

// Searches fragment by fragment
class VdbSearch :: FragmentMatchIterator : public VdbSearch :: MatchIterator 
{
public:
    FragmentMatchIterator ( SearchBlock* p_block, const std::string& p_accession )
    :   MatchIterator ( p_block, p_accession ),
        m_coll ( ncbi :: NGS :: openReadCollection ( p_accession ) ), 
        m_readIt ( m_coll . getReads ( Read :: all ) )
    {
        m_readIt . nextRead ();
    }
    virtual ~FragmentMatchIterator ()
    {
    }
    
    virtual bool NextMatch ( std::string& p_fragmentId )
    {
        do
        {
            while ( m_readIt . nextFragment () )
            {   // report one match per fragment
                StringRef bases = m_readIt . getFragmentBases ();
                if ( m_searchBlock -> FirstMatch ( bases . data (), bases . size () ) )
                {
                    p_fragmentId = m_readIt . getFragmentId () . toString ();
                    return true;
                }
            }
        }
        while ( m_readIt . nextRead () );
        return false;
    }

private: 
    ngs::ReadCollection m_coll;
    ngs::ReadIterator   m_readIt;
};

// Searches blob by blob
class VdbSearch :: BlobMatchIterator : public VdbSearch :: MatchIterator 
{
public:
    BlobMatchIterator ( SearchBlock* p_block, const std::string& p_accession )
    :   MatchIterator ( p_block, p_accession ),
        m_coll ( ncbi :: NGS_VDB :: openVdbReadCollection ( p_accession ) ), 
        m_blobIt ( m_coll . getFragmentBlobs() ),
        m_startInBlob ( 0 )    
    {
        m_blobIt . nextBlob ();        
    }
    virtual ~BlobMatchIterator ()
    {
    }
    
    virtual bool NextMatch ( std::string& p_fragmentId )
    {
        do
        {
            if ( VdbSearch :: logResults )
            {
                cout << "m_startInBlob=" << m_startInBlob << endl;
            }
            
            uint64_t hitStart;
            uint64_t hitEnd;
            while ( m_searchBlock -> FirstMatch ( m_blobIt . Data () + m_startInBlob, m_blobIt . Size () - m_startInBlob, hitStart, hitEnd  ) )
            {
                if ( VdbSearch :: logResults )
                {
                    cout << "hitStart=" << hitStart << " hitEnd=" << hitEnd << endl;
                }
                
                // convert to offsets from the strat of the blob
                hitStart += m_startInBlob;
                hitEnd += hitEnd;
                
                uint64_t fragEnd;
                bool biological;
                m_blobIt . GetFragmentInfo ( hitStart, p_fragmentId, fragEnd, biological );
                if ( VdbSearch :: logResults )
                {
                    cout << "fragId=" << p_fragmentId << " fragEnd=" << fragEnd << " biological=" << ( biological ? "true" : "false" ) << endl;
                }
                // TODO: if not biological, resume search from fragEnd
                if ( biological )
                {
                    if ( hitEnd < fragEnd || // inside a fragment, report and move to the next fragment
                        m_searchBlock -> FirstMatch ( m_blobIt . Data () + hitStart, fragEnd - hitStart  ) ) // result crosses fragment boundary, retry within the fragment
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
                }
                // false hit, move on to the next fragment
                m_startInBlob = fragEnd;
            }
            m_startInBlob = 0;
        }
        while ( m_blobIt . nextBlob () );
        return false;
    }
    
private: 
    ngs::vdb::VdbReadCollection      m_coll;
    ngs::vdb::FragmentBlobIterator   m_blobIt;
    uint64_t                         m_startInBlob;
};

//////////////////// VdbSearch

bool VdbSearch :: logResults = false;
bool VdbSearch :: useBlobSearch = false;

void 
VdbSearch :: CheckArguments ( bool p_isExpression, unsigned int p_minScorePct) throw ( invalid_argument )
{
    if ( p_isExpression && m_algorithm != NucStrstr ) 
    {
        throw invalid_argument ( "query expressions are only supported for NucStrstr" );
    }
    if ( p_minScorePct != 100 )
    {
        switch ( m_algorithm )
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

VdbSearch :: VdbSearch ( Algorithm p_algorithm, const std::string& p_query, bool p_isExpression, unsigned int p_minScorePct, unsigned int p_threads )  
    throw ( invalid_argument )
:   m_algorithm ( p_algorithm ),
    m_query ( p_query ),
    m_isExpression ( p_isExpression ),
    m_minScorePct ( p_minScorePct ),
    m_threads ( p_threads ),
    m_output ( 0 ),
    m_searchBlock ( 0 )
{
    CheckArguments ( p_isExpression, p_minScorePct );
}

VdbSearch :: VdbSearch ( const string& p_algorithm, const std::string& p_query, bool p_isExpression, unsigned int p_minScorePct, unsigned int p_threads )  
    throw ( invalid_argument )
:   m_query ( p_query ),
    m_isExpression ( p_isExpression ),
    m_minScorePct ( p_minScorePct ),
    m_threads ( p_threads ),
    m_output ( 0 ),
    m_searchBlock ( 0 )
{
    if ( ! SetAlgorithm ( p_algorithm ) )
    {
        throw invalid_argument ( string ( "unrecognized algorithm: " ) + p_algorithm );
    }
    CheckArguments ( p_isExpression, p_minScorePct );
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
    delete m_searchBlock;
    delete m_output;    
}

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

bool 
VdbSearch :: SetAlgorithm ( const std :: string& p_algStr )
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


void 
VdbSearch :: AddAccession ( const string& p_accession ) throw ( ErrorMsg )
{
    SearchBlock* sb = SearchBlockFactory ( m_query, m_isExpression, m_algorithm, m_minScorePct );
    if ( useBlobSearch )
    {
        m_searches . push ( new BlobMatchIterator ( sb, p_accession ) );
    }
    else
    {
        m_searches . push ( new FragmentMatchIterator ( sb, p_accession ) );
    }
}

rc_t CC VdbSearch :: SearchThread ( const KThread *self, void *data )
{
    assert ( data );
    SearchThreadBlock& sb = * reinterpret_cast < SearchThreadBlock* > ( data );
    assert ( sb . m_lock );
    
    while ( true )
    {
        KLockAcquire ( sb . m_lock );
        if ( sb . m_search . empty () )
        {
            KLockUnlock ( sb . m_lock );
            break;
        }

        MatchIterator* it = sb . m_search . front ();
        sb . m_search . pop ();
        
        KLockUnlock ( sb . m_lock );
        
        string fragmentId;
        while ( it -> NextMatch ( fragmentId ) )
        {
            sb . m_output . Push ( it -> AccessionName (), fragmentId );
        }
        delete it;
    }
    
    sb . m_output . ProducerDone();
    return 0;    
}

bool 
VdbSearch :: NextMatch ( string& p_accession, string& p_fragmentId ) throw ( ErrorMsg )
{
    if ( m_threads > 0 )
    {   
        if ( m_output == 0 ) // first call to NextMatch() - set up
        {
            size_t threadNum = min<size_t> ( m_threads, m_searches . size () );
            m_output = new OutputQueue ( threadNum );
            m_searchBlock = new SearchThreadBlock ( m_searches, *m_output ); 
            for ( unsigned  int i = 0 ; i != threadNum; ++i )
            {
                KThread* t;
                rc_t rc = KThreadMake ( &t, SearchThread, m_searchBlock );
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
            if ( m_searches . front () -> NextMatch ( p_fragmentId ) )
            {
                p_accession = m_searches. front () -> AccessionName ();
                return true;
            }
            delete m_searches . front ();
            m_searches . pop ();
        }
        
        return false;
    }
}

//////////////////// SearchBlock factory

VdbSearch :: SearchBlock* 
VdbSearch :: SearchBlockFactory ( const string& p_query, bool p_isExpression, Algorithm p_algorithm, unsigned int p_minScorePct )
{
    switch ( p_algorithm )
    {
        case VdbSearch :: FgrepDumb: 
        case VdbSearch :: FgrepBoyerMoore:
        case VdbSearch :: FgrepAho:
            return new FgrepSearch ( p_query, p_algorithm );
            
        case VdbSearch :: AgrepDP:
        case VdbSearch :: AgrepWuManber:
        case VdbSearch :: AgrepMyers:
        case VdbSearch :: AgrepMyersUnltd:
            return new AgrepSearch ( p_query, p_algorithm, p_minScorePct );
        
        case VdbSearch :: NucStrstr:
            return new NucStrstrSearch ( p_query, p_isExpression );
            
        case VdbSearch :: SmithWaterman:
            return new SmithWatermanSearch ( p_query, p_minScorePct );
            
        default:
            throw ( ErrorMsg ( "SearchBlockFactory: unsupported algorithm" ) );
    }
}
