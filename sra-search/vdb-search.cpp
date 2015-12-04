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
#include <search/grep.h>
#include <search/nucstrstr.h>
#include <../libs/search/search-priv.h>

#include <atomic32.h>

#include <klib/time.h>

#include <kproc/thread.h>
#include <kproc/lock.h>

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


//////////////////// VdbSearch

bool VdbSearch :: logResults = false;

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
    ALG ( FgrepDumb ),
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
    m_searches . push ( new MatchIterator ( SearchBlockFactory ( m_query, m_isExpression, m_algorithm, m_minScorePct ), 
                                            p_accession ) );
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

//////////////////// VdbSearch :: MatchIterator
VdbSearch :: MatchIterator :: MatchIterator ( SearchBlock* p_block, const std::string& p_accession )
:   m_coll ( ncbi :: NGS :: openReadCollection ( p_accession ) ), 
    m_searchBlock ( p_block ),
    m_readIt ( m_coll . getReads ( Read :: all ) )
{
    m_readIt . nextRead ();
}

VdbSearch :: MatchIterator :: ~MatchIterator ()
{
    delete m_searchBlock;
}

bool 
VdbSearch :: MatchIterator :: NextMatch ( std::string& p_fragmentId )
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

std::string 
VdbSearch :: MatchIterator :: AccessionName () const
{
    return m_coll . getName ();
}

//////////////////// SearchBlock subclasses

class FgrepSearch : public VdbSearch :: SearchBlock
{   
public:
    FgrepSearch ( const string& p_query, VdbSearch :: Algorithm p_algorithm )
    {
        m_query[0] = p_query . c_str(); // this object will not outlive it master who owns the query string
        
        rc_t rc = 0;
        switch ( p_algorithm )
        {
        case VdbSearch :: FgrepDumb:
            rc = FgrepMake ( & m_fgrep, FGREP_MODE_ACGT | FGREP_ALG_DUMB, m_query, 1 );
            break;
        case VdbSearch :: FgrepBoyerMoore:
            rc = FgrepMake ( & m_fgrep, FGREP_MODE_ACGT | FGREP_ALG_BOYERMOORE, m_query, 1 );
            break;
        case VdbSearch :: FgrepAho:
            rc = FgrepMake ( & m_fgrep, FGREP_MODE_ACGT | FGREP_ALG_AHOCORASICK, m_query, 1 );
            break;
        default:
            throw ( ErrorMsg ( "FgrepSearch: unsupported algorithm" ) );
        }
        if ( rc != 0 )
        {
            throw ( ErrorMsg ( "FgrepMake failed" ) );
        }
    }
    virtual ~FgrepSearch ()
    {
        FgrepFree ( m_fgrep );
    }    
    
    virtual bool FirstMatch ( const char* p_bases, size_t p_size )
    {
        FgrepMatch matchinfo;
        return  FgrepFindFirst ( m_fgrep, p_bases, p_size, & matchinfo ) != 0;
    }

private:
    struct Fgrep*   m_fgrep;
    const char*     m_query[1];
};

class AgrepSearch : public VdbSearch :: SearchBlock
{   
public:
    AgrepSearch ( const string& p_query, VdbSearch :: Algorithm p_algorithm, uint8_t p_minScorePct )
    : m_minScorePct ( p_minScorePct )
    {
        m_query = p_query . c_str(); // this object will not outlive it master who owns the query string
        
        rc_t rc = 0;
        switch ( p_algorithm )
        {
        case VdbSearch :: AgrepDP:
            rc = AgrepMake ( & m_agrep, AGREP_MODE_ASCII | AGREP_ALG_DP, m_query );
            break;
        case VdbSearch :: AgrepWuManber:
            rc = AgrepMake ( & m_agrep, AGREP_MODE_ASCII | AGREP_ALG_WUMANBER, m_query );
            break;
        case VdbSearch :: AgrepMyers:
            rc = AgrepMake ( & m_agrep, AGREP_MODE_ASCII | AGREP_ALG_MYERS, m_query );
            break;
        case VdbSearch :: AgrepMyersUnltd:
            rc = AgrepMake ( & m_agrep, AGREP_MODE_ASCII | AGREP_ALG_MYERS_UNLTD, m_query );
            break;
        default:
            throw ( ErrorMsg ( "AgrepSearch: unsupported algorithm" ) );
        }
        if ( rc != 0 )
        {
            throw ( ErrorMsg ( "AgrepMake failed" ) );
        }
    }
    virtual ~AgrepSearch ()
    {
        AgrepWhack ( m_agrep );
    }    
    
    virtual bool FirstMatch ( const char* p_bases, size_t p_size )
    {
        AgrepMatch matchinfo;
        bool ret = AgrepFindFirst ( m_agrep, 
                                    strlen ( m_query ) * ( 100 - m_minScorePct ) / 100, // 0 = perfect match
                                    p_bases, 
                                    p_size, 
                                    & matchinfo ) != 0; 
        if ( ret && VdbSearch :: logResults )
        {
            cout << "Pattern='" << m_query << "' " 
                 << "Score=" << matchinfo . score 
                 << " Substr(" << matchinfo . position 
                 << ", " << matchinfo . length << ") = '" 
                 << string ( p_bases + matchinfo . position, matchinfo . length ) << "'" 
                 << endl;
        }
        return ret;
    }

private:
    struct Agrep*   m_agrep;
    const char*     m_query;
    uint8_t         m_minScorePct;
};


class NucStrstrBlock : public VdbSearch :: SearchBlock
{   
public:
    NucStrstrBlock ( const string& p_query, bool p_positional = false )
    {
        
        m_query = p_query . c_str(); // this object will not outlive it master who owns the query string

        rc_t rc = NucStrstrMake ( &m_nss, p_positional ? 1 : 0, m_query, p_query . size () );
        if ( rc != 0 )
        {
            throw ( ErrorMsg ( "NucStrstrMake failed" ) );
        }
    }
    virtual ~NucStrstrBlock ()
    {
        NucStrstrWhack ( m_nss );
    }    
    
    virtual bool FirstMatch ( const char* p_bases, size_t p_size )
    {
        // convert p_bases to 2na packed since nucstrstr works with that format only
        const size_t bufSize = p_size / 4 + 1 + 16; // NucStrstrSearch expects the buffer to be at least 16 bytes longer than the sequence
        unsigned char* buf2na = new unsigned char [ bufSize ];
        ConvertAsciiTo2NAPacked ( p_bases, p_size, buf2na, bufSize );        
            
        unsigned int selflen;
        bool ret = NucStrstrSearch ( m_nss, reinterpret_cast < const void * > ( buf2na ), 0, p_size, & selflen ) != 0;
        delete [] buf2na;
        return ret;
    }
    
    static void ConvertAsciiTo2NAPacked ( const char* pszRead, size_t nReadLen, unsigned char* pBuf2NA, size_t nBuf2NASize )
    {
        static unsigned char map [ 1 << ( sizeof ( char ) * 8 ) ];
        map['A'] = map['a'] = 0;
        map['C'] = map['c'] = 1;
        map['G'] = map['g'] = 2;
        map['T'] = map['t'] = 3;
        
        static size_t shiftLeft [ 4 ] = { 6, 4, 2, 0 };

        fill ( pBuf2NA, pBuf2NA + nBuf2NASize, 0 );
        
        for ( size_t iChar = 0; iChar < nReadLen; ++iChar )
        {
            size_t iByte = iChar / 4;
            if ( iByte > nBuf2NASize )
            {
                assert ( false );
                break;
            }

            pBuf2NA[iByte] |= map [ size_t ( pszRead [ iChar ] ) ] << shiftLeft [ iChar % 4 ];
        }
    }

private:
    NucStrstr*   m_nss;
    const char*  m_query;
};

class SmithWatermanSearch : public VdbSearch :: SearchBlock
{   
public:
    SmithWatermanSearch ( const string& p_query, uint8_t p_minScorePct )
    :   m_query ( p_query . c_str() ),
        m_querySize ( p_query . size() ),
        m_matrix ( 0 ),
        m_matrixSize ( 0 ),
        m_minScorePct ( p_minScorePct )        
    {
    }
    virtual ~SmithWatermanSearch ()
    {
        delete [] m_matrix;
    }    
    
    virtual bool FirstMatch ( const char* p_bases, size_t p_size )
    {
        const size_t Rows = p_size + 1;
        const size_t Cols = m_querySize + 1;
        const size_t MatrixSize = Rows * Cols;
        if ( MatrixSize > m_matrixSize )
        {
            delete [] m_matrix;
            m_matrix = new int [ MatrixSize ];
            m_matrixSize = MatrixSize;
        }
        int maxScore = -1;
        rc_t rc = calculate_similarity_matrix ( p_bases, p_size, m_query, m_querySize, m_matrix, true, & maxScore ); 
        if ( rc != 0 )
        {
            throw ( ErrorMsg ( "calculate_similarity_matrix failed" ) );
        }
        unsigned int scoreThreshold = ( m_querySize * 2 ) * m_minScorePct / 100; // m_querySize * 2 == exact match
        bool ret = maxScore >= scoreThreshold;
        if ( ret && VdbSearch :: logResults )
        {
            cout << "Pattern='" << m_query << "' " 
                 << "Score=" << maxScore 
                 << " scoreThreshold=" << scoreThreshold 
                 << endl;
        }
        return ret;
    }

private:
    const char*     m_query;
    const size_t    m_querySize;
    int*            m_matrix;
    size_t          m_matrixSize;
    uint8_t         m_minScorePct;
};

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
            return new NucStrstrBlock ( p_query, p_isExpression );
            
        case VdbSearch :: SmithWaterman:
            return new SmithWatermanSearch ( p_query, p_minScorePct );
            
        default:
            throw ( ErrorMsg ( "SearchBlockFactory: unsupported algorithm" ) );
    }
}
