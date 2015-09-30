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

using namespace std;
using namespace ngs;

bool VdbSearch :: logResults = false;

VdbSearch :: VdbSearch ()
:   m_algorithm ( Default )
{
}

VdbSearch :: ~VdbSearch ()
{
    while ( ! m_searches . empty () )
    {
        delete m_searches . front ();
        m_searches . pop ();
    }
}

void 
VdbSearch :: SetQuery ( const string& p_query )
{
    m_query = p_query;
}

void 
VdbSearch :: SetAlgorithm ( Algorithm p_algorithm )
{
    m_algorithm = p_algorithm;
}

bool 
VdbSearch :: SetAlgorithm ( const std :: string& p_algStr )
{
    if      ( p_algStr == "FgrepDumb" )         SetAlgorithm ( VdbSearch :: FgrepDumb );
    else if ( p_algStr == "FgrepBoyerMoore" )   SetAlgorithm ( VdbSearch :: FgrepBoyerMoore );
    else if ( p_algStr == "FgrepAho" )          SetAlgorithm ( VdbSearch :: FgrepAho );
    else if ( p_algStr == "AgrepDP" )           SetAlgorithm ( VdbSearch :: AgrepDP );
    else if ( p_algStr == "AgrepWuManber" )     SetAlgorithm ( VdbSearch :: AgrepWuManber );
    else if ( p_algStr == "AgrepMyers" )        SetAlgorithm ( VdbSearch :: AgrepMyers );
    else if ( p_algStr == "AgrepMyersUnltd" )   SetAlgorithm ( VdbSearch :: AgrepMyersUnltd );
    else if ( p_algStr == "NucStrstr" )         SetAlgorithm ( VdbSearch :: NucStrstr );
    else return false;
    
    return true;
}


void 
VdbSearch :: AddAccession ( const string& p_accession ) throw ( ErrorMsg )
{
    m_searches . push ( new MatchIterator ( SearchBlockFactory ( m_query, m_algorithm ), p_accession ) );
}

bool 
VdbSearch :: NextMatch ( string& p_accession, string& p_fragmentId ) throw ( ErrorMsg )
{
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

class FgrepSearch : public SraSearch :: SearchBlock
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

class AgrepSearch : public SraSearch :: SearchBlock
{   
public:
    AgrepSearch ( const string& p_query, VdbSearch :: Algorithm p_algorithm )
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
        bool ret = AgrepFindFirst ( m_agrep, 0, p_bases, p_size, & matchinfo ) != 0; // 0 = perfect match
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
};


class NucStrstrBlock : public SraSearch :: SearchBlock
{   
public:
    NucStrstrBlock ( const string& p_query )
    {
        
        m_query = p_query . c_str(); // this object will not outlive it master who owns the query string

        rc_t rc = NucStrstrMake ( &m_nss, 0, m_query, p_query . size () );
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

//////////////////// SearchBlock factory

SraSearch :: SearchBlock* 
VdbSearch :: SearchBlockFactory ( const string& p_query, Algorithm p_algorithm )
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
            return new AgrepSearch ( p_query, p_algorithm );
        
        case VdbSearch :: NucStrstr:
            return new NucStrstrBlock ( p_query );
            
        default:
            throw ( ErrorMsg ( "SearchBlockFactory: unsupported algorithm" ) );
    }
}
