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

#include "searchblock.hpp"

#include <search/grep.h>
#include <search/nucstrstr.h>
#include <search/smith-waterman.h>

#include <cstring>

using namespace std;
using namespace ngs;

//////////////////// VdbSearch :: SearchBlock subclasses

FgrepSearch :: FgrepSearch ( const string& p_query, VdbSearch :: Algorithm p_algorithm ) throw ( ngs :: ErrorMsg )
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
        throw ( ErrorMsg ( "FgrepMake() failed" ) );
    }
}

FgrepSearch :: ~FgrepSearch ()
{
    FgrepFree ( m_fgrep );
}    
    
bool
FgrepSearch :: FirstMatch ( const char* p_bases, size_t p_size, uint64_t& p_hitStart, uint64_t& p_hitEnd ) throw ( ngs :: ErrorMsg )
{
    FgrepMatch matchinfo;
    bool ret = FgrepFindFirst ( m_fgrep, p_bases, p_size, & matchinfo ) != 0;
    if ( ret )
    {
        p_hitStart = matchinfo . position;
        p_hitEnd = p_hitStart + matchinfo . length;
        if ( VdbSearch :: logResults )
        {
            cout << "Pattern='" << m_query[0] << "' " 
                 << " Substr(" << matchinfo . position 
                 << ", " << matchinfo . length << ") = '" 
                 << string ( p_bases + matchinfo . position, matchinfo . length ) << "'" 
                 << endl;
        }
    }
    return ret;
}

AgrepSearch :: AgrepSearch ( const string& p_query, VdbSearch :: Algorithm p_algorithm, uint8_t p_minScorePct ) throw ( ngs :: ErrorMsg )
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

AgrepSearch ::  ~AgrepSearch ()
{
    AgrepWhack ( m_agrep );
}    
    
bool 
AgrepSearch :: FirstMatch ( const char* p_bases, size_t p_size, uint64_t& p_hitStart, uint64_t& p_hitEnd ) throw ( ngs :: ErrorMsg )
{
    AgrepMatch matchinfo;
    bool ret = AgrepFindFirst ( m_agrep, 
                                strlen ( m_query ) * ( 100 - m_minScorePct ) / 100, // 0 = perfect match
                                p_bases, 
                                p_size, 
                                & matchinfo ) != 0; 
    if ( ret )
    {
        p_hitStart = matchinfo . position;
        p_hitEnd = p_hitStart + matchinfo . length;       
        if ( VdbSearch :: logResults )
        {
            cout << "Pattern='" << m_query << "' " 
                 << "Score=" << matchinfo . score 
                 << " Substr(" << matchinfo . position 
                 << ", " << matchinfo . length << ") = '" 
                 << string ( p_bases + matchinfo . position, matchinfo . length ) << "'" 
                 << endl;
        }
    }
    return ret;
}

NucStrstrSearch :: NucStrstrSearch ( const string& p_query, bool p_positional, bool p_useBlobSearch ) throw ( ngs :: ErrorMsg )
:   m_positional ( p_positional || p_useBlobSearch ), // always use positional mode when searching blob-by-blob since it is the only one reporting position of the match
    m_query ( p_query )
{
    rc_t rc = NucStrstrMake ( &m_nss, m_positional ? 1 : 0, m_query . c_str (), m_query . size () ); 
    if ( rc != 0 )
    {
        throw ( ErrorMsg ( "NucStrstrMake() failed" ) );
    }
}

NucStrstrSearch ::  ~NucStrstrSearch ()
{
    NucStrstrWhack ( m_nss );
}

bool 
NucStrstrSearch :: FirstMatch ( const char* p_bases, size_t p_size ) throw ( ngs :: ErrorMsg )
{
    // convert p_bases to 2na packed since nucstrstr works with that format only
    const size_t bufSize = p_size / 4 + 1 + 16; // NucStrstrSearch expects the buffer to be at least 16 bytes longer than the sequence
    unsigned char* buf2na = new unsigned char [ bufSize ];
    ConvertAsciiTo2NAPacked ( p_bases, p_size, buf2na, bufSize );        
        
    unsigned int selflen;
    bool ret = ::NucStrstrSearch ( m_nss, reinterpret_cast < const void * > ( buf2na ), 0, p_size, & selflen ) != 0;
    delete [] buf2na;
    return ret;
}
    
bool 
NucStrstrSearch :: FirstMatch ( const char* p_bases, size_t p_size, uint64_t& p_hitStart, uint64_t& p_hitEnd ) throw ( ngs :: ErrorMsg )
{
    if ( ! m_positional )
    {   // Should not land here since NucStrstrSearch::CanUseBlobs() returns false.
        throw ( ErrorMsg ( "NucStrstrSearch: non-positional search in a blob is not supported" ) );
    }
    
    // convert p_bases to 2na packed since nucstrstr works with that format only
    const size_t bufSize = p_size / 4 + 1 + 16; // NucStrstrSearch expects the buffer to be at least 16 bytes longer than the sequence
    unsigned char* buf2na = new unsigned char [ bufSize ];
    ConvertAsciiTo2NAPacked ( p_bases, p_size, buf2na, bufSize );        
        
    unsigned int selflen;
    int pos = ::NucStrstrSearch ( m_nss, reinterpret_cast < const void * > ( buf2na ), 0, p_size, & selflen );
    bool ret = pos > 0;
    if ( ret )
    {
        p_hitStart = pos - 1;
        p_hitEnd = p_hitStart + selflen;       
    }
    delete [] buf2na;
    return ret;
}
    
void
NucStrstrSearch :: ConvertAsciiTo2NAPacked ( const char* pszRead, size_t nReadLen, unsigned char* pBuf2NA, size_t nBuf2NASize )
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

SmithWatermanSearch :: SmithWatermanSearch ( const string& p_query, uint8_t p_minScorePct )
:   m_query ( p_query . c_str() ),
    m_querySize ( p_query . size() ),
    m_matrixSize ( 0 ),
    m_minScorePct ( p_minScorePct )        
{
    rc_t rc = SmithWatermanMake ( &m_sw, p_query . c_str() );
    if ( rc != 0 )
    {
        throw ( ErrorMsg ( "SmithWatermanMake() failed" ) );
    }
}

SmithWatermanSearch :: ~SmithWatermanSearch ()
{
    SmithWatermanWhack ( m_sw );
}    

bool
SmithWatermanSearch :: CanUseBlobs () const 
{ 
    // As blob size grows, SW quickly becomes unusable (VDB-3019)
    return false; 
}

bool
SmithWatermanSearch :: FirstMatch ( const char* p_bases, size_t p_size, uint64_t& p_hitStart, uint64_t& p_hitEnd ) throw ( ngs :: ErrorMsg )
{
    SmithWatermanMatch matchinfo;
    unsigned int scoreThreshold = ( m_querySize * 2 ) * m_minScorePct / 100; // m_querySize * 2 == exact match
    rc_t rc = SmithWatermanFindFirst ( m_sw, scoreThreshold, p_bases, p_size, & matchinfo );     
    if ( rc == 0 )
    {
        if ( VdbSearch :: logResults )
        {
            cout << "Pattern='" << m_query << "' " 
                    << "Score=" << matchinfo . score 
                    << " scoreThreshold=" << scoreThreshold 
                    << endl;
        }
        p_hitStart = matchinfo . position;
        p_hitEnd = p_hitStart + matchinfo . length;
        return true;        
    }
    else if ( GetRCObject ( rc ) == (RCObject)rcQuery  && GetRCState ( rc ) == rcNotFound )
    {
        return false;
    }
    throw ( ErrorMsg ( "SmithWatermanFindFirst() failed" ) );
}
