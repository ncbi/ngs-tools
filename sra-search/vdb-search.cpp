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

using namespace std;
using namespace ngs;

class VdbSearch :: Search 
{
public:
    Search( const std::string& p_query, VdbSearch :: Algorithm, const std::string& accession );
    ~Search ();
    
    bool NextMatch ( std::string& fragmentId );
    
    std::string AccessionName () const;
    
private: 
    ngs::ReadCollection m_coll;
    ngs::ReadIterator   m_readIt;
    
    // for now, Fgrep-oriented
    struct Fgrep*   m_fgrep;
    const char*     m_query[1];
};

VdbSearch :: VdbSearch ( const string& p_query )
:   m_query ( p_query ),
    m_algorithm ( Default )
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
VdbSearch :: SetAlgorithm ( Algorithm p_algorithm )
{
    m_algorithm = p_algorithm;
}

bool 
VdbSearch :: SetAlgorithm ( const std :: string& p_algStr )
{
    if ( p_algStr == "FgrepDumb" )
    {
        SetAlgorithm ( VdbSearch :: FgrepDumb );
    }
    else if ( p_algStr == "FgrepBoyerMoore" )
    {
        SetAlgorithm ( VdbSearch :: FgrepBoyerMoore );
    }
    else
    {
        return false;
    }
    return true;
}


void 
VdbSearch :: AddAccession ( const string& p_accession ) throw ( ErrorMsg )
{
    m_searches . push ( new Search ( m_query, m_algorithm, p_accession ) );
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


//////////////////// VdbSearch :: Search
VdbSearch :: Search :: Search( const string& p_query, VdbSearch :: Algorithm p_algorithm, const string& p_accession )
:   m_coll ( ncbi :: NGS :: openReadCollection ( p_accession ) ), 
    m_readIt ( m_coll . getReads ( Read :: all ) )
{
    m_query[0] = p_query . c_str(); // this object will not outlive it master who owns the query string
    
    switch ( p_algorithm )
    {
    case VdbSearch :: FgrepDumb:
        {
            rc_t rc = FgrepMake ( & m_fgrep, FGREP_MODE_ASCII | FGREP_ALG_DUMB, m_query, 1 );
            if ( rc != 0 )
            {
                throw ( ErrorMsg ( "FgrepMake failed" ) );
            }
        }
        break;
    case VdbSearch :: FgrepBoyerMoore:
        {
            rc_t rc = FgrepMake ( & m_fgrep, FGREP_MODE_ASCII | FGREP_ALG_BOYERMOORE, m_query, 1 );
            if ( rc != 0 )
            {
                throw ( ErrorMsg ( "FgrepMake failed" ) );
            }
        }
        break;
    default:
        throw ( ErrorMsg ( "unsupported VdbSearch :: Algorithm" ) );
    }
}

VdbSearch :: Search :: ~Search ()
{
    FgrepFree ( m_fgrep );
}    
        
bool 
VdbSearch :: Search :: NextMatch ( string& p_fragmentId )
{
    while ( m_readIt . nextRead () )
    {
        while ( m_readIt . nextFragment () )
        {
            StringRef bases = m_readIt . getFragmentBases ();
            FgrepMatch matchinfo;
            if ( FgrepFindFirst ( m_fgrep, bases . data (), bases . size (), & matchinfo ) != 0 )
            {
                p_fragmentId = m_readIt . getFragmentId () . toString ();
                return true;
            }
        }
    }
    return false;
}
        
string 
VdbSearch :: Search :: AccessionName () const
{
    return m_coll . getName ();
}

