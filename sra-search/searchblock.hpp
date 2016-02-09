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

#ifndef _hpp_searchblock_
#define _hpp_searchblock_

#include "vdb-search.hpp"

struct Fgrep;
struct Agrep;
union NucStrstr;
struct SmithWaterman;

//////////////////// VdbSearch :: SearchBlock subclasses

class FgrepSearch : public VdbSearch :: SearchBlock
{   
public:
    FgrepSearch ( const std::string& p_query, VdbSearch :: Algorithm p_algorithm ) throw ( ngs :: ErrorMsg );
    virtual ~FgrepSearch (); 
    
    virtual bool FirstMatch ( const char* p_bases, size_t p_size, uint64_t& p_hitStart, uint64_t& p_hitEnd ) throw ( ngs :: ErrorMsg );

private:
    struct Fgrep*   m_fgrep;
    const char*     m_query[1];
};

class AgrepSearch : public VdbSearch :: SearchBlock
{   
public:
    AgrepSearch ( const std::string& p_query, VdbSearch :: Algorithm p_algorithm, uint8_t p_minScorePct ) throw ( ngs :: ErrorMsg );
    virtual ~AgrepSearch ();
    
    virtual bool FirstMatch ( const char* p_bases, size_t p_size, uint64_t& p_hitStart, uint64_t& p_hitEnd ) throw ( ngs :: ErrorMsg );

private:
    struct Agrep*   m_agrep;
    const char*     m_query;
    uint8_t         m_minScorePct;
};

class NucStrstrSearch : public VdbSearch :: SearchBlock
{   
public:
    NucStrstrSearch ( const std::string& p_query, bool p_positional = false ) throw ( ngs :: ErrorMsg );
    virtual ~NucStrstrSearch ();
    
    virtual bool FirstMatch ( const char* p_bases, uint64_t p_size ) throw ( ngs :: ErrorMsg );  
    virtual bool FirstMatch ( const char* p_bases, size_t p_size, uint64_t& p_hitStart, uint64_t& p_hitEnd ) throw ( ngs :: ErrorMsg ); // will throw if not positional
    
private:
    static void ConvertAsciiTo2NAPacked ( const char* pszRead, size_t nReadLen, unsigned char* pBuf2NA, size_t nBuf2NASize );
    
    bool                m_positional;
    std::string         m_query;
    union NucStrstr*    m_nss;
};

class SmithWatermanSearch : public VdbSearch :: SearchBlock
{   
public:
    SmithWatermanSearch ( const std::string& p_query, uint8_t p_minScorePct );
    virtual ~SmithWatermanSearch ();
    
    virtual bool FirstMatch ( const char* p_bases, size_t p_size, uint64_t& p_hitStart, uint64_t& p_hitEnd ) throw ( ngs :: ErrorMsg );  

private:
    const char*             m_query;
    size_t                  m_querySize;
    size_t                  m_matrixSize;
    uint8_t                 m_minScorePct;
    struct SmithWaterman*   m_sw;
};

#endif
