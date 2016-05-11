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

#ifndef _hpp_vdb_search_
#define _hpp_vdb_search_

#include <string>
#include <vector>
#include <queue>
#include <stdexcept>

#include <klib/rc.h>

#include <ngs/ReadCollection.hpp>

struct KThread;

class VdbSearch
{
public:
    // base class of a hierarchy implementing various search algorithms 
    class SearchBlock
    {   
    public:
        virtual ~SearchBlock () {}
        
        virtual bool FirstMatch ( const char* p_bases, size_t p_size ) throw ( ngs :: ErrorMsg )
        {
            uint64_t hitStart;
            uint64_t hitEnd;
            return FirstMatch ( p_bases, p_size, hitStart, hitEnd );
        }
        virtual bool FirstMatch ( const char* p_bases, size_t p_size, uint64_t& hitStart, uint64_t& hitEnd ) throw ( ngs :: ErrorMsg ) = 0;
    };
    
public:
    static bool logResults;

public:
    // Search algorithms supported by this class
    typedef enum
    {
        FgrepDumb = 0, 
        FgrepBoyerMoore,
        FgrepAho,
        AgrepDP,
        AgrepWuManber,
        AgrepMyers,
        AgrepMyersUnltd,
        NucStrstr,
        SmithWaterman,
    } Algorithm;

    typedef std :: vector < std :: string >  SupportedAlgorithms;

    // enum Algorithm values correspond to indexes in the container returned here
    static SupportedAlgorithms GetSupportedAlgorithms (); 
    
public:
    VdbSearch ( Algorithm, 
                const std::string& query, 
                bool isExpression, 
                bool p_useBlobSearch, 
                unsigned int p_minScorePct = 100, 
                unsigned int p_threads = 0 ) throw ( std :: invalid_argument );
    VdbSearch ( const std :: string& algorithm, 
                const std::string& query, 
                bool isExpression, 
                bool p_useBlobSearch, 
                unsigned int p_minScorePct = 100, 
                unsigned int p_threads = 0 ) throw ( std :: invalid_argument );
    
    ~VdbSearch ();
    
    Algorithm GetAlgorithm () const { return m_algorithm; }
    
    void AddAccession ( const std::string& ) throw ( ngs :: ErrorMsg );
    
    bool NextMatch ( std::string& accession, std::string& fragmentId ) throw ( ngs :: ErrorMsg );
    
private:

    // create search blocks corresponding to current configuration
    class SearchBlockFactory
    {
    public:
        SearchBlockFactory ( const std :: string&   p_query, 
                             bool                   p_isExpression, 
                             bool                   p_useBlobSearch, 
                             Algorithm              p_algorithm, 
                             unsigned int           p_minScorePct );
                             
        bool CanUseBlobs () const;
        void SetAllgorithm ( Algorithm p_algorithm) { m_algorithm = p_algorithm; }
                             
        SearchBlock* MakeSearchBlock () const;
        
    private:
        const std :: string&    m_query; 
        bool                    m_isExpression;
        bool                    m_useBlobSearch; 
        Algorithm               m_algorithm; 
        unsigned int            m_minScorePct;
    };
    
private:

    // A buffer that can be searched for one or more matches, bound to an accession and an engine-side algorithm;
    // Subclasses implement fragment-based (1 hit) or blob-based (multiple hits) search 
    class SearchBuffer
    {
    public:
        SearchBuffer ( SearchBlock*, const std::string& accession );
        virtual ~SearchBuffer();
        
        virtual bool NextMatch ( std::string& fragmentId ) = 0;
        
        std::string AccessionName () const { return m_accession; }
        
    protected: 
        SearchBlock* m_searchBlock;  // owned here
        std::string m_accession; 
    };
    
    class FragmentSearchBuffer;
    class BlobSearchBuffer; 

private:

    // a VDB-agnostic iterator. Subclasses implement different styles of retrieval (fragment/blob based, SRA/WGS schema)
    class MatchIterator 
    {
    public:
        MatchIterator ( SearchBlockFactory&, const std::string& accession );
        virtual ~MatchIterator ();
        
        virtual SearchBuffer* NextBuffer () = 0;
        
    protected:
        const SearchBlockFactory&   m_factory; 
        std::string                 m_accession;
    };
    
    class FragmentMatchIterator; // fragment-based search
    class BlobMatchIterator; // blob-based search
    
private:

    typedef std::queue < MatchIterator* > SearchQueue;
    typedef std :: vector < KThread* > ThreadPool;
    
    struct SearchThreadBlock;
    class OutputQueue;
    
    static rc_t CC SearchAccPerThread ( const struct KThread *self, void *data );
    static rc_t CC SearchBlobPerThread ( const struct KThread *self, void *data );

private:

    void CheckArguments ( bool isExpression, unsigned int minScorePct) throw ( std :: invalid_argument );
    
    bool SetAlgorithm ( const std :: string& p_algStr );

private:
    std::string     m_query;
    bool            m_isExpression;
    bool            m_useBlobSearch;  
    Algorithm       m_algorithm;
    unsigned int    m_minScorePct;
    unsigned int    m_threads;
    bool            m_blobPerThread;  
    
    SearchBlockFactory m_sbFactory;    

    SearchBuffer*   m_buf;

    SearchQueue     m_searches;
        
    // used with multi-threading        
    OutputQueue*                m_output;  
    struct SearchThreadBlock*   m_searchBlock;
    ThreadPool                  m_threadPool;
};

#endif
