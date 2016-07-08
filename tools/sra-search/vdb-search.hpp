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

#include "searchblock.hpp"

struct KThread;
class MatchIterator;
class SearchBuffer;

class VdbSearch
{
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

    struct Settings
    {
        Algorithm                   m_algorithm;    // default FgrepDumb
        std::string                 m_query;
        std::vector < std::string > m_accessions;
        bool                        m_isExpression;     // default false
        unsigned int                m_minScorePct;      // default 100
        unsigned int                m_threads;          // default 2
        bool                        m_useBlobSearch;    // default true
        bool                        m_referenceDriven;  // default false

        Settings ();
        bool SetAlgorithm ( const std :: string& algorithm );
    };

public:
    VdbSearch ( const Settings& settings ) throw ( std :: invalid_argument );
    ~VdbSearch ();

    bool NextMatch ( std::string& accession, std::string& fragmentId );

private:

    typedef std::queue < MatchIterator* > SearchQueue;
    typedef std :: vector < KThread* > ThreadPool;

    struct SearchThreadBlock;
    class OutputQueue;

    // create search blocks corresponding to current configuration
    class SearchBlockFactory : public SearchBlock :: Factory
    {
    public:
        SearchBlockFactory ( const Settings& p_settings );

        virtual SearchBlock* MakeSearchBlock () const;

    private:
        const Settings& m_settings; // not a copy, since the settings may be changed post-creation
    };

    static rc_t CC SearchAccPerThread ( const struct KThread *self, void *data );
    static rc_t CC SearchBlobPerThread ( const struct KThread *self, void *data );

private:
    Settings m_settings;

    SearchBlockFactory m_sbFactory;

    SearchBuffer*   m_buf;

    SearchQueue     m_searches;

    // used with multi-threading
    OutputQueue*                m_output;
    struct SearchThreadBlock*   m_searchBlock;
    ThreadPool                  m_threadPool;
};

#endif
