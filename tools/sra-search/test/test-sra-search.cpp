/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author'm_s official duties as a United States Government employee and
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
#include "searchblock.hpp"
#include "searchbuffer.hpp"
#include "referencematchiterator.hpp"

#include <set>

#include <ktst/unit_test.hpp>

#include <vdb/manager.h> // VDBManager
#include <vdb/database.h>
#include <vdb/table.h>
#include <vdb/cursor.h>
#include <vdb/blob.h>
#include <vdb/vdb-priv.h>
#include <sra/sraschema.h> // VDBManagerMakeSRASchema
#include <vdb/schema.h> /* VSchemaRelease */

#include <../libs/vdb/blob-priv.h>

#include <search/grep.h>

using namespace std;
using namespace ncbi::NK;
using namespace ngs;

TEST_SUITE(SraSearchTestSuite);
#if 0
// SearchBlock

TEST_CASE ( FgrepDumb )
{
    FgrepSearch sb ( "CTA", FgrepSearch :: FgrepDumb );
    uint64_t hitStart = 0;
    uint64_t hitEnd = 0;
    const string Bases = "ACTGACTAGTCA";
    REQUIRE ( sb.FirstMatch ( Bases.c_str(), Bases.size(), hitStart, hitEnd ) );
    REQUIRE_EQ ( (uint64_t)5, hitStart );
    REQUIRE_EQ ( (uint64_t)8, hitEnd );
}

TEST_CASE ( SearchFgrepBoyerMoore )
{
    FgrepSearch sb ( "CTA", FgrepSearch :: FgrepBoyerMoore );
    uint64_t hitStart = 0;
    uint64_t hitEnd = 0;
    const string Bases = "ACTGACTAGTCA";
    REQUIRE ( sb.FirstMatch ( Bases.c_str(), Bases.size(), hitStart, hitEnd ) );
    REQUIRE_EQ ( (uint64_t)5, hitStart );
    REQUIRE_EQ ( (uint64_t)8, hitEnd );
}

TEST_CASE ( SearchFgrepAho )
{
    FgrepSearch sb ( "CTA", FgrepSearch :: FgrepAho );
    uint64_t hitStart = 0;
    uint64_t hitEnd = 0;
    const string Bases = "ACTGACTAGTCA";
    REQUIRE ( sb.FirstMatch ( Bases.c_str(), Bases.size(), hitStart, hitEnd ) );
    REQUIRE_EQ ( (uint64_t)5, hitStart );
    REQUIRE_EQ ( (uint64_t)8, hitEnd );
}

TEST_CASE ( SearchAgrepDP )
{
    AgrepSearch sb ( "CTA", AgrepSearch :: AgrepDP, 100 );
    uint64_t hitStart = 0;
    uint64_t hitEnd = 0;
    const string Bases = "ACTGACTAGTCA";
    REQUIRE ( sb.FirstMatch ( Bases.c_str(), Bases.size(), hitStart, hitEnd ) );
    REQUIRE_EQ ( (uint64_t)5, hitStart );
    REQUIRE_EQ ( (uint64_t)8, hitEnd );
}

TEST_CASE ( SearchNucStrstr_NoExpr_NoCoords )
{
    NucStrstrSearch sb ( "CTA", false );
    const string Bases = "ACTGACTAGTCA";
    REQUIRE ( sb.FirstMatch ( Bases.c_str(), Bases.size() ) );
}

TEST_CASE ( SearchNucStrstr_NoExpr_Coords_NotSupported )
{
    NucStrstrSearch sb ( "CTA", false );
    uint64_t hitStart = 0;
    uint64_t hitEnd = 0;
    const string Bases = "ACTGACTAGTCA";
    REQUIRE_THROW ( sb.FirstMatch ( Bases.c_str(), Bases.size(), hitStart, hitEnd ) ); // not supported
}

TEST_CASE ( SearchNucStrstr_Expr_Coords )
{
    NucStrstrSearch sb ( "CTA", true );
    uint64_t hitStart = 0;
    uint64_t hitEnd = 0;
    const string Bases = "ACTGACTAGTCA";
    REQUIRE ( sb.FirstMatch ( Bases.c_str(), Bases.size(), hitStart, hitEnd ) );
    REQUIRE_EQ ( (uint64_t)5, hitStart );
    REQUIRE_EQ ( (uint64_t)8, hitEnd );
}

TEST_CASE ( SearchSmithWaterman_Coords_NotSupported )
{
    SmithWatermanSearch sb ( "CTA", 100 );
    uint64_t hitStart = 0;
    uint64_t hitEnd = 0;
    const string Bases = "ACTGACTAGTCA";
    REQUIRE ( sb.FirstMatch ( Bases.c_str(), Bases.size(), hitStart, hitEnd ) );
    REQUIRE_EQ ( (uint64_t)5, hitStart );
    REQUIRE_EQ ( (uint64_t)8, hitEnd );
}
#endif
// VdbSearch

class VdbSearchFixture
{
public:
    VdbSearchFixture ()
    : m_s ( 0 )
    {
        VdbSearch :: logResults = false;
    }
    ~VdbSearchFixture ()
    {
        delete m_s;
    }

    void LogResults()
    {
        VdbSearch :: logResults = true;
    }

    void Setup ( const string& p_query, VdbSearch :: Algorithm p_algorithm, const string& p_accession, bool p_expression = false, unsigned int p_threads = 0, bool p_blobBased = false )
    {
        delete m_s;
        VdbSearch :: Settings s;
        s . m_algorithm = p_algorithm;
        s . m_query = p_query;
        s . m_accessions . push_back ( p_accession );
        s . m_isExpression = p_expression;
        s . m_threads = p_threads;
        s . m_useBlobSearch = p_blobBased;
        m_s = new VdbSearch ( s );
    }
    void Setup ( const string& p_query, VdbSearch :: Algorithm p_algorithm, const vector < string > & p_accessions, unsigned int p_threads = 0 )
    {
        delete m_s;
        VdbSearch :: Settings s;
        s . m_algorithm = p_algorithm;
        s . m_query = p_query;
        s . m_accessions = p_accessions;
        s . m_threads = p_threads;
        m_s = new VdbSearch ( s );
    }
    void SetupWithScore ( const string& p_query, VdbSearch :: Algorithm p_algorithm, const string& p_accession, unsigned int p_minScore, bool p_blobBased = false )
    {
        delete m_s;
        VdbSearch :: Settings s;
        s . m_algorithm = p_algorithm;
        s . m_query = p_query;
        s . m_accessions . push_back ( p_accession );
        s . m_minScorePct =  p_minScore;
        s . m_useBlobSearch = p_blobBased;
        m_s = new VdbSearch ( s );
    }
    void SetupWithReference ( const string& p_query, VdbSearch :: Algorithm p_algorithm, const string& p_accession, const string& p_referenceName = string() )
    {
        delete m_s;
        VdbSearch :: Settings s;
        s . m_algorithm = p_algorithm;
        s . m_query = p_query;
        s . m_accessions . push_back ( p_accession );
        s . m_referenceDriven = true;
        if ( ! p_referenceName . empty () )
        {
            s . m_references . push_back ( p_referenceName );
        }
        m_s = new VdbSearch ( s );
    }

    const string& NextFragmentId ()
    {
        if ( ! m_s -> NextMatch ( m_accession, m_fragment ) )
        {
            throw logic_error ( "VdbSearchFixture::NextFragmentId : NextMatch() failed" );
        }
        return m_fragment;
    }

    VdbSearch* m_s;
    string m_accession;
    string m_fragment;
};
#if 0
//TODO: bad query string
//TODO: bad accession name

FIXTURE_TEST_CASE ( Create_Destroy, VdbSearchFixture )
{
    VdbSearch :: Settings s;
    s . m_query = "ACGT";
    m_s = new VdbSearch ( s );
}

FIXTURE_TEST_CASE ( SupportedAlgorithms, VdbSearchFixture )
{
    VdbSearch :: SupportedAlgorithms algs = m_s -> GetSupportedAlgorithms ();
    REQUIRE_LT ( ( size_t ) 0, algs . size () );
}

FIXTURE_TEST_CASE ( SingleAccession_FirstMatches, VdbSearchFixture )
{
    const string Accession = "SRR000001";
    Setup ( "A", VdbSearch :: FgrepDumb, Accession ); // will hit (almost) every fragment

    REQUIRE_EQ ( string ( "SRR000001.FR0.1" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR1.2" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR1.3" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( SingleAccession_FirstMatches_BlobBased_WGS, VdbSearchFixture )
{
    const string Accession = "ALWZ01";
    Setup ( "A", VdbSearch :: FgrepDumb, Accession, false, 0, true ); // will hit (almost) every fragment

    REQUIRE_EQ ( string ( "ALWZ01.FR0.1" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "ALWZ01.FR0.2" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( FgrepDumb_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    Setup ( "ATTAGC", VdbSearch :: FgrepDumb, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( FgrepBoyerMoore_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    Setup ( "ATTAGC", VdbSearch :: FgrepBoyerMoore, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( FgrepAho_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    Setup ( "ATTAGC", VdbSearch :: FgrepAho, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepDP_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{   /* VDB-2681: AgrepDP algorithm is broken */
    Setup ( "ATTAGC", VdbSearch :: AgrepDP, "SRR000001" );
    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepWuManber_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    Setup ( "ATTAGC", VdbSearch :: AgrepWuManber, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepMyers_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    Setup ( "ATTAGC", VdbSearch :: AgrepMyers, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepMyersUnltd_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    Setup ( "ATTAGC", VdbSearch :: AgrepMyersUnltd, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( NucStrstr_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    Setup ( "ATTAGC", VdbSearch :: NucStrstr, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( SmithWaterman_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    Setup ( "ATTAGC", VdbSearch :: SmithWaterman, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

#if REPORT_MULTIPLE_HITS_IN_ONE_FRAGMENT
FIXTURE_TEST_CASE ( SingleAccession_HitsInsideOneFragment, VdbSearchFixture )
{
    Setup ( "AT", VdbSearch :: FgrepDumb, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.1" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.1" ), NextFragmentId () );
    //TODO: verify positions when supported
}
#endif

FIXTURE_TEST_CASE ( MultipleAccessions, VdbSearchFixture )
{
    const string Sra1 = "SRR600096";
    const string Sra2 = "SRR000001";
    vector<string> acc;
    acc. push_back(Sra1);
    acc. push_back(Sra2);
    Setup ( "ACGTACG", VdbSearch :: NucStrstr, acc );

    REQUIRE_EQ ( Sra1 + ".FR1.5", NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR0.26",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR0.717",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR0.951",  NextFragmentId () );
}

FIXTURE_TEST_CASE ( NucStrstr_Expression, VdbSearchFixture )
{
    Setup ( "AAAAAAACCCCCCC||ATTAGC", VdbSearch :: NucStrstr, "SRR000001", true );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( Expression_OnlyForNucStrstr, VdbSearchFixture )
{
    REQUIRE_THROW ( Setup ( "AAAAAAA||ATTAGC", VdbSearch :: SmithWaterman, "SRR000001", true ) );
}

// Imperfect matches
FIXTURE_TEST_CASE ( FgrepDumb_ImperfectMatch_Unsupported, VdbSearchFixture )
{
    REQUIRE_THROW ( SetupWithScore ( "ATTAGCATTAGC", VdbSearch :: FgrepDumb, "SRR000001", 90 ) );
}
FIXTURE_TEST_CASE ( FgrepBoyerMoore_ImperfectMatch_Unsupported, VdbSearchFixture )
{
    REQUIRE_THROW ( SetupWithScore ( "ATTAGCATTAGC", VdbSearch :: FgrepBoyerMoore, "SRR000001", 90 ) );
}
FIXTURE_TEST_CASE ( FgrepAho_ImperfectMatch_Unsupported, VdbSearchFixture )
{
    REQUIRE_THROW ( SetupWithScore ( "ATTAGCATTAGC", VdbSearch :: FgrepAho, "SRR000001", 90 ) );
}

FIXTURE_TEST_CASE ( AgrepDP_ImperfectMatch, VdbSearchFixture )
{
    SetupWithScore ( "ATTAGCATTAGC", VdbSearch :: AgrepDP, "SRR000001", 90 );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3608" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepWuManber_ImperfectMatch, VdbSearchFixture )
{
    SetupWithScore ( "ATTAGCATTAGC", VdbSearch :: AgrepWuManber, "SRR000001", 90 );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3608" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepMyers_ImperfectMatch, VdbSearchFixture )
{
    SetupWithScore ( "ATTAGCATTAGC", VdbSearch :: AgrepMyers, "SRR000001", 90 );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3608" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepMyersUnltd_ImperfectMatch, VdbSearchFixture )
{
    SetupWithScore ( "ATTAGCATTAGC", VdbSearch :: AgrepMyersUnltd, "SRR000001", 90 );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3608" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( SmithWaterman_ImperfectMatch, VdbSearchFixture )
{   // SW scoring function is different from Agrep's, so the results are slightly different
// Lately the SW scoring functuion seems to have changed, and the results are very, very different
    SetupWithScore ( "ATTAGCATTAGC", VdbSearch :: SmithWaterman, "SRR000001", 90 );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.183" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
}

///////// Multi threading
FIXTURE_TEST_CASE ( Threads_RandomCrash, VdbSearchFixture )
{
    Setup ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, "SRR000001", false, 4, true ); // 4 blob-based threads on one run

    unsigned int count = 0;
    while (  m_s -> NextMatch ( m_accession, m_fragment ) )  // used to have a random crash inside VDB
    {
        ++count;
    }
    REQUIRE_EQ ( 12u, count );
}

#if TOO_SLOW_FOR_A_UNIT_TEST
FIXTURE_TEST_CASE ( MultipleAccessions_Threaded_Unsorted, VdbSearchFixture )
{
    const string Sra1 = "SRR600094";
    const string Sra2 = "SRR600095";
    Setup ( "ACGTAGGGTCC", VdbSearch :: NucStrstr, Sra2, false, 2 );
    m_s -> AddAccession ( Sra1 );

    set <string> frags;
    while ( m_s -> NextMatch ( m_accession, m_fragment ) )
    {
        frags.insert(m_fragment);
    }

    set <string> :: const_iterator it = frags . begin ();

    // this is a straight up alphanumerical sort, so the order here is not quite the same as one would expect (accession/read#/frag#)
    REQUIRE_EQ ( Sra1 + ".FR0.101990",  *it++ );
    REQUIRE_EQ ( Sra1 + ".FR0.101991",  *it++ );
    REQUIRE_EQ ( Sra1 + ".FR0.1053648", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR0.1053650", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR0.1053651", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR0.1053652", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR0.1561682", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR0.1667877", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR0.2625526", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR0.2805749", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR0.324216",  *it++ );
    REQUIRE_EQ ( Sra1 + ".FR1.101989",  *it++ );
    REQUIRE_EQ ( Sra1 + ".FR1.1053649", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR1.1053653", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR1.1561683", *it++ );
    REQUIRE_EQ ( Sra1 + ".FR1.2625553", *it++ );
    REQUIRE_EQ ( Sra2 + ".FR0.1746431", *it++ );
    REQUIRE_EQ ( Sra2 + ".FR1.1034389", *it++ );
    REQUIRE_EQ ( Sra2 + ".FR1.1746425", *it++ );
    REQUIRE_EQ ( Sra2 + ".FR1.1746434", *it++ );
    REQUIRE_EQ ( Sra2 + ".FR1.694078",  *it++ );
    REQUIRE_EQ ( Sra2 + ".FR1.69793",   *it++ );

    REQUIRE ( it == frags . end () );
}
#endif

FIXTURE_TEST_CASE ( SingleAccession_Threaded_OnBlobs, VdbSearchFixture )
{
    const string Sra1 = "SRR600094";
    Setup ( "ACGTAGGGTCC", VdbSearch :: NucStrstr, Sra1, false, 2, true );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) );
    REQUIRE_EQ ( Sra1 + ".FR1.101989",  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) );
    REQUIRE_EQ ( Sra1 + ".FR0.101990",  m_fragment );
    // now, let all the threads finish
    while ( m_s -> NextMatch ( m_accession, m_fragment ) )
    {
    }
}

//TODO: stop multi-threaded search before the end

// Reference-driven mode

class DummySearchBlockFactory : public SearchBlock :: Factory
{
public:
    virtual SearchBlock* MakeSearchBlock () const
    {
        return new FgrepSearch ( "", FgrepSearch :: FgrepDumb );
    }
};

FIXTURE_TEST_CASE ( ReferenceMatchIterator_Construct, VdbSearchFixture )
{
    DummySearchBlockFactory factory;
    ReferenceMatchIterator it ( factory, "SRR600094" );
    SearchBuffer* buf = it . NextBuffer ();
    REQUIRE_NOT_NULL ( buf );
    delete buf;
}

FIXTURE_TEST_CASE ( ReferenceMatchIterator_AccessionName, VdbSearchFixture )
{
    string accName ( "SRR600094" );
    DummySearchBlockFactory factory;
    ReferenceMatchIterator it ( factory, accName );
    SearchBuffer* buf = it . NextBuffer ();
    REQUIRE_EQ ( accName, buf -> AccessionName () );
    delete buf;
}

FIXTURE_TEST_CASE ( ReferenceMatchIterator_BufferId, VdbSearchFixture )
{
    string accName ( "SRR833251" );
    DummySearchBlockFactory factory;
    ReferenceMatchIterator it ( factory, accName );

    SearchBuffer* buf = it . NextBuffer ();
    REQUIRE_EQ ( string("gi|169794206|ref|NC_010410.1|"), buf -> BufferId () );
    delete buf;
}

#endif

FIXTURE_TEST_CASE ( ReferenceDriven_AllReferences, VdbSearchFixture )
{
    SetupWithReference ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, "SRR600094" );
/*
SRR600094.FR1.101989
SRR600094.FR0.101990
SRR600094.FR0.101991
SRR600094.FR0.324216    Not reported in reference mode b/c matches are in clipped portions of the read
SRR600094.FR0.1053648
SRR600094.FR1.1053649
SRR600094.FR0.1053650
SRR600094.FR0.1053651
SRR600094.FR0.1053652
SRR600094.FR1.1053653
SRR600094.FR0.1561682
SRR600094.FR1.1561683
SRR600094.FR0.1667877   Not reported in reference mode b/c matches are in the clipped portion of the read
SRR600094.FR0.2625526
SRR600094.FR1.2625553
SRR600094.FR0.2805749
*/

    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.101989" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.101990" ),  m_fragment );
            REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.101989" ),  m_fragment ); //TODO: remove duplicates
            REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.101990" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.101991" ),  m_fragment );
            REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.101991" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.1053649" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1053650" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1053648" ),  m_fragment );
        REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.1053649" ),  m_fragment );
        REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1053650" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1053651" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1053652" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.1053653" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1561682" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.1561683" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.2625526" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.2625553" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.2805749" ),  m_fragment );

    REQUIRE ( ! m_s -> NextMatch ( m_accession, m_fragment ) );
}

//TODO: reference mode on a non-CSRA object
//TODO: no double-reporting of fragments (may be related to how far we advance the refrence position after a match has been found and processed)
//TODO: circular references

int
main( int argc, char *argv [] )
{
    return SraSearchTestSuite(argc, argv);
}
