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

    void Setup ( const string& p_query, VdbSearch :: Algorithm p_algorithm, const string& p_accession = string() )
    {
        m_settings . m_algorithm = p_algorithm;
        m_settings . m_query = p_query;
        if ( ! p_accession . empty () )
        {
            m_settings . m_accessions . push_back ( p_accession );
        }
        delete m_s;
        m_s = new VdbSearch ( m_settings );
    }

    void SetupSingleThread ( const string& p_query, VdbSearch :: Algorithm p_algorithm, const string& p_accession = string() )
    {
        m_settings . m_threads = 0;
        Setup ( p_query, p_algorithm, p_accession );
    }

    void SetupMultiThread ( const string& p_query, VdbSearch :: Algorithm p_algorithm, unsigned int p_threads, bool p_blobs, const string& p_accession = string() )
    {
        m_settings . m_threads = p_threads;
        m_settings . m_useBlobSearch = p_blobs;
        Setup ( p_query, p_algorithm, p_accession );
    }

    const string& NextFragmentId ()
    {
        if ( ! m_s -> NextMatch ( m_accession, m_fragment ) )
        {
            throw logic_error ( "VdbSearchFixture::NextFragmentId : NextMatch() failed" );
        }
        return m_fragment;
    }

    VdbSearch :: Settings m_settings;
    VdbSearch* m_s;
    string m_accession;
    string m_fragment;
};

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
    SetupSingleThread ( "A", VdbSearch :: FgrepDumb, Accession ); // will hit (almost) every fragment

    REQUIRE_EQ ( string ( "SRR000001.FR0.1" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR1.2" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR1.3" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( SingleAccession_MaxMatches, VdbSearchFixture )
{
    const string Accession = "SRR000001";
    m_settings . m_maxMatches = 3;
    SetupSingleThread ( "A", VdbSearch :: FgrepDumb, Accession ); // will hit (almost) every fragment, only return the first 3

    REQUIRE_EQ ( string ( "SRR000001.FR0.1" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR1.2" ), NextFragmentId () );
    REQUIRE ( ! m_s -> NextMatch ( m_accession, m_fragment ) );
}

FIXTURE_TEST_CASE ( FgrepDumb_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    SetupSingleThread ( "ATTAGC", VdbSearch :: FgrepDumb, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( FgrepBoyerMoore_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    SetupSingleThread ( "ATTAGC", VdbSearch :: FgrepBoyerMoore, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( FgrepAho_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    SetupSingleThread ( "ATTAGC", VdbSearch :: FgrepAho, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepDP_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{   /* VDB-2681: AgrepDP algorithm is broken */
    SetupSingleThread ( "ATTAGC", VdbSearch :: AgrepDP, "SRR000001" );
    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepWuManber_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    SetupSingleThread ( "ATTAGC", VdbSearch :: AgrepWuManber, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepMyers_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    SetupSingleThread ( "ATTAGC", VdbSearch :: AgrepMyers, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepMyersUnltd_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    SetupSingleThread ( "ATTAGC", VdbSearch :: AgrepMyersUnltd, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( NucStrstr_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    SetupSingleThread ( "ATTAGC", VdbSearch :: NucStrstr, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( SmithWaterman_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{
    SetupSingleThread ( "ATTAGC", VdbSearch :: SmithWaterman, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( MultipleAccessions, VdbSearchFixture )
{
    const string Sra1 = "SRR600096";
    const string Sra2 = "SRR000001";
    m_settings . m_accessions . push_back(Sra1);
    m_settings . m_accessions . push_back(Sra2);
    SetupSingleThread ( "ACGTACG", VdbSearch :: NucStrstr );

    REQUIRE_EQ ( Sra1 + ".FR1.5", NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR0.26",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR0.717",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR0.951",  NextFragmentId () );
}

FIXTURE_TEST_CASE ( NucStrstr_Expression, VdbSearchFixture )
{
    m_settings . m_isExpression = true;
    SetupSingleThread ( "AAAAAAACCCCCCC||ATTAGC", VdbSearch :: NucStrstr, "SRR000001" );

    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( Expression_OnlyForNucStrstr, VdbSearchFixture )
{
    m_settings . m_isExpression = true;
    REQUIRE_THROW ( SetupSingleThread ( "AAAAAAA||ATTAGC", VdbSearch :: SmithWaterman, "SRR000001" ) );
}
// Imperfect matches
FIXTURE_TEST_CASE ( FgrepDumb_ImperfectMatch_Unsupported, VdbSearchFixture )
{
    m_settings . m_minScorePct =  90;
    REQUIRE_THROW ( SetupSingleThread ( "ATTAGCATTAGC", VdbSearch :: FgrepDumb, "SRR000001" ) );
}
FIXTURE_TEST_CASE ( FgrepBoyerMoore_ImperfectMatch_Unsupported, VdbSearchFixture )
{
    m_settings . m_minScorePct =  90;
    REQUIRE_THROW ( SetupSingleThread ( "ATTAGCATTAGC", VdbSearch :: FgrepBoyerMoore, "SRR000001" ) );
}
FIXTURE_TEST_CASE ( FgrepAho_ImperfectMatch_Unsupported, VdbSearchFixture )
{
    m_settings . m_minScorePct =  90;
    REQUIRE_THROW ( SetupSingleThread ( "ATTAGCATTAGC", VdbSearch :: FgrepAho, "SRR000001" ) );
}

FIXTURE_TEST_CASE ( AgrepDP_ImperfectMatch, VdbSearchFixture )
{
    m_settings . m_minScorePct =  90;
    SetupSingleThread ( "ATTAGCATTAGC", VdbSearch :: AgrepDP, "SRR000001" );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3608" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepWuManber_ImperfectMatch, VdbSearchFixture )
{
    m_settings . m_minScorePct =  90;
    SetupSingleThread ( "ATTAGCATTAGC", VdbSearch :: AgrepWuManber, "SRR000001" );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3608" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepMyers_ImperfectMatch, VdbSearchFixture )
{
    m_settings . m_minScorePct =  90;
    SetupSingleThread ( "ATTAGCATTAGC", VdbSearch :: AgrepMyers, "SRR000001" );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3608" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( AgrepMyersUnltd_ImperfectMatch, VdbSearchFixture )
{
    m_settings . m_minScorePct =  90;
    SetupSingleThread ( "ATTAGCATTAGC", VdbSearch :: AgrepMyersUnltd, "SRR000001" );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3608" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( SmithWaterman_ImperfectMatch, VdbSearchFixture )
{   // SW scoring function is different from Agrep's, so the results are slightly different
    m_settings . m_minScorePct =  90;
    SetupSingleThread ( "ATTAGCATTAGC", VdbSearch :: SmithWaterman, "SRR000001" );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.183" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
}

///////// Multi threading

FIXTURE_TEST_CASE ( SingleAccession_FirstMatches_BlobBased_WGS, VdbSearchFixture )
{
    const string Accession = "ALWZ01";
    SetupMultiThread ( "A", VdbSearch :: FgrepDumb, 2, true, Accession ); // will hit (almost) every fragment
    // observe some results coming back and stop by destroying the search object
    for ( unsigned int i = 0 ; i < 20; ++i )
    {
        REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) );
    }
}

FIXTURE_TEST_CASE ( Threads_RandomCrash, VdbSearchFixture )
{
    SetupMultiThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, 4, true, "SRR000001" ); // 4 blob-based threads on one run

    unsigned int count = 0;
    while (  m_s -> NextMatch ( m_accession, m_fragment ) )  // used to have a random crash inside VDB
    {
        ++count;
    }
    REQUIRE_EQ ( 12u, count );
}

FIXTURE_TEST_CASE ( SingleAccession_Threaded_OnBlobs, VdbSearchFixture )
{
    const string Sra1 = "SRR600094";
    SetupMultiThread ( "ACGTAGGGTCC", VdbSearch :: NucStrstr, 2, false, Sra1 );

    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) );
    CHECK_EQ ( Sra1 + ".FR1.101989",  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) );
    CHECK_EQ ( Sra1 + ".FR0.101990",  m_fragment );
}

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
    REQUIRE_NOT_NULL ( buf );
    REQUIRE_EQ ( accName, buf -> AccessionName () );
    delete buf;
}

FIXTURE_TEST_CASE ( ReferenceMatchIterator_BufferId, VdbSearchFixture )
{
    string accName ( "SRR833251" );
    DummySearchBlockFactory factory;
    ReferenceMatchIterator it ( factory, accName );

    SearchBuffer* buf = it . NextBuffer ();
    REQUIRE_NOT_NULL ( buf );
    REQUIRE_EQ ( string("gi|169794206|ref|NC_010410.1|"), buf -> BufferId () );
    delete buf;
}

FIXTURE_TEST_CASE ( ReferenceDriven_ReferenceNotFound, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    m_settings . m_references . push_back ( ReferenceSpec ( "NOT_ME_GUV" ) );
    SetupSingleThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, "SRR600094" );
    REQUIRE ( ! m_s -> NextMatch ( m_accession, m_fragment ) );
}

FIXTURE_TEST_CASE ( ReferenceDriven_NotCSRA, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    SetupSingleThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, "SRR000001" );
    // No references or alignments in the archive, this becomes a scan of SEQUENCE table
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR000001.FR0.28322" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR000001.FR0.65088" ),  m_fragment );
    // etc
}

FIXTURE_TEST_CASE ( ReferenceDriven_SingleReference_SingleAccession, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    m_settings . m_useBlobSearch  = false;
    m_settings . m_references . push_back ( ReferenceSpec ( "NC_000007.13" ) );
    SetupSingleThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, "SRR600094" );

    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.1053649" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1053650" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1053648" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1053651" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1053652" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.1053653" ),  m_fragment );
    REQUIRE ( ! m_s -> NextMatch ( m_accession, m_fragment ) );
}

FIXTURE_TEST_CASE ( ReferenceDriven_AllReferences_NoDuplicates, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    m_settings . m_useBlobSearch  = false;
    SetupSingleThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, "SRR600094" );

/*
SRR600094.FR1.101989
SRR600094.FR0.101990
SRR600094.FR0.101991
SRR600094.FR0.324216    Not reported in reference mode b/c matches are in clipped portions of the read
SRR600094.FR1.1053649
etc
*/
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.101989" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.101990" ),  m_fragment );
            // REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.101989" ),  m_fragment ); // used to be duplicates
            // REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.101990" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.101991" ),  m_fragment );
}

FIXTURE_TEST_CASE ( ReferenceDriven_AllReferences_NoDuplicates_Blobs, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    m_settings . m_useBlobSearch  = true;
    SetupSingleThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, "SRR600094" );

    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.101989" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.101990" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.101991" ),  m_fragment );
}

FIXTURE_TEST_CASE ( ReferenceDriven_NoBlobs_SingleSlice_SingleAccession, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    m_settings . m_useBlobSearch  = false;
    m_settings . m_references . push_back ( ReferenceSpec ( "NC_000007.13", 81000000, 105000000 ) );
    SetupSingleThread ( "ACGTAGGGTC", VdbSearch :: FgrepDumb, "SRR600094" );

    // Match on NC_000007.13 at 104,782,835-104,782,845
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1125868" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1125870" ),  m_fragment );
    // Match on NC_000007.13 at 81,579,623-81,579,633 (reverse)
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.1094914" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.1094915" ),  m_fragment );

    REQUIRE ( ! m_s -> NextMatch ( m_accession, m_fragment ) );
}

FIXTURE_TEST_CASE ( ReferenceDriven_Blobs_SingleSlice_SingleAccession, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    m_settings . m_useBlobSearch  = true;
    m_settings . m_references . push_back ( ReferenceSpec ( "NC_000007.13", 81575001, 105000000 ) );
    SetupSingleThread ( "ACGTAGGGTC", VdbSearch :: FgrepDumb, "SRR600094" );

    // Match on NC_000007.13 at 81,579,623-81,579,633 (reverse)
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.1094914" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR1.1094915" ),  m_fragment );
    // Match on NC_000007.13 at 104,782,835-104,782,845
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1125868" ),  m_fragment );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) ); REQUIRE_EQ ( string ( "SRR600094.FR0.1125870" ),  m_fragment );

    REQUIRE ( ! m_s -> NextMatch ( m_accession, m_fragment ) );
}

//TODO: reference-driven, specify a single reference slice, match against different accessions
//TODO: reference-driven, specify multiple reference slices, single accession
//TODO: reference-driven, specify multiple reference slices, match against different accessions


// An NGS-based way to find alignments on a blob boundary:
// FIXTURE_TEST_CASE(temporary, CSRA1_Fixture)
// {
//     // we should filter out secondary alignments with SEQ_SPOT_ID == 0
//     ngs :: ReadCollection run = ncbi :: NGS :: openReadCollection ( "SRR600094" );
//     ngs :: ReferenceIterator refIt = run . getReferences ();
//     while (refIt . nextReference ())
//     {
//         ngs :: AlignmentIterator alIt = refIt . getAlignments(ngs::Alignment::all);
//         while (alIt.nextAlignment())
//         {
//             int64_t begin = alIt.getAlignmentPosition();
//             int64_t end = begin + alIt.getAlignmentLength();
//             if (begin / 5000 != end / 5000)
//             {
//                 cout << alIt.getAlignmentId() << ": " << refIt . getCanonicalName() << " " << begin << "-" << end << "'" << alIt . getReferenceBases() . toString () << "'" << endl;
//                 exit(0);
//             }
//         }
//     }
// }

FIXTURE_TEST_CASE ( ReferenceDriven_MatchAcrossBlobBoundary, VdbSearchFixture )
{
    const string Query = "TTGAAGAGATCCGACATCA";
    m_settings . m_referenceDriven = true;
    m_settings . m_useBlobSearch  = true;
    m_settings . m_accessions . push_back("SRR600094");
    m_settings . m_references . push_back ( ReferenceSpec ( "NC_000001.10", 14936, 15011 ) ); // crosses the blob boundary (5000)

    SetupSingleThread ( Query, VdbSearch :: FgrepDumb );
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) );
    REQUIRE_EQ ( string("SRR600094.FR0.1647758"), m_fragment ); // aligns into the above region, across the blob boundary
}

//TODO: reference-driven, across the "end" of a circular reference

//TODO: reference-driven, specify multiple reference names, single accession
//TODO: reference-driven, specify multiple reference names, match against different accessions

int
main( int argc, char *argv [] )
{
    return SraSearchTestSuite(argc, argv);
}
