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
#include "ngs-vdb.hpp"

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

// SearchBlock

TEST_CASE ( FgrepDumb )
{
    FgrepSearch sb ( "CTA", VdbSearch :: FgrepDumb );
    uint64_t hitStart = 0;
    uint64_t hitEnd = 0;
    const string Bases = "ACTGACTAGTCA";
    REQUIRE ( sb.FirstMatch ( Bases.c_str(), Bases.size(), hitStart, hitEnd ) );
    REQUIRE_EQ ( (uint64_t)5, hitStart );
    REQUIRE_EQ ( (uint64_t)8, hitEnd );
}

TEST_CASE ( SearchFgrepBoyerMoore )
{
    FgrepSearch sb ( "CTA", VdbSearch :: FgrepBoyerMoore );
    uint64_t hitStart = 0;
    uint64_t hitEnd = 0;
    const string Bases = "ACTGACTAGTCA";
    REQUIRE ( sb.FirstMatch ( Bases.c_str(), Bases.size(), hitStart, hitEnd ) );
    REQUIRE_EQ ( (uint64_t)5, hitStart );
    REQUIRE_EQ ( (uint64_t)8, hitEnd );
}

TEST_CASE ( SearchFgrepAho )
{
    FgrepSearch sb ( "CTA", VdbSearch :: FgrepAho );
    uint64_t hitStart = 0;
    uint64_t hitEnd = 0;
    const string Bases = "ACTGACTAGTCA";
    REQUIRE ( sb.FirstMatch ( Bases.c_str(), Bases.size(), hitStart, hitEnd ) );
    REQUIRE_EQ ( (uint64_t)5, hitStart );
    REQUIRE_EQ ( (uint64_t)8, hitEnd );
}

TEST_CASE ( SearchAgrepDP )
{
    AgrepSearch sb ( "CTA", VdbSearch :: AgrepDP, 100 );
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
        m_s = 0;
        
        VdbSearch :: useBlobSearch = p_blobBased;
        
        m_s = new VdbSearch ( p_algorithm, p_query, p_expression, 100, p_threads );
        m_s -> AddAccession ( p_accession );
    }
    void SetupWithScore ( const string& p_query, VdbSearch :: Algorithm p_algorithm, const string& p_accession, unsigned int p_minScore, bool p_blobBased = false )
    {
        delete m_s;
        m_s = 0;
        
        VdbSearch :: useBlobSearch = p_blobBased;
        
        m_s = new VdbSearch ( p_algorithm, p_query, false, p_minScore  );
        m_s -> AddAccession ( p_accession );
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

//TODO: bad query string
//TODO: bad accession name

FIXTURE_TEST_CASE ( Create_Destroy, VdbSearchFixture )
{
    m_s = new VdbSearch ( VdbSearch :: FgrepDumb, "ACGT", false );
    REQUIRE_EQ ( VdbSearch :: FgrepDumb, m_s -> GetAlgorithm () ); 
    delete m_s;
    m_s = 0;
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
#if 0
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
    Setup ( "ACGTACG", VdbSearch :: NucStrstr, Sra1 );
    m_s -> AddAccession ( Sra2 );
    
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

/* TODO
FIXTURE_TEST_CASE ( MultipleAccessions_Threaded_Sorted, VdbSearchFixture )
{
    
    const string Sra1 = "SRR600094";
    const string Sra2 = "SRR600095";
    Setup ( "ACGTAGGGTCC", VdbSearch :: NucStrstr, Sra2, false, 2 );
    m_s -> AddAccession ( Sra1 );
    // sorted in order accession, read#, frag#
    // accessions sort in the order they were added (in this case sra2, sra1)
    REQUIRE_EQ ( Sra2 + ".FR1.69793",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR1.694078",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR1.1034389",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR1.1746425",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR0.1746431",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR1.1746434",   NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR1.101989", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.101990", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.101991", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.324216", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.1053648", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR1.1053649", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.1053650", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.1053651", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.1053652", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR1.1053653", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.1561682", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR1.1561683", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.1667877", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.2625526", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR1.2625553", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.2805749", NextFragmentId () );
    

    REQUIRE ( ! m_s -> NextMatch ( m_accession, m_fragment ) );
}
*/

//TODO: stop multi-threaded search before the end 

/////////////////// Search directly in Blobs ( this should go to ncbi-vdb/test/search )

class VdbBlobFixture
{
public:
    VdbBlobFixture()
    : mgr(0), curs(0), col_idx(~0)
    {
        if ( VDBManagerMakeRead(&mgr, NULL) != 0 )
            throw logic_error ( "VdbBlobFixture: VDBManagerMakeRead failed" );
    }
    
    ~VdbBlobFixture()
    {
        if ( mgr && VDBManagerRelease ( mgr ) != 0 )
            throw logic_error ( "~VdbBlobFixture: VDBManagerRelease failed" );
        if ( curs && VCursorRelease ( curs ) != 0 )
            throw logic_error ( "~VdbBlobFixture: VCursorRelease failed" );
    }
    
    rc_t Setup( const char * acc )
    {
        const VDatabase *db = NULL;
        rc_t rc = VDBManagerOpenDBRead ( mgr, &db, NULL, acc );
        rc_t rc2;

        const VTable *tbl = NULL;
        if (rc == 0) 
        {
            rc = VDatabaseOpenTableRead ( db, &tbl, "SEQUENCE" );
        }
        else 
        {
            rc = VDBManagerOpenTableRead ( mgr, &tbl, NULL, acc );
            if (rc != 0) 
            {
                VSchema *schema = NULL;
                rc = VDBManagerMakeSRASchema(mgr, &schema);
                if ( rc != 0 )
                {
                    return rc;
                }
                    
                rc = VDBManagerOpenTableRead ( mgr, &tbl, schema, acc );
                
                rc2 = VSchemaRelease ( schema );
                if ( rc == 0 )
                    rc = rc2;
            }
        }

        if ( rc == 0 )
        {
            rc = VTableCreateCursorRead(tbl, &curs);
            if ( rc == 0 )
            {
                col_idx = ~0;
                rc = VCursorAddColumn ( curs, &col_idx, "READ" );
                if ( rc == 0 )
                {
                    rc = VCursorOpen(curs);
                }
            }
        }
        
        rc2 = VTableRelease ( tbl );
        if ( rc == 0 )
            rc = rc2;
        
        if ( db != NULL )
        {
            rc2 = VDatabaseRelease ( db );
            if ( rc == 0 )
               rc = rc2;
        }
            
        return rc;
    }
    
    const VDBManager * mgr;
    const VCursor * curs;
    uint32_t col_idx;
};


FIXTURE_TEST_CASE ( BlobSearch_FgrepDumb, VdbBlobFixture )
{
    const string Accession = "SRR000001";
    const char* query[] = { "TCAGA" };
    
    REQUIRE_RC ( Setup ( Accession.c_str() ) );
    REQUIRE_RC ( VCursorOpen (curs ) );
    
    REQUIRE_RC ( VCursorSetRowId (curs, 1 ) );
    REQUIRE_RC ( VCursorOpenRow (curs ) );
    
    struct VBlob *blob;
    REQUIRE_RC ( VCursorGetBlob ( curs, (const VBlob**)&blob, col_idx ) );
    
    int64_t first;
    uint64_t count;
    REQUIRE_RC ( VBlobIdRange ( blob, &first, &count ) );
    REQUIRE_EQ ( (int64_t)1, first );    
    REQUIRE_EQ ( (uint64_t)4, count );

    struct Fgrep* fgrep;
    REQUIRE_RC ( FgrepMake ( & fgrep, FGREP_MODE_ASCII | FGREP_ALG_DUMB, &query[0], 1 ) );
    
    size_t blobLength = BlobBufferBytes ( blob );
    FgrepMatch matchinfo;
    
    int32_t curStart = 0; 
    REQUIRE_EQ ( (uint32_t)1, FgrepFindFirst ( fgrep, ((const char*)blob->data.base) + curStart, blobLength - curStart, & matchinfo ) );

    REQUIRE_EQ ( (int32_t)0,                matchinfo.position );
    REQUIRE_EQ ( (int32_t)strlen(query[0]), matchinfo.length );

    curStart += matchinfo.position + matchinfo.length; 
    REQUIRE_EQ ( (uint32_t)1, FgrepFindFirst ( fgrep, ((const char*)blob->data.base) + curStart, blobLength - curStart, & matchinfo ) );
    
    REQUIRE_EQ ( (int32_t)792,              matchinfo.position );
    REQUIRE_EQ ( string ( query[0] ), string ( ((const char*)blob->data.base) + curStart + matchinfo.position, matchinfo.length ) );
    
    curStart += matchinfo.position + matchinfo.length;
    REQUIRE_EQ ( (uint32_t)0, FgrepFindFirst ( fgrep, ((const char*)blob->data.base) + curStart, blobLength - curStart, & matchinfo ) );
    
    //TODO: convert positions into rowIds
    
    FgrepFree ( fgrep );

    REQUIRE_RC ( VCursorCloseRow (curs ) );
    REQUIRE_RC ( VBlobRelease ( blob ) );
}

//////////////////////////////////////////// ngs-vdb ( migrate these tests to wherever ngs-vdb code moves, eventually )

#include <ngs-vdb.hpp>

using namespace ngs::vdb;

TEST_CASE ( VdbReadCollection_WGS_Construct )
{
    VdbReadCollection coll = ncbi :: NGS_VDB :: openVdbReadCollection ( "ALWZ01" );
    FragmentBlobIterator fragIt = coll.getFragmentBlobs();    
    REQUIRE ( fragIt . nextBlob () );
    REQUIRE_EQ ( (size_t)948, fragIt . Size () );
    const string BeginsWith = "GCCTCTCTCTC"; 
    REQUIRE_EQ ( BeginsWith, string ( fragIt . Data (), BeginsWith.size() ) );
}

TEST_CASE ( VdbReadCollection_WGS_RowInfo )
{
    VdbReadCollection coll = ncbi :: NGS_VDB :: openVdbReadCollection ( "ALWZ01" );
    FragmentBlobIterator fragIt = coll.getFragmentBlobs();    
    REQUIRE ( fragIt . nextBlob () );
    
    uint64_t offset = 0;
    int64_t rowId = 0;
    uint64_t nextRowStart = 0;
    fragIt . GetRowInfo ( offset, rowId, nextRowStart );
    REQUIRE_EQ ( (int64_t)1, rowId );
    REQUIRE_EQ ( (uint64_t)227, nextRowStart );
    
    offset = nextRowStart - 1;
    fragIt . GetRowInfo ( offset, rowId, nextRowStart );
    REQUIRE_EQ ( (int64_t)1, rowId );
    REQUIRE_EQ ( (uint64_t)227, nextRowStart );
    
    offset = nextRowStart;
    fragIt . GetRowInfo ( offset, rowId, nextRowStart );
    REQUIRE_EQ ( (int64_t)2, rowId );
    REQUIRE_EQ ( (uint64_t)227+245, nextRowStart );
}    

//TODO: accessing blob without a call to nextBlob()
#endif
int
main( int argc, char *argv [] )
{
    return SraSearchTestSuite(argc, argv);
}
