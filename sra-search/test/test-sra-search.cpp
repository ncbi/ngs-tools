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

//////////////////////////////////////////// ngs-vdb ( migrate these tests to wherever ngs-vdb code moves, eventually )

#include <ngs-vdb.hpp>

using namespace ngs::vdb;

static const string AccessionWGS            = "ALWZ01";
static const string AccessionSRA_table      = "SRR000123";
static const string AccessionSRA_database   = "SRR600096";
static const string AccessionCSRA1          = "SRR1063272";

class BlobIteratorHelper 
{
public:    
    BlobIteratorHelper ( const string& p_accession )
    :   coll ( ncbi :: NGS_VDB :: openVdbReadCollection ( p_accession ) ),
        fragIt ( coll.getFragmentBlobs() )
    {
        if ( ! fragIt . nextBlob () )
        {
            throw logic_error ( "BlobIteratorHelper : nextBlob() failed" );
        }
    }
    
    bool VerifyFragInfo ( uint64_t p_offset, const string& p_fragId,  uint64_t p_nextFragmentStart, bool p_biological = true )
    {
        string fragId;
        uint64_t nextFragmentStart = 0;
        bool biological;
        fragIt . GetFragmentInfo ( p_offset, fragId, nextFragmentStart, biological );
        if ( p_fragId != fragId )
        {
            cout << "fragId expected = " << p_fragId << ", actual = " << fragId << endl;
            return false;
        }
        if ( p_nextFragmentStart != nextFragmentStart )
        {
            cout << "nextFragmentStart expected = " << p_nextFragmentStart << ", actual = " << nextFragmentStart << endl;
            return false;
        }
        if ( p_biological != biological )
        {
            cout << "biological expected = " << p_biological << ", actual = " << biological << endl;
            return false;
        }
        return true;
    }

    VdbReadCollection coll;
    FragmentBlobIterator fragIt;    
};

TEST_CASE ( VdbReadCollection_WGS_Construct )
{
    BlobIteratorHelper bit ( AccessionWGS );
    REQUIRE_EQ ( (size_t)948, bit . fragIt . Size () );
    const string BeginsWith = "GCCTCTCTCTC"; 
    REQUIRE_EQ ( BeginsWith, string ( bit . fragIt . Data (), BeginsWith . size() ) );
}

TEST_CASE ( VdbReadCollection_WGS_RowInfo )
{
    BlobIteratorHelper bit ( AccessionWGS );
    REQUIRE ( bit . VerifyFragInfo ( 0,   AccessionWGS + ".FR0.1", 227 ) );
    REQUIRE ( bit . VerifyFragInfo ( 226, AccessionWGS + ".FR0.1", 227 ) );
    REQUIRE ( bit . VerifyFragInfo ( 227, AccessionWGS + ".FR0.2", 227+245 ) );
}    

TEST_CASE ( VdbReadCollection_SRAtable_Construct )
{
    BlobIteratorHelper bit ( AccessionSRA_table );
    REQUIRE_EQ ( (size_t)858, bit . fragIt . Size () );
    const string BeginsWith = "TCAGTTTCG"; 
    REQUIRE_EQ ( BeginsWith, string ( bit . fragIt . Data (), BeginsWith.size() ) );
}

TEST_CASE ( VdbReadCollection_SRAtable_RowInfo )
{
    BlobIteratorHelper bit ( AccessionSRA_table );
    REQUIRE ( bit . VerifyFragInfo ( 0, "", 4, false ) );
    REQUIRE ( bit . VerifyFragInfo ( 4,       AccessionSRA_table + ".FR0.1", 157, true ) );
    REQUIRE ( bit . VerifyFragInfo ( 156,     AccessionSRA_table + ".FR0.1", 157, true ) );
    REQUIRE ( bit . VerifyFragInfo ( 157, "", 157+4, false ) );
    REQUIRE ( bit . VerifyFragInfo ( 157+4,   AccessionSRA_table + ".FR0.2", 157+4+276, true ) );
}

TEST_CASE ( VdbReadCollection_SRAdatabase_Construct )
{
    BlobIteratorHelper bit ( AccessionSRA_database );
    REQUIRE_EQ ( (size_t)5600, bit . fragIt . Size () );
    const string BeginsWith = "TACGGAGGGGGCTA"; 
    REQUIRE_EQ ( BeginsWith, string ( bit . fragIt . Data (), BeginsWith.size() ) );
}

TEST_CASE ( VdbReadCollection_SRAdatabase_RowInfo )
{
    BlobIteratorHelper bit ( AccessionSRA_database );
    REQUIRE ( bit . VerifyFragInfo ( 0,   AccessionSRA_database + ".FR0.1", 175 ) );
    REQUIRE ( bit . VerifyFragInfo ( 174, AccessionSRA_database + ".FR0.1", 175 ) );
    REQUIRE ( bit . VerifyFragInfo ( 175, AccessionSRA_database + ".FR1.1", 350 ) ); 
}

TEST_CASE ( VdbReadCollection_CSRA1_Construct )
{
    BlobIteratorHelper bit ( AccessionCSRA1 );
    REQUIRE_EQ ( (size_t)808, bit . fragIt . Size () );
    const string BeginsWith = "ACTCGACATTCTGCC"; 
    REQUIRE_EQ ( BeginsWith, string ( bit . fragIt . Data (), BeginsWith.size() ) );
}

TEST_CASE ( VdbReadCollection_CSRA1_RowInfo )
{
    BlobIteratorHelper bit ( AccessionCSRA1 );
    REQUIRE ( bit . VerifyFragInfo ( 0,   AccessionCSRA1 + ".FR0.1", 101 ) );
    REQUIRE ( bit . VerifyFragInfo ( 100, AccessionCSRA1 + ".FR0.1", 101 ) );
    REQUIRE ( bit . VerifyFragInfo ( 101, AccessionCSRA1 + ".FR1.1", 202 ) ); 
    REQUIRE ( bit . VerifyFragInfo ( 202, AccessionCSRA1 + ".FR0.2", 202+101 ) ); 
}

//TODO: accessing blob without a call to nextBlob()

int
main( int argc, char *argv [] )
{
    return SraSearchTestSuite(argc, argv);
}
