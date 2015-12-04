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

#include "vdb-search.cpp"

#include <set>

#include <ktst/unit_test.hpp> 

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
    
    void Setup ( const string& p_query, VdbSearch :: Algorithm p_algorithm, const string& p_accession, bool p_expression = false, unsigned int p_threads = 0 )
    {
        delete m_s;
        m_s = 0;
        
        m_s = new VdbSearch ( p_algorithm, p_query, p_expression, 100, p_threads );
        m_s -> AddAccession ( p_accession );
    }
    void SetupWithScore ( const string& p_query, VdbSearch :: Algorithm p_algorithm, const string& p_accession, unsigned int p_minScore )
    {
        delete m_s;
        m_s = 0;
        
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
    
    REQUIRE ( m_s -> NextMatch ( m_accession, m_fragment ) );
    REQUIRE_EQ ( Accession, m_accession );
    REQUIRE_EQ ( string ( "SRR000001.FR0.1" ), m_fragment );
    
    // multiple fragments per read
    REQUIRE_EQ ( Accession, m_accession );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2" ), NextFragmentId () );
    REQUIRE_EQ ( Accession, m_accession );
    REQUIRE_EQ ( string ( "SRR000001.FR1.2" ), NextFragmentId () );
    
    REQUIRE_EQ ( Accession, m_accession );
    REQUIRE_EQ ( string ( "SRR000001.FR0.3" ), NextFragmentId () );
    REQUIRE_EQ ( Accession, m_accession );
    REQUIRE_EQ ( string ( "SRR000001.FR1.3" ), NextFragmentId () );
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
    SetupWithScore ( "ATTAGCATTAGC", VdbSearch :: SmithWaterman, "SRR000001", 90 );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.183" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( MultipleAccessions_Threaded_Unsorted, VdbSearchFixture )
{
    
    const string Sra1 = "SRR600094";
    const string Sra2 = "SRR600095";
    Setup ( "ACGTAGGGTCC", VdbSearch :: NucStrstr, Sra2, false, 2 );
    m_s -> AddAccession ( Sra1 );
    
    set <string> frags;
    while (  m_s -> NextMatch ( m_accession, m_fragment ) )
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

//TODO: bad query string
///TODO: bad accession name

int
main( int argc, char *argv [] )
{
    return SraSearchTestSuite(argc, argv);
}
