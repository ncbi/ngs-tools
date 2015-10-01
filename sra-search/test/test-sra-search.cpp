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

#include <ktst/unit_test.hpp> 

using namespace std;
using namespace ncbi::NK;
using namespace ngs;

TEST_SUITE(SraSearchTestSuite);

class VdbSearchFixture
{
public:
    VdbSearchFixture ()
    {
        VdbSearch :: logResults = false;
    }
    
    void LogResults() 
    { 
        VdbSearch :: logResults = true; 
    }
    
    void Setup ( const string& p_query, VdbSearch :: Algorithm p_algorithm, const string& p_accession)
    {
        m_s . SetQuery ( p_query );
        m_s . SetAlgorithm ( p_algorithm );
        m_s . AddAccession ( p_accession );
    }
    
    const string& NextFragmentId () 
    {
        if ( ! m_s . NextMatch ( m_accession, m_fragment ) )
        {
            throw logic_error ( "VdbSearchFixture::NextFragmentId : NextMatch() failed" );
        }
        return m_fragment;
    }

    VdbSearch m_s;
    string m_accession;
    string m_fragment;
};

FIXTURE_TEST_CASE ( Create_Destroy, VdbSearchFixture )
{
    m_s . SetQuery ( "ACGT" );
}

FIXTURE_TEST_CASE ( Algorithm_Default, VdbSearchFixture )
{
    m_s . SetQuery ( "ACGT" );
    REQUIRE_EQ ( VdbSearch :: Default, m_s . GetAlgorithm () ); 
}

FIXTURE_TEST_CASE ( SetAlgorithm_Bad, VdbSearchFixture )
{
    m_s . SetQuery ( "ACGT" );
    REQUIRE ( ! m_s . SetAlgorithm ( "bad" ) );
}

// algorithm selection based on a string (for command line)
FIXTURE_TEST_CASE ( SetAlgorithmString_FgrepDumb, VdbSearchFixture )
{
    REQUIRE ( m_s . SetAlgorithm ( "FgrepDumb" ) );
    REQUIRE_EQ ( VdbSearch :: FgrepDumb, m_s . GetAlgorithm () ); 
}
FIXTURE_TEST_CASE ( SetAlgorithmString_FgrepBoyerMoore, VdbSearchFixture )
{
    REQUIRE ( m_s . SetAlgorithm ( "FgrepBoyerMoore" ) );
    REQUIRE_EQ ( VdbSearch :: FgrepBoyerMoore, m_s . GetAlgorithm () ); 
}
FIXTURE_TEST_CASE ( SetAlgorithmString_FgrepAho, VdbSearchFixture )
{
    REQUIRE ( m_s . SetAlgorithm ( "FgrepAho" ) );
    REQUIRE_EQ ( VdbSearch :: FgrepAho, m_s . GetAlgorithm () ); 
}
FIXTURE_TEST_CASE ( SetAlgorithmString_AgrepDP_Disabled, VdbSearchFixture )
{
    REQUIRE ( ! m_s . SetAlgorithm ( "AgrepDP" ) ); // using AgrepDB from command line is disabled
}
FIXTURE_TEST_CASE ( SetAlgorithmString_AgrepWuManber, VdbSearchFixture )
{
    REQUIRE ( m_s . SetAlgorithm ( "AgrepWuManber" ) );
    REQUIRE_EQ ( VdbSearch :: AgrepWuManber, m_s . GetAlgorithm () ); 
}
FIXTURE_TEST_CASE ( SetAlgorithmString_AgrepMyers, VdbSearchFixture )
{
    REQUIRE ( m_s . SetAlgorithm ( "AgrepMyers" ) );
    REQUIRE_EQ ( VdbSearch :: AgrepMyers, m_s . GetAlgorithm () ); 
}
FIXTURE_TEST_CASE ( SetAlgorithmString_AgrepMyersUnltd, VdbSearchFixture )
{
    REQUIRE ( m_s . SetAlgorithm ( "AgrepMyersUnltd" ) );
    REQUIRE_EQ ( VdbSearch :: AgrepMyersUnltd, m_s . GetAlgorithm () ); 
}
FIXTURE_TEST_CASE ( SetAlgorithmString_NucStrstr, VdbSearchFixture )
{
    REQUIRE ( m_s . SetAlgorithm ( "NucStrstr" ) );
    REQUIRE_EQ ( VdbSearch :: NucStrstr, m_s . GetAlgorithm () ); 
}
//

FIXTURE_TEST_CASE ( SingleAccession_FirstMatches, VdbSearchFixture )
{
    const string Accession = "SRR000001";
    Setup ( "A", VdbSearch :: FgrepDumb, Accession ); // will hit (almost) every fragment
    
    REQUIRE ( m_s . NextMatch ( m_accession, m_fragment ) );
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

#if SHOW_UNIMPLEMENTED
FIXTURE_TEST_CASE ( AgrepDP_SingleAccession_HitsAcrossFragments, VdbSearchFixture )
{   /* VDB-2681: AgrepDP algorithm is broken */
    Setup ( "ATTAGC", VdbSearch :: AgrepDP, "SRR000001" );
LogResults();    
    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
}
#endif

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
    m_s . SetQuery ( "ACGTACG" );
    m_s . SetAlgorithm ( VdbSearch :: NucStrstr );
    m_s . AddAccession ( Sra1 );
    m_s . AddAccession ( Sra2 );
    
    REQUIRE_EQ ( Sra1 + ".FR1.5", NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR0.26",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR0.717",   NextFragmentId () );
    REQUIRE_EQ ( Sra2 + ".FR0.951",  NextFragmentId () );
}

//TODO: bad query string
//TODO: bad m_accession name
//TODO: multiple accessions

int
main( int argc, char *argv [] )
{
    return SraSearchTestSuite(argc, argv);
}
