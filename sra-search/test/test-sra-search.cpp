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

#include "vdb-search.cpp"

#include <ktst/unit_test.hpp> 

using namespace std;
using namespace ncbi::NK;
using namespace ngs;

TEST_SUITE(SraSearchTestSuite);

TEST_CASE ( Create_Destroy )
{
    VdbSearch s ( "ACGT" );
}

TEST_CASE ( Algorithm_Default )
{
    VdbSearch s ( "ACGT" );
    REQUIRE_EQ ( VdbSearch :: Default, s . GetAlgorithm () ); 
}

TEST_CASE ( SetAlgorithm_Bad )
{
    VdbSearch s ( "ACGT" );
    REQUIRE ( ! s . SetAlgorithm ( "bad" ) );
}

TEST_CASE ( SetAlgorithmString_FgrepDumb )
{
    VdbSearch s ( "ACGT" );
    REQUIRE ( s . SetAlgorithm ( "FgrepDumb" ) );
    REQUIRE_EQ ( VdbSearch :: FgrepDumb, s . GetAlgorithm () ); 
}
TEST_CASE ( SetAlgorithmString_FgrepBoyerMoore )
{
    VdbSearch s ( "ACGT" );
    REQUIRE ( s . SetAlgorithm ( "FgrepBoyerMoore" ) );
    REQUIRE_EQ ( VdbSearch :: FgrepBoyerMoore, s . GetAlgorithm () ); 
}

TEST_CASE ( SingleAccession_FirstMatch )
{
    const string Accession = "SRR000001";
    const string Query = "AC";

    VdbSearch s ( Query );
    s . SetAlgorithm ( VdbSearch :: FgrepDumb );
    s . AddAccession ( Accession );
    string accession;
    string fragmentId;
    REQUIRE ( s . NextMatch ( accession, fragmentId ) );
    REQUIRE_EQ ( Accession, accession );
    REQUIRE_EQ ( string ( "SRR000001.FR0.1" ), fragmentId );
}

TEST_CASE ( FgrepDumb_SingleAccession_MultipleHitsAcrossFragments )
{
    const string Accession = "SRR000001";
    const string Query = "ATTAGC";

    VdbSearch s ( Query );
    s . SetAlgorithm ( VdbSearch :: FgrepDumb );
    s . AddAccession ( Accession );
    
    string accession;
    string fragmentId;
    REQUIRE ( s . NextMatch ( accession, fragmentId ) );
    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), fragmentId );
    
    REQUIRE ( s . NextMatch ( accession, fragmentId ) );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), fragmentId );
    
    REQUIRE ( s . NextMatch ( accession, fragmentId ) );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), fragmentId );
    
    //TODO: verify match positions when supported
}

TEST_CASE ( FgrepBoyerMoore_SingleAccession_MultipleHitsAcrossFragments )
{
    const string Accession = "SRR000001";
    const string Query = "ATTAGC";

    VdbSearch s ( Query );
    s . SetAlgorithm ( VdbSearch :: FgrepBoyerMoore );
    s . AddAccession ( Accession );
    
    string accession;
    string fragmentId;
    REQUIRE ( s . NextMatch ( accession, fragmentId ) );
    REQUIRE_EQ ( string ( "SRR000001.FR0.23" ), fragmentId );
    
    REQUIRE ( s . NextMatch ( accession, fragmentId ) );
    REQUIRE_EQ ( string ( "SRR000001.FR0.36" ), fragmentId );
    
    REQUIRE ( s . NextMatch ( accession, fragmentId ) );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), fragmentId );
    
    //TODO: verify match positions when supported
}

#if REPORT_MULTIPLE_HITS_IN_ONE_FRAGMENT
TEST_CASE ( SingleAccession_MultipleHitsInsideOneFragment )
{
    const string Accession = "SRR000001";
    const string Query = "AT";

    VdbSearch s ( Query );
    s . SetAlgorithm ( VdbSearch :: FgrepDumb );
    s . AddAccession ( Accession );
    
    string accession;
    string fragmentId;
    REQUIRE ( s . NextMatch ( accession, fragmentId ) );
    REQUIRE_EQ ( string ( "SRR000001.FR0.1" ), fragmentId );
    
    REQUIRE ( s . NextMatch ( accession, fragmentId ) );
    REQUIRE_EQ ( string ( "SRR000001.FR0.1" ), fragmentId );
    //TODO: verify positions when supported
}
#endif

//TODO: bad query string
//TODO: bad accession name
//TODO: multiple accessions

int
main( int argc, char *argv [] )
{
    return SraSearchTestSuite(argc, argv);
}
