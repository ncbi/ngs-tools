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

// this file/project is for the slow-running test cases of the sra-search unit test suite

// the faster test cases and main()
#include "test-sra-search.cpp"

FIXTURE_TEST_CASE ( SmithWaterman_ImperfectMatch, VdbSearchFixture )
{   // SW scoring function is different from Agrep's, so the results are slightly different
    m_settings . m_minScorePct =  90;
    SetupSingleThread ( "ATTAGCATTAGC", VdbSearch :: SmithWaterman, "SRR000001" );
    REQUIRE_EQ ( string ( "SRR000001.FR0.141" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.183" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR000001.FR0.2944" ), NextFragmentId () );
}

///////// Multi threading

FIXTURE_TEST_CASE ( SingleAccession_Threaded, VdbSearchFixture )
{
    const string Sra1 = "SRR600094";
    SetupMultiThread ( "ACGTAGGGTCC", VdbSearch :: NucStrstr, 2, Sra1 );

    REQUIRE_EQ ( Sra1 + ".FR1.101989", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.101990", NextFragmentId () );
}

FIXTURE_TEST_CASE ( SingleAccession_Threaded_Blobs, VdbSearchFixture )
{
    const string Sra1 = "SRR600094";
    SetupMultiThread ( "ACGTAGGGTCC", VdbSearch :: NucStrstr, 2, Sra1 );
    m_settings . m_useBlobSearch = true;

    REQUIRE_EQ ( Sra1 + ".FR1.101989", NextFragmentId () );
    REQUIRE_EQ ( Sra1 + ".FR0.101990", NextFragmentId () );
}

FIXTURE_TEST_CASE ( Threads_RandomCrash, VdbSearchFixture )
{   // used to crash randomly
    SetupMultiThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, 4, "SRR000001" ); // 4 blob-based threads on one run
    m_settings . m_useBlobSearch = true;

    unsigned int count = 0;
    while ( NextMatch () )  // used to have a random crash inside VDB
    {
        ++count;
    }
    REQUIRE_EQ ( 12u, count );
}

// Reference-driven mode

FIXTURE_TEST_CASE ( ReferenceDriven_SingleReference_SingleAccession, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    m_settings . m_useBlobSearch  = false;
    m_settings . m_references . push_back ( ReferenceSpec ( "NC_000007.13" ) );
    SetupSingleThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, "SRR600094" );

    REQUIRE_EQ ( string ( "SRR600094.FR1.1053649" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053650" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053648" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053651" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053652" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR1.1053653" ), NextFragmentId () );
    REQUIRE ( ! NextMatch () );
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
    REQUIRE_EQ ( string ( "SRR600094.FR1.101989" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.101990" ), NextFragmentId () );
            // REQUIRE_EQ ( string ( "SRR600094.FR1.101989" ), NextFragmentId () ); // used to be duplicates
            // REQUIRE_EQ ( string ( "SRR600094.FR0.101990" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.101991" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( ReferenceDriven_AllReferences_NoDuplicates_Blobs, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    m_settings . m_useBlobSearch  = true;
    SetupSingleThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, "SRR600094" );

    REQUIRE_EQ ( string ( "SRR600094.FR1.101989" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.101990" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.101991" ), NextFragmentId () );
}

FIXTURE_TEST_CASE ( ReferenceDriven_MultipleReferences_SingleAccession, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    m_settings . m_useBlobSearch  = true;
    m_settings . m_references . push_back ( ReferenceSpec ( "NC_000007.13" ) );
    m_settings . m_references . push_back ( ReferenceSpec ( "NC_000001.10" ) );
    SetupSingleThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb, "SRR600094" );

    // on NC_000007.13
    REQUIRE_EQ ( string ( "SRR600094.FR1.1053649" ),  NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053650" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053648" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053651" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053652" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR1.1053653" ), NextFragmentId () );
    // NC_000001.10
    REQUIRE_EQ ( string ( "SRR600094.FR1.101989" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.101990" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.101991" ), NextFragmentId () );
    // there are matches on other references, not reported here
    REQUIRE ( ! NextMatch () );
}

FIXTURE_TEST_CASE ( ReferenceDriven_MultipleReferences_MultipleAccessions, VdbSearchFixture )
{
    m_settings . m_referenceDriven = true;
    m_settings . m_useBlobSearch  = true;
    m_settings . m_accessions . push_back ( "SRR600095" );
    m_settings . m_accessions . push_back ( "SRR600094" );
    m_settings . m_references . push_back ( ReferenceSpec ( "NC_000007.13" ) );
    m_settings . m_references . push_back ( ReferenceSpec ( "NC_000001.10" ) );
    SetupSingleThread ( "ACGTAGGGTCC", VdbSearch :: FgrepDumb );

    // SRR600095, on NC_000007.13
    REQUIRE_EQ ( string ( "SRR600095.FR1.694078" ), NextFragmentId () );
    // SRR600094, NC_000001.10
    REQUIRE_EQ ( string ( "SRR600095.FR1.69793" ), NextFragmentId () );
    // SRR600094, on NC_000007.13
    REQUIRE_EQ ( string ( "SRR600094.FR1.1053649" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053650" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053648" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053651" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.1053652" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR1.1053653" ), NextFragmentId () );
    // SRR600094, NC_000001.10
    REQUIRE_EQ ( string ( "SRR600094.FR1.101989" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.101990" ), NextFragmentId () );
    REQUIRE_EQ ( string ( "SRR600094.FR0.101991" ), NextFragmentId () );
    // there are matches on other references, not reported here
    REQUIRE ( ! NextMatch () );
}



