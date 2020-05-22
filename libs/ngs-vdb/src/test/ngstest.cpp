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

/**
* Unit tests for low-level NGS functions
*/

// suppress macro max from windows.h
#define NOMINMAX

#define __mod__ "NGS_TEST"
#define __file__  __FILE__
#define __fext__  "cpp"

#include "ngs_c_fixture.hpp"

#include <kfc/xc.h>

#include <SRA_ReadGroupInfo.h>
#include <SRA_Statistics.h>
#include <NGS_Cursor.h>
#include <SRA_Read.h>

#include <NGS_Id.h>

#include <klib/namelist.h>

#include <kfg/kfg-priv.h>
#include <kfg/repository.h>

#include <vdb/table.h>
#include <vdb/database.h>
#include <vdb/blob.h>

#include <stdexcept>
#include <cstring>
#include <limits>
#include <cmath>

using namespace std;
using namespace ncbi::NK;

TEST_SUITE(NgsTestSuite);

const char* SRA_Accession = "SRR000001";
const char* SRA_Accession_WithReadGroups = "SRR006061";
const char* SRADB_Accession_WithBamHeader = "SRR600096";

class ReadGroupInfo_Fixture : public NGS_C_Fixture
{
public:
    ReadGroupInfo_Fixture()
    : m_tbl(0), m_rgi(0)
    {
    }
    ~ReadGroupInfo_Fixture()
    {
    }

    void MakeSRA( const char* acc )
    {
        m_tbl = openTable( acc );
        m_rgi = SRA_ReadGroupInfoMake ( m_ctx, m_tbl );
    }
    void MakeSRADB( const char* acc )
    {
        const VDatabase* db = openDB( acc );
        if ( VDatabaseOpenTableRead ( db, & m_tbl, "SEQUENCE" ) != 0 )
            throw logic_error ("ReadGroupInfo_Fixture::MakeSRADB VDatabaseOpenTableRead failed");
        VDatabaseRelease ( db );

        m_rgi = SRA_ReadGroupInfoMake ( m_ctx, m_tbl );
    }

    virtual void Release()
    {
        if (m_ctx != 0)
        {
            if (m_rgi != 0)
                SRA_ReadGroupInfoRelease ( m_rgi, m_ctx );
            if ( m_tbl != 0 )
                VTableRelease ( m_tbl );
        }
        NGS_C_Fixture :: Release ();
    }

    const VTable*               m_tbl;
    const SRA_ReadGroupInfo*    m_rgi;
};


FIXTURE_TEST_CASE ( ReadGroupInfo_Make, ReadGroupInfo_Fixture )
{
    ENTRY;
    MakeSRA ( SRA_Accession_WithReadGroups );

    REQUIRE ( ! FAILED () );
    REQUIRE_NOT_NULL ( m_rgi );

    EXIT;
}

FIXTURE_TEST_CASE ( ReadGroupInfo_Count, ReadGroupInfo_Fixture )
{
    ENTRY;
    MakeSRA ( SRA_Accession_WithReadGroups );

    REQUIRE_EQ ( (uint32_t)144, m_rgi -> count );

    EXIT;
}

FIXTURE_TEST_CASE ( ReadGroupInfo_Access, ReadGroupInfo_Fixture )
{
    ENTRY;
    MakeSRA ( SRA_Accession_WithReadGroups );
    REQUIRE ( ! FAILED () );

    REQUIRE_NOT_NULL ( m_rgi -> groups [ 2 ] . name );
    REQUIRE_EQ ( string ( "S104_V2" ), toString ( m_rgi -> groups [ 2 ] . name,    ctx ) );

    REQUIRE_NULL ( m_rgi -> groups [ 2 ] . bam_LB );
    REQUIRE_NULL ( m_rgi -> groups [ 2 ] . bam_SM );

    REQUIRE_EQ ( (uint64_t)3263,        m_rgi -> groups [ 2 ] . min_row );
    REQUIRE_EQ ( (uint64_t)6140,        m_rgi -> groups [ 2 ] . max_row );
    REQUIRE_EQ ( (uint64_t)2878,        m_rgi -> groups [ 2 ] . row_count );
    REQUIRE_EQ ( (uint64_t)759518,      m_rgi -> groups [ 2 ] . base_count );
    REQUIRE_EQ ( (uint64_t)653032,      m_rgi -> groups [ 2 ] . bio_base_count );

    EXIT;
}
#if SHOW_UNIMPLEMENTED
FIXTURE_TEST_CASE ( ReadGroupInfo_BamHeader, ReadGroupInfo_Fixture )
{
    ENTRY;
    MakeSRADB ( SRADB_Accession_WithBamHeader );
    REQUIRE ( ! FAILED () );

    REQUIRE_NOT_NULL ( m_rgi -> groups [ 0 ] . name );
    REQUIRE_EQ ( string ( "A1DLC.1" ), toString ( m_rgi -> groups [ 0 ] . name, ctx ) );

    REQUIRE_NOT_NULL ( m_rgi -> groups [ 0 ] . bam_LB );
    REQUIRE_EQ ( string ( "Solexa-112136" ), toString ( m_rgi -> groups [ 0 ] . bam_LB,  ctx ) );

    REQUIRE_NOT_NULL ( m_rgi -> groups [ 0 ] . bam_SM );
    REQUIRE_EQ ( string ( "12341_SN_05_1" ), toString ( m_rgi -> groups [ 0 ] . bam_SM,  ctx ) );

    EXIT;
}
#endif
FIXTURE_TEST_CASE ( ReadGroupInfo_Find_Found, ReadGroupInfo_Fixture )
{
    ENTRY;
    MakeSRA ( SRA_Accession_WithReadGroups );

    NGS_String * s = NGS_StringMake ( ctx, "S104_V2", strlen ( "S104_V2" ) );
    REQUIRE_EQ ( (uint32_t)2, SRA_ReadGroupInfoFind ( m_rgi, ctx, NGS_StringData(s, ctx), NGS_StringSize(s, ctx) ) );
    REQUIRE ( ! FAILED () );
    NGS_StringRelease ( s, ctx );

    EXIT;
}


#if 0
FIXTURE_TEST_CASE ( ReadGroupInfo_PrintAll, ReadGroupInfo_Fixture )
{
    ENTRY;
    MakeSRA ( SRA_Accession_WithReadGroups );
//    MakeSRA ( "SRR000001" );

    for ( uint32_t  i = 0; i < m_rgi -> count; ++i )
    {
        const NGS_String * name = m_rgi -> groups [ i ] . name;
        REQUIRE_NOT_NULL ( name );
        cout << NGS_StringData ( name, ctx ) << " "
             << m_rgi -> groups [ i ] . min_row         << " "
             << m_rgi -> groups [ i ] . max_row         << " "
             << m_rgi -> groups [ i ] . row_count       << " "
             << m_rgi -> groups [ i ] . base_count      << " "
             << m_rgi -> groups [ i ] . bio_base_count  << " "
             << endl;
    }

    EXIT;
}
#endif

// NGS_Id

class Id_Fixture : public NGS_C_Fixture
{
public:
    Id_Fixture()
    : run ( 0 ), id ( 0 )
    {
    }
    ~Id_Fixture()
    {
    }

    void Release()
    {
        if (m_ctx != 0)
        {
            if (run != 0)
                NGS_StringRelease ( run, m_ctx );
            if (id != 0)
                NGS_StringRelease ( id, m_ctx );
        }
        NGS_C_Fixture :: Release ();
    }

    void MakeId ( enum NGS_Object object, int64_t rowId )
    {
        run  = NGS_StringMake ( m_ctx, "run", strlen ( "run" ) );
        id = NGS_IdMake ( m_ctx, run, object, rowId );
    }
    void MakeFragmentId ( bool alignment, int64_t rowId, uint32_t frag_num )
    {
        run  = NGS_StringMake ( m_ctx, "run", strlen ( "run" ) );
        id = NGS_IdMakeFragment ( m_ctx, run, alignment, rowId, frag_num);
    }

    NGS_String * run;
    NGS_String * id;
};


//TODO: error cases

FIXTURE_TEST_CASE(NGS_IdMake_Read, Id_Fixture)
{
    ENTRY;

    MakeId ( NGSObject_Read, 12345678 );
    REQUIRE_EQ ( string ( "run.R.12345678" ), string ( NGS_StringData ( id, ctx ), NGS_StringSize ( id, ctx ) ) );

    EXIT;
}

FIXTURE_TEST_CASE(NGS_IdMake_Primary, Id_Fixture)
{
    ENTRY;

    MakeId ( NGSObject_PrimaryAlignment, 12345678 );
    REQUIRE_EQ ( string ( "run.PA.12345678" ), string ( NGS_StringData ( id, ctx ), NGS_StringSize ( id, ctx ) ) );

    EXIT;
}

FIXTURE_TEST_CASE(NGS_IdMake_Secondary, Id_Fixture)
{
    ENTRY;

    MakeId ( NGSObject_SecondaryAlignment, 12345678 );
    REQUIRE_EQ ( string ( "run.SA.12345678" ), string ( NGS_StringData ( id, ctx ), NGS_StringSize ( id, ctx ) ) );

    EXIT;
}

FIXTURE_TEST_CASE(NGS_IdMake_ReadFragment, Id_Fixture)
{
    ENTRY;

    MakeFragmentId ( false, 12345678, 1 );
    REQUIRE_EQ ( string ( "run.FR1.12345678" ), string ( NGS_StringData ( id, ctx ), NGS_StringSize ( id, ctx ) ) );

    EXIT;
}

FIXTURE_TEST_CASE(NGS_IdMake_AlignmentFragment, Id_Fixture)
{
    ENTRY;

    MakeFragmentId ( true, 12345678, 2 );
    REQUIRE_EQ ( string ( "run.FA2.12345678" ), string ( NGS_StringData ( id, ctx ), NGS_StringSize ( id, ctx ) ) );

    EXIT;
}

FIXTURE_TEST_CASE(NGS_Id_Parse_Read, Id_Fixture)
{
    ENTRY;

    MakeId ( NGSObject_Read, 12345678 );
    struct NGS_Id parsed = NGS_IdParse ( NGS_StringData(id, ctx), NGS_StringSize(id, ctx), ctx );
    REQUIRE_EQ ( string ( "run" ), string ( parsed . run . addr, parsed . run . len ) );
    REQUIRE_EQ ( (int32_t)NGSObject_Read, parsed . object );
    REQUIRE_EQ ( (int64_t)12345678, parsed . rowId );
    REQUIRE_EQ ( (uint32_t)0, parsed . fragId );

    EXIT;
}

FIXTURE_TEST_CASE(NGS_Id_Parse_Primary, Id_Fixture)
{
    ENTRY;

    MakeId ( NGSObject_PrimaryAlignment, 12345678 );
    struct NGS_Id parsed = NGS_IdParse ( NGS_StringData(id, ctx), NGS_StringSize(id, ctx), ctx );
    REQUIRE_EQ ( string ( "run" ), string ( parsed . run . addr, parsed . run . len ) );
    REQUIRE_EQ ( (int32_t)NGSObject_PrimaryAlignment, parsed . object );
    REQUIRE_EQ ( (int64_t)12345678, parsed . rowId );
    REQUIRE_EQ ( (uint32_t)0, parsed . fragId );

    EXIT;
}

FIXTURE_TEST_CASE(NGS_Id_Parse_Secondary, Id_Fixture)
{
    ENTRY;

    MakeId ( NGSObject_SecondaryAlignment, 12345678 );
    struct NGS_Id parsed = NGS_IdParse ( NGS_StringData(id, ctx), NGS_StringSize(id, ctx), ctx );
    REQUIRE_EQ ( string ( "run" ), string ( parsed . run . addr, parsed . run . len ) );
    REQUIRE_EQ ( (int32_t)NGSObject_SecondaryAlignment, parsed . object );
    REQUIRE_EQ ( (int64_t)12345678, parsed . rowId );
    REQUIRE_EQ ( (uint32_t)0, parsed . fragId );

    EXIT;
}

FIXTURE_TEST_CASE(NGS_Id_Parse_ReadFragment, Id_Fixture)
{
    ENTRY;

    MakeFragmentId ( false, 12345678, 1 );
    struct NGS_Id parsed = NGS_IdParse ( NGS_StringData(id, ctx), NGS_StringSize(id, ctx), ctx );
    REQUIRE_EQ ( string ( "run" ), string ( parsed . run . addr, parsed . run . len ) );
    REQUIRE_EQ ( (int32_t)NGSObject_ReadFragment, parsed . object );
    REQUIRE_EQ ( (int64_t)12345678, parsed . rowId );
    REQUIRE_EQ ( (uint32_t)1, parsed . fragId );

    EXIT;
}

FIXTURE_TEST_CASE(NGS_Id_Parse_AlignmentFragment, Id_Fixture)
{
    ENTRY;

    MakeFragmentId ( true, 12345678, 2 );
    struct NGS_Id parsed = NGS_IdParse ( NGS_StringData(id, ctx), NGS_StringSize(id, ctx), ctx );
    REQUIRE_EQ ( string ( "run" ), string ( parsed . run . addr, parsed . run . len ) );
    REQUIRE_EQ ( (int32_t)NGSObject_AlignmentFragment, parsed . object );
    REQUIRE_EQ ( (int64_t)12345678, parsed . rowId );
    REQUIRE_EQ ( (uint32_t)2, parsed . fragId );

    EXIT;
}

// NGS_Statistics
//TODO: decide whether to allow overwriting
//TODO: type conversions

TEST_CASE(NGS_Statistics_Make)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_AddU64)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddU64 ( stats, ctx, "path", 1 );
	REQUIRE_EQ ( ( uint32_t ) NGS_StatisticValueType_UInt64, NGS_StatisticsGetValueType ( stats, ctx, "path" ) );
    REQUIRE_EQ ( (uint64_t)1, NGS_StatisticsGetAsU64 ( stats, ctx, "path" ) );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_AddI64)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddI64 ( stats, ctx, "path", -12 );
	REQUIRE_EQ ( ( uint32_t ) NGS_StatisticValueType_Int64, NGS_StatisticsGetValueType ( stats, ctx, "path" ) );
    REQUIRE_EQ ( (int64_t)-12, NGS_StatisticsGetAsI64 ( stats, ctx, "path" ) );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_AddString)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    const char* cstr = "blah";
    NGS_String * str = NGS_StringMake ( ctx, cstr, strlen ( cstr ) );
    NGS_StatisticsAddString ( stats, ctx, "path", str );
	REQUIRE_EQ ( ( uint32_t ) NGS_StatisticValueType_String, NGS_StatisticsGetValueType ( stats, ctx, "path" ) );
    REQUIRE_EQ ( string ( cstr ), toString ( NGS_StatisticsGetAsString ( stats, ctx, "path" ), ctx, true ) );

    NGS_StringRelease ( str, ctx );
    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_AddDouble)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddDouble ( stats, ctx, "path", 3.1415926 );
	REQUIRE_EQ ( ( uint32_t ) NGS_StatisticValueType_Real, NGS_StatisticsGetValueType ( stats, ctx, "path" ) );
    REQUIRE_EQ ( 3.1415926, NGS_StatisticsGetAsDouble ( stats, ctx, "path" ) );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_AddNotANumber)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddDouble ( stats, ctx, "path", std::numeric_limits<double>::quiet_NaN() );
    REQUIRE_FAILED ();

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_OverwriteU64)
{   // currently, an attempt to overwrite a path throws
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddU64 ( stats, ctx, "path", 1 );
    NGS_StatisticsAddU64 ( stats, ctx, "path", 2 );
    REQUIRE_FAILED ();

    REQUIRE_EQ ( (uint64_t)1, NGS_StatisticsGetAsU64 ( stats, ctx, "path" ) );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_FindFound)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddU64 ( stats, ctx, "path1", 1 );
    NGS_StatisticsAddU64 ( stats, ctx, "path2", 2 );
    NGS_StatisticsAddU64 ( stats, ctx, "path3", 3 );

    REQUIRE_EQ ( (uint64_t)1, NGS_StatisticsGetAsU64 ( stats, ctx, "path1" ) );
    REQUIRE_EQ ( (uint64_t)2, NGS_StatisticsGetAsU64 ( stats, ctx, "path2" ) );
    REQUIRE_EQ ( (uint64_t)3, NGS_StatisticsGetAsU64 ( stats, ctx, "path3" ) );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}


TEST_CASE(NGS_Statistics_FindNotFound)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddU64 ( stats, ctx, "path1", 1 );
    NGS_StatisticsAddU64 ( stats, ctx, "path2", 2 );
    NGS_StatisticsAddU64 ( stats, ctx, "path3", 3 );

    NGS_StatisticsGetAsU64 ( stats, ctx, "path4" );
    REQUIRE_FAILED ();

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}


TEST_CASE(NGS_Statistics_Iterate)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddU64 ( stats, ctx, "path3", 3 );
    NGS_StatisticsAddU64 ( stats, ctx, "path1", 1 );
    NGS_StatisticsAddU64 ( stats, ctx, "path2", 2 );

    const char* path;
    REQUIRE ( NGS_StatisticsNextPath ( stats, ctx, "", & path ) );
    REQUIRE_EQ ( string ( "path1" ), string ( path ) );

    REQUIRE ( NGS_StatisticsNextPath ( stats, ctx, "path1", & path ) );
    REQUIRE_EQ ( string ( "path2" ), string ( path ) );

    REQUIRE ( NGS_StatisticsNextPath ( stats, ctx, "path2", & path ) );
    REQUIRE_EQ ( string ( "path3" ), string ( path ) );

    REQUIRE ( ! NGS_StatisticsNextPath ( stats, ctx, "path3", & path ) );
    REQUIRE_NULL ( path );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionU64)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddU64 ( stats, ctx, "path", 1 );
    REQUIRE_EQ ( (int64_t)1, NGS_StatisticsGetAsI64 ( stats, ctx, "path" ) );
    REQUIRE_EQ ( 1.0, NGS_StatisticsGetAsDouble( stats, ctx, "path" ) );
    REQUIRE_EQ ( string ( "1" ), toString ( NGS_StatisticsGetAsString( stats, ctx, "path" ), ctx , true ) );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionU64_Error)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    // MAX_U64 throws when reading as I64
    NGS_StatisticsAddU64 ( stats, ctx, "path", std::numeric_limits<uint64_t>::max() );
    NGS_StatisticsGetAsI64 ( stats, ctx, "path" );
    REQUIRE_FAILED ();

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionI64)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddI64 ( stats, ctx, "path", 1 );
    REQUIRE_EQ ( (uint64_t)1, NGS_StatisticsGetAsU64 ( stats, ctx, "path" ) );
    REQUIRE_EQ ( 1.0, NGS_StatisticsGetAsDouble( stats, ctx, "path" ) );
    REQUIRE_EQ ( string ( "1" ), toString ( NGS_StatisticsGetAsString( stats, ctx, "path" ), ctx , true ) );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionI64_Error)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    // negatives throw when reading as U64
    NGS_StatisticsAddI64 ( stats, ctx, "path", -1 );
    NGS_StatisticsGetAsU64 ( stats, ctx, "path" );
    REQUIRE_FAILED ();

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionReal)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddDouble ( stats, ctx, "path", 3.14 );

    // GetAsU64 truncates
    REQUIRE_EQ ( (uint64_t)3, NGS_StatisticsGetAsU64 ( stats, ctx, "path" ) );
    // GetAsI64 truncates
    REQUIRE_EQ ( (int64_t)3, NGS_StatisticsGetAsI64 ( stats, ctx, "path" ) );
    // GetAsString converts with a default precision */
    REQUIRE_EQ ( string ( "3.140000" ), toString ( NGS_StatisticsGetAsString( stats, ctx, "path" ), ctx , true ) );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionReal_Negative)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    NGS_StatisticsAddDouble ( stats, ctx, "path", -1.1 );
    // GetAsU64 throws
    NGS_StatisticsGetAsU64( stats, ctx, "path" );
    REQUIRE_FAILED ();
    // GetAsI64 truncates
    REQUIRE_EQ ( (int64_t)-1, NGS_StatisticsGetAsI64 ( stats, ctx, "path" ) );
    // GetAsString converts with a default precision */
    REQUIRE_EQ ( string ( "-1.100000" ), toString ( NGS_StatisticsGetAsString( stats, ctx, "path" ), ctx , true ) );

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionReal_ErrorSize)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    // throws when too big for a 64 bit number
    NGS_StatisticsAddDouble ( stats, ctx, "path", std::numeric_limits<double>::max() );
    NGS_StatisticsGetAsU64( stats, ctx, "path" );
    REQUIRE_FAILED ();
    NGS_StatisticsGetAsI64( stats, ctx, "path" );
    REQUIRE_FAILED ();

    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionString)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    const char* cstr = "   \t3.14"; /* leading space is ok */
    NGS_String * str = NGS_StringMake ( ctx, cstr, strlen ( cstr ) );
    NGS_StatisticsAddString ( stats, ctx, "path", str );

    // GetAsU64 truncates
    REQUIRE_EQ ( (uint64_t)3, NGS_StatisticsGetAsU64 ( stats, ctx, "path" ) );
    // GetAsI64 truncates
    REQUIRE_EQ ( (int64_t)3, NGS_StatisticsGetAsI64 ( stats, ctx, "path" ) );
    REQUIRE_EQ ( 3.14, NGS_StatisticsGetAsDouble( stats, ctx, "path" ) );

//TODO: more conversions to real
// "  +3.14"
// "  -3.14"
// "  -3E2"
// "  -3e2"
// "  -3.14e2"
// "  -0xF.0"
// "  -0xF.0P2" binary exponent P/p might not work on MSVS (MSVS doc on strtod mentions D/d without explanation)
// "  -0xFp2"
// "  -0xFp-2"
// "  INf"
// "  iNFiNiTY"
// "  NAN"
// "  nAn"
// "  naN(blah)"

    NGS_StringRelease ( str, ctx );
    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionString_BigUInt)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    const char* cstr = "18446744073709551615"; // std::numeric_limits<uint64_t>::max()
    NGS_String * str = NGS_StringMake ( ctx, cstr, strlen ( cstr ) );
    NGS_StatisticsAddString ( stats, ctx, "path", str );

    REQUIRE_EQ ( std::numeric_limits<uint64_t>::max(), NGS_StatisticsGetAsU64 ( stats, ctx, "path" ) );

    NGS_StringRelease ( str, ctx );
    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionString_BigUInt_Error)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    const char* cstr = "18446744073709551616"; // std::numeric_limits<uint64_t>::max() + 1
    NGS_String * str = NGS_StringMake ( ctx, cstr, strlen ( cstr ) );
    NGS_StatisticsAddString ( stats, ctx, "path", str );

    NGS_StatisticsGetAsU64 ( stats, ctx, "path" );
    REQUIRE_FAILED ();

    NGS_StringRelease ( str, ctx );
    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionString_BigInt)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    const char* cstr = "9223372036854775807"; // std::numeric_limits<int64_t>::max()
    NGS_String * str = NGS_StringMake ( ctx, cstr, strlen ( cstr ) );
    NGS_StatisticsAddString ( stats, ctx, "path", str );

    REQUIRE_EQ ( std::numeric_limits<int64_t>::max(), NGS_StatisticsGetAsI64 ( stats, ctx, "path" ) );

    NGS_StringRelease ( str, ctx );
    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionString_BigInt_Error)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    const char* cstr = "9223372036854775808"; // std::numeric_limits<int64_t>::max() + 1
    NGS_String * str = NGS_StringMake ( ctx, cstr, strlen ( cstr ) );
    NGS_StatisticsAddString ( stats, ctx, "path", str );

    NGS_StatisticsGetAsI64 ( stats, ctx, "path" );
    REQUIRE_FAILED ();

    NGS_StringRelease ( str, ctx );
    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

TEST_CASE(NGS_Statistics_ConversionString_TrailingSpace)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_Statistics * stats = SRA_StatisticsMake ( ctx );

    const char* cstr = "   \t3.14 \t\n  "; /* trailing space is an error*/
    NGS_String * str = NGS_StringMake ( ctx, cstr, strlen ( cstr ) );
    NGS_StatisticsAddString ( stats, ctx, "path", str );

    NGS_StatisticsGetAsDouble( stats, ctx, "path" );
    REQUIRE_FAILED ();

    NGS_StringRelease ( str, ctx );
    NGS_StatisticsRelease ( stats, ctx );
    REQUIRE ( ! FAILED () );
}

//////////////////////////////////////////// Errors opening read collection

#define BAD_ACCESSION "that refuses to open"
TEST_CASE(NGS_FailedToOpen)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );

    NGS_ReadCollectionMake ( ctx, BAD_ACCESSION);
    REQUIRE ( ctx_xc_isa ( ctx, xcTableOpenFailed ) );

    KConfig* kfg;
    REQUIRE_RC ( KConfigMakeLocal ( &kfg, NULL ) );
    const KRepositoryMgr* repoMgr;
    REQUIRE_RC ( KConfigMakeRepositoryMgrRead ( kfg, &repoMgr ) );
    if ( KRepositoryMgrHasRemoteAccess ( repoMgr ) )
    {
        string startsWith = "Cannot open accession '" BAD_ACCESSION "', rc = RC";
        REQUIRE_EQ ( startsWith, string ( WHAT () ) . substr ( 0, startsWith . size () ) );
    }
    else
    {
        string startsWith = "Cannot open accession '" BAD_ACCESSION "', rc = RC";
        string what = WHAT ();
        REQUIRE_EQ ( startsWith, what . substr ( 0, startsWith . size () ) );
        string endsWith = "Note: remote access is disabled in the configuration";
        REQUIRE_EQ ( endsWith, what . substr ( what . size() - endsWith . size () ) );
    }
    REQUIRE_FAILED ();
    REQUIRE_RC ( KRepositoryMgrRelease ( repoMgr ) );
    REQUIRE_RC ( KConfigRelease ( kfg ) );
}

//////////////////////////////////////////// Open read collection using system path

#ifdef WIN32
#define SYSPATH_ACCESSION ".\\data\\SysPathTest"
#else
#define SYSPATH_ACCESSION "./data/SysPathTest"
#endif
TEST_CASE(NGS_OpenBySysPath)
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcRow, rcAccessing );
    NGS_ReadCollection * readColl = NGS_ReadCollectionMake( ctx, SYSPATH_ACCESSION );
    REQUIRE ( ! FAILED () );
    REQUIRE_NOT_NULL ( readColl );
    NGS_ReadCollectionRelease ( readColl, ctx );
    REQUIRE ( ! FAILED () );
}


//////////////////////////////////////////// NGS_Cursor

class NGSCursorFixture : public NGS_C_Fixture
{
public:
    NGSCursorFixture()
    : m_cursor ( 0 )
    {
    }

    ~NGSCursorFixture()
    {
    }

    virtual void Release()
    {
        if (m_ctx != 0)
        {
            NGS_CursorRelease ( m_cursor, m_ctx );
        }
        NGS_C_Fixture :: Release ();
    }

    void MakeCursor ( const string& p_acc, const string& p_table, const char * col_specs[], uint32_t num_cols )
    {
        const VDatabase *db;
        THROW_ON_RC ( VDBManagerOpenDBRead ( m_ctx -> rsrc -> vdb, & db, NULL, "%s", p_acc . c_str () ) );

        VTable* tbl;
        THROW_ON_RC ( VDatabaseOpenTableRead ( db, (const VTable**)&tbl, p_table . c_str () ) );

        m_cursor = NGS_CursorMake ( m_ctx, tbl, col_specs, num_cols ); // this will add the first column to the cursor

        THROW_ON_RC ( VTableRelease ( tbl ) );
        THROW_ON_RC ( VDatabaseRelease ( db ) );
    }

    const NGS_Cursor* m_cursor;
};

FIXTURE_TEST_CASE ( NGS_CursorMakeDB_fails, NGSCursorFixture )
{
    ENTRY;
    const VDatabase *db;
    THROW_ON_RC ( VDBManagerOpenDBRead ( m_ctx -> rsrc -> vdb, & db, NULL, "%s", SRADB_Accession_WithBamHeader ) );

    NGS_String* runName = NGS_StringMake ( ctx, SRADB_Accession_WithBamHeader, string_size ( SRADB_Accession_WithBamHeader ) );
    m_cursor = NGS_CursorMakeDb ( m_ctx, db, runName, "bad table", sequence_col_specs, seq_NUM_COLS );
    REQUIRE ( FAILED () );
    REQUIRE_NULL ( m_cursor );
    // make sure RC is reported
    REQUIRE ( string ( WHAT () ) . find ( "rc =" ) != string :: npos);
    CLEAR ();

    NGS_StringRelease ( runName, ctx );

    REQUIRE_RC ( VDatabaseRelease ( db ) );
    EXIT;
}


FIXTURE_TEST_CASE ( NGS_Cursor_GetColumnIndex_adds_column, NGSCursorFixture )
{
    ENTRY;
    MakeCursor ( SRADB_Accession_WithBamHeader, "SEQUENCE", sequence_col_specs, seq_NUM_COLS ); // this will add the first column (READ) to the cursor
    REQUIRE ( ! FAILED () );
    REQUIRE_NOT_NULL ( m_cursor );

    REQUIRE_NE ( (uint32_t)0, NGS_CursorGetColumnIndex ( m_cursor, ctx, seq_READ_LEN ) ); // this should add the column we are requesting to the cursor

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_Cursor_Leak_when_Make_fails, NGSCursorFixture )
{   // use valgrind to detect the absence of a leak
    ENTRY;

    const char * bogus_col_specs [] = { "not a column at all!" };
    MakeCursor ( SRADB_Accession_WithBamHeader, "SEQUENCE", bogus_col_specs, 1 );
    REQUIRE ( FAILED () );
    REQUIRE_NULL ( m_cursor );
    CLEAR();

    EXIT;
}

FIXTURE_TEST_CASE ( NGS_Cursor_GetVBlob, NGSCursorFixture )
{
    ENTRY;
    MakeCursor ( SRADB_Accession_WithBamHeader, "SEQUENCE", sequence_col_specs, seq_NUM_COLS ); // this will add the first column (READ) to the cursor
    REQUIRE ( ! FAILED () );
    REQUIRE_NOT_NULL ( m_cursor );

    const int64_t rowId = 10;
    const struct VBlob* blob = NGS_CursorGetVBlob ( m_cursor, ctx, rowId, seq_READ );
    REQUIRE ( ! FAILED () );
    REQUIRE_NOT_NULL ( blob );

    int64_t first;
    uint64_t count;
    REQUIRE_RC ( VBlobIdRange ( blob, &first, &count ) );
    REQUIRE_LE ( first, rowId );
    REQUIRE_GT ( (int64_t)( first + count ), rowId );

    REQUIRE_RC ( VBlobRelease ( blob ) );

    EXIT;
}

//////////////////////////////////////////// Main

extern "C"
{

#include <kapp/args.h>

ver_t CC KAppVersion ( void )
{
    return 0x1000000;
}
rc_t CC UsageSummary (const char * progname)
{
    return 0;
}

rc_t CC Usage ( const Args * args )
{
    return 0;
}

const char UsageDefaultName[] = "test-ngs";

rc_t CC KMain ( int argc, char *argv [] )
{
    KConfigDisableUserSettings();
    int ret = NgsTestSuite(argc, argv);
    return (rc_t)ret;
}

}

