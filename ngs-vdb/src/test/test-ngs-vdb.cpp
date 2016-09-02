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

#include <ktst/unit_test.hpp>

#include <ngs-vdb/inc/NGS-VDB.hpp>

#define __mod__     "TEST_NGS_VDB"
#define __file__    "test-ngs-vdb"
#define __fext__    "cpp"
#include <kfc/ctx.h>

#include <kfc/rsrc.h>
#include <kfc/except.h>

#include <../libs/ngs/NGS_ReadCollection.h>
#include <../libs/ngs/NGS_FragmentBlobIterator.h>
#include <../libs/ngs/NGS_FragmentBlob.h>
#include <../libs/ngs/NGS_Reference.h>
#include <../libs/ngs/NGS_ReferenceBlobIterator.h>
#include <../libs/ngs/NGS_ReferenceBlob.h>

using namespace std;
using namespace ncbi::NK;
using namespace ngs;
using namespace ncbi :: ngs :: vdb;

const String SRA_Accession = "SRR000001";
const String SRA_Accession_Prime = "SRR000002";
const String CSRA1_Accession = "SRR1063272";
const String CSRA1_Accession_Prime = "SRR1063273";
const String CSRA1_Accession_WithRepeats = "SRR600094";

TEST_SUITE(NgsVdbTestSuite);

#define ENTRY \
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing ); \
    m_ctx = ctx; \

#define EXIT \
    REQUIRE ( ! FAILED () ); \
    Release()

#define REQUIRE_FAILED() ( REQUIRE ( FAILED () ), CLEAR() )

class KfcFixture
{
public:
    KfcFixture()
    :   m_ctx(0)
    {
    }
    virtual ~KfcFixture()
    {
        Release();
    }

    virtual void Release()
    {
        m_ctx = 0;
    }

    const KCtx* m_ctx;  // points into the test case's local memory
};


class FragmentBlobFixture : public KfcFixture
{
public:
    FragmentBlobFixture()
    :   m_readColl ( 0 ),
        m_iter ( 0 )
    {
    }

    virtual void Release()
    {
        if (m_ctx != 0)
        {
            if ( m_iter != 0 )
            {
                NGS_FragmentBlobIteratorRelease ( m_iter, m_ctx );
            }
            if ( m_readColl != 0 )
            {
                NGS_ReadCollectionRelease ( m_readColl, m_ctx );
            }
            m_ctx = 0; // a pointer into the caller's local memory
        }
    }

    void MakeIterator ( ctx_t ctx, const string& p_acc  )
    {
        FUNC_ENTRY ( ctx, rcSRA, rcArc, rcAccessing );
        ON_FAIL ( m_readColl = NGS_ReadCollectionMake ( ctx, p_acc . c_str () ) )
        {
            throw logic_error ("NGS_ReadCollectionMake failed");
        }
        ON_FAIL ( m_iter = NGS_ReadCollectionGetFragmentBlobs ( m_readColl, ctx ) )
        {
            throw logic_error ("NGS_ReadCollectionGetFragmentBlobs failed");
        }
    }

    NGS_ReadCollection * m_readColl;
    NGS_FragmentBlobIterator* m_iter;
};


/// FragmentBlob

FIXTURE_TEST_CASE ( FragmentBlob_Create_Size, FragmentBlobFixture )
{
    ENTRY;
    MakeIterator ( ctx, SRA_Accession);

    TRY ( NGS_FragmentBlob* ref = NGS_FragmentBlobIteratorNext ( m_iter, ctx ) )
    {
        FragmentBlob b ( ref );
        REQUIRE_EQ ( (uint64_t)1080, b . Size() );
    }

    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlob_Assign, FragmentBlobFixture )
{
    ENTRY;
    MakeIterator ( ctx, SRA_Accession );

    FragmentBlob b1 ( NGS_FragmentBlobIteratorNext ( m_iter, ctx ) );
    FragmentBlob b2 ( NGS_FragmentBlobIteratorNext ( m_iter, ctx ) );
    b1 = b2;
    REQUIRE_EQ ( (uint64_t)913, b1 . Size() );

    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlob_Copy, FragmentBlobFixture )
{
    ENTRY;
    MakeIterator ( ctx, SRA_Accession );

    FragmentBlob b1 ( NGS_FragmentBlobIteratorNext ( m_iter, ctx ) );
    FragmentBlob b2 ( b1 );
    REQUIRE_EQ ( (uint64_t)1080, b2 . Size() );

    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlob_Data, FragmentBlobFixture )
{
    ENTRY;
    MakeIterator ( ctx, SRA_Accession );

    FragmentBlob b ( NGS_FragmentBlobIteratorNext ( m_iter, ctx ) );
    REQUIRE_EQ ( string ("TCAGAT"), string ( (const char*)b . Data(), 6 ) );

    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlob_GetFragmentInfo_Biological, FragmentBlobFixture )
{
    ENTRY;
    MakeIterator ( ctx, SRA_Accession );

    FragmentBlob b ( NGS_FragmentBlobIteratorNext ( m_iter, ctx ) );
    std::string fragId;
    uint64_t startInBlob = 0;
    uint64_t lengthInBases = 0;
    bool biological = false;
    b . GetFragmentInfo ( 300, fragId, startInBlob, lengthInBases, biological );
    REQUIRE_EQ ( SRA_Accession+".FR0.2", fragId );
    REQUIRE_EQ ( (uint64_t)288, startInBlob );
    REQUIRE_EQ ( (uint64_t)115, lengthInBases );
    REQUIRE    ( biological );

    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlob_GetRowRange, FragmentBlobFixture )
{
    ENTRY;
    MakeIterator ( ctx, SRA_Accession );

    FragmentBlob b ( NGS_FragmentBlobIteratorNext ( m_iter, ctx ) );
    int64_t first=0;
    uint64_t count=0;
    b . GetRowRange ( first, count );
    REQUIRE_EQ ( (int64_t)1, first );
    REQUIRE_EQ ( (uint64_t)4, count );

    b = NGS_FragmentBlobIteratorNext ( m_iter, ctx );
    b . GetRowRange ( first, count );
    REQUIRE_EQ ( (int64_t)5, first );
    REQUIRE_EQ ( (uint64_t)4, count );

    EXIT;
}

/// FragmentBlobIterator

FIXTURE_TEST_CASE ( FragmentBlobIterator_Create, KfcFixture )
{
    ENTRY;
    TRY ( NGS_ReadCollection * readColl = NGS_ReadCollectionMake ( ctx, SRA_Accession . c_str () ) )
    {
        TRY ( FragmentBlobIteratorRef ref = NGS_ReadCollectionGetFragmentBlobs ( readColl, ctx ) )
        {
            FragmentBlobIterator blobIt ( ref );
            //TODO: Verify
        }
        NGS_ReadCollectionRelease ( readColl, ctx );
    }
    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlobIterator_Assign, KfcFixture )
{
    ENTRY;
    TRY ( NGS_ReadCollection * readColl_1 = NGS_ReadCollectionMake ( ctx, SRA_Accession . c_str () ) )
    {
        TRY ( FragmentBlobIteratorRef ref1 = NGS_ReadCollectionGetFragmentBlobs ( readColl_1, ctx ) )
        {
            FragmentBlobIterator blobIt1 ( ref1 );

            TRY ( NGS_ReadCollection * readColl_2 = NGS_ReadCollectionMake ( ctx, SRA_Accession_Prime . c_str () ) )
            {
                TRY ( FragmentBlobIteratorRef ref2 = NGS_ReadCollectionGetFragmentBlobs ( readColl_2, ctx ) )
                {
                    FragmentBlobIterator blobIt2 ( ref2 );

                    blobIt2 = blobIt1;
                    //TODO: Verify
                }
                NGS_ReadCollectionRelease ( readColl_2, ctx );
            }
            NGS_ReadCollectionRelease ( readColl_1, ctx );
        }
    }
    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlobIterator_Copy, KfcFixture )
{
    ENTRY;
    TRY ( NGS_ReadCollection * readColl = NGS_ReadCollectionMake ( ctx, SRA_Accession . c_str () ) )
    {
        TRY ( FragmentBlobIteratorRef ref1 = NGS_ReadCollectionGetFragmentBlobs ( readColl, ctx ) )
        {
            FragmentBlobIterator blobIt1 ( ref1 );
            FragmentBlobIterator blobIt2 ( blobIt1 );
            //TODO: Verify
        }
        NGS_ReadCollectionRelease ( readColl, ctx );
    }
    EXIT;
}

TEST_CASE ( FragmentBlobIterator_HasMore )
{
    VdbReadCollection coll = NGS_VDB :: openVdbReadCollection ( SRA_Accession . c_str () );
    FragmentBlobIterator blobIt = coll . getFragmentBlobs ();
    REQUIRE ( blobIt . hasMore () );
}

TEST_CASE ( FragmentBlobIterator_nextBlob )
{
    VdbReadCollection coll = NGS_VDB :: openVdbReadCollection ( SRA_Accession . c_str () );
    FragmentBlobIterator blobIt = coll . getFragmentBlobs ();
    FragmentBlob blob = blobIt . nextBlob ();
    //TODO: Verify
}

/// VdbReadCollection

TEST_CASE ( VdbReadCollection_CreateFromReadCollection )
{
    VdbReadCollection coll ( ncbi :: NGS :: openReadCollection ( SRA_Accession ) );
    REQUIRE_EQ ( SRA_Accession, coll.toReadCollection().getName() );
}

TEST_CASE ( VdbReadCollection_CopyConstruct )
{
    VdbReadCollection coll1 ( NGS_VDB :: openVdbReadCollection ( SRA_Accession ) );
    VdbReadCollection coll2 ( coll1 );
    REQUIRE_EQ ( SRA_Accession, coll2.toReadCollection().getName() );
}

TEST_CASE ( VdbReadCollection_Assign )
{
    VdbReadCollection coll1 ( NGS_VDB :: openVdbReadCollection ( SRA_Accession ) );
    VdbReadCollection coll2 ( NGS_VDB :: openVdbReadCollection ( SRA_Accession_Prime ) );
    coll1 = coll2;
    REQUIRE_EQ ( SRA_Accession_Prime, coll1.toReadCollection().getName() );
}

TEST_CASE ( VdbReadCollection_GetFragmentBlobs )
{
    VdbReadCollection coll ( NGS_VDB :: openVdbReadCollection ( SRA_Accession ) );
    FragmentBlobIterator blobIt = coll . getFragmentBlobs ();
    //TODO: Verify
}

/// ReferenceBlob

class ReferenceBlobFixture : public KfcFixture
{
public:
    ReferenceBlobFixture()
    :   m_readColl ( 0 ),
        m_iter ( 0 )
    {
    }

    virtual void Release()
    {
        if (m_ctx != 0)
        {
            if ( m_iter != 0 )
            {
                NGS_ReferenceBlobIteratorRelease ( m_iter, m_ctx );
            }
            if ( m_readColl != 0 )
            {
                NGS_ReadCollectionRelease ( m_readColl, m_ctx );
            }
            m_ctx = 0; // a pointer into the caller's local memory
        }
    }

    void MakeBlobIterator ( ctx_t ctx, const string& p_acc, const string& p_refName )
    {
        FUNC_ENTRY ( ctx, rcSRA, rcArc, rcAccessing );
        ON_FAIL ( m_readColl = NGS_ReadCollectionMake ( ctx, p_acc . c_str () ) )
        {
            throw logic_error ("NGS_ReadCollectionMake failed");
        }
        TRY ( NGS_Reference* ref = NGS_ReadCollectionGetReference ( m_readColl, ctx, p_refName . c_str () ) )
        {
            ON_FAIL ( m_iter = NGS_ReferenceGetBlobs ( ref, ctx ) )
            {
                throw logic_error ("NGS_ReferenceGetBlobs failed");
            }
            ON_FAIL ( NGS_ReferenceRelease ( ref, ctx ) )
            {
                throw logic_error ("NGS_ReferenceRelease failed");
            }
        }
        CATCH_ALL()
        {
            throw logic_error ("NGS_ReadCollectionGetReference failed");
        }
    }

    ReferenceBlob GetBlobByNumber ( ctx_t ctx, const string& p_acc, const string& p_refName, unsigned int p_index = 1 )
    {
        FUNC_ENTRY ( ctx, rcSRA, rcArc, rcAccessing );
        MakeBlobIterator ( ctx, p_acc, p_refName );
        NGS_ReferenceBlob* blob = 0;
        for ( unsigned int i = 0 ; i < p_index; ++i )
        {
            ON_FAIL ( blob = NGS_ReferenceBlobIteratorNext ( m_iter, ctx ) )
            {
                throw logic_error ("NGS_ReferenceBlobIteratorNext failed");
            }
        }
        return ReferenceBlob ( blob );
    }

    ReferenceBlob GetBlobByRowId ( ctx_t ctx, const string& p_acc, const string& p_refName, unsigned int p_rowId )
    {
        FUNC_ENTRY ( ctx, rcSRA, rcArc, rcAccessing );
        MakeBlobIterator ( ctx, p_acc, p_refName );
        NGS_ReferenceBlob* blob = 0;
        do
        {
            ON_FAIL ( blob = NGS_ReferenceBlobIteratorNext ( m_iter, ctx ) )
            {
                throw logic_error ("NGS_ReferenceBlobIteratorNext failed");
            }
            int64_t firstRow;
            uint64_t rowCount;
            NGS_ReferenceBlobRowRange ( blob, ctx, &firstRow, &rowCount );
            if ( p_rowId >= firstRow && p_rowId < firstRow + rowCount )
            {
                return ReferenceBlob ( blob );
            }

        }
        while ( blob != 0 );
        throw logic_error ("GetBlobByRowId: row not found");
    }

    NGS_ReadCollection*         m_readColl;
    NGS_ReferenceBlobIterator*  m_iter;
};

FIXTURE_TEST_CASE ( ReferenceBlob_Create_Size, ReferenceBlobFixture )
{
    ENTRY;
    REQUIRE_EQ ( (uint64_t)20000, GetBlobByNumber ( ctx, CSRA1_Accession, "supercont2.1" ) . Size () );
    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlob_Create_Data, ReferenceBlobFixture )
{
    ENTRY;
    ReferenceBlob b = GetBlobByNumber ( ctx, CSRA1_Accession, "supercont2.1" );
    REQUIRE_EQ ( string("GAATTCTAAA"), string ( b . Data (), 10 ) );
    REQUIRE_EQ ( string("TGTCACATCA"), string ( b . Data () + b . Size () - 10, 10 ) );
    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlob_Assign, ReferenceBlobFixture )
{
    ENTRY;

    ReferenceBlob b1 ( GetBlobByNumber ( ctx, CSRA1_Accession, "supercont2.1", 1 ) );
    ReferenceBlob b2 ( GetBlobByNumber ( ctx, CSRA1_Accession, "supercont2.1", 2 ) );
    REQUIRE_NE ( string ( b1 . Data (), b1 . Size () ), string ( b2 . Data (), b2 . Size () ) );
    b1 = b2;
    REQUIRE_EQ ( string ( b1 . Data (), b1 . Size () ), string ( b2 . Data (), b2 . Size () ) );

    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlob_Copy, ReferenceBlobFixture )
{
    ENTRY;
    ReferenceBlob b1 ( GetBlobByNumber ( ctx, CSRA1_Accession, "supercont2.1", 1 ) );
    ReferenceBlob b2 ( b1 );
    REQUIRE_EQ ( string ( b1 . Data (), b1 . Size () ), string ( b2 . Data (), b2 . Size () ) );
    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlob_GetRowRange, ReferenceBlobFixture )
{
    ENTRY;

    ReferenceBlob b ( GetBlobByNumber ( ctx, CSRA1_Accession, "supercont2.1", 1 ) );
    int64_t first=0;
    uint64_t count=0;
    b . GetRowRange ( first, count );
    REQUIRE_EQ ( (int64_t)1, first );
    REQUIRE_EQ ( (uint64_t)4, count );

    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlob_FindRepeat_NotFound, ReferenceBlobFixture )
{
    ENTRY;

    ReferenceBlob b ( GetBlobByNumber ( ctx, CSRA1_Accession, "supercont2.1", 1 ) );
    uint64_t inReference = 0;
    uint64_t nextInBlob = 0;
    uint32_t repeatCount = 0;
    uint64_t increment = 999;
    REQUIRE ( ! b . FindRepeat ( 0, nextInBlob, inReference, repeatCount, increment ) );

    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlob_FindRepeat_Found, ReferenceBlobFixture )
{
    ENTRY;

    const int64_t repeatedRowId = 96;
    ReferenceBlob b ( GetBlobByRowId ( ctx, CSRA1_Accession_WithRepeats, "NC_000001.10", repeatedRowId ) );   /* this blob consists of 9 repeated all-N rows */
    int64_t first=0;
    uint64_t count=0;
    b . GetRowRange ( first, count );

    uint64_t nextInBlob = 0;
    uint64_t inReference = 0;
    uint32_t repeatCount = 0;
    uint64_t increment = 0;
    REQUIRE ( b . FindRepeat ( 0, nextInBlob, inReference, repeatCount, increment ) );
    REQUIRE_EQ ( (uint64_t)( repeatedRowId - first ) * 5000 , nextInBlob );
    REQUIRE_EQ ( (uint64_t)( repeatedRowId - 1 ) * 5000, inReference );
    REQUIRE_EQ ( (uint32_t)9, repeatCount );
    REQUIRE_EQ ( (uint64_t)5000, increment );

    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlob_GetResolveOffset, ReferenceBlobFixture )
{
    ENTRY;

    ReferenceBlob b ( GetBlobByNumber ( ctx, CSRA1_Accession, "supercont2.1", 2 ) );
    uint64_t inReference = 0;
    uint32_t repeatCount = 0;
    uint64_t increment = 999;
    b . ResolveOffset ( 1000, inReference, repeatCount, increment );
    REQUIRE_EQ ( (uint64_t)21000, inReference );
    REQUIRE_EQ ( (uint32_t)1, repeatCount );
    REQUIRE_EQ ( (uint64_t)0, increment ); // 0 if no repeat

    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlob_GetResolveOffset_WithRepeat, ReferenceBlobFixture )
{
    ENTRY;

    const int64_t repeatedRowId = 96;
    ReferenceBlob b ( GetBlobByRowId ( ctx, CSRA1_Accession_WithRepeats, "NC_000001.10", repeatedRowId ) );   /* this blob consists of 9 repeated all-N rows */
    int64_t first=0;
    uint64_t count=0;
    b . GetRowRange ( first, count );

    uint64_t inReference = 0;
    uint32_t repeatCount = 0;
    uint64_t increment = 0;
    b . ResolveOffset ( ( repeatedRowId - first ) * 5000 + 1000, inReference, repeatCount, increment );
    REQUIRE_EQ ( (uint64_t)( repeatedRowId - 1 ) * 5000 + 1000, inReference );
    REQUIRE_EQ ( (uint32_t)9, repeatCount );
    REQUIRE_EQ ( (uint64_t)5000, increment );

    EXIT;
}

/// ReferenceBlobIterator

FIXTURE_TEST_CASE ( ReferenceBlobIterator_Create_hasMore, ReferenceBlobFixture )
{
    ENTRY;
    MakeBlobIterator ( ctx, CSRA1_Accession, "supercont2.1" );
    ReferenceBlobIterator iter ( m_iter );
    REQUIRE ( iter. hasMore () );
    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlobIterator_Next, ReferenceBlobFixture )
{
    ENTRY;
    MakeBlobIterator ( ctx, CSRA1_Accession, "supercont2.1" );
    ReferenceBlobIterator iter ( m_iter );
    ReferenceBlob b = iter . nextBlob ();
    REQUIRE_EQ ( string("GAATTCTAAA"), string ( b . Data (), 10 ) );
    b = iter . nextBlob ();
    REQUIRE_EQ ( string("CTTCGAGCAC"), string ( b . Data (), 10 ) );
    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlobIterator_Assign, ReferenceBlobFixture )
{
    ENTRY;
    MakeBlobIterator ( ctx, CSRA1_Accession, "supercont2.1" );
    ReferenceBlobIterator iter1 ( m_iter );
    ReferenceBlobIterator iter2 ( m_iter );
    iter1 . nextBlob ();
    iter2 = iter1;
    //TODO: verify
    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlobIterator_Copy, ReferenceBlobFixture )
{
    ENTRY;
    MakeBlobIterator ( ctx, CSRA1_Accession, "supercont2.1" );
    ReferenceBlobIterator iter1 ( m_iter );
    ReferenceBlobIterator iter2 ( iter1 );
    //TODO: verify
    EXIT;
}

FIXTURE_TEST_CASE ( ReferenceBlobIterator_FullScan, ReferenceBlobFixture )
{
    ENTRY;
    MakeBlobIterator ( ctx, CSRA1_Accession, "supercont2.1" );
    ReferenceBlobIterator iter ( m_iter );
    size_t count = 0;
    while ( iter . hasMore () )
    {
        ReferenceBlob b = iter . nextBlob ();
        ++count;
    }
    REQUIRE_EQ ( (size_t)12, count );
    EXIT;
}

int
main( int argc, char *argv [] )
{
    return NgsVdbTestSuite(argc, argv);
}
