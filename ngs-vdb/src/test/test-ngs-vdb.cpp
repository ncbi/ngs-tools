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

using namespace std;
using namespace ncbi::NK;
using namespace ngs;
using namespace ncbi :: ngs :: vdb;

const String SRA_Accession = "SRR000001";
const String SRA_Accession_Prime = "SRR000002";

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
    :   m_ctx(0),
        m_readColl(0),
        m_iter(0)
    {
    }
    ~KfcFixture()
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

    void MakeIterator ( const string& p_acc, ctx_t ctx )
    {
        FUNC_ENTRY ( ctx, rcSRA, rcArc, rcAccessing );
        ON_FAIL ( m_readColl = NGS_ReadCollectionMake ( ctx, p_acc . c_str () ) )
        {
            throw ("NGS_ReadCollectionMake failed");
        }
        ON_FAIL ( m_iter = NGS_ReadCollectionGetFragmentBlobs ( m_readColl, ctx ) )
        {
            throw ("NGS_ReadCollectionGetFragmentBlobs failed");
        }
    }

    const KCtx* m_ctx;  // points into the test case's local memory
    NGS_ReadCollection * m_readColl;
    NGS_FragmentBlobIterator* m_iter;
};


/// FragmentBlob

FIXTURE_TEST_CASE ( FragmentBlob_Create_Size, KfcFixture )
{
    ENTRY;
    MakeIterator ( SRA_Accession, ctx );

    TRY ( NGS_FragmentBlob* ref = NGS_FragmentBlobIteratorNext ( m_iter, ctx ) )
    {
        REQUIRE_EQ ( (uint64_t)1080, FragmentBlob ( ref ) . Size() );
    }

    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlob_Assign, KfcFixture )
{
    ENTRY;
    MakeIterator ( SRA_Accession, ctx );

    FragmentBlob b1 ( NGS_FragmentBlobIteratorNext ( m_iter, ctx ) );
    FragmentBlob b2 ( NGS_FragmentBlobIteratorNext ( m_iter, ctx ) );
    b1 = b2;
    REQUIRE_EQ ( (uint64_t)913, b1 . Size() );

    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlob_Copy, KfcFixture )
{
    ENTRY;
    MakeIterator ( SRA_Accession, ctx );

    FragmentBlob b1 ( NGS_FragmentBlobIteratorNext ( m_iter, ctx ) );
    FragmentBlob b2 ( b1 );
    REQUIRE_EQ ( (uint64_t)1080, b2 . Size() );

    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlob_Data, KfcFixture )
{
    ENTRY;
    MakeIterator ( SRA_Accession, ctx );

    FragmentBlob b ( NGS_FragmentBlobIteratorNext ( m_iter, ctx ) );
    REQUIRE_EQ ( string ("TCAGAT"), string ( (const char*)b . Data(), 6 ) );

    EXIT;
}

FIXTURE_TEST_CASE ( FragmentBlob_GetFragmentInfo_Biological, KfcFixture )
{
    ENTRY;
    MakeIterator ( SRA_Accession, ctx );

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

FIXTURE_TEST_CASE ( FragmentBlob_GetRowRange, KfcFixture )
{
    ENTRY;
    MakeIterator ( SRA_Accession, ctx );

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

int
main( int argc, char *argv [] )
{
    return NgsVdbTestSuite(argc, argv);
}
