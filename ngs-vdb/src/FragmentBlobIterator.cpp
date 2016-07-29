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

#include <ngs-vdb/inc/FragmentBlobIterator.hpp>

#define __mod__     "NGS_VDB"
#define __file__    "FragmentBlobIterator"
#define __fext__    "cpp"
#include <kfc/ctx.h>

#include <kfc/except.h>

#include <ngs/itf/ErrBlock.hpp>

#include <../libs/ngs/NGS_FragmentBlobIterator.h>
#include <../libs/ngs/NGS_ErrBlock.h>

using namespace ncbi :: ngs :: vdb;

FragmentBlobIterator :: FragmentBlobIterator ( FragmentBlobIteratorRef ref ) throw ()
: self ( ref )
{
}

FragmentBlobIterator &
FragmentBlobIterator :: operator = ( const FragmentBlobIterator & obj ) throw ( :: ngs :: ErrorMsg )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    TRY ( NGS_FragmentBlobIteratorRelease ( self, ctx) )
    {
        self = NGS_FragmentBlobIteratorDuplicate ( obj . self, ctx);
    }
    if ( FAILED () )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
    return *this;
}

FragmentBlobIterator :: FragmentBlobIterator ( const FragmentBlobIterator & obj ) throw ( :: ngs :: ErrorMsg )
: self ( 0 )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    TRY ( NGS_FragmentBlobIteratorRelease ( self, ctx) )
    {
        self = NGS_FragmentBlobIteratorDuplicate ( obj . self, ctx);
    }
    if ( FAILED () )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
}

FragmentBlobIterator :: ~ FragmentBlobIterator () throw ()
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    NGS_FragmentBlobIteratorRelease ( self, ctx);
}

bool
FragmentBlobIterator :: hasMore() const throw ( :: ngs :: ErrorMsg )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    bool ret = false;
    ON_FAIL ( ret = NGS_FragmentBlobIteratorHasMore ( self, ctx ) )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
    return ret;
}

FragmentBlob
FragmentBlobIterator :: nextBlob() throw ( :: ngs :: ErrorMsg )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    NGS_FragmentBlob* blob;
    ON_FAIL ( blob = NGS_FragmentBlobIteratorNext ( self, ctx ) )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
    if ( blob == 0 )
    {
        throw :: ngs :: ErrorMsg( "No more blobs" );
    }
    return FragmentBlob ( blob );
}

