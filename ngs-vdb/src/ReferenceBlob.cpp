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

#include <ngs-vdb/inc/ReferenceBlob.hpp>

#define __mod__     "NGS_VDB"
#define __file__    "ReferenceBlob"
#define __fext__    "cpp"
#include <kfc/ctx.h>

#include <kfc/except.h>

#include <ngs/itf/ErrBlock.hpp>

#include <../libs/ngs/NGS_ReferenceBlob.h>
#include <../libs/ngs/NGS_ErrBlock.h>
#include <../libs/ngs/NGS_String.h>

using namespace ncbi :: ngs :: vdb;

ReferenceBlob :: ReferenceBlob ( ReferenceBlobRef ref ) throw ()
: self ( 0 )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    self = NGS_ReferenceBlobDuplicate ( ref, ctx);
}

ReferenceBlob &
ReferenceBlob :: operator = ( const ReferenceBlob & obj ) throw ( :: ngs :: ErrorMsg )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    TRY ( NGS_ReferenceBlobRelease ( self, ctx) )
    {
        self = NGS_ReferenceBlobDuplicate ( obj . self, ctx);
    }
    if ( FAILED () )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
    return *this;
}

ReferenceBlob :: ReferenceBlob ( const ReferenceBlob & obj ) throw ( :: ngs :: ErrorMsg )
: self ( 0 )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    TRY ( NGS_ReferenceBlobRelease ( self, ctx) )
    {
        self = NGS_ReferenceBlobDuplicate ( obj . self, ctx);
    }
    if ( FAILED () )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
}

ReferenceBlob :: ~ ReferenceBlob () throw ()
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    NGS_ReferenceBlobRelease ( self, ctx);
}

const char*
ReferenceBlob :: Data() const throw ()
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    return (const char*) NGS_ReferenceBlobData ( self, ctx );
}

uint64_t
ReferenceBlob :: Size() const throw ()
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    uint64_t ret = 0;
    ON_FAIL ( ret = NGS_ReferenceBlobSize ( self, ctx ) )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
    return ret;
}

void
ReferenceBlob :: GetRowRange ( int64_t& p_first, uint64_t& p_count ) const throw ( :: ngs :: ErrorMsg )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    ON_FAIL ( NGS_ReferenceBlobRowRange ( self, ctx,  & p_first, & p_count ) )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
}

void
ReferenceBlob :: ResolveOffset ( uint64_t p_inBlob, uint64_t& p_inReference, uint32_t& p_repeatCount, uint64_t& p_increment ) const
                    throw ( :: ngs :: ErrorMsg )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    ON_FAIL ( NGS_ReferenceBlobResolveOffset ( self, ctx, p_inBlob, & p_inReference, & p_repeatCount, & p_increment ) )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
}

bool
ReferenceBlob :: FindRepeat ( uint64_t p_startInBlob, uint64_t& p_nextInBlob, uint64_t& p_inReference, uint32_t& p_repeatCount, uint64_t& p_increment ) const
                    throw ( :: ngs :: ErrorMsg )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    bool ret = false;
    ON_FAIL ( ret = NGS_ReferenceBlobFindRepeat ( self, ctx, p_startInBlob, & p_nextInBlob, & p_inReference, & p_repeatCount, & p_increment ) )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
    return ret;
}
