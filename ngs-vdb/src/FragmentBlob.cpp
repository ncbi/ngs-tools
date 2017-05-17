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

#include <ngs-vdb/inc/FragmentBlob.hpp>

#define __mod__     "NGS_VDB"
#define __file__    "FragmentBlob"
#define __fext__    "cpp"
#include <kfc/ctx.h>

#include <kfc/except.h>

#include <ngs/itf/ErrBlock.hpp>

#include <../libs/ngs/NGS_FragmentBlob.h>
#include <../libs/ngs/NGS_Id.h>
#include <../libs/ngs/NGS_ErrBlock.h>
#include <../libs/ngs/NGS_String.h>
#include <../libs/ngs/NGS_Id.h>

using namespace ncbi :: ngs :: vdb;

FragmentBlob :: FragmentBlob ( FragmentBlobRef ref ) throw ()
: self ( 0 )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    self = NGS_FragmentBlobDuplicate ( ref, ctx);
}

FragmentBlob &
FragmentBlob :: operator = ( const FragmentBlob & obj ) throw ( :: ngs :: ErrorMsg )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    TRY ( NGS_FragmentBlobRelease ( self, ctx) )
    {
        self = NGS_FragmentBlobDuplicate ( obj . self, ctx);
    }
    if ( FAILED () )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
    return *this;
}

FragmentBlob :: FragmentBlob ( const FragmentBlob & obj ) throw ( :: ngs :: ErrorMsg )
: self ( 0 )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    TRY ( NGS_FragmentBlobRelease ( self, ctx) )
    {
        self = NGS_FragmentBlobDuplicate ( obj . self, ctx);
    }
    if ( FAILED () )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
}

FragmentBlob :: ~ FragmentBlob () throw ()
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    NGS_FragmentBlobRelease ( self, ctx);
}

const char*
FragmentBlob :: Data() const throw ()
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    return (const char*) NGS_FragmentBlobData ( self, ctx );
}

uint64_t
FragmentBlob :: Size() const throw ()
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    uint64_t ret = 0;
    ON_FAIL ( ret = NGS_FragmentBlobSize ( self, ctx ) )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
    return ret;
}

void
FragmentBlob :: GetFragmentInfo ( uint64_t p_offset, std::string & p_fragId, uint64_t& p_startInBlob, uint64_t& p_lengthInBases, bool& p_biological ) const
    throw ( :: ngs :: ErrorMsg )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    int64_t rowId;
    uint64_t fragStart;
    uint64_t baseCount;
    int32_t fragNum;
    TRY ( NGS_FragmentBlobInfoByOffset ( self, ctx, p_offset, & rowId, & fragStart, & baseCount, & fragNum ) )
    {
        if ( fragNum >= 0 )
        {
            TRY ( const NGS_String * run = NGS_FragmentBlobRun ( self, ctx ) )
            {
                TRY ( const NGS_String* readId = NGS_IdMakeFragment ( ctx, run, false, rowId, fragNum ) )
                {
                    TRY ( p_fragId = std::string ( NGS_StringData ( readId, ctx ), NGS_StringSize ( readId, ctx ) ) );
                    {
                        p_biological = true;
                    }
                    NGS_StringRelease ( readId, ctx );
                }
            }
        }
        else
        {
            p_biological = false;
            p_fragId . empty ();
        }
    }
    if ( FAILED () )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
    else
    {
        p_startInBlob = fragStart;
        p_lengthInBases = baseCount;
    }
}

void
FragmentBlob :: GetRowRange ( int64_t& first, uint64_t& count ) const throw ( :: ngs :: ErrorMsg )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcArc, rcAccessing );
    ON_FAIL ( NGS_FragmentBlobRowRange ( self, ctx,  &first, &count ) )
    {
        :: ngs :: ErrBlock err;
        NGS_ErrBlockThrow ( &err, ctx );
        err.Throw();
    }
}
