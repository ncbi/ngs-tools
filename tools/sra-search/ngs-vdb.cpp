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

#include "ngs-vdb.hpp"

#include <ngs/itf/ErrBlock.hpp>

using namespace ngs;
using namespace ngs::vdb;
using namespace std;

#include "NGS_VDB.h"

//////// class FragmentBlob

FragmentBlob :: FragmentBlob ( const FragmentBlob & obj ) throw ( ErrorMsg )
:   m_self ( obj . m_self ),
    m_coll ( obj . m_coll )
{
    NGS_VDB_BlobAddRef ( m_self );
}

FragmentBlob :: FragmentBlob ( BlobRef ref ) throw ()
{
    throw ( ErrorMsg ( " FragmentBlob :: FragmentBlob ( BlobRef ref ) not implemented " ) );
}

FragmentBlob &
FragmentBlob :: operator = ( const FragmentBlob & obj ) throw ( ErrorMsg )
{
    throw ( ErrorMsg ( " FragmentBlob :: operator = not implemented " ) );
}

FragmentBlob :: ~ FragmentBlob () throw ()
{
    NGS_VDB_BlobRelease ( m_self );
}

const char*
FragmentBlob :: Data() const throw ()
{
    if ( m_self == 0 )
    {
        throw ( ErrorMsg ( " FragmentBlob :: Data(NULL) " ) );
    }
    return (const char*)NGS_VDB_BlobData ( m_self );
}

uint64_t
FragmentBlob :: Size() const throw ()
{
    if ( m_self == 0 )
    {
        throw ( ErrorMsg ( " FragmentBlob :: Size(NULL) " ) );
    }
    return NGS_VDB_BlobSize ( m_self  );
}

void
FragmentBlob :: GetFragmentInfo ( uint64_t p_offset, string& p_fragId, uint64_t& p_nextFragStart, bool& p_biological ) const throw ( ErrorMsg )
{
    if ( m_self == 0 )
    {
        throw ( ErrorMsg ( " FragmentBlob :: GetRowInfo(NULL) " ) );
    }
    char* fragId;
    ErrBlock err;
    NGS_VDB_BlobRowInfo ( m_coll, m_self, p_offset, &fragId, &p_nextFragStart, &p_biological, & err );
    err . Check ();
    if ( fragId != 0 )
    {
        p_fragId = fragId;
        free ( fragId );
    }
    else
    {
        p_fragId = string ();
    }
}

FragmentBlob :: FragmentBlob ( struct NGS_VDB_ReadCollection* p_coll )
:   m_self ( 0 ),
    m_coll ( p_coll )
{
}

//////// class FragmentBlobIterator

FragmentBlobIterator :: FragmentBlobIterator ( BlobRef ref ) throw ( )
:   FragmentBlob ( ref )
{
    throw ( ErrorMsg ( " FragmentBlobIterator :: FragmentBlobIterator ( BlobRef ref ) not implemented " ) );
}

FragmentBlobIterator :: FragmentBlobIterator ( const FragmentBlobIterator & obj ) throw ( ErrorMsg )
:   FragmentBlob ( obj )
{
    throw ( ErrorMsg ( " FragmentBlobIterator :: FragmentBlobIterator ( const FragmentBlobIterator & obj ) not implemented " ) );
}

FragmentBlobIterator :: ~ FragmentBlobIterator () throw ()
{
}

FragmentBlobIterator &
FragmentBlobIterator :: operator = ( const FragmentBlobIterator & obj ) throw ( ErrorMsg )
{
    throw ( ErrorMsg ( " FragmentBlobIterator :: operator = not implemented " ) );
}

bool
FragmentBlobIterator :: nextBlob () throw ( ErrorMsg )
{
    ErrBlock err;
    struct VBlob* next = NGS_VDB_ReadCollectionNextBlob ( m_coll, m_self, &err );
    err . Check ();
    NGS_VDB_BlobRelease ( m_self );
    m_self = next;
    return m_self != 0;
}

FragmentBlobIterator :: FragmentBlobIterator ( struct NGS_VDB_ReadCollection* p_coll )
: FragmentBlob ( p_coll )
{
}

//////// class VdbReadCollection

VdbReadCollection :: VdbReadCollection ( const VdbReadCollection & obj ) throw ()
{
    throw ( ErrorMsg ( " VdbReadCollection :: VdbReadCollection ( const VdbReadCollection & obj ) not implemented " ) );
}

VdbReadCollection :: VdbReadCollection ( ReadCollectionRef ref ) throw ()
{
    throw ( ErrorMsg ( " VdbReadCollection :: VdbReadCollection ( ReadCollectionRef ref ) not implemented " ) );
}

VdbReadCollection :: VdbReadCollection ( const String & spec ) throw ()
{
    ErrBlock err;
    m_coll = NGS_VDB_ReadCollectionMake ( spec.c_str(), & err );
    err . Check ();
}

VdbReadCollection :: ~ VdbReadCollection () throw ()
{
    ErrBlock err;
    NGS_VDB_ReadCollectionRelease ( m_coll, & err );
}

VdbReadCollection &
VdbReadCollection :: operator = ( const VdbReadCollection & obj ) throw ()
{
    throw ( ErrorMsg ( " VdbReadCollection :: operator = not implemented " ) );
}

FragmentBlobIterator
VdbReadCollection :: getFragmentBlobs () const throw ( ErrorMsg )
{
    return FragmentBlobIterator ( m_coll );
}

//////// class  ncbi ::  NGS_VDB

VdbReadCollection
ncbi ::  NGS_VDB :: openVdbReadCollection ( const String & spec ) throw ( ErrorMsg )
{
    return VdbReadCollection ( spec );
}

