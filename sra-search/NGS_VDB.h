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

#ifndef _h_ngs_vdb_
#define _h_ngs_vdb_

#include <stdint.h>

struct VBlob;
typedef struct NGS_VDB_ReadCollection NGS_VDB_ReadCollection;

#ifdef __cplusplus
extern "C" {
#endif

/* TODO: error reporting */

NGS_VDB_ReadCollection *    NGS_VDB_ReadCollectionMake ( const char * spec );
void                        NGS_VDB_ReadCollectionRelease ( NGS_VDB_ReadCollection * self );
struct VBlob*               NGS_VDB_ReadCollectionNextBlob ( NGS_VDB_ReadCollection * self, struct VBlob* );
const char*                 NGS_VDB_ReadCollectionRowIdToFragmentId ( NGS_VDB_ReadCollection * self, int64_t rowId );

const void*     NGS_VDB_BlobData ( const struct VBlob* );
uint64_t        NGS_VDB_BlobSize ( const struct VBlob* );
void            NGS_VDB_BlobRowInfo ( const struct VBlob*,  uint64_t offset, int64_t* rowId, uint64_t* nextRowStart );
void            NGS_VDB_BlobRelease ( struct VBlob* );

#ifdef __cplusplus
}
#endif

#endif 
