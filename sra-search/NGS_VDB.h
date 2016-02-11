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
#include <stdbool.h>

struct VBlob;
struct NGS_ErrBlock_v1;
typedef struct NGS_ErrBlock_v1 NGS_VDB_ErrBlock; 

typedef struct NGS_VDB_ReadCollection NGS_VDB_ReadCollection;

#ifdef __cplusplus
extern "C" {
#endif

/* TODO: error reporting */

NGS_VDB_ReadCollection *    NGS_VDB_ReadCollectionMake ( const char * spec, NGS_VDB_ErrBlock * err );
void                        NGS_VDB_ReadCollectionRelease ( NGS_VDB_ReadCollection * self, NGS_VDB_ErrBlock * err  );
struct VBlob*               NGS_VDB_ReadCollectionNextBlob ( NGS_VDB_ReadCollection * self, struct VBlob*, NGS_VDB_ErrBlock * err  );

/* sets *fragId to NULL for technical fragments */
void                        NGS_VDB_BlobRowInfo ( NGS_VDB_ReadCollection * self, const struct VBlob*,  uint64_t offset, const char** fragId, uint64_t* nextFragStart, bool* biological, NGS_VDB_ErrBlock * err  );

const void*     NGS_VDB_BlobData ( const struct VBlob* );
uint64_t        NGS_VDB_BlobSize ( const struct VBlob* );
void            NGS_VDB_BlobRelease ( struct VBlob* );

#ifdef __cplusplus
}
#endif

#endif 
