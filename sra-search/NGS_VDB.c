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

#include "NGS_VDB.h"

#define __mod__ "sra-search"
#define __file__ "NGS_VDB"
#define __fext__ "c"
#include <kfc/ctx.h>

#include <kfc/xc.h>
#include <kfc/except.h>
#include <kfc/rsrc.h>

#include <klib/printf.h>

#include <vdb/database.h>
#include <vdb/cursor.h>
#include <vdb/table.h>
#include <vdb/blob.h>

extern "C"
{
    #include <../libs/vdb/blob-priv.h>
    #include <../libs/vdb/page-map.h>
}

struct NGS_VDB_ReadCollection {
    char*               name;
    const VDatabase*    db;
    VCursor*            curs;    
    uint32_t            read_col_idx;    
};

NGS_VDB_ReadCollection * 
NGS_VDB_ReadCollectionMake ( const char * spec )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );

    if ( spec == NULL )
        USER_ERROR ( xcParamNull, "NULL read-collection specification string" );
    else if ( spec [ 0 ] == 0 )
        USER_ERROR ( xcStringEmpty, "empty read-collection specification string" );
    else
    {
        rc_t rc;
        const VDatabase *db;

        /* the first order of work is to determine what type of object is in "spec" */
        const VDBManager * mgr = ctx -> rsrc -> vdb;
        assert ( mgr != NULL );

        /* try as VDB database */
        rc = VDBManagerOpenDBRead ( mgr, & db, NULL, "%s", spec );
        if ( rc == 0 ) 
        {   // WGS
            NGS_VDB_ReadCollection * ret = (NGS_VDB_ReadCollection *) malloc ( sizeof ( NGS_VDB_ReadCollection ) );
            if ( ret != NULL )
            {
                ret -> name = string_dup_measure ( spec, NULL );
                if ( ret -> name != NULL )
                {
                    ret -> db = db;
                    ret -> curs = NULL;
                    return ret;
                }
                SYSTEM_ERROR ( xcNoMemory, "allocating NGS_VDB_ReadCollection -> name ( '%s' )", spec );
                free ( ret );
            }
            VDatabaseRelease ( db );
            SYSTEM_ERROR ( xcNoMemory, "allocating NGS_VDB_ReadCollection ( '%s' )", spec );
        }
    }
    return NULL;
}

void 
NGS_VDB_ReadCollectionRelease ( NGS_VDB_ReadCollection * self )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );
    if ( self->curs != NULL )
    {
        VCursorRelease ( self->curs );
    }
    if ( self->db != NULL )
    {
        VDatabaseRelease ( self->db );
    }
    free ( self -> name );
    free ( self );
}

struct VBlob* 
NGS_VDB_ReadCollectionNextBlob ( NGS_VDB_ReadCollection * self, struct VBlob* p_blob )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );
    if ( p_blob == NULL )
    {   /* open cursor, retrieve the first blob */
        VTable* tbl;
        rc_t rc = VDatabaseOpenTableRead ( self->db, (const VTable**)&tbl, "SEQUENCE" );
        if ( rc == 0 )
        {
            rc = VTableCreateCursorRead ( tbl, (const VCursor**)&self->curs );
            if ( rc == 0 )
            {
                rc = VCursorAddColumn ( self->curs, &self->read_col_idx, "READ" );
                if ( rc == 0 )
                {
                    rc = VCursorOpen ( self->curs );
                    if ( rc == 0 )
                    {
                        rc = VCursorSetRowId ( self->curs, 1 );
                        if ( rc == 0 )
                        {
                            rc = VCursorOpenRow ( self->curs );
                            if ( rc == 0 )
                            {
                                struct VBlob *blob;
                                rc = VCursorGetBlob ( self->curs, (const VBlob**)&blob, self->read_col_idx );
                                if ( rc == 0 )
                                {
                                    VTableRelease ( tbl );
                                    return blob;
                                }
                                else
                                {
                                    INTERNAL_ERROR ( xcUnexpected, "VCursorGetBlob() rc = %R", rc );            
                                }
                            }                            
                            else
                            {
                                INTERNAL_ERROR ( xcUnexpected, "VCursorOpenRow() rc = %R", rc );            
                            }
                        }
                        else
                        {
                            INTERNAL_ERROR ( xcUnexpected, "VCursorSetRowId() rc = %R", rc );            
                        }
                    }
                    else
                    {
                        INTERNAL_ERROR ( xcUnexpected, "VCursorOpen() rc = %R", rc );            
                    }
                }
                else
                {
                    INTERNAL_ERROR ( xcUnexpected, "VCursorAddColumn() rc = %R", rc );            
                }
                VCursorRelease ( self->curs );
                self->curs = NULL;
            }
            else
            {
                INTERNAL_ERROR ( xcUnexpected, "VTableCreateCursorRead() rc = %R", rc );            
            }
            
            VTableRelease ( tbl );
        }
        else
        {
            INTERNAL_ERROR ( xcUnexpected, "VDatabaseOpenTableRead(SEQUENCE) rc = %R", rc );            
        }
        
    }
    return NULL;    
}

const char*
NGS_VDB_ReadCollectionRowIdToFragmentId ( NGS_VDB_ReadCollection * self, int64_t rowId )
{   // WGS only for now
    static char buf[1024]; 
    string_printf ( buf, sizeof ( buf ), NULL, "%s.FR0.%li", self -> name, rowId );
    return buf;
}


const void*     
NGS_VDB_BlobData ( struct VBlob* p_blob )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );
    return p_blob->data.base;
}

uint64_t        
NGS_VDB_BlobSize ( struct VBlob* p_blob )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );
    return KDataBufferBytes ( &p_blob->data );
}

void            
NGS_VDB_BlobRowInfo ( struct VBlob* p_blob,  uint64_t p_offset, int64_t* p_rowId, uint64_t* p_nextRowStart )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );
    
    *p_rowId = 0;
    *p_nextRowStart = 0;
    
    int64_t first;
    uint64_t count;
    rc_t rc = VBlobIdRange ( p_blob, &first, &count );
    if ( rc == 0  )
    {
        PageMapIterator pmIt;
        rc = PageMapNewIterator ( (const PageMap*)p_blob->pm, &pmIt, 0,count );
        if ( rc == 0 )
        {
            row_count_t rowInBlob = 0;
            do 
            {
                elem_count_t length = PageMapIteratorDataLength ( &pmIt );     
                elem_count_t offset = PageMapIteratorDataOffset ( &pmIt );   
                row_count_t  repeat = PageMapIteratorRepeatCount ( &pmIt );
                
                if ( p_offset < offset + length * repeat )
                {
                    *p_rowId = first + rowInBlob + ( p_offset - offset ) / length; 
                    *p_nextRowStart = offset + length; 
                    return;
                }       
                ++rowInBlob;
            }
            while ( PageMapIteratorNext ( &pmIt ) ); 
        }
        else
        {
            INTERNAL_ERROR ( xcUnexpected, "PageMapNewIterator() rc = %R", rc );            
        }
    }
    else
    {
        INTERNAL_ERROR ( xcUnexpected, "VBlobIdRange() rc = %R", rc );            
    }
    
    *p_rowId = 0;
    *p_nextRowStart = 0;
    // TODO: how to report an error?         
}    

void 
NGS_VDB_BlobRelease ( struct VBlob* p_blob )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );
    VBlobRelease ( p_blob );    
}
