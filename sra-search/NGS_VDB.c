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
#include <klib/rc.h>

#include <kfg/kfg-priv.h>
#include <kfg/repository.h>

#include <vdb/database.h>
#include <vdb/cursor.h>
#include <vdb/table.h>
#include <vdb/blob.h>

#include <sra/sraschema.h>

#include <stdio.h>

#include <../libs/vdb/blob-priv.h>
#include <../libs/vdb/page-map.h>

#include <../libs/ngs/NGS_Id.h>
#include <../libs/ngs/NGS_String.h>
#include <../libs/ngs/NGS_Cursor.h>
#include <../libs/ngs/NGS_ErrBlock.h>

#include <../libs/ngs/SRA_Read.h>

struct NGS_VDB_ReadCollection {
    NGS_String*         name;
    const NGS_Cursor*   curs;    
    NGS_String*         last_frag_id;
};

static 
VTable*
GetTable ( ctx_t ctx, const char * spec )
{
    FUNC_ENTRY ( ctx, rcSRA, rcDatabase, rcConstructing );
    
    rc_t rc;
    const VDatabase *db;
    VTable* ret;

    const struct VDBManager * mgr = ctx -> rsrc -> vdb;
    assert ( mgr != NULL );

    /* try as VDB database */
    rc = VDBManagerOpenDBRead ( mgr, & db, NULL, "%s", spec );
    if ( rc == 0 )
    {
        rc_t rc = VDatabaseOpenTableRead ( db, (const VTable**)&ret, "SEQUENCE" );
        if ( rc == 0 )
        {
            VDatabaseRelease ( db );
            return ret;
        }
        INTERNAL_ERROR ( xcUnimplemented, "Cannot open accession '%s' as an SRA database.", spec );
    }
    else 
    {   /* try as VDB table */
        VSchema *sra_schema;
        rc = VDBManagerMakeSRASchema ( mgr, & sra_schema );
        if ( rc != 0 )
            INTERNAL_ERROR ( xcUnexpected, "failed to make default SRA schema: rc = %R", rc );
        else
        {
            rc = VDBManagerOpenTableRead ( mgr, (const VTable**)&ret, sra_schema, "%s", spec );
            VSchemaRelease ( sra_schema );

            if ( rc == 0 )
            {   /* VDB-2641: examine the schema name to make sure this is an SRA table */
                char ts_buff[1024];
                rc = VTableTypespec ( ret, ts_buff, sizeof ( ts_buff ) );
                if ( rc != 0 )
                {
                    INTERNAL_ERROR ( xcUnexpected, "VTableTypespec failed: rc = %R", rc );
                }
                else
                {
                    const char SRA_PREFIX[] = "NCBI:SRA:";
                    size_t pref_size = sizeof ( SRA_PREFIX ) - 1;
                    if ( string_match ( SRA_PREFIX, pref_size, ts_buff, string_size ( ts_buff ), ( uint32_t ) pref_size, NULL ) == pref_size )
                    {
                        return ret;
                    }
                    INTERNAL_ERROR ( xcUnimplemented, "Cannot open accession '%s' as an SRA table.", spec );
                }
                VTableRelease ( ret );
            }
            else
            {
                KConfig* kfg = NULL;
                const KRepositoryMgr* repoMgr = NULL;
                if ( KConfigMakeLocal ( & kfg, NULL ) != 0 || 
                        KConfigMakeRepositoryMgrRead ( kfg, & repoMgr ) != 0 ||
                        KRepositoryMgrHasRemoteAccess ( repoMgr ) )
                {
                    INTERNAL_ERROR ( xcUnimplemented, "Cannot open accession '%s'.", spec );
                }
                else
                {
                    INTERNAL_ERROR ( xcUnimplemented, "Cannot open accession '%s'. Note: remote access is disabled in the configuration.", spec );
                }
                KRepositoryMgrRelease ( repoMgr );
                KConfigRelease ( kfg );
            }
        }
    }
    return NULL;
} 

NGS_VDB_ReadCollection * 
NGS_VDB_ReadCollectionMake ( const char * spec, NGS_VDB_ErrBlock * err  )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );

    if ( spec == NULL )
        USER_ERROR ( xcParamNull, "NULL read-collection specification string" );
    else if ( spec [ 0 ] == 0 )
        USER_ERROR ( xcStringEmpty, "empty read-collection specification string" );
    else
    {
        TRY ( VTable* tbl = GetTable ( ctx, spec ) )
        {   /* open cursor */
            TRY ( const NGS_Cursor* curs = NGS_CursorMake ( ctx, tbl, sequence_col_specs, seq_NUM_COLS ) )
            {
                NGS_VDB_ReadCollection * ret = (NGS_VDB_ReadCollection *) malloc ( sizeof ( NGS_VDB_ReadCollection ) );
                VTableRelease ( tbl );
                if ( ret != NULL )
                {
                    ret -> name = NGS_StringMakeCopy ( ctx, spec, string_size ( spec ) );
                    if ( ret -> name != NULL )
                    {
                        ret -> curs = curs;
                        ret -> last_frag_id = NULL;
                        return ret;
                    }
                    SYSTEM_ERROR ( xcNoMemory, "allocating NGS_VDB_ReadCollection -> name ( '%s' )", spec );
                    free ( ret );
                }
                else
                {
                    SYSTEM_ERROR ( xcNoMemory, "allocating NGS_VDB_ReadCollection ( '%s' )", spec );
                }
                NGS_CursorRelease ( curs, ctx );
            }
            VTableRelease ( tbl );
        }
    }
    NGS_ErrBlockThrow ( err, ctx );
    return NULL;
}

void 
NGS_VDB_ReadCollectionRelease ( NGS_VDB_ReadCollection * self, NGS_VDB_ErrBlock * err  )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );
    if ( self -> last_frag_id != NULL )
    {
        NGS_StringRelease ( self -> last_frag_id, ctx );
    }
    NGS_StringRelease ( self -> name, ctx );
    NGS_CursorRelease ( self -> curs, ctx );
    free ( self );
    NGS_ErrBlockThrow ( err, ctx );
}

static
struct VBlob*
NGS_VDB_ReadCollectionGetBlob ( ctx_t ctx, NGS_VDB_ReadCollection * self, int64_t p_rowId )
{
    FUNC_ENTRY ( ctx, rcSRA, rcDatabase, rcConstructing );
    
    TRY ( const struct VCursor* vcurs = NGS_CursorGetVCursor ( self -> curs ) )
    {
        rc_t rc = VCursorSetRowId ( vcurs, p_rowId );
        if ( rc == 0 )
        {
            rc = VCursorOpenRow ( vcurs );
            if ( rc == 0 )
            {
                struct VBlob *ret = NULL;
                rc = VCursorGetBlob ( vcurs, (const VBlob**)&ret, NGS_CursorGetColumnIndex ( self -> curs, ctx, seq_READ ) );
                if ( rc == 0 || 
                    GetRCObject ( rc ) == rcRow && GetRCState( rc ) == rcNotFound )
                {
                    rc = VCursorCloseRow ( vcurs );
                    if ( rc == 0 )
                    {
                        return ret;                
                    }
                }
                else
                {
                    VCursorCloseRow ( vcurs );
                }
                INTERNAL_ERROR ( xcUnexpected, "VCursorGetBlob() rc = %R", rc );            
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
    return NULL;
}

struct VBlob* 
NGS_VDB_ReadCollectionNextBlob ( NGS_VDB_ReadCollection * self, struct VBlob* p_blob, NGS_VDB_ErrBlock * err )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );
    struct VBlob* ret = NULL;
    if ( p_blob == NULL )
    {
        ret = NGS_VDB_ReadCollectionGetBlob ( ctx, self, 1 );
    }
    else
    {
        int64_t first;
        uint64_t count;
        rc_t rc = VBlobIdRange ( p_blob, &first, &count );
        if ( rc == 0  )
        {
            ret = NGS_VDB_ReadCollectionGetBlob ( ctx, self, first + count );
        }
        else
        {
            INTERNAL_ERROR ( xcUnexpected, "VBlobIdRange() rc = %R", rc );            
        }
    }
    NGS_ErrBlockThrow ( err, ctx );
    return ret;    
}

const void*     
NGS_VDB_BlobData ( const struct VBlob* p_blob )
{
    return p_blob->data.base;
}

uint64_t        
NGS_VDB_BlobSize ( const struct VBlob* p_blob )
{
    return KDataBufferBytes ( &p_blob->data );
}

static
void 
GetFragment ( ctx_t ctx, NGS_VDB_ReadCollection * self, int64_t p_rowId, uint64_t p_offsetInRow, uint32_t* p_bioFragNum, uint32_t* p_nextFragStart, bool* p_biological )
{   
    FUNC_ENTRY ( ctx, rcSRA, rcDatabase, rcConstructing );
    uint32_t elem_bits; 
    const void *base;
    uint32_t boff; 
    uint32_t row_len;
    TRY ( NGS_CursorCellDataDirect ( self -> curs, 
                                     ctx,
                                     p_rowId,
                                     seq_READ_LEN, 
                                     & elem_bits, 
                                     & base,
                                     & boff, 
                                     & row_len ) )
    {
        uint32_t i = 0 ; 
        uint64_t offset = 0;
        uint32_t bioFragNum = 0;
        assert ( base != NULL ); 
        assert ( elem_bits % 8 == 0 );
        assert ( boff == 0 );
        
        while ( i < row_len )
        {
            uint64_t frag_length;
            switch ( elem_bits )
            {
                case 64:
                {
                    frag_length = ( (const uint64_t*)base ) [ i ];
                    break;
                }
                case 32:
                {
                    frag_length = ( (const uint32_t*)base ) [ i ];
                    break;
                }
                case 16:
                {
                    frag_length = ( (const uint16_t*)base ) [ i ];
                    break;
                }
                case 8:
                {
                    frag_length = ( (const uint8_t*)base ) [ i ];
                    break;
                }
                default:
                {
                    INTERNAL_ERROR ( xcUnexpected, "Unexpected elem_bits: %u", elem_bits );
                    return;
                }
            }
            
            {
                uint32_t frag_type_elem_bits; 
                const void *frag_type_base;
                uint32_t frag_type_boff; 
                uint32_t frag_type_row_len;
                TRY ( NGS_CursorCellDataDirect ( self -> curs, 
                                                ctx,
                                                p_rowId,
                                                seq_READ_TYPE, 
                                                & frag_type_elem_bits, 
                                                & frag_type_base,
                                                & frag_type_boff, 
                                                & frag_type_row_len ) )
                {
                    const uint8_t* frag_types = (const uint8_t*)frag_type_base;
                    bool isBiological; 
                    assert ( frag_type_row_len == row_len );
                    assert ( frag_type_base != NULL ); 
                    assert ( frag_type_elem_bits == 8 );
                    assert ( frag_type_boff == 0 );
                    isBiological = frag_types [ i ] & READ_TYPE_BIOLOGICAL;
                    if ( p_offsetInRow < offset + frag_length )
                    {
                        *p_nextFragStart = ( uint32_t ) ( offset + frag_length );
                        *p_biological = isBiological;
                        *p_bioFragNum = *p_biological ? bioFragNum : 0;
                        return;
                    }
                    
                    if ( isBiological )
                    {
                        ++ bioFragNum;
                    }
                }
            }
            offset += frag_length;
            ++i;
        }
        INTERNAL_ERROR ( xcUnexpected, "Invalid READ_LEN, SEQUENCE.rowdId = %li, offset = %lu", p_rowId, p_offsetInRow );
    } 
}

void            
NGS_VDB_BlobRowInfo ( NGS_VDB_ReadCollection * self, const struct VBlob* p_blob, uint64_t p_offset, const char** p_fragId, uint64_t* p_nextFragStart, bool* p_biological, NGS_VDB_ErrBlock * err  )
{
    HYBRID_FUNC_ENTRY ( rcSRA, rcDatabase, rcConstructing );
    int64_t first;
    uint64_t count;
    rc_t rc = VBlobIdRange ( p_blob, &first, &count );
    if ( rc == 0  )
    {
    
        PageMapIterator pmIt;
        rc = PageMapNewIterator ( (const PageMap*)p_blob->pm, &pmIt, 0, count );
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
                    int64_t rowId = first + rowInBlob + ( p_offset - offset ) / length;
                    uint32_t bioFragNum;
                    uint32_t nextFragStart;
                    bool biological;
                    TRY ( GetFragment ( ctx, self, rowId, p_offset - offset, & bioFragNum, & nextFragStart, & biological ) )
                    {
                        if ( self -> last_frag_id != NULL )
                        {
                            NGS_StringRelease ( self -> last_frag_id, ctx );
                            self -> last_frag_id = NULL;
                        }
                        if ( biological )
                        {
                            ON_FAIL ( self -> last_frag_id = NGS_IdMakeFragment ( ctx, self -> name, false, rowId, bioFragNum ) ) 
                            {
                                break;
                            }
                            *p_fragId = NGS_StringData ( self -> last_frag_id , ctx );                
                        }
                        else
                        {
                            *p_fragId = NULL;
                        }
                        /* recalculate nextFragStart from relative to fragment to relative to blob */
                        *p_nextFragStart =  offset + nextFragStart; 
                        *p_biological = biological;
                        return;
                    } 
                    break;            
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
    
    NGS_ErrBlockThrow ( err, ctx );
}    

void 
NGS_VDB_BlobRelease ( struct VBlob* p_blob )
{
    VBlobRelease ( p_blob );    
}
