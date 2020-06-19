/*======================================================================
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
* ====================================================================/======
*
*/

#include <sstream>
#include <iostream>

#include <klib/log.h>

#include "tl_traceinfo.hpp"
#include "tl_names.hpp"


using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceInfo
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_TraceInfo :: TL_TraceInfo ()
{
}   /* TL_TraceInfo :: TL_TraceInfo () */

TL_TraceInfo :: ~TL_TraceInfo ()
{
}   /* TL_TraceInfo :: ~TL_TraceInfo () */

size_t
TL_TraceInfo :: FieldQty () const
{
    throw TL_Exception ( "Unimplemented methdo TL_TraceInfo :: FieldQty () called" );
}   /* TL_TraceInfo :: FieldQty () */

const TL_TraceField & 
TL_TraceInfo :: Field ( size_t FieldIndex ) const
{
    throw TL_Exception ( "Unimplemented methdo TL_TraceInfo :: Field () called" );
}   /* TL_TraceInfo :: Field () */

size_t
TL_TraceInfo :: RowQty () const
{
    throw TL_Exception ( "Unimplemented methdo TL_TraceInfo :: RowQty () called" );
}   /* TL_TraceInfo :: RowQty () */

const TL_SVec & 
TL_TraceInfo :: Row ( size_t RowIndex ) const
{
    throw TL_Exception ( "Unimplemented methdo TL_TraceInfo :: Row () called" );
}   /* TL_TraceInfo :: Row () */

bool
TL_TraceInfo :: HasField ( const string & FieldName) const
{
    throw TL_Exception ( "Unimplemented methdo TL_TraceInfo :: HasField () called" );
}   /* TL_TraceInfo :: HasField () */

const string &
TL_TraceInfo :: Value ( size_t RowIndex, const string & FieldName) const
{
    throw TL_Exception ( "Unimplemented methdo TL_TraceInfo :: Value () called" );
}   /* TL_TraceInfo :: Value () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_SimpleTraceInfo
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_SimpleTraceInfo :: TL_SimpleTraceInfo ()
:   TL_TraceInfo ()
,   _fields ()
,   _field_index ()
,   _rows ()
{
}   /* TL_SimpleTraceInfo :: TL_SimpleTraceInfo () */

TL_SimpleTraceInfo :: ~TL_SimpleTraceInfo ()
{
    Reset ();
}   /* TL_SimpleTraceInfo :: ~TL_SimpleTraceInfo () */

void
TL_SimpleTraceInfo :: Reset ()
{
    _fields . clear ();
    _field_index . clear ();
    _rows . clear ();
}   /* TL_SimpleTraceInfo :: Reset () */

void
TL_SimpleTraceInfo :: AddField ( const string & FieldName )
{
    const TL_TraceField & Fld = TL_TraceFields :: Find ( FieldName );
    if ( Fld . IsEmpty () ) {
        throw TL_Exception ( string ( "Can not find definition for a field \"" ) + FieldName + "\"" );
    }

    _field_index [ FieldName ] = _fields . size ();
    _fields . insert ( _fields . end (), Fld );
}   /* TL_SimpleTraceInfo :: Reset () */

size_t
TL_SimpleTraceInfo :: FieldQty () const
{
    return _fields . size ();
}   /* TL_SimpleTraceInfo :: Reset () */

const TL_TraceField &
TL_SimpleTraceInfo :: Field ( size_t FieldIndex ) const
{
    if ( _fields . size () <= FieldIndex ) {
        throw TL_Exception ( "TL_SimpleTraceInfo :: Field () : invalid field index" );
    }

    return _fields [ FieldIndex ]; 
}   /* TL_SimpleTraceInfo :: Reset () */

void
TL_SimpleTraceInfo :: AddRow ( const TL_SVec & Row )
{
    if ( _fields . size () == 0 ) {
        throw TL_Exception ( "TL_SimpleTraceInfo :: AddRow () : set fields first" );
    }

    if ( _fields . size () != Row . size () ) {
        throw TL_Exception ( "TL_SimpleTraceInfo :: AddRow () : invalid row size" );
    }

    _rows . insert ( _rows . end (), Row );
}   /* TL_SimpleTraceInfo :: Reset () */

size_t
TL_SimpleTraceInfo :: RowQty () const
{
    return _rows . size ();
}   /* TL_SimpleTraceInfo :: Reset () */

const TL_SVec &
TL_SimpleTraceInfo :: Row ( size_t RowIndex ) const
{
    if ( _rows . size () <= RowIndex ) {
        throw TL_Exception ( "TL_SimpleTraceInfo :: Row () : invalid row index" );
    }

    return _rows [ RowIndex ];
}   /* TL_SimpleTraceInfo :: Reset () */

bool
TL_SimpleTraceInfo :: HasField ( const string & FieldName ) const
{
    return _field_index . find ( FieldName ) != _field_index . end ();
}   /* TL_SimpleTraceInfo :: HasField () */

const string &
TL_SimpleTraceInfo :: Value (
                            size_t RowIndex,
                            const string & FieldName
) const
{
    TL_SIMap :: const_iterator It = _field_index . find ( FieldName );
    if ( It == _field_index . end () ) {
        // Agree, it was wrong idea :D
        // throw TL_Exception ( string ( "TL_SimpleTraceInfo :: value () : invalid filed name '" ) + FieldName + "'" );

        return TL_StringU :: EmptyString ();
    }

    return Row ( RowIndex ) [ It -> second ];
}   /* TL_SimpleTraceInfo :: Value () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_SimpleTraceInfoLoader
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_SimpleTraceInfoLoader :: TL_SimpleTraceInfoLoader ()
{
}   /* TL_SimpleTraceInfoLoader :: TL_SimpleTraceInfoLoader () */

TL_SimpleTraceInfoLoader :: ~TL_SimpleTraceInfoLoader ()
{
}   /* TL_SimpleTraceInfoLoader :: ~TL_SimpleTraceInfoLoader () */

void
TL_SimpleTraceInfoLoader :: Load (
                            TL_SimpleTraceInfo & TraceInfo,
                            const string & Path
)
{
    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TRACE INFO] [LOADING] for [$(file)]",
            "severity=debug,file=$s",
            Path . c_str ()
            )
            );

    TL_OVec Vec;
    TL_OwpVector :: Load ( Path, Vec );

    Load ( TraceInfo, Vec );
}   /* TL_SimpleTraceInfoLoader :: Load () */

void
TL_SimpleTraceInfoLoader :: Load (
                            TL_SimpleTraceInfo & TraceInfo,
                            const TL_OVec & Vec
)
{
    _read_fields ( TraceInfo, Vec );
    _make_rows ( TraceInfo, Vec );
}   /* TL_SimpleTraceInfoLoader :: Load () */

void
TL_SimpleTraceInfoLoader :: _read_fields (
                                TL_SimpleTraceInfo & TraceInfo,
                                const TL_OVec & Vec
)
{
    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TRACE INFO] [READING] fields",
            "severity=debug"
            )
            );
        /* This is a list of all fields. Sorry using brutal forse
         * to make that set
         */

    TL_CISSet Keys;
    for (
        TL_OVec :: const_iterator It = Vec . begin ();
        It != Vec . end ();
        It ++
    ) {
        TL_SVec V;

        if ( It -> ListKeys ( V ) ) {
            for (
                TL_SVec :: const_iterator Ti = V . begin ();
                Ti != V . end ();
                Ti ++
            ) {
                Keys . insert ( * Ti );
            }
        }
    }

        /* Getting all fields ready
         */
    TraceInfo . Reset ();
    for (
        TL_SSet :: const_iterator It = Keys . begin ();
        It != Keys . end ();
        It ++
    ) {
        TraceInfo . AddField ( * It );
    }

    PLOGMSG(
            klogInfo,
            (klogInfo,
            "[TRACE INFO] [$(records)] records [$(fields)] fields",
            "severity=info,records=%d,fields=%d",
            Vec . size (),
            TraceInfo . FieldQty ()
            )
            );
}   /* TL_SimpleTraceInfoLoader :: _read_fields () */

void
TL_SimpleTraceInfoLoader :: _make_rows (
                                TL_SimpleTraceInfo & TraceInfo,
                                const TL_OVec & Vec
)
{
    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TRACE INFO] [PREPARING METADATA]",
            "severity=debug"
            )
            );

    for (
        TL_OVec :: const_iterator It = Vec . begin ();
        It != Vec . end ();
        It ++
    ) {
        TL_SVec Row;
        for ( size_t llp = 0; llp < TraceInfo . FieldQty (); llp ++ ) {

            string Key = TraceInfo . Field ( llp ) . Name ();

            Row . insert ( Row . end (), It -> Value ( Key ) );
        }

        TraceInfo . AddRow ( Row );
    }

}   /* TL_SimpleTraceInfoLoader :: _make_rows () */

