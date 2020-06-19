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

#include <string>
#include <iostream>
#include <sstream>

#include "tl_tracefields.hpp"
#include "tl_tracefields_init.hpp"

#include "tl_util.hpp"
#include "tl_log.hpp"


using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceField
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_TraceField :: TL_TraceField ()
{
    _reset ();
}   /* TL_TraceField :: TL_TraceField () */

TL_TraceField :: TL_TraceField ( const TL_TraceField & Field )
{
    _reset ();
    _set ( Field . _owp );
}   /* TL_TraceField :: TL_TraceField () */

TL_TraceField :: ~TL_TraceField ()
{
    _reset ();
}   /* TL_TraceField :: ~TL_TraceField () */

TL_TraceField &
TL_TraceField :: operator = ( const TL_TraceField & Field )
{
    if ( this != & Field ) {
        * this = Field . _owp;
    }
    return * this;
}   /* TL_TraceField :: TL_TraceField () */

TL_TraceField &
TL_TraceField :: operator = ( const TL_Owp & Desc )
{
    _set ( Desc );

    return * this;
}   /* TL_TraceField :: TL_TraceField () */

void
TL_TraceField :: _reset ()
{
    _name = "";
    _mandatory = false;
    _log_pos = false;
    _description = "";
    _common = false;
    _can_be_ignored = false;
    _can_be_synonym = false;
    _rfc = false;
    _table_name = "";
    _field_prim = "";
    _field_norm = "";
    _static = false;
    _searchable = false;
    _for_dump = false;
    _field_pos = 0;
    _vdb_name = "";
    _vdb_size = sizeof ( char );
    _vdb_type = "ascii";
    _vdb_type_cast = "";
    _vdb_compression = "";
    _vdb_default_value = "";
    _comment . clear ();
    _deprecated = false;

    _owp . Clear ();
}   /* TL_TraceField :: _reset () */

#define TL_CHECK_SET_BOOL(Fld,Name,Owp) { string Val = Owp . Value ( Name##_TAG ); Fld = Val . empty () ? Name##_DEF : ( Val == TL_BOOL_TRUE ); }

#define TL_CHECK_SET_STRING(Fld,Name,Owp) { string Val = Owp . Value ( Name##_TAG ); Fld = Val . empty () ? Name##_DEF : Val; }

#define TL_CHECK_SET_INT(Fld,Name,Owp) { string Val = Owp . Value ( Name##_TAG ); Fld = Val . empty () ? Name##_DEF : atoi ( Val . c_str () ); }


void
TL_TraceField :: _set ( const TL_Owp & Owp )
{
    _owp . Clear ();

    if ( Owp . IsEmpty () ) {
        _reset ();
    }
    else {
        _owp = Owp;

        _name = _owp . Name ();

        string Value;

        TL_CHECK_SET_BOOL ( _mandatory, TL_MANDATORY, _owp );
        TL_CHECK_SET_BOOL ( _log_pos, TL_LOG_POS, _owp );
        TL_CHECK_SET_STRING ( _description, TL_DESCRIPTION, _owp );
        TL_CHECK_SET_BOOL ( _common, TL_COMMON, _owp );
        TL_CHECK_SET_BOOL ( _can_be_ignored, TL_CB_IGNORED, _owp );
        TL_CHECK_SET_BOOL ( _can_be_synonym, TL_CB_SYNONYM, _owp );
        TL_CHECK_SET_BOOL ( _rfc, TL_RFC, _owp );
        TL_CHECK_SET_STRING ( _table_name, TL_TABLE_NAME, _owp );
        TL_CHECK_SET_STRING ( _field_prim, TL_FIELD_PRIM, _owp );
        TL_CHECK_SET_STRING ( _field_norm, TL_FIELD_NORM, _owp );
        TL_CHECK_SET_BOOL ( _static, TL_STATIC, _owp );
        TL_CHECK_SET_BOOL ( _searchable, TL_SEARCHABLE, _owp );
        TL_CHECK_SET_BOOL ( _for_dump, TL_FOR_DUMP, _owp );
        TL_CHECK_SET_INT ( _field_pos, TL_FIELD_POS, _owp );
        TL_CHECK_SET_STRING ( _vdb_name, TL_VDB_NAME, _owp );
        TL_CHECK_SET_INT ( _vdb_size, TL_VDB_SIZE, _owp );
        TL_CHECK_SET_STRING ( _vdb_type, TL_VDB_TYPE, _owp );
        TL_CHECK_SET_STRING ( _vdb_type_cast, TL_VDB_TYPE_CAST, _owp );
        TL_CHECK_SET_STRING ( _vdb_compression, TL_VDB_COMPRESSION, _owp );
        TL_CHECK_SET_STRING ( _vdb_default_value, TL_VDB_DEF_VALUE, _owp );
        TL_CHECK_SET_STRING ( _comment, TL_COMMENT, _owp );
        TL_CHECK_SET_BOOL ( _deprecated, TL_DEPRECATED, _owp );

        _validate ();
    }
}   /* TL_TraceField :: _set () */

void
TL_TraceField :: _validate () const
{
    /* JOJOBA */
}   /* TL_TraceField :: _validate () */

#define TL_READ_SET_BOOL(Value,Name,Owp) if ( Value != Name##_DEF ) Owp . SetValue ( Name##_TAG, ( Value ? TL_BOOL_TRUE : TL_BOOL_FALSE ) );

#define TL_READ_SET_STRING(Value,Name,Owp) if ( Value != Name##_DEF ) Owp . SetValue ( Name##_TAG, Value );

#define TL_READ_SET_INT(Value,Name,Owp) if ( Value != Name##_DEF ) { stringstream Str; Str << Value; Owp . SetValue ( Name##_TAG, Str . str () ); }

TL_Owp
TL_TraceField :: Compile () const
{
    TL_Owp Owp ( Name () );

    TL_READ_SET_BOOL ( _mandatory, TL_MANDATORY, Owp );
    TL_READ_SET_BOOL ( _log_pos, TL_LOG_POS, Owp );
    TL_READ_SET_STRING ( _description, TL_DESCRIPTION, Owp );
    TL_READ_SET_BOOL ( _common, TL_COMMON, Owp );
    TL_READ_SET_BOOL ( _can_be_ignored, TL_CB_IGNORED, Owp );
    TL_READ_SET_BOOL ( _can_be_synonym, TL_CB_SYNONYM, Owp );
    TL_READ_SET_BOOL ( _rfc, TL_RFC, Owp );
    TL_READ_SET_STRING ( _table_name, TL_TABLE_NAME, Owp );
    TL_READ_SET_STRING ( _field_prim, TL_FIELD_PRIM, Owp );
    TL_READ_SET_STRING ( _field_norm, TL_FIELD_NORM, Owp );
    TL_READ_SET_BOOL ( _static, TL_STATIC, Owp );
    TL_READ_SET_BOOL ( _searchable, TL_SEARCHABLE, Owp );
    TL_READ_SET_BOOL ( _for_dump, TL_FOR_DUMP, Owp );
    TL_READ_SET_INT ( _field_pos, TL_FIELD_POS, Owp );
    TL_READ_SET_STRING ( _vdb_name, TL_VDB_NAME, Owp );
    TL_READ_SET_INT ( _vdb_size, TL_VDB_SIZE, Owp );
    TL_READ_SET_STRING ( _vdb_type, TL_VDB_TYPE, Owp );
    TL_READ_SET_STRING ( _vdb_type_cast, TL_VDB_TYPE_CAST, Owp );
    TL_READ_SET_STRING ( _vdb_compression, TL_VDB_COMPRESSION, Owp );
    TL_READ_SET_STRING ( _vdb_default_value, TL_VDB_DEF_VALUE, Owp );
    TL_READ_SET_STRING ( _comment, TL_COMMENT, Owp );
    TL_READ_SET_BOOL ( _deprecated, TL_DEPRECATED, Owp );


    return Owp;
}   /* TL_TraceField :: Compile () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceFieldsRep
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class TL_TraceFieldsRep {
public:
    TL_TraceFieldsRep ();
    ~TL_TraceFieldsRep ();

    TL_TFMapI Begin ();
    TL_TFMapI End ();

    bool Has ( const string & FieldName );
    bool Has ( const TL_TraceField & Field );

    const TL_TraceField & Find ( const string & Name );

private:
    void _InitMap ();
    void _DisposeMap ();

    TL_TFMap _map;
    bool _map_inited;
};

static const TL_TraceField _sEmptyTraceField;

static TL_TraceFieldsRep _sTL_TraceFieldsRep;

TL_TraceFieldsRep :: TL_TraceFieldsRep ()
:   _map ()
,   _map_inited ( false )
{
}   /* TL_TraceFieldsRep :: TL_TraceFieldsRep () */

TL_TraceFieldsRep :: ~TL_TraceFieldsRep ()
{
    TL_TRY {
        _DisposeMap ();
    }
    TL_CATCH_R
}   /* TL_TraceFieldsRep :: ~TL_TraceFieldsRep () */

void
TL_TraceFieldsRep :: _InitMap ()
{
    if ( ! _map_inited ) {
        TL_OVec Descriptions;
        TL_GetFieldDescriptions ( Descriptions );

        for (
            TL_OVec :: const_iterator It = Descriptions . begin ();
            It != Descriptions . end ();
            It ++
        ) {
/* JOJOBA :: think about Name of field name
 */
            _map [ It -> Name () ] = * It;
        }

        _map_inited = true;
    }

}   /* TL_TraceFieldsRep :: _InitMap () */

void
TL_TraceFieldsRep :: _DisposeMap ()
{
    _map_inited = false;
    _map . clear ();
}   /* TL_TraceFieldsRep :: _DisposeMap () */

bool
TL_TraceFieldsRep :: Has ( const string & Name )
{
    _InitMap ();

    return _map . find ( Name ) != _map . end ();
}   /* TL_TraceFieldsRep :: Has () */

bool
TL_TraceFieldsRep :: Has ( const TL_TraceField & Field )
{
    _InitMap ();

    return Has ( Field . Name () );
}   /* TL_TraceFieldsRep :: Has () */

const TL_TraceField &
TL_TraceFieldsRep :: Find ( const string & Name )
{
    _InitMap ();

    TL_TFMap :: const_iterator Pos = _map . find ( Name );
    if ( Pos != _map . end () ) {
        return Pos -> second;
    }

    return _sEmptyTraceField;
}   /* TL_TraceFieldsRep :: Find () */

TL_TFMapI
TL_TraceFieldsRep :: Begin ()
{
    _InitMap ();

    return _map . begin ();
}   /* TL_TraceFieldsRep :: Begin () */

TL_TFMapI
TL_TraceFieldsRep :: End ()
{
    _InitMap ();

    return _map . end ();
}   /* TL_TraceFieldsRep :: End () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceFields
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_TraceFields :: TL_TraceFields ()
{
}   /* TL_TraceFields :: TL_TraceFields () */

TL_TFMapI
TL_TraceFields :: Begin ()
{
    return _sTL_TraceFieldsRep . Begin ();
}   /* TL_TraceFields :: Begin () */

TL_TFMapI
TL_TraceFields :: End ()
{
    return _sTL_TraceFieldsRep . End ();
}   /* TL_TraceFields :: End () */

bool
TL_TraceFields :: Has ( const string & FieldName )
{
    return _sTL_TraceFieldsRep . Has ( FieldName );
}   /* TL_TraceFields :: Has () */

bool
TL_TraceFields :: Has ( const TL_TraceField & Field )
{
    return _sTL_TraceFieldsRep . Has ( Field );
}   /* TL_TraceFields :: Has () */

const TL_TraceField &
TL_TraceFields :: Find ( const string & Name )
{
    return _sTL_TraceFieldsRep . Find ( Name );
}   /* TL_TraceFields :: Find () */

static
void
_dumpHead ()
{
    cout << "//" << endl;
    cout << "// This file contains description of Trace database fields" << endl;
    cout << "// Format is very simple:" << endl;
    cout << "//" << endl;
    cout << "// FIELD = CANONIC_FIELD_NAME" << endl;
    cout << "// PROPERTY = PROPERTY_VALUE" << endl;
    cout << "// ..." << endl;
    cout << "// COMMENT = Comment line one" << endl;
    cout << "// COMMENT = Comment line two" << endl;
    cout << "// COMMENT = Comment line three" << endl;
    cout << "// <empty line>" << endl;
    cout << "//" << endl;
    cout << "// Directives always starts with '#' character" << endl;
    cout << "// There is only one directive for now:" << endl;
    cout << "//" << endl;
    cout << "// #include filename" << endl;
    cout << "//" << endl;
    cout << "// That directive will load additional description fields from file" << endl;
    cout << "//" << endl;
    cout << "// If line starts with double slash '//', it considered as comment" << endl;
    cout << "//" << endl;
    cout << endl;
}   /* _dumpHead () */

static
void
_dumpComments ()
{
    cout << endl;
    cout << "field = $COMMENTS$" << endl;
    cout << "field_comment = column name in the RFC" << endl;
    cout << "name_comment = column name in VDB" << endl;
    cout << "mandatory_comment = mandatory field or not (former 'presence')" << endl;
    cout << "log_pos_comment = if to log, which position to put into" << endl;
    cout << "description_comment = meaning of the field" << endl;
    cout << "common_comment = can be specified under common section" << endl;
    cout << "can_be_ignored_comment = can be ignored during parsing" << endl;
    cout << "rfc_comment = field listed in RFC document" << endl;
    cout << "table_name_comment = table wich normalizes it" << endl;
    cout << "field_prim_comment = primary field name in Trace table" << endl;
    cout << "field_norm_comment = normailzed field (field with the actual value)" << endl;
    cout << "static_comment = what vocabulary to use ? static" << endl;
    cout << "searchable_comment = allow to search by this field" << endl;
    cout << endl;
}   /* _dumpComments () */

static
void
_dumpField ( const TL_TraceField & Field )
{
    cout << "field = " << TL_StringU :: ToUpper ( Field . Name () ) << endl;

    cout << "// original Trace fields" << endl;
    cout << "name = " << Field . VdbName () << endl;
    cout << "mandatory = " << ( Field . IsMandatory () ? "true" : "false" ) << endl;
    cout << "log_pos = " << Field . LogPos () << endl;
    cout << "description = " << Field . Description () << endl;
    cout << "common = " << ( Field . IsCommon () ? "true" : "false" ) << endl;
    cout << "can_be_ignored = " << ( Field . CanBeIgnored () ? "true" : "false" ) << endl;
    cout << "can_be_synonym = " << ( Field . CanBeSynonym () ? "true" : "false" ) << endl;
    cout << "rfc = " << ( Field . IsRfc () ? "true" : "false" ) << endl;

    cout << "// obsolete original Trace fields" << endl;
    cout << "table_name = " << Field . TableName () << endl;
    cout << "field_prim = " << Field . FieldPrim () << endl;
    cout << "field_norm = " << Field . FieldNorm () << endl;
    cout << "static = " << ( Field . IsStatic () ? "true" : "false" ) << endl;
    cout << "searchable = " << ( Field . IsSearchable () ? "true" : "false" ) << endl;
    cout << "for_dump = " << ( Field . IsForDump () ? "true" : "false" ) << endl;
    cout << "field_pos = " << Field . FieldPos () << endl;

    cout << "// VDB realted fields" << endl;
    cout << "vdb_name = " << Field . VdbName () << endl;
    cout << "vdb_size = " << Field . VdbSize () << endl;
    cout << "vdb_type = " << Field . VdbType () << endl;
    cout << "vdb_type_cast = " << Field . VdbTypeCast () << endl;
    cout << "vdb_compression = " << Field . VdbCompression () << endl;
    cout << "vdb_defauld_value = " << Field . VdbDefaultValue () << endl;
    cout << "comments = " << Field . Comment () << endl;

    cout << "// accessibility" << endl;
    cout << "deprecated = " << ( Field . IsDeprecated () ? "true" : "false" ) << endl;

    cout << endl;
}   /* _dumpField () */

void
TL_TraceFields :: Dump ()
{
    _dumpHead ();

    _dumpComments ();

    for ( TL_TFMapI Pos = Begin (); Pos != End (); Pos ++ ) {
        _dumpField ( Pos -> second );
    }
}   /* TL_TraceFields :: Dump () */
