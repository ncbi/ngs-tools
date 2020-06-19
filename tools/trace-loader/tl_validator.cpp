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

#include <math.h>

#include <klib/log.h>

#include "tl_validator.hpp"
#include "tl_traceinfo.hpp"
#include "tl_tracedata.hpp"
#include "tl_names.hpp"
#include "tl_vrules.hpp"

#include <stdio.h>


using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ Weird stuff
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

#define _INVALID_ROW "INVALID_ROW"

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _TraceInfoValidator
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _TraceInfoValidator {
public:
    typedef enum {
                eInvalid = -1,
                eReady = 0,
                eLoaded,
                eValid
            } eState;
public:
    _TraceInfoValidator (
                        TL_TraceConfig * Config,
                        size_t FailureRatioPerCent = 0
                        );
    ~_TraceInfoValidator ();

    bool Validate ();

    void Export ( TL_OVec & Ret );

    void Export ( const string & Path = "" );

private:
    string _normalize_key ( const string & Key ) const;

    void _load_info ();
    void _validate_info ();

    void _validate_fields ();
    void _check_files ();
    void _misc_checks ();
    void _check_externals ();

    void _invalidate_row (
                        size_t RowIdx,
                        const string & Message
                        );
    bool _valid_row ( size_t RowIdx ) const;
    bool _check_field_is_set (
                        size_t RowIdx,
                        const string & FieldName
                        );
    bool _check_file (
                    TL_TraceData * TraceData,
                    size_t RowIdx,
                    const string & FieldName
                    );
    void _bad_file_complain (
                        size_t RowIdx,
                        const string & FieldName,
                        const string & FileName,
                        const string & Action
                        );

private:
    eState _state;

    TL_TraceConfig *  _config;

    size_t _failure_ratio;

    TL_Owp _common_fields;
    TL_OVec _trace_info;
};

_TraceInfoValidator :: _TraceInfoValidator (
                                            TL_TraceConfig * Config,
                                            size_t FailureRatioPerCent
)
:   _state ( eReady )
,   _config ( Config )
,   _failure_ratio ( FailureRatioPerCent )
,   _common_fields ()
,   _trace_info ()
{
}   /* _TraceInfoValidator :: _TraceInfoValidator () */

_TraceInfoValidator :: ~_TraceInfoValidator ()
{
    _common_fields . Clear ();
    _trace_info . clear ();
    _failure_ratio = 0;
    _config = NULL;
    _state = eReady;
}   /* _TraceInfoValidator :: ~_TraceInfoValidator () */

bool
_TraceInfoValidator :: Validate ()
{
    _load_info ();
    _validate_info ();

    return _state == eValid;
}   /* _TraceInfoValidator :: Validate () */

void
_TraceInfoValidator :: Export ( TL_OVec & Ret )
{
    Ret . empty ();

    Validate ();

    if ( _state == eValid ) {
        for ( size_t llp = 0; llp < _trace_info . size (); llp ++ ) {
             if ( _valid_row ( llp ) ) {
                Ret . insert ( Ret . end (), _trace_info [ llp ] );
            }
        }
    }
    else {
        throw TL_Exception ( "_TraceInfoValidator :: Export (): Can not export invalid trace info" );
    }
}   /* _TraceInfoValidator :: Export () */

void
_TraceInfoValidator :: Export ( const string & Path )
{
    string OutPath = Path . empty ()
                                    ? _config -> TraceInfoOwpPath ()
                                    : Path
                                    ;
    if ( OutPath . empty () ) {
        throw TL_Exception ( "_TraceInfoValidator :: Export (): Can not find appropriate file to store data" );
    }

    if ( TL_Io :: Exists ( OutPath ) ) {
        TL_Io :: Remove ( OutPath );
    }

    TL_OVec Vec;
    Export ( Vec );

        /*  First we will write object
         */
    TL_OwpVector :: Store ( OutPath, Vec );

        /*  Second we will read it and compare
         *  Of course, immediate reading after writing is not warranted
         */
    TL_OVec Vec1;
    TL_OwpVector :: Load ( OutPath, Vec1 );

    if ( Vec . size () != Vec1 . size () ) {
        throw TL_Exception ( "_TraceInfoValidator :: Export (): Something weird happened during writing" );
    }

    for ( size_t llp = 0; llp < Vec . size (); llp ++ ) {
        if ( Vec [ llp ] != Vec1 [ llp ] ) {
            throw TL_Exception ( "_TraceInfoValidator :: Export (): Something weird happened during writing: exported and readed objects are different" );
        }
    }
}   /* _TraceInfoValidator :: Export () */

string
_TraceInfoValidator :: _normalize_key ( const string & Key ) const
{
    return TL_StringU :: ToPrintable ( TL_StringU :: ToLower ( Key ) );
}   /* _TraceInfoValidator :: _normalize_key () */

    /*  Kinda long method, but who cares :D
     */
void
_TraceInfoValidator :: _load_info ()
{
        /* Attempt of loading was made already */
    if ( _state != eReady ) {
        return;
    }

        /*  Set invalid state
         */
    _state = eInvalid;

        /* Here we are reading common fields */
    _common_fields . Clear ();

    TL_LineReader Reader ( _config -> TraceInfoPath () );
    string Line;

    while ( Reader . NextLine () ) {
        Line = Reader . Line ();
        if ( Line . empty () ) {
            continue;
        }

        if ( ! TL_StringU :: IsPrint ( Line ) ) {
            throw TL_Exception ( string ( "ERROR: Line contains non printable characters at line " ) + TL_StringU :: ToStr < uint32_t > ( Reader . LineNo () ) );
        }

        string Key, Value;
        if ( ! TL_StringU :: Split ( Line, '=', Key, Value ) ) {
            break;
        }

        Key = _normalize_key ( Key );

            /*  Because of those fields are common, we do not 
             *  think much about duplicae fields, but ... report
             */
        if ( ! _common_fields . Value ( Key ) . empty () ) {
            PLOGMSG(
                    klogWarn,
                    (klogWarn,
                    "WARN: Ignoring duplicate declaration of common field '$(key)' ($(value))",
                    "severity=warn,key=%s,value=%s",
                    Key . c_str (),
                    Value . c_str ()
                    )
                    );
        }

        _common_fields . SetValue ( Key, Value );
    }

        /*  The next line, after common fields is <tab> separated
         *  fields definition line
         */
    TL_SVec Vec;
    size_t Q = TL_StringU :: Tokenize ( Line, '\t', Vec );
    if ( Q == 0 ) {
        throw TL_Exception ( "ERROR: Invalid TraceIfo file formant" );
    }

        /*  Trying to define list of fields to initialize
         */
    TL_SVec fVec;
    TL_SSet fSet;

    for ( size_t llp = 0; llp < Q; llp ++ ) {
        string Key = _normalize_key ( Vec [ llp ] );
        string rKey = Key;

            /* WARNING: There is very strange case, cuz field was
             *          in both sections : common and non-common.
             *          So, we will drop all duplicated fields
             */
        if ( fSet . find ( Key ) != fSet . end ()
            || _common_fields . Has ( Key ) ) {
            PLOGMSG(
                    klogWarn,
                    (klogWarn,
                    "WARN: Ignoring duplicate declaration of field '$(key)'",
                    "severity=warn,key=%s",
                    Key . c_str ()
                    )
                    );

            rKey = "";
        }

        fVec . insert ( fVec . end (), rKey );
        fSet . insert ( Key );
    }

        /*  Here we are reading non-common meta data
         */
    _trace_info . clear ();
    size_t RowNo = 0;

    fSet . clear ();    // will reuse to check duplicate rows

    while ( Reader . NextLine () ) {
        Line = Reader . Line ();
        if ( Line . empty () ) {
            continue;
        }

        if ( ! TL_StringU :: IsPrint ( Line ) ) {
            throw TL_Exception ( string ( "ERROR: Line contains non printable characters at line " ) + TL_StringU :: ToStr < uint32_t > ( Reader . LineNo () ) );
        }

        TL_SVec tVec;
        size_t Q = TL_StringU :: Tokenize ( Line, '\t', tVec, true );
        if ( Q == 0 ) {
            throw TL_Exception ( string ( "ERROR: Invalid TraceIfo file formant at line " ) + TL_StringU :: ToStr < uint32_t > ( Reader . LineNo () ) );
        }

            /*  Trailing \t - sometimes happens
             */
        if ( fVec . size () == ( Q - 1 )
            && TL_StringU :: EndsWith ( Line, "\t" ) ) {
            Q --;
        }

            /*  Sometimes last column value is not assigned
             */
        if ( fVec . size () == ( Q + 1 )
            && TL_StringU :: EndsWith ( Line, "\t" ) ) {
            Q ++;
            tVec . insert ( tVec . end (), "" );
        }

        if ( fVec . size () != Q ) {
            throw TL_Exception ( string ( "ERROR: Number of declared columns does not match provided number for row at line " ) + TL_StringU :: ToStr < uint32_t > ( Reader . LineNo () ) );
        }

        TL_Owp Owp = _common_fields;

        RowNo ++;
        Owp . SetName ( TL_StringU :: ToStr < uint32_t > ( RowNo ) );

        string NormString;

        for ( size_t llp = 0; llp < Q; llp ++ ) {
            string Key = fVec [ llp ];
            string Value = tVec [ llp ];
            if ( ! Key . empty () ) {
                Owp . SetValue ( Key, Value );
                NormString += _normalize_key ( Value );
            }
        }

        if ( fSet . find ( NormString ) != fSet . end () ) {
            PLOGMSG(
                    klogWarn,
                    (klogWarn,
                    "WARN: skipping duplicate row No $(row) ($(line))",
                    "severity=warn,row=%d,line=%s",
                    Reader . LineNo (),
                    Line . c_str ()
                    )
                    );
        }
        else {
            _trace_info . insert ( _trace_info . end (), Owp );
            fSet . insert ( NormString );
        }
    }

        /*  Here all fields are loaded, and we do some initial checks.
         */
    if ( _trace_info . size () == 0 ) {
        throw TL_Exception ( "FATAL: There is nothing to load, bailing out" );
    }

        /*  Reporting parsed lines
         */
    PLOGMSG(
            klogInfo,
            (klogInfo,
            "[TR_INF] parsed $(size) lines",
            "severity=info,size=%d",
            _trace_info . size ()
            )
            );

        /*  Here we do tweak TraceConfig for failure threshold
         */
    _config -> ReportModule () .SetItemCount ( _trace_info . size () );


    _state = eLoaded;
}   /* _TraceInfoValidator :: _load_info () */

void
_TraceInfoValidator :: _validate_info ()
{
        /* Attempt of validation was made already */
    if ( _state != eLoaded ) {
        return;
    }

        /*  Set invalid state
         */
    _state = eInvalid;

    PLOGMSG(
            klogInfo,
            (klogInfo,
            "[TR_INF] start trace info validation of $(size) records",
            "severity=info,size=%d",
            _trace_info . size ()
            )
            );

        /*  Here is almost original Volodya plan Validation of all
         *  necessary fields present. Crashable.
         */
    _validate_fields ();

        /*  Checking external rules
         */
    _check_externals ();

        /*  Some miscaltionous checks
         */
    _misc_checks ();

        /*  Checking files
         */
    _check_files ();

        /*  Here all fields are loaded, and we do some initial checks.
         */
    size_t ValidQty = 0;
    for ( size_t llp = 0; llp < _trace_info . size (); llp ++ ) {
         if ( _valid_row ( llp ) ) {
            ValidQty ++;
        }
    }
    if ( ValidQty == 0 ) {
        throw TL_Exception ( "FATAL: AFTER FIELD VALIDATION : there is nothing to load, bailing out" );
    }

    PLOGMSG(
            klogInfo,
            (klogInfo,
            "[TR_INF] validated $(size) rows of $(qty)",
            "severity=info,size=%d,qty=%d",
            ValidQty,
            _trace_info . size ()
            )
            );

        /*  Set valid state
         */
    _state = eValid;
}   /* _TraceInfoValidator :: _validate_info () */

void
_TraceInfoValidator :: _validate_fields ()
{
    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TR_INF] field validation ... ",
            "severity=debug"
            )
            );

    /*  Fields validation.
     *  1) we should check that all mandatory fields presents
     *  2) we should check that all fiedls declared exists
     */
        /*  Checking that all mandatory fidsa present
         */
    for (
        TL_TFMapI Pos = TL_TraceFields :: Begin();
        Pos != TL_TraceFields :: End ();
        Pos ++
    ) {
        string FieldName = Pos -> second . Name ();

            /*  Common fields check ... Not sure if we need it
             *  fall outimediately
             */
        if ( _common_fields . Has ( FieldName )
            && ! Pos -> second . IsCommon ()
        ) {
            throw TL_Exception ( string ( "ERROR: Incorrect common field: '" ) + FieldName + "'" );
        }

            /*  Required fields ... opposite to 'canBeIgnored'
             */
        if ( ! _trace_info [ 0 ] . Has ( FieldName )
            && ! Pos -> second . CanBeIgnored ()
            && Pos -> second . IsMandatory ()
        ) {
            throw TL_Exception ( string ( "ERROR: Missing required field: '" ) + FieldName + "'" );
        }
    }

        /*  Checking that all declared fields do exists for every row
         *  quite strange check, cuz if invalid fiedl presents, it
         *  presents in any row ...
         */
    for ( size_t llp = 0; llp < _trace_info . size (); llp ++ ) {
            /*  row could be already invalid
             */
        if ( ! _valid_row ( llp ) ) {
            continue;
        }

        TL_SVec Vec;
        _trace_info [ llp ] . ListKeys ( Vec );

        for (
            TL_SVec :: const_iterator It = Vec . begin ();
            It != Vec . end ();
            It ++
        ) {
                /*  Field does not exists, and row should be removed
                 */
            if ( ! TL_TraceFields :: Has ( * It ) ) {
                _invalidate_row ( llp, string ( "ERROR: Missing required field: '" ) + * It + "'" );

                break;
            }
        }
    }

    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TR_INF] PASSED ( field validation )",
            "severity=debug"
            )
            );
}   /* _TraceInfoValidator :: _validate_fields () */

void
_TraceInfoValidator :: _check_externals ()
{
    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TR_INF] perfoming external checks ... ",
            "severity=debug"
            )
            );

    TL_VRules _rules;
    _rules . Init ( _config -> ValidationRules () );

    for ( size_t llp = 0; llp < _trace_info . size (); llp ++ ) {
            /*  row could be already invalid
             */
        if ( ! _valid_row ( llp ) ) {
            continue;
        }

        string Reason;
        if ( ! _rules . Validate ( _trace_info [ llp ], Reason ) ) {
            _invalidate_row ( llp, Reason );
        }

    }

    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TR_INF] PASSED ( external checks )",
            "severity=debug"
            )
            );

}   /* _TraceInfoValidator :: _check_externals () */

void
_TraceInfoValidator :: _misc_checks ()
{
        /*  There are miscaltinuus checks which present in old 
         *  tracedb.cpp code ... will do comments
         *  Those are record by record checks
         */

    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TR_INF] perfoming misc checks ... ",
            "severity=debug"
            )
            );

    TL_SSet TraceNames;

    for ( size_t llp = 0; llp < _trace_info . size (); llp ++ ) {
            /*  row could be already invalid
             */
        if ( ! _valid_row ( llp ) ) {
            continue;
        }

            /*   checking that every record has "trace_name" defined
             *   tracedb.cpp:535,551,570
             */
        if ( ! _check_field_is_set ( llp, _TL_TRACE_NAME ) ) continue;

            /*   checking that "trace_name" is valid
             *   tracedb.cpp:425
             */
        string tName = _trace_info [ llp ] . Value ( _TL_TRACE_NAME );
        const char * InvChars = "~`!@$%^&*()+={}[]:\";'<>?,/\\|";
        if ( strpbrk ( tName . c_str (), InvChars ) != NULL ) {
            _invalidate_row ( llp, string ( "Field '" ) + _TL_TRACE_NAME + "' has invalid value '" + tName + "'" );
            continue;
        }

            /*  trace_name should be unique
             *   tracedb.cpp:493
             */
        if ( TraceNames . find ( tName ) != TraceNames . end () ) {
            _invalidate_row ( llp, string ( "Field '" ) + _TL_TRACE_NAME + "' has duplicate name '" + tName + "'" );
            continue;
        }

            /*   checking that every record has "species_code"
             *   tracedb.cpp:582
             */
        string sPec = _trace_info [ llp ] . Value ( _TL_SPECIES_CODE );
        if ( sPec . empty () ) {
            _invalidate_row ( llp, string ( "Field '" ) + _TL_SPECIES_CODE + "' is not set" );
            continue;
        }
        else {
                /* That one is at : tracedb.cpp:440
                 */
            if ( sPec == "454" ) {
                _invalidate_row ( llp, string ( "Field '" ) + _TL_SPECIES_CODE + "' can not be '454'" );
                continue;
            }
        }

            /*   checking that every record has "trace_type_code"
             *   tracedb.cpp:591
             */
        if ( ! _check_field_is_set ( llp, _TL_TRACE_TYPE_CODE ) ) continue;

            /*   checking that every record has "trace_format"
             *   tracedb.cpp:591
             */
        if ( ! _check_field_is_set ( llp, _TL_TRACE_FORMAT ) ) continue;

            /*   checking that every record has "source_type"
             *   tracedb.cpp:591
             */
        if ( ! _check_field_is_set ( llp, _TL_SOURCE_TYPE ) ) continue;

            /*  if there is base_file defined, there should be peak_file
             *  and qual_file difined.
             *  tracedb.cpp:881
             */
        string bFile = _trace_info [ llp ] . Value ( _TL_BASE_FILE );
        if ( ! bFile . empty () ) {
            if ( ! _check_field_is_set ( llp, _TL_PEAK_FILE ) ) continue;
            if ( ! _check_field_is_set ( llp, _TL_QUAL_FILE ) ) continue;
        }

            /*  That isners here, to be sure
             */
        TraceNames . insert ( tName );
    }

    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TR_INF] PASSED ( misc checks )",
            "severity=debug"
            )
            );
}   /* _TraceInfoValidator :: _misc_checks () */

void
_TraceInfoValidator :: _check_files ()
{
    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TR_INF] checking files ... ",
            "severity=debug"
            )
            );

    TL_TraceData TraceData;
    TraceData . Init ( _config );

    for ( size_t llp = 0; llp < _trace_info . size (); llp ++ ) {
            /*  row could be already invalid
             */
        if ( ! _valid_row ( llp ) ) {
            continue;
        }

        if ( ! _check_file ( & TraceData, llp, _TL_TRACE_FILE ) ) {
            continue;
        }

        if ( ! _check_file ( & TraceData, llp, _TL_BASE_FILE ) ) {
            continue;
        }

        if ( ! _check_file ( & TraceData, llp, _TL_PEAK_FILE ) ) {
            continue;
        }

        if ( ! _check_file ( & TraceData, llp, _TL_QUAL_FILE ) ) {
            continue;
        }
    }

    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[TR_INF] PASSED ( checking files )",
            "severity=debug"
            )
            );
}   /* _TraceInfoValidator :: _check_files () */

void
_TraceInfoValidator :: _invalidate_row (
                                        size_t RowIdx,
                                        const string & Message
)
{
    if ( _config -> ReportModule () . ReportErrCheckIfShouldFail () ) {
            /*  Here we should fail
             */
        PLOGMSG(
                klogErr,
                (klogErr,
                "$(reason)",
                "severity=error,reason=%s",
                Message . c_str ()
                )
                );
        throw TL_Exception ( Message );
    }

    PLOGMSG(
            klogWarn,
            (klogWarn,
            "$(message)",
            "severity=warning,message=%s",
            Message . c_str ()
            )
            );

    if ( 0 <= RowIdx && RowIdx < _trace_info . size () ) {
        _trace_info [ RowIdx ] . SetValue ( _INVALID_ROW, Message );
    }
}   /* _TraceInfoValidator :: _invalidate_row () */

bool
_TraceInfoValidator :: _valid_row ( size_t RowIdx ) const
{
    if ( 0 <= RowIdx && RowIdx < _trace_info . size () ) {
        return ! _trace_info [ RowIdx ] . Has ( _INVALID_ROW );
    }

    return false;
}   /* _TraceInfoValidator :: _valid_row () */

bool
_TraceInfoValidator :: _check_field_is_set (
                                            size_t RowIdx,
                                            const string & FieldName
)
{
    if ( 0 < RowIdx && RowIdx < _trace_info . size () ) {
        string Val = _trace_info [ RowIdx ] . Value ( FieldName );
        if ( ! Val . empty () ) {
            return true;
        }

        _invalidate_row (
                    RowIdx,
                    string ( "Field '" ) + FieldName + "' is not set"
                    );
    }

    return false;
}   /* _TraceInfoValidator :: _check_field_is_set () */

bool
_TraceInfoValidator :: _check_file (
                                    TL_TraceData * TraceData,
                                    size_t RowIdx,
                                    const string & FieldName
)
{
    string fName = _trace_info [ RowIdx ] . Value ( FieldName );
    if ( fName . empty () ) {
        return true;
    }
    PLOGMSG(
            klogDebug,
            (klogDebug,
            "Checking file [$(file_name)]",
            "severity=debug,file_name=%s",
            fName . c_str ()
            )
            );

    if ( ! TraceData -> FileExists ( fName ) ) {
        _bad_file_complain (
                        RowIdx,
                        FieldName,
                        fName,
                        "File does not exists"
                        );
        return false;
    }

    uint64_t FileSize = TraceData -> FileSize ( fName );
    if ( FileSize == 0 ) {
        _bad_file_complain (
                        RowIdx,
                        FieldName,
                        fName,
                        "File has zero length"
                        );
        return false;
    }

    return true;
}   /* _TraceInfoValidator :: _check_file () */

void
_TraceInfoValidator :: _bad_file_complain (
                                        size_t RowIdx,
                                        const string & FieldName,
                                        const string & FileName,
                                        const string & Action
)
{
    stringstream Str;

    if ( FileName . empty () ) {
        Str << Action << " [" << FieldName << "]";
    }
    else {
        Str << Action << " [" << FieldName << "] associated with [" << FileName << "]";
    }

    _invalidate_row ( RowIdx, Str . str () );
}   /* _TraceInfoValidator :: _bad_file_complain () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceInfoValidator
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_TraceInfoValidator :: TL_TraceInfoValidator ()
:   _validator ( NULL )
{
}   /* TL_TraceInfoValidator :: TL_TraceInfoValidator () */

TL_TraceInfoValidator :: ~TL_TraceInfoValidator ()
{
    Dispose ();
}   /* TL_TraceInfoValidator :: ~TL_TraceInfoValidator () */

void
TL_TraceInfoValidator :: Init (
                            TL_TraceConfig * Config,
                            size_t FailureRatioPerCent
)
{
    Dispose ();

    _TraceInfoValidator * T = new _TraceInfoValidator (
                                                    Config,
                                                    FailureRatioPerCent
                                                    );

    _validator = T;
}   /* TL_TraceInfoValidator :: Init () */

void
TL_TraceInfoValidator :: Dispose ()
{
    if ( _validator != NULL ) {
        _TraceInfoValidator * T = ( _TraceInfoValidator * ) _validator;

            /* Once is enough */
        _validator = NULL;

        delete T;
    }
}   /* TL_TraceInfoValidator :: Dispose () */

bool
TL_TraceInfoValidator :: Validate ()
{
    if ( _validator == NULL ) {
        throw TL_Exception ( "TL_TraceInfoValidator :: Validate (): uninitialized validator object" );
    }

    return ( ( _TraceInfoValidator * ) _validator ) -> Validate ();
}   /* TL_TraceInfoValidator :: Validate () */

void
TL_TraceInfoValidator :: Export ( TL_OVec & Ret )
{
    if ( _validator == NULL ) {
        throw TL_Exception ( "TL_TraceInfoValidator :: Export (): uninitialized validator object" );
    }

    return ( ( _TraceInfoValidator * ) _validator ) -> Export ( Ret );
}   /* TL_TraceInfoValidator :: Export () */

void
TL_TraceInfoValidator :: Export ( const string & Path )
{
    if ( _validator == NULL ) {
        throw TL_Exception ( "TL_TraceInfoValidator :: Export (): uninitialized validator object" );
    }

    return ( ( _TraceInfoValidator * ) _validator ) -> Export ( Path );
}   /* TL_TraceInfoValidator :: Export () */
