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

#include <iostream>

#include <unistd.h>
#include <math.h>

#include "tl_util.hpp"
#include "tl_log.hpp"
#include "tl_traceconfig.hpp"


using namespace std;
using namespace _tl_;

#define _DEFAULT_ERR_COUNT   100
#define _DEFAULT_ERR_PERCENT 5

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceConfig
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_TraceConfig :: TL_TraceConfig ()
:   _trace_info_path ( "" )
,   _trace_info_owp_path ( "" )
,   _trace_schema ( "" )
,   _trace_output ( "" )
,   _validation_rules ( "" )
{
    _report . SetContext ( _trace_info_path );
    _report . SetMeanings ( "traces" );
}   /* TL_TraceConfig :: TL_TraceConfig () */

TL_TraceConfig :: ~TL_TraceConfig ()
{
    TL_TRY {
        Reset ();
    }
    TL_CATCH_R
}   /* TL_TraceConfig :: ~TL_TraceConfig () */

void
TL_TraceConfig :: Reset ()
{
    _report . SetContext ();
    _report . SetMeanings ();

    _validation_rules . clear ();
    _trace_output . clear ();
    _trace_schema . clear ();
    _trace_info_owp_path . clear ();
    _trace_info_path . clear ();
}   /* TL_TraceConfig :: Reset () */

static
string
_replaceSlashes ( const string & Path )
{
        /* Quite lazy, but works ...
         */
    string RetVal = Path;
    for ( uint32_t llp = 0; llp < RetVal . length (); llp ++ ) {
        if ( RetVal [ llp ] == '\\' ) {
            RetVal [ llp ] = '/';
        }
    }
    return RetVal;
}   /* _replaceSlashes () */

static
string
_normalizePath ( const string & Path )
{
    string Ret = _replaceSlashes ( Path );

        /* Trimming trailing slashes
         */
    while ( ! Ret . empty () && Ret [ Ret . length () - 1 ] == '/' ) {
        Ret = Ret . substr ( 0, Ret . length () - 1 );
    }

    return Ret;
}   /* _normalizePath () */

static
string
_dirNameAsShellDoIt ( const string & Path )
{
    string Ret = _normalizePath ( Path );

        /* I sez like Shell do it 
         */
    if ( Ret . empty () || Ret == "." || Ret == ".." ) {
        return ".";
    }

    size_t Pos = Ret . rfind ( '/' );
    if ( Pos != string :: npos ) {
        Ret = Ret . substr ( 0, Pos );
    }

    return Ret;
}   /* _dirNameAsShellDoIt () */

string
TL_TraceConfig :: NormalizePath ( const string & Path ) const
{
    if ( Path . empty () ) {
        return Path;
    }

    string Ret = _replaceSlashes ( Path );

    if ( Ret [ 0 ] != '/' ) {
        Ret = _dirNameAsShellDoIt ( _trace_info_path )
                                    + "/" + _replaceSlashes ( Path );
    }

    return Ret;
}   /* TL_TraceConfig :: NormalizePath () */

void
TL_TraceConfig :: SetTraceInfoPath ( const std::string & InfoPath )
{
    _trace_info_path . clear ();
    _trace_info_owp_path . clear ();

    if ( InfoPath . empty () ) {
        throw TL_Exception ( "TraceConfig: empty TraceInfoPath passed" );
    }

    string Path = _normalizePath ( InfoPath );

    if ( ! TL_Io :: Exists ( Path ) ) {
        throw TL_Exception ( string ( "Can not stat TraceInfoPath \"" ) + InfoPath + "\"" );
    }

    _trace_info_path = Path;
    _trace_info_owp_path = _trace_info_path + ".owp";

    _report . SetContext ( _trace_info_path );
}   /* TL_TraceConfig :: SetTraceInfoPath () */

void
TL_TraceConfig :: SetTraceSchemaPath ( const std::string & SchemaPath )
{
    _trace_schema . clear ();

    if ( SchemaPath . empty () ) {
        throw TL_Exception ( "TraceConfig: empty TraceSchemaPath passed" );
    }

    string Path = _normalizePath ( SchemaPath );

    if ( ! TL_Io :: IsFile ( Path ) ) {
        throw TL_Exception ( string ( "Can not stat TraceSchemaPath \"" ) + SchemaPath + "\"" );
    }

    _trace_schema = Path;
}   /* TL_TraceConfig :: SetTraceSchemaPath () */

void
TL_TraceConfig :: SetTraceOutputPath ( const std::string & OutputPath )
{
    _trace_output = _normalizePath ( OutputPath );
}   /* TL_TraceConfig :: SetTraceOutputPath () */

void
TL_TraceConfig :: SetTraceValidationRulesPath (
                                const std::string & ValidationRulesPath
)
{
    _validation_rules . clear ();

    if ( ValidationRulesPath . empty () ) {
        throw TL_Exception ( "TraceConfig: empty TraceValidationRules passed" );
    }

    string Path = _normalizePath ( ValidationRulesPath );

        /* Ha-ha-ha ... stupid check */
    if ( TL_Io :: IsDir ( Path ) ) {
        throw TL_Exception ( string ( "Can not stat TraceValidationRules \"" ) + ValidationRulesPath + "\"" );
    }

    _validation_rules = Path;
}   /* TL_TraceConfig :: SetTraceValidationRulesPath () */

void
TL_TraceConfig :: SaySomething () const
{
#ifdef JOJOBA
    cout << "TRACE INFO [" << _trace_info_path << "]" << endl;
    cout << " TRACE OWP [" << _trace_info_owp_path << "]" << endl;
    cout << "    SCHEMA [" << _trace_schema << "]" << endl;
    if ( _trace_output . empty () ) {
        cout << "       OUT [STDOUT]" << endl;
    }
    else {
        cout << "       OUT [" << _trace_output << "]" << endl;
    }
    cout << "     RULEZ [" << _validation_rules << "]" << endl;
    cout << "MX ERR CNT [" << _max_err_count << "]" << endl;
    cout << "MX ERR PCT [" << _max_err_percent << "]" << endl;
#endif /* JOJOBA */
}   /* TL_TraceConfig :: SaySomething () */
