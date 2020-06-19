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

#include <exception>
#include <string>

#include <kapp/main.h>
#include <kapp/args.h>
#include <kapp/log-xml.h>

#include <klib/out.h>
#include <klib/rc.h>

#include "tl_exception.hpp"
#include "tl_util.hpp"
#include "tl_owp.hpp"
#include "tl_traceinfo.hpp"
#include "tl_traceconfig.hpp"
#include "tl_trace_dp.hpp"
#include "tl_trace_fwa.hpp"
#include "tl_tracedata.hpp"
#include "tl_tracefields_init.hpp"
#include "tl_validator.hpp"
#include "tl_log.hpp"

#include <general-writer.hpp>

#include <stdio.h>

using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_* 
 *
 * Lyrics
 *
 * That program accepts these arguments:
 *
 *     --trace-info <path>
 *     --schema <path>
 *     --val-rules <path>
 *     --output <path>
 *     -E|--max-err-count <number>
 *     --max-err-percent <number>
 *
 * As many other standard options ...
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_* 
 *
 * Arguementes !?
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

/*  Output Path: could be real path or it will pritn to stdout
 */
static const char * _OutputPathName = "output";
static const char * _OutputPathAlias = NULL;
static const char * _OutputPathUsage [] = {
    "Output path for General Writer. Optional. If none: prints to STDOUT. Please don't make print him to STDOUT, unless it is piped to GeneralLoader", NULL
};

/*  TraceInfo Path: mandatory
 */
static const char * _TraceInfoName = "trace-info";
static const char * _TraceInfoAlias = NULL;
static const char * _TraceInfoUsage [] = {
    "Path to trace info file. Mandatory", NULL
};

/*  TraceSchema Path: mandatory
 */
static const char * _TraceSchemaName = "schema";
static const char * _TraceSchemaAlias = NULL;
static const char * _TraceSchemaUsage [] = {
    "Schema file to use. Mandatory.", NULL
};

/*  Path to a file with validation rules : mandatory
 */
static const char * _ValRulesName = "val-rules";
static const char * _ValRulesAlias = NULL;
static const char * _ValRulesUsage [] = {
    "Path to a file with validation rules. Optional.", NULL
};

/*  Max errors count : optional
 */
static const char * _ErrCountName = "max-err-count";
static const char * _ErrCountAlias = "E";
static const char * _ErrCountUsage [] = {
    "Set maximum amount of errors from TRACE submission to ignore. Optional. Default value 100. If not set, suppressed by \"max-err-percent\".", NULL
};

/*  Max errors percent, if set - suppresses error count : optional
 */
static const char * _ErrPercentName = "max-err-percent";
static const char * _ErrPercentAlias = NULL;
static const char * _ErrPercentUsage [] = {
    "Set maximum percent of errors from amount of TRACE submission to ignore, but not more than 100. Optional. Default value is 5%. Suppresses \"max-err-count\".", NULL
};


OptDef _TraceLoaderOptions [] =
{
        /* OutputPath */
    {   _OutputPathName,     /* long name */
        _OutputPathAlias,    /* single character names */
        NULL,               /* function which generate help string */
        _OutputPathUsage,    /* help string(s) */
        1,                  /* max allowed values QTY, 0 - unlimited */
        true,               /* option requires argument */
        false                /* option is mandatory */
    },
        /* TraceInfo */
    {   _TraceInfoName,     /* long name */
        _TraceInfoAlias,    /* single character names */
        NULL,               /* function which generate help string */
        _TraceInfoUsage,    /* help string(s) */
        1,                  /* max allowed values QTY, 0 - unlimited */
        true,               /* option requires argument */
        true                /* option is mandatory */
    },
        /* TraceSchema */
    {   _TraceSchemaName,   /* long name */
        _TraceSchemaAlias,  /* single character names */
        NULL,               /* function which generate help string */
        _TraceSchemaUsage,  /* help string(s) */
        1,                  /* max allowed values QTY, 0 - unlimited */
        true,               /* option requires argument */
        true                /* option is mandatory */
    },
        /* Validation Rules */
    {   _ValRulesName,      /* long name */
        _ValRulesAlias,     /* single character names */
        NULL,               /* function which generate help string */
        _ValRulesUsage,     /* help string(s) */
        1,                  /* max allowed values QTY, 0 - unlimited */
        true,               /* option requires argument */
        true                /* option is mandatory */
    },
        /* Validation Rules */
    {   _ErrCountName,      /* long name */
        _ErrCountAlias,     /* single character names */
        NULL,               /* function which generate help string */
        _ErrCountUsage,     /* help string(s) */
        1,                  /* max allowed values QTY, 0 - unlimited */
        true,               /* option requires argument */
        false               /* option is mandatory */
    },
        /* Validation Rules */
    {   _ErrPercentName,    /* long name */
        _ErrPercentAlias,   /* single character names */
        NULL,               /* function which generate help string */
        _ErrPercentUsage,   /* help string(s) */
        1,                  /* max allowed values QTY, 0 - unlimited */
        true,               /* option requires argument */
        false               /* option is mandatory */
    }
};

const char * _TraceLoaderOptionsParams [] =
{
    /* order here is same as in OptDef array above!!! */
    "path",
    "path",
    "path",
    "path",
    "path",
    "number",
    "number",
    "",
};

string
_getOptionValue (
                Args * TheArgs,
                const string & OptionName,
                bool Mondatory
)
{
    string RetVal = "";
    uint32_t OptQty = 0;
    rc_t RCt = ArgsOptionCount ( TheArgs, OptionName . c_str (), & OptQty );
    if ( RCt == 0 ) {
        if ( Mondatory ) {
            if ( OptQty == 0 ) {
                RCt = RC ( rcApp, rcArgv, rcAccessing, rcParam, rcInsufficient );
            }
        }

        if ( 1 < OptQty ) {
            RCt = RC ( rcApp, rcArgv, rcAccessing, rcParam, rcExcessive );
        }

        if ( RCt == 0 ) {
            if ( OptQty == 1 ) {
                const void * Value;
                RCt = ArgsOptionValue ( TheArgs, OptionName . c_str (), 0, & Value );
                if ( RCt == 0 ) {
                    RetVal = string ( ( char * ) Value );
                }
            }
        }
    }
    else {
        if ( Mondatory ) {
            RCt = RC ( rcApp, rcArgv, rcAccessing, rcParam, rcInsufficient );
        }
    }

    if ( RCt != 0 ) {
        throw TL_Exception ( string ( "Can not read value for option \"" ) + OptionName + "\"", RCt );
    }

    return RetVal;
}   /* _getOptionValue () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_* 
 *
 * Usages and helps and others
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

ver_t CC
KAppVersion ()
{
    return 0x01000000;
}   /* KAppVersion () */

rc_t
UsageSummary ( const char * ProgName )
{
    return KOutMsg (
                "Usage:\n"
                "\t%s [options] \n"
                "\n"
                "Summary:\n"
                "\tSubmits trace data to GenralWriter\n"
                "\n"
                , ProgName
                );
}   /* UsageSummary () */

char const UsageDefaultName[] = "trace-loader";

rc_t CC
Usage ( const Args * TheArgs )
{
    rc_t RCt;
    const char * ProgName = UsageDefaultName;
    const char * FullPath = UsageDefaultName;

    RCt = 0;

    if ( TheArgs == NULL ) {
        RCt  = RC ( rcApp, rcArgv, rcAccessing, rcSelf, rcNull );
    }
    else    {
        RCt = ArgsProgram ( TheArgs, &FullPath, &ProgName );
    }

    if ( RCt != 0 )
        ProgName = FullPath = UsageDefaultName;

    UsageSummary ( ProgName );

    size_t ArgsQty = sizeof ( _TraceLoaderOptions )
                                / sizeof ( _TraceLoaderOptions [ 0 ] );
    for ( size_t llp = 0; llp < ArgsQty; llp ++ ) {
        if( _TraceLoaderOptions [ llp ] . required
            && _TraceLoaderOptions [ llp ] . help [ 0 ] != NULL ) {
            HelpOptionLine (
                            _TraceLoaderOptions [ llp ] . aliases,
                            _TraceLoaderOptions [ llp ] . name,
                            _TraceLoaderOptionsParams [ llp ],
                            _TraceLoaderOptions [ llp ] . help
                            );
        }
    }

    OUTMSG(("\nOptions:\n"));
    for ( size_t llp = 0; llp < ArgsQty; llp ++ ) {
        if( ! _TraceLoaderOptions [ llp ] . required
            && _TraceLoaderOptions [ llp ] . help [ 0 ] != NULL ) {
            HelpOptionLine (
                            _TraceLoaderOptions [ llp ] . aliases,
                            _TraceLoaderOptions [ llp ] . name, \
                            _TraceLoaderOptionsParams [ llp ],
                            _TraceLoaderOptions [ llp ] . help
                            );
        }
        }
    XMLLogger_Usage();

    HelpOptionsStandard ();
    HelpVersion ( FullPath, KAppVersion() );

    return RCt;
}   /* Usage () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_* 
 *
 * Main
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

void _checkInitConfig ( Args * TheArgs, TL_TraceConfig & Config );
void _loadTrace ( TL_TraceConfig * Config );

rc_t CC
KMain ( int argc, char ** argv )
{
    rc_t RCt = 0;
    Args * TheArgs = NULL;
    uint32_t pCount = 0;
    const XMLLogger * XmlLogger = NULL;

    TL_TraceConfig Config;  // Because it is here

    try {
        RCt = ArgsMakeAndHandle (
                                & TheArgs,
                                argc,
                                argv,
                                2,
                                _TraceLoaderOptions,
                                sizeof ( _TraceLoaderOptions )
                                                    / sizeof ( OptDef ),
                                XMLLogger_Args,
                                XMLLogger_ArgsQty
                                );

        if ( RCt == 0 ) {
            RCt = XMLLogger_Make ( & XmlLogger, NULL, TheArgs );
            if ( RCt == 0 ) {
                RCt = ArgsParamCount ( TheArgs, & pCount );
                if ( RCt == 0 ) {
                    if ( pCount != 0 ) {
                        RCt = RC ( rcApp, rcArgv, rcAccessing, rcParam, rcExcessive );

                        MiniUsage ( TheArgs );
                    }
                    else {

                        _checkInitConfig ( TheArgs, Config );

// JOJOBA Config . SaySomething ();

                        _loadTrace ( & Config );


                        Config . ReportModule () . ReportSuccess ();
                    }
                }
                XMLLogger_Release ( XmlLogger );
            }

            ArgsWhack ( TheArgs );

        }
    }
    catch ( TL_Exception & E ) {
        PLOGERR (
            klogErr,
            ( klogErr,
            E . rc (),
            "Failed because : $(reason)",
            "severity=error,reason=%s",
            E . what ()
            )
            );
        Config . ReportModule () . ReportFailure ();

        RCt = E . rc () == 0 ? 1 : E . rc ();
    }
    catch ( exception & E ) {
        PLOGMSG (
            klogErr,
            ( klogErr,
            "Failed because : $(reason)",
            "severity=error,reason=%s",
            E . what ()
            )
            );
        Config . ReportModule () . ReportFailure ();
        RCt = 1;
    }
    catch ( ... ) {
        PLOGMSG (
            klogErr,
            ( klogErr,
            "Failed because : $(reason)",
            "severity=error,reason=%s",
            "reason unknown"
            )
            );
        Config . ReportModule () . ReportFailure ();
        RCt = 1;
    }

    return RCt;
}   /* main () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_* 
 *
 * Some stuff
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
void
_checkInitConfig ( Args * TheArgs, TL_TraceConfig & Config )
{
    string OptionValue = _getOptionValue (
                                    TheArgs,
                                    _OutputPathName,
                                    false
                                    );
    Config . SetTraceOutputPath ( OptionValue );

    OptionValue = _getOptionValue ( TheArgs, _TraceInfoName, true );
    Config . SetTraceInfoPath ( OptionValue );

    OptionValue = _getOptionValue ( TheArgs, _TraceSchemaName, true );
    Config . SetTraceSchemaPath ( OptionValue );

    OptionValue = _getOptionValue ( TheArgs, _ValRulesName, true );
    Config . SetTraceValidationRulesPath ( OptionValue );

    OptionValue = _getOptionValue ( TheArgs, _ErrCountName, false );
    if ( ! OptionValue . empty () ) {
        Config . ReportModule () . SetMaxErrCount ( OptionValue );
    }

    OptionValue = _getOptionValue ( TheArgs, _ErrPercentName, false );
    if ( ! OptionValue . empty () ) {
        Config . ReportModule () . SetMaxErrPercent ( OptionValue );
    }

}   /* _checkInitConfig () */

void
_loadTrace ( TL_TraceConfig * Config )
{
    PLOGMSG(
        klogInfo,
        (klogInfo,
        "Loading traces",
        "severity=info,file=%s",
        Config -> TraceInfoPath () . c_str ()
        )
        );

    TL_TraceInfoValidator Validator;
    Validator . Init ( Config );
    Validator . Validate ();

    TL_OVec Vec;
    Validator . Export ( Vec );

    TL_SimpleTraceInfo TraceInfo;
    TL_SimpleTraceInfoLoader () . Load ( TraceInfo, Vec );

    TL_TraceData TraceData;
    TraceData . Init ( Config );

    PLOGMSG(
        klogInfo,
        (klogInfo,
        "[Loading Batch]",
        "severity=info"
        )
        );

    TL_TraceFieldAdapter FA ( Config -> TraceOutput () );
    FA . Init ( Config -> TraceSchema () );
    FA . OpenWriter ( TraceInfo );

    for ( size_t llp = 0; llp < TraceInfo . RowQty (); llp ++ ) {
        try {
            TL_Trace_DP DP;
            DP . Set ( TraceData, TraceInfo, llp );
            FA . WriteTraceData ( DP, TraceInfo, llp );

            Config -> ReportModule () . ReportProgress ();
            Config -> ReportModule () . AddToReport (
                            TraceInfo . Value ( llp, _TL_TRACE_FILE )
                            );
        }
        catch ( exception & E ) {
            if ( Config -> ReportModule () . ReportErrCheckIfShouldFail () ) {
                throw;
            }

            PLOGMSG(
                    klogWarn,
                    (klogWarn,
                    "$(message)",
                    "severity=warning,message=%s",
                    E . what ()
                    )
                    );
        }
        catch ( ... ) {
            if ( Config -> ReportModule () . ReportErrCheckIfShouldFail () ) {
                throw;
            }

            PLOGMSG(
                    klogWarn,
                    (klogWarn,
                    "unknown exception handled",
                    "severity=warning"
                    )
                    );
        }
    }

    TraceData . Dispose ();
    TraceInfo . Reset ();
    Validator . Dispose ();
}   /* _loadTrace () */
