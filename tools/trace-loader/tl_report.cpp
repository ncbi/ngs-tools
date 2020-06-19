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
#include "tl_report.hpp"


using namespace std;
using namespace _tl_;

#define _DEFAULT_ERR_COUNT   100
#define _DEFAULT_ERR_PERCENT 5

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_Report
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_Report :: TL_Report ()
:   _max_err_count ( _DEFAULT_ERR_COUNT )
,   _max_err_percent ( _DEFAULT_ERR_PERCENT )
,   _max_err_count_set ( false )
,   _fail_threshold ( 0 )
,   _processed_count ( 0 )
,   _failed_count ( 0 )
,   _item_count ( 0 )
,   _last_gauge_pos ( 0 )
,   _to_report ()
,   _context ( "" )
,   _meanings ( "" )
{
}   /* TL_Report :: TL_Report () */

TL_Report :: ~TL_Report ()
{
    TL_TRY {
        Reset ();
    }
    TL_CATCH_R
}   /* TL_Report :: ~TL_Report () */

void
TL_Report :: Reset ()
{
    _context . clear ();
    _to_report . clear ();
    _max_err_count_set = false;
    _fail_threshold = 0;
    _processed_count = 0;
    _failed_count = 0;
    _item_count = 0;
    _last_gauge_pos = 0;
    _max_err_percent = _DEFAULT_ERR_PERCENT;
    _max_err_count = _DEFAULT_ERR_COUNT;
}   /* TL_Report :: Reset () */

void
TL_Report :: SetMaxErrCount ( const std::string & ErrCount )
{
        /*  Resetting that
         */
    _fail_threshold = 0;
    _processed_count = 0;
    _failed_count = 0;
    _item_count = 0;
    _last_gauge_pos = 0;

    _max_err_count = ErrCount . empty ()
                    ? _DEFAULT_ERR_COUNT
                    : TL_StringU :: FromStr < uint32_t > ( ErrCount )
                    ;
    _max_err_count_set = true;
}   /* TL_Report :: SetMaxErrCount () */

void
TL_Report :: SetMaxErrPercent ( const std::string & ErrPercent )
{
        /*  Resetting that
         */
    _fail_threshold = 0;
    _processed_count = 0;
    _failed_count = 0;
    _item_count = 0;
    _last_gauge_pos = 0;

    _max_err_percent = ErrPercent . empty ()
                    ?  _DEFAULT_ERR_PERCENT
                    : TL_StringU :: FromStr < uint32_t > ( ErrPercent )
                    ;
    if ( 100 < _max_err_percent ) {
        _max_err_percent = 100;
    }

    _max_err_count_set = false;
}   /* TL_Report :: SetMaxErrCount () */

void
TL_Report :: SetItemCount ( size_t ItemCount )
{
    _fail_threshold = 0;
    _processed_count = 0;
    _failed_count = 0;
    _last_gauge_pos = 0;
    _to_report . clear ();

    if ( ItemCount == 0 ) {
        throw TL_Exception ( "TraceConfig: ZERO SubmissionCount set" );
    }

    _item_count = ItemCount;

    if ( _max_err_count_set ) {
        _fail_threshold = _max_err_count;
    }
    else {
            /*  Here we are :LOL:
             */
        _fail_threshold = ( int ) floor ( ( ( ( float ) _item_count / ( float ) 100 ) * ( float ) _max_err_percent ) + 0.5f );

        if ( _fail_threshold == 0 ) {
            _fail_threshold = 1;
        }

        if ( 100 < _fail_threshold ) {
            _fail_threshold = 100;
        }
    }

    if ( _item_count < _fail_threshold ) {
        _fail_threshold = _item_count;
    }

    PLOGMSG(
            klogInfo,
            (klogInfo,
            "[TR_CFG] setting failure threshold $(fail_threshold)",
            "severity=info,fail_threshold=%u",
            _fail_threshold
            )
            );
}   /* TL_Report :: SetItemCount () */

bool
TL_Report :: ReportErrCheckIfShouldFail ()
{
        /*  Updateing stats 
         */
    _failed_count ++;

    if ( _fail_threshold <= _failed_count ) {
        PLOGMSG(
                klogWarn,
                (klogWarn,
                "There were too many errors ( all ($(qty)), but failure threshold ($(thres))",
                "severity=warnning,qty=%u,thres=%u",
                _failed_count,
                _fail_threshold
                )
                );
        return true;
    }

    return false;
}   /* TL_Report :: ReportErrCheckIfShouldFail () */

void
TL_Report :: ReportProgress ()
{
        /*  Updateing stats 
         */
    _processed_count ++;

        /*  Updateing gauge
         */
    size_t NewGauge =
        ( int ) floor ( ( ( ( float ) _processed_count / ( float ) _item_count ) * ( float ) 100 ) + 0.5f );

    if ( _last_gauge_pos < NewGauge ) {
        _last_gauge_pos = NewGauge;

        PLOGMSG(
                klogInfo,
                (klogInfo,
                "processed $(percent)% ($(total) $(meanings), $(errors) failed)",
                "severity=status,percent=%u,total=%u,meanings=%s,errors=%u",
                _last_gauge_pos,
                _processed_count,
                ( _meanings . empty () ? "traces" : _meanings . c_str () ),
                _failed_count
                )
                );
    }
}   /* TL_Report :: ReportProgress () */

void
TL_Report :: AddToReport ( const string & Path )
{
    _to_report . insert ( _to_report . end (), Path );
}   /* TL_Report :: AddToReport () */

void
TL_Report :: ReportSuccess ()
{
    string FormatStr = string ( "severity=total," ) + _meanings + "=%u,errors=%u,status=success";
        /*  First is a status
         */
    PLOGMSG(
            klogInfo,
            (klogInfo,
            "loaded",
            FormatStr . c_str (),
            _processed_count,
            _failed_count
            )
            );

        /*  Second are produced file(s)
         */
    for ( 
        TL_SVec :: const_iterator It = _to_report . begin ();
        It != _to_report . end ();
        It ++
    ) {
        PLOGMSG(
                klogInfo,
                (klogInfo,
                "loaded",
                "severity=result,file=%s,spot_count=%d",
                It -> c_str (),
                1
                )
                );
    }
}   /* TL_Report :: ReportSuccess () */

void
TL_Report :: ReportFailure ()
{
        /*  First is a status
         */
    PLOGMSG(
            klogInfo,
            (klogInfo,
            "program failed",
            "severity=total,file=%s,errors=%u,status=failure",
            ( _context . empty () ? "not set" : _context . c_str () ),
            _failed_count
            )
            );

}   /* TL_Report :: ReportFailure () */
