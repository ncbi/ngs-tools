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

#include <klib/rc.h>
#include <klib/namelist.h>

#include "tl_tracedata.hpp"
#include "tl_names.hpp"
#include "tl_tracedata_reader.hpp"


using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceDataReader
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

TL_TraceDataReader :: TL_TraceDataReader ()
{
}   /* TL_TraceDataReader :: TL_TraceDataReader () */

TL_TraceDataReader :: ~TL_TraceDataReader ()
{
}   /* TL_TraceDataReader :: ~TL_TraceDataReader () */

void
TL_TraceDataReader :: Read ( const TL_TraceFile & File )
{
    throw TL_Exception ( "TL_TraceDataReader :: Read () : kinda do not call unimplemented methods :LOL:" );
}   /* TL_TraceDataReader :: Read () */

void
TL_TraceDataReader :: Export ( TL_Traces & Traces ) const
{
    throw TL_Exception ( "TL_TraceDataReader :: Write () : kinda do not call unimplemented methods :LOL:" );
}   /* TL_TraceDataReader :: Export () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceDataReaderWithExport
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

TL_TraceDataReaderWithExport :: TL_TraceDataReaderWithExport ()
:   _samples_A ()
,   _samples_C ()
,   _samples_G ()
,   _samples_T ()
,   _max_sample_value ( 0 )
,   _bases ()
,   _peak_index ()
,   _prob_A ()
,   _prob_C ()
,   _prob_G ()
,   _prob_T ()
,   _valid_scores ( false )
,   _prob_sub ()
,   _prob_ins ()
,   _prob_del ()
,   _extra_probs ( false )
,   _comments ( "" )
,   _private_data ()
,   _clip_quality_right ( 0 )
,   _clip_quality_left ( 0 )
,   _clip_adapter_right ( 0 )
,   _clip_adapter_left ( 0 )
{
}   /* TL_TraceDataReaderWithExport :: TL_TraceDataReaderWithExport () */

TL_TraceDataReaderWithExport :: ~TL_TraceDataReaderWithExport ()
{
}   /* TL_TraceDataReaderWithExport :: ~TL_TraceDataReaderWithExport () */

void
TL_TraceDataReaderWithExport :: Export ( TL_Traces & Traces ) const
{
        /* Samples first
         */
    Traces . Samples () . SampleA () = _samples_A;
    Traces . Samples () . SampleC () = _samples_C;
    Traces . Samples () . SampleG () = _samples_G;
    Traces . Samples () . SampleT () = _samples_T;

        /* Bases next
         */
    Traces . Bases () . Bases () = _bases;
    Traces . Bases () . PeakIndex () = _peak_index;
    if ( _valid_scores ) {
        Traces . Bases () . ProbA () = _prob_A;
        Traces . Bases () . ProbC () = _prob_C;
        Traces . Bases () . ProbG () = _prob_G;
        Traces . Bases () . ProbT () = _prob_T;
        Traces . Bases () . SetValidScores ( _valid_scores );
    }
    if ( _extra_probs ) {
        Traces . Bases () . ProbSub () = _prob_sub;
        Traces . Bases () . ProbIns () = _prob_ins;
        Traces . Bases () . ProbDel () = _prob_del;
        Traces . Bases () . SetExtraProbs ( _extra_probs );
    }

        /* Other misc
         */
    if ( _comments . size () != 0 ) {
        Traces . SetComments ( _comments );
    }

    if ( _private_data . size () ) {
        Traces . SetPrivateData ( _private_data );
    }

    if ( _clip_quality_left != 0 )
        Traces . SetClipQualityLeft ( _clip_quality_left );
    if ( _clip_quality_right != 0 )
        Traces . SetClipQualityRight ( _clip_quality_right );
    if ( _clip_adapter_left != 0 )
        Traces . SetClipAdapterLeft ( _clip_adapter_left );
    if ( _clip_adapter_right != 0 )
        Traces . SetClipAdapterRight ( _clip_adapter_right );
}   /* TL_TraceDataReaderWithExport :: Export () */

void
TL_TraceDataReaderWithExport :: Reset ()
{
    _samples_A . reset ();
    _samples_C . reset ();
    _samples_G . reset ();
    _samples_T . reset ();
    _max_sample_value = 0;
    _bases . reset ();
    _peak_index . reset ();
    _prob_A . reset ();
    _prob_C . reset ();
    _prob_G . reset ();
    _prob_T . reset ();
    _valid_scores = false;
    _prob_sub . reset ();
    _prob_ins . reset ();
    _prob_del . reset ();
    _extra_probs = false;
    _comments . clear ();
    _private_data . reset ();
    _clip_quality_right = 0;
    _clip_quality_left = 0;
    _clip_adapter_right = 0;
    _clip_adapter_left = 0;
}   /* TL_TraceDataReaderWithExport :: Reset () */
