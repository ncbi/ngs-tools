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
#include <sstream>
#include <set>

#include <sra/types.h>
#include <insdc/sra.h>

#include "tl_util.hpp"
#include "tl_names.hpp"

#include "tl_trace_dp.hpp"


using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceDataProvider - or field which comes with schema, and ...
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_TraceDataProvider :: TL_TraceDataProvider ()
{
}   /* TL_TraceDataProvider :: TL_TraceDataProvider () */

TL_TraceDataProvider :: ~TL_TraceDataProvider ()
{
}   /* TL_TraceDataProvider :: ~TL_TraceDataProvider () */

const std::string &
TL_TraceDataProvider :: Name () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: Name ()\"" );
}   /* TL_TraceDataProvider :: Name () */

const std::string &
TL_TraceDataProvider :: ProgramID () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: ProgramID ()\"" );
}   /* TL_TraceDataProvider :: ProgramID () */

const char_a_t &
TL_TraceDataProvider :: Bases () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: Bases ()\"" );
}   /* TL_TraceDataProvider :: Bases () */

const uint32_a_t &
TL_TraceDataProvider :: PeakIndex () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: PeakIndex ()\"" );
}   /* TL_TraceDataProvider :: PeakIndex () */

const uchar_a_t &
TL_TraceDataProvider :: ProbScores () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: ProbScores ()\"" );
}   /* TL_TraceDataProvider :: ProbScores () */

uint16_t
TL_TraceDataProvider :: Bases20 () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: Bases20 ()\"" );
}   /* TL_TraceDataProvider :: Bases20 () */

uint16_t
TL_TraceDataProvider :: Bases40 () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: Bases40 ()\"" );
}   /* TL_TraceDataProvider :: Bases40 () */

uint16_t
TL_TraceDataProvider :: Bases60 () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: Bases60 ()\"" );
}   /* TL_TraceDataProvider :: Bases60 () */

const uint64_a_t &
TL_TraceDataProvider :: SampleCombined () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: SampleCombined ()\"" );
}   /* TL_TraceDataProvider :: SampleCombined () */

uint16_t
TL_TraceDataProvider :: MaxTraceVal () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: MaxTraceVal ()\"" );
}   /* TL_TraceDataProvider :: MaxTraceVal () */

uint32_t
TL_TraceDataProvider :: ClipQualityRight () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: ClipQualityRight ()\"" );
}   /* TL_TraceDataProvider :: ClipQualityRight () */

uint32_t
TL_TraceDataProvider :: ClipQualityLeft () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: ClipQualityLeft ()\"" );
}   /* TL_TraceDataProvider :: ClipQualityLeft () */

const string &
TL_TraceDataProvider :: Comments () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: Comments ()\"" );
}   /* TL_TraceDataProvider :: Comments () */

const char_a_t &
TL_TraceDataProvider :: FlowChars () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: FlowChars ()\"" );
}   /* TL_TraceDataProvider :: FlowChars () */

const char_a_t &
TL_TraceDataProvider :: KeySequence () const
{
    throw TL_Exception ( "Unimplemented method \"TL_TraceDataProvider :: KeySequence ()\"" );
}   /* TL_TraceDataProvider :: KeySequence () */

