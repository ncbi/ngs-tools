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

#ifndef _tl_trace_dp_
#define _tl_trace_dp_

#include <string>

#include "tl_types.hpp"


namespace _tl_ {

    /*  Simple abstract class with pure virtual methods
     *  It will be feeded to General Writer to retrieve data to write
     */
class TL_TraceDataProvider {
public:
    TL_TraceDataProvider ();
    virtual ~TL_TraceDataProvider ();

        /*  Name : TRACE_NAME
         */
    virtual const std::string & Name () const = 0;

        /*  ProgramId : PROGRAM_ID
         */
    virtual const std::string & ProgramID () const = 0;

        /*  Bases : "READ" and "READ_LEN"
         */
    virtual const char_a_t & Bases () const = 0;

        /*  Peaks : "POSITION"
         */
    virtual const uint32_a_t & PeakIndex () const = 0;

        /*  ProbScores : "QUALITY" and "BASES_(20/40/60)"
         */
    virtual const uchar_a_t & ProbScores () const = 0;
    virtual uint16_t Bases20 () const = 0;
    virtual uint16_t Bases40 () const = 0;
    virtual uint16_t Bases60 () const = 0;

        /*  Traces : "SIGNAL" and "TRACE_MAX_VALUE"
         */
    virtual const uint64_a_t & SampleCombined () const = 0;
    virtual uint16_t MaxTraceVal () const = 0;

        /*  Clips, only two : QUALITY_LEFT and QUALITY_RIGHT
         */
    virtual uint32_t ClipQualityRight () const = 0;
    virtual uint32_t ClipQualityLeft () const = 0;

        /*  FlowChars and KeySequence : FLOW_CHARS and KEY_SEQUENCE
         */
    virtual const char_a_t & FlowChars () const = 0;
    virtual const char_a_t & KeySequence () const = 0;

        /*  Comment are not part of Trace table
         */
    virtual const std::string & Comments () const = 0;

        /*  ExtendedData are stored everywhere
         */
    virtual const std::string & ExtendedData () const = 0;
};

}   /* namespace _tl_ */

#endif /* _tl_trace_dp_ */

