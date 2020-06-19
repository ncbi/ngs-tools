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

#ifndef _tl_tracedata_reader_
#define _tl_tracedata_reader_

#include <string>
#include <map>
#include <list>

#include "tl_util.hpp"
#include "tl_tracedata.hpp"

namespace _tl_ {

/*******************************************************************
 *  That class will give us a joy of reading Trace Data
 *  Something as general as a dawg
 *******************************************************************/
class TL_TraceDataReader {
public:
    TL_TraceDataReader ();
    virtual ~TL_TraceDataReader ();

    virtual void Read ( const TL_TraceFile & File ) = 0;
    virtual void Export ( TL_Traces & Traces ) const = 0;


};

/*******************************************************************
 *  That has unified export stuff :D
 *  Basically, that is almost full set of fields we need
 *  to store to VDB
 *******************************************************************/
class TL_TraceDataReaderWithExport {
public:
    TL_TraceDataReaderWithExport ();
    ~TL_TraceDataReaderWithExport ();

    void Export ( TL_Traces & Traces ) const;

    void Reset ();

protected:

        /* Samples */
    uint16_a_t _samples_A;
    uint16_a_t _samples_C;
    uint16_a_t _samples_G;
    uint16_a_t _samples_T;
    uint16_t _max_sample_value;

        /* Bases, peakindex and qual scores */
    char_a_t _bases;
    uint32_a_t _peak_index;
    uchar_a_t _prob_A;
    uchar_a_t _prob_C;
    uchar_a_t _prob_G;
    uchar_a_t _prob_T;
    bool _valid_scores;

    uchar_a_t _prob_sub;            /* SCF */
    uchar_a_t _prob_ins;            /* SCF */
    uchar_a_t _prob_del;            /* SCF */
    bool _extra_probs;              /* SCF */

        /* Misc and other */
    std::string _comments;          /* ZTR, ABI, SCF, SFF */
    uchar_a_t _private_data;        /* SCF */

    uint32_t _clip_quality_right;   /* ZTR, SFF */
    uint32_t _clip_quality_left;    /* ZTR, SFF */

    uint32_t _clip_adapter_right;   /* come with TRACEINFO.tbl */
    uint32_t _clip_adapter_left;    /* come with TRACEINFO.tbl */
};

}   /* namespace _tl_ */

#endif /* _tl_tracedata_reader_ */

