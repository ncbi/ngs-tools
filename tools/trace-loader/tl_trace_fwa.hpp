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

#ifndef _tl_trace_fwa_
#define _tl_trace_fwa_

#include <string>
#include <map>
#include <vector>

#include "tl_types.hpp"
#include "tl_owp.hpp"
#include "tl_exception.hpp"

#include <general-writer.hpp>

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * That one needed to configure GeneralWriter and provide correct prin
 * ting to the writer :D
 *
 * All methods are static
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

namespace ncbi {

    class GeneralWriter;

};

namespace _tl_ {

class TL_TraceField;
class TL_TraceInfo;
class TL_StreamId;

    /*  Simple abstract class with pure virtual methods
     */
class TL_TraceFieldDataProvider {
public:
    TL_TraceFieldDataProvider ();
    virtual ~TL_TraceFieldDataProvider ();

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

};

class TL_TraceFieldAdapter {
public:
    typedef std::map < std::string, TL_StreamId * > SMap;

public:
    TL_TraceFieldAdapter ( const std::string & TraceOut );
    ~TL_TraceFieldAdapter ();

    void Init ( const std::string & Schema );

    void AddColumn ( const TL_TraceField & Field );

    bool HasColumn ( const std::string & ColumnName ) const;

    void OpenWriter ( const TL_TraceInfo & Info );

    void WriteTraceData (
                    const TL_TraceDataProvider & DataProvider,
                    const TL_TraceInfo & Info,
                    size_t RowIdx
                    );

protected:
    void _Reset ();

    void _AddMandatoryColumns ();
    void _AddOtherColumns ( const TL_TraceInfo & Info );
    void _SetColumnDefaults ();

    TL_StreamId * GetStreamId (
                        const std::string & Name,
                        bool Throw = true
                        );

    void _FillStream ( const std::string & StreamName, size_t Qty );

    void _WrtieDerivatives (
                        const TL_TraceInfo & TraceInfo,
                        size_t RowIdx
                        );

private:
    const std::string _output;
    ncbi :: GeneralWriter * _writer;

    int _table_id;

    SMap _stream_ids;

    size_t _filler_capacity;
    char * _filler;
};


}   /* namespace _tl_ */

#endif /* _tl_trace_fwa_ */

