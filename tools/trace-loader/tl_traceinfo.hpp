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

#ifndef _tl_traceinfo_
#define _tl_traceinfo_

#include <string>
#include <list>

#include "tl_exception.hpp"
#include "tl_util.hpp"
#include "tl_types.hpp"
#include "tl_tracefields.hpp"

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 * 
 * Lyrics:
 *
 * That file contains description of TL_TraceInfo and TL_SimpleTraceInfo
 * classes. They used by WriterAdapter to read and write actual data.
 *
 * TL_TraceInfo is an abstract class which contains description of
 * valid TraceSubmission. It includes:
 *     Fields and their descriptions for submission
 *     Meta information for submitted traces
 *     Methods to access Field descriptions and Meta data
 *
 * The column description is vector of TL_TraceFields
 * The metadata is accessed by rows, and each row is a vector or
 * strings, which are values.
 *
 * TL_SimpleTraceInfo - class which methods to fill trace info
 * 
 * NOTE : we suppose that all info is valid, so there is no any checks
 * 
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

namespace _tl_ {

/*)))
 ///   TL_TraceInfo - basic class
(((*/
class TL_TraceInfo {
public:
    TL_TraceInfo ();
    virtual ~TL_TraceInfo ();

        /*  Field descriptions return vector of TL_TraceFields
         */
    virtual size_t FieldQty () const = 0;
    virtual const TL_TraceField & Field ( size_t FieldIndex ) const = 0;

        /*  Accessing Metadata by row number, for this class
         *  it is always zero rows and all them empty
         */
    virtual size_t RowQty () const = 0;
    virtual const TL_SVec & Row ( size_t RowIndex ) const = 0;

        /*  Checking if field with name exists
         */
    virtual bool HasField ( const std::string & FieldName ) const = 0;

        /*  Retrieving metadata by field name and row number
         */
    virtual const std::string & Value (
                            size_t RowIndex,
                            const std::string & FieldName
                            ) const = 0;
};  /* TL_TraceInfo */

/*)))
 ///   TL_SimpleTraceInfo - simple class
(((*/
class TL_SimpleTraceInfo : public TL_TraceInfo {
public:
    TL_SimpleTraceInfo ();
    ~TL_SimpleTraceInfo ();

        /*  NOTE: Call it if You are planning to reuse object
         */
    void Reset ();

        /*  NOTE: set fields before adding rows, or there will be
         *        exception throwed
         */
    void AddField ( const std::string & FieldName );
    size_t FieldQty () const;
    const TL_TraceField & Field ( size_t FieldIndex ) const;

        /*  NOTE: set fields before adding rows, or there will be
         *        exception throwed
         */
    void AddRow ( const TL_SVec & Row );
    size_t RowQty () const;
    const TL_SVec & Row ( size_t RowIndex ) const;

    bool HasField ( const std::string & FieldName ) const;

    const std::string & Value (
                            size_t RowIndex,
                            const std::string & FieldName
                            ) const;

private:
        /*  Apparently map < string, int > would be enough, but ...
         *  may be I will put here fake fields later :D
         */
    TL_TFVec _fields;
    TL_CISIMap _field_index;
    TL_SVVec _rows;
};

/*)))
 ///   TL_SimpleTraceInfoLoader - simple OWP style data loader
(((*/
class TL_SimpleTraceInfoLoader {
public:
    TL_SimpleTraceInfoLoader ();
    ~TL_SimpleTraceInfoLoader ();

    void Load ( TL_SimpleTraceInfo & Info, const std::string & Path );

    void Load ( TL_SimpleTraceInfo & Info, const TL_OVec & Vec );

private:
    void _read_fields (
                    TL_SimpleTraceInfo & TraceInfo,
                    const TL_OVec & Vec
                    );
    void _make_rows (
                    TL_SimpleTraceInfo & TraceInfo,
                    const TL_OVec & Vec
                    );
};

}   /* namespace _tl_ */

#endif /* _tl_traceinfo_ */

