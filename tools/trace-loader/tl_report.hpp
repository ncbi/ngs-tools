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

#ifndef _tl_report_
#define _tl_report_

#include <string>
#include <map>
#include <list>

#include "tl_exception.hpp"

namespace _tl_ {

/*))
  ||  That class is for error handling, etc
  ||
  ((*/
class TL_Report {
public:
    TL_Report ();
    ~TL_Report ();

    void Reset ();

        /*  ERROR HANDLING: That prolly should not be here :OLOLO:
         */
        /*  Calling one of those methods will reset all error related
         *  variables ... particularly, all errors amount computation 
         *  will be set to zero
         */
    void SetMaxErrCount ( const std::string & ErrCount = "" );
    void SetMaxErrPercent ( const std::string & ErrPercent = "" );

    inline size_t MaxErrCount () const
                                { return _max_err_count; };
    inline size_t MaxErrPercent () const
                                { return _max_err_percent; };


        /*  Calling this will recalculate error threshold and reset
         *  error count
         */
    void SetItemCount ( size_t ItemCount );
    inline size_t ItemCount ( size_t ItemCount ) const
                                { return _item_count; };

        /*  Here we are dropping milestones
         *  We should fail out if it returns True
         */
    bool ReportErrCheckIfShouldFail ();

        /*  Here we will paing gauge
         */
    void ReportProgress ();

        /*  Here we add file and spot to report
         */
    void AddToReport ( const std::string & Path );

        /*  Here we will report total
         */
    void ReportSuccess ();
    void ReportFailure ();

        /*  Context
         */
    inline const std::string & Context () const
                                            { return _context; };
    inline void SetContext ( const std::string & Cnt = "" )
                                            { _context = Cnt; };
    inline const std::string & Meanings () const
                                            { return _meanings; };
    inline void SetMeanings ( const std::string & Mng = "" )
                                            { _meanings = Mng; };

private:
    size_t _max_err_count;
    size_t _max_err_percent;
    bool _max_err_count_set;
    size_t _fail_threshold;
    size_t _processed_count;
    size_t _failed_count;
    size_t _item_count;
    size_t _last_gauge_pos;

    TL_SVec _to_report;

    std::string _context;
    std::string _meanings;
};

}   /* namespace _tl_ */

#endif /* _tl_report_ */

