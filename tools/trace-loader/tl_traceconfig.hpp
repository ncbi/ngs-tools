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

#ifndef _tl_traceconfig_
#define _tl_traceconfig_

#include <string>
#include <map>
#include <list>

#include "tl_exception.hpp"
#include "tl_report.hpp"

namespace _tl_ {

/*))
  ||  it it too late to develop that class :LOL:
  ||
  ((*/
class TL_TraceConfig {
public:
    TL_TraceConfig ();
    ~TL_TraceConfig ();

    void Reset ();

        /* Resolving Path */
    std::string NormalizePath ( const std::string & Path ) const;

        /* Data accessing methods */
    inline const std::string & TraceInfoPath () const
                                    { return _trace_info_path; };

    inline const std::string & TraceInfoOwpPath () const
                                    { return _trace_info_owp_path; };

    inline const std::string & TraceSchema () const
                                { return _trace_schema; };

    inline const std::string & TraceOutput () const
                                { return _trace_output; };


    inline const std::string & ValidationRules () const
                                { return _validation_rules; };

        /* Data setting methods */
    void SetTraceInfoPath ( const std::string & InfoPath );

    void SetTraceSchemaPath ( const std::string & SchemaPath );

    void SetTraceOutputPath ( const std::string & OutputPath );

    void SetTraceValidationRulesPath (
                            const std::string & ValidationRulesPath
                            );

        /*  ERROR HANDLING: and other reporting
         */
    inline TL_Report & ReportModule ()
                                { return _report; };

        /* Very personal method */
    void SaySomething () const;

private:
    std::string _trace_info_path;
    std::string _trace_info_owp_path;

    std::string _trace_schema;

    std::string _trace_output;

    std::string _validation_rules;

    TL_Report _report;
};

}   /* namespace _tl_ */

#endif /* _tl_traceconfig_ */

