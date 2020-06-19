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

#ifndef _tl_validator_
#define _tl_validator_

#include <string>
#include <list>

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * LYrIcs:
 *
 *    That class loads trace info file and does it's validation
 *
 *       a) Syntax validation
 *       b) Content validation
 *
 *    Attempt to made external validator.
 *
 *    As original validator, it will crash in Syntax errors, and will
 *    fall after some amount of errors in content validation.
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

#include "tl_traceinfo.hpp"
#include "tl_traceconfig.hpp"
#include "tl_owp.hpp"

namespace _tl_ {

class TL_TraceConfig;

class TL_TraceInfoValidator {
public:
    TL_TraceInfoValidator ();
    ~TL_TraceInfoValidator ();

    void Init ( 
                TL_TraceConfig * Config,
                size_t FailureRatioPerCent = 5
                );
    void Dispose ();

        /*  Just validates TraceInfo
         */
    bool Validate ();

        /*  Returns valid TraceInfo OWP
         */
    void Export ( TL_OVec & ValidInfo );

        /*  Creates and stores TraceInfo OWP to some file
         *  If Path is empty, the path from config will be taken
         */
    void Export ( const std::string & Path = "" );

private:
    void * _validator;

    TL_TraceConfig * _config;
};

}   /* namespace _tl_ */

#endif /* _tl_validator_ */

