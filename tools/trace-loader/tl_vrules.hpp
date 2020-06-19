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

#ifndef _tl_vrules_
#define _tl_vrules_

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * LYrIcs:
 *
 *    That class contains simple configurable external validator
 *    Class receives as an imput file with rules in Owp format.
 *    Description of that file look in tl_vrules.cpp
 *
 *    It made to reduce coding ...
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

#include "tl_traceinfo.hpp"
#include "tl_traceconfig.hpp"
#include "tl_owp.hpp"

namespace _tl_ {

class TL_TraceInfo;

class TL_VRules {
public:
    TL_VRules ();
    ~TL_VRules ();

    void Init ( const std::string & Path );
    void Dispose ();

        /*  if validation failed, Message could contain reason why
         */
    bool Validate ( const TL_Owp & Row, std::string & Message );
};

}   /* namespace _tl_ */

#endif /* _tl_vrules_ */

