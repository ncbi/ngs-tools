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

#ifndef _tl_exception_
#define _tl_exception_

#include <exception>
#include <string>

#include <kfc/defs.h>   /* definition of rc_t */
#include <kfc/rc.h>     /* RCModule/Target/Context/Object */

namespace _tl_ {

class TL_Exception : public std::exception {

public:
    TL_Exception ( rc_t rc );
    TL_Exception ( const std::string & Msg, rc_t rc = 0 );

    virtual ~TL_Exception () throw ();

    const char *what () const throw ();

    inline const rc_t rc () const { return _rc; };

    const enum RCModule  module () const;
    const enum RCTarget  target () const;
    const enum RCContext context () const;
    const enum RCObject  object () const;
    const enum RCState   state () const;

private:
    rc_t _rc;

    std::string _message;
};

}   /* namespace _tl_ */

#endif /* _tl_exception_ */

