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

#include <exception>
#include <string>

#include <klib/printf.h>
#include <klib/rc.h>

#include "tl_exception.hpp"

using namespace std;
using namespace _tl_;

TL_Exception :: TL_Exception ( rc_t rc )
:   _rc ( rc )
,   _message ( "" )
{
    char BB [ 4096 ];
    size_t NumWrit;

    rc_t RCt = string_printf ( BB, sizeof ( BB ), & NumWrit, "%R", rc );
    if ( RCt != 0 ) {
        RCt = string_printf (
                            BB,
                            sizeof ( BB ),
                            & NumWrit,
                            "Unknown error reported ( RC = %d )",
                            rc
                            );
    }

    if ( RCt != 0 ) {
        _message = "Help me! I feel bad";
    }
    else {
        _message = BB;
    }
}   /* TL_Exception :: TL_Exception () */

TL_Exception :: TL_Exception ( const string & Message, rc_t rc )
:   _rc ( rc )
,   _message ( Message )
{
    /* JOJOBA to think what to do with RC */
}   /* TL_Exception :: TL_Exception () */

TL_Exception :: ~TL_Exception () throw ()
{
    _rc = 0;
}   /* TL_Exception :: ~TL_Exception () */

const char *
TL_Exception :: what () const throw ()
{
    return _message . c_str ();
}   /* TL_Exception :: what () */

const enum RCModule
TL_Exception :: module () const
{
    return GetRCModule ( _rc );
}   /* TL_Exception :: module () */

const enum RCTarget
TL_Exception :: target () const
{
    return GetRCTarget ( _rc );
}   /* TL_Exception :: target () */

const enum RCContext
TL_Exception :: context () const
{
    return GetRCContext ( _rc );
}   /* TL_Exception :: context () */

const enum RCObject
TL_Exception :: object () const
{
    return GetRCObject ( _rc );
}   /* TL_Exception :: object () */

const enum RCState
TL_Exception :: state () const
{
    return GetRCState ( _rc );
}   /* TL_Exception :: state () */

