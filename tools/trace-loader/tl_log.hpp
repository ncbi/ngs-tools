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
 *    Couple macroces to prevent much copy/pasta
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

#include <klib/log.h>
#include "tl_exception.hpp"

    /*  Just a try
     */
#define TL_TRY     try

    /*  Catches, and reports
     */
#define TL_CATCH_R  \
            catch ( TL_Exception & E ) {    \
                PLOGERR (    \
                        klogErr,    \
                        ( klogErr, E . rc (), "[Exception handled] $(msg)", "msg=%s", E . what () )    \
                        );    \
            }    \
            catch ( exception & E ) {    \
                PLOGMSG (    \
                        klogErr,    \
                        ( klogErr, "[Exception handled] $(msg)", "msg=%s", E . what () )    \
                        );    \
            }    \
            catch ( int & i ) {    \
                PLOGMSG (    \
                        klogErr,    \
                        ( klogErr, "[Exception handled] error code: $(err_c)", "err_c=%d", i )    \
                        );    \
            }    \
            catch ( ... ) {    \
                PLOGMSG (    \
                        klogErr,    \
                        ( klogErr, "[Exception handled] unknown", "severity=error" )    \
                        );    \
            }

    /*  Catches, reports and trhow
     */
#define TL_CATCH_RT  \
            catch ( TL_Exception & E ) {    \
                PLOGERR (    \
                        klogErr,    \
                        ( klogErr, E . rc (), "[Exception handled] $(msg)", "msg=%s", E . what () )    \
                        );    \
                throw;      \
            }    \
            catch ( exception & E ) {    \
                PLOGMSG (    \
                        klogErr,    \
                        ( klogErr, "[Exception handled] $(msg)", "msg=%s", E . what () )    \
                        );    \
                throw;      \
            }    \
            catch ( int & i ) {    \
                PLOGMSG (    \
                        klogErr,    \
                        ( klogErr, "[Exception handled] error code: $(err_c)", "err_c=%d", i )    \
                        );    \
                throw;      \
            }    \
            catch ( ... ) {    \
                PLOGMSG (    \
                        klogErr,    \
                        ( klogErr, "[Exception handled] unknown", "severity=error" )    \
                        );    \
                throw;      \
            }

    /*  Catches, reports and returns RC code
     */
#define TL_CATCH_RR(RCt)  \
            catch ( TL_Exception & E ) {    \
                RCt = E . rc ();    \
                PLOGERR (    \
                        klogErr,    \
                        ( klogErr, RCt, "[Exception handled] $(msg)", "msg=%s", E . what () )    \
                        );    \
            }    \
            catch ( exception & E ) {    \
                RCt = 1;    \
                PLOGERR (    \
                        klogErr,    \
                        ( klogErr, RCt, "[Exception handled] $(msg)", "msg=%s", E . what () )    \
                        );    \
            }    \
            catch ( int & i ) {    \
                RCt = i;    \
                PLOGERR (    \
                        klogErr,    \
                        ( klogErr, 1, "[Exception handled] error code: $(err_c)", "err_c=%d", i )    \
                        );    \
            }    \
            catch ( ... ) {    \
                RCt = 1;    \
                PLOGERR (    \
                        klogErr,    \
                        ( klogErr, RCt, "[Exception handled] unknown", "severity=error" )    \
                        );    \
            }

    /*  Catches, reports, dispose and trhow
     */
#define TL_CATCH_RDT(DispFn)  \
            catch ( TL_Exception & E ) {    \
                PLOGERR (    \
                        klogErr,    \
                        ( klogErr, E . rc (), "[Exception handled] $(msg)", "msg=%s", E . what () )    \
                        );    \
                DispFn ();  \
                throw;      \
            }    \
            catch ( exception & E ) {    \
                PLOGMSG (    \
                        klogErr,    \
                        ( klogErr, "[Exception handled] $(msg)", "msg=%s", E . what () )    \
                        );    \
                DispFn ();  \
                throw;      \
            }    \
            catch ( int & i ) {    \
                PLOGMSG (    \
                        klogErr,    \
                        ( klogErr, "[Exception handled] error code: $(err_c)", "err_c=%d", i )    \
                        );    \
                DispFn ();  \
                throw;      \
            }    \
            catch ( ... ) {    \
                PLOGMSG (    \
                        klogErr,    \
                        ( klogErr, "[Exception handled] unknown", "severity=error" )    \
                        );    \
                DispFn ();  \
                throw;      \
            }

#endif /* _tl_vrules_ */

