// $Id: trace_database.hpp,v 2.41 2009/07/29 20:05:30 yaschenk Exp $
// $Date: 2009/07/29 20:05:30 $
//
#ifndef _tl_names_
#define _tl_names_

#include "tl_tracefields_init.hpp"

namespace _tl_ {

#define _TL_NULL_ENTRY                   "null"

#define _TL_SFF                          "sff"

#define _TL_CLIP_ADAPTER_LEFT            "clip_adapter_left"
#define _TL_CLIP_ADAPTER_RIGHT           "clip_adapter_right"
#define _TL_NAME                         "NAME"
#define _TL_PLATFORM                     "PLATFORM"
//-- deprecated #define _TL_READ_FILTER                  "READ_FILTER"
//-- deprecated #define _TL_SPOT_FILTER                  "SPOT_FILTER"
//-- deprecated #define _TL_TRACE_FILTER                 "TRACE_FILTER"
#define _TL_READ_TYPE                    "READ_TYPE"
#define _TL_POSITION                     "POSITION"
#define _TL_QUALITY                      "QUALITY"
// #define _TL_QUALITY                      "ORIGINAL_QUALITY"
#define _TL_READ                         "READ"
#define _TL_READ_LEN                     "READ_LEN"
#define _TL_SIGNAL                       "SIGNAL"
#define _TL_SIGNAL_LEN                   "SIGNAL_LEN"

#define _TL_SUBMISSION_NEW               "new"
#define _TL_SUBMISSION_UPDATE            "update"

#define _TL_SUBMISSION_WITHDRAW          "withdraw"

#define _TL_SUBMISSION_UPDATE_INFO       "updateinfo"
#define _TL_SUBMISSION_UPDATE_INFO1      "update info"

#define _TL_SUBMISSION_NEW_TRACELESS     "new traceless"
#define _TL_SUBMISSION_NEW_TRACELESS1    "newtraceless"
#define _TL_SUBMISSION_NEW_TRACELESS2    "traceless"

#define _TL_SUBMISSION_UPDATE_TRACELESS  "update traceless"
#define _TL_SUBMISSION_UPDATE_TRACELESS1 "updatetraceless"

#define _TL_STRATEGY_POOL_CLONE          "poolclone"
#define _TL_STRATEGY_POOL_CLONE1         "pool_clone"
#define _TL_STRATEGY_POOL_CLONE2         "pool clone"
#define _TL_STRATEGY_POOL_CLONE3         "poolclone"
#define _TL_TRACE_TYPE_SHOTGUN           "shotgun"

#define _TL_454                          "454"
#define _TL_HOMO_SAPIENS_                "homo sapiens"

#define _TL_NREADS                       "NREADS"
#define _TL_READ_START                   "READ_START"
#define _TL_OUT_NAME                     "_out_name"

} /* namespace _tl_ */

#endif // _tl_names_
