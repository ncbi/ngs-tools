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

#ifndef _tl_types_
#define _tl_types_

#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <algorithm>

namespace _tl_ {

typedef std::map < std::string, std::string > TL_SSMap;
typedef std::map < std::string, int > TL_SIMap;
typedef std::map < int, std::string > TL_ISMap;
typedef std::map < int, int > TL_IIMap;
typedef std::list < std::string > TL_SList;
typedef std::vector < std::string > TL_SVec;
typedef std::vector < TL_SVec > TL_SVVec;    // will think about that
typedef std::vector < int > TL_IVec;
typedef std::set < std::string > TL_SSet;

/*  Case Insensitive comparator, string to string */
struct TL_SSCompare {
    bool
    operator () (
                const std :: string & Left,
                const std :: string & Right
    ) const
    {
        std :: string Lstr = Left;
        std :: transform ( Lstr . begin (), Lstr . end (), Lstr . begin (), :: tolower );
        std :: string Rstr = Right;
        std :: transform ( Rstr . begin (), Rstr . end (), Rstr . begin (), :: tolower );
        return Lstr < Rstr;
    };
};

/*  Case Insencitive containers ... STRING2STRING
*/
typedef std::map < std::string, std::string, TL_SSCompare > TL_CISSMap;
typedef std::map < std::string, int, TL_SSCompare > TL_CISIMap;
typedef std::set < std::string, TL_SSCompare > TL_CISSet;


}   /* namespace _tl_ */

#define TL_BOOL_TRUE           "true"
#define TL_BOOL_FALSE          "false"

#define TL_NAME_TAG            "name"
#define TL_NAME_DEF            ""

#define TL_MANDATORY_TAG       "mandatory"
#define TL_MANDATORY_DEF       false

#define TL_LOG_POS_TAG         "log_pos"
#define TL_LOG_POS_DEF         false

#define TL_DESCRIPTION_TAG     "description"
#define TL_DESCRIPTION_DEF     ""

#define TL_COMMON_TAG          "common"
#define TL_COMMON_DEF          false

#define TL_CB_IGNORED_TAG      "can_be_ignored"
#define TL_CB_IGNORED_DEF      false

#define TL_CB_SYNONYM_TAG      "can_be_synonym"
#define TL_CB_SYNONYM_DEF      false

#define TL_RFC_TAG             "rfc"
#define TL_RFC_DEF             false

#define TL_TABLE_NAME_TAG      "table_name"
#define TL_TABLE_NAME_DEF      ""

#define TL_FIELD_PRIM_TAG      "field_prim"
#define TL_FIELD_PRIM_DEF      ""

#define TL_FIELD_NORM_TAG      "field_norm"
#define TL_FIELD_NORM_DEF      ""

#define TL_STATIC_TAG          "static"
#define TL_STATIC_DEF          false

#define TL_SEARCHABLE_TAG      "searchable"
#define TL_SEARCHABLE_DEF      false

#define TL_FOR_DUMP_TAG        "for_dump"
#define TL_FOR_DUMP_DEF        false

#define TL_FIELD_POS_TAG       "fiedl_pos"
#define TL_FIELD_POS_DEF       0

#define TL_VDB_NAME_TAG        "vdb_name"
#define TL_VDB_NAME_DEF        ""

#define TL_VDB_SIZE_TAG        "vdb_size"
#define TL_VDB_SIZE_DEF        sizeof ( char )

#define TL_VDB_TYPE_TAG        "vdb_type"
#define TL_VDB_TYPE_DEF        "ascii"

#define TL_VDB_TYPE_CAST_TAG   "vdb_type_cast"
#define TL_VDB_TYPE_CAST_DEF   ""

#define TL_VDB_COMPRESSION_TAG "vdb_compression"
#define TL_VDB_COMPRESSION_DEF "zip_encoding"

#define TL_VDB_DEF_VALUE_TAG   "vdb_default_value"
#define TL_VDB_DEF_VALUE_DEF   ""

#define TL_COMMENT_TAG         "comment"
#define TL_COMMENT_DEF         ""

#define TL_DEPRECATED_TAG      "deprecated"
#define TL_DEPRECATED_DEF      false

#endif /* _tl_types_ */

