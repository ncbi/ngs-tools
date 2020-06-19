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

#ifndef _tl_tracefields_
#define _tl_tracefields_

#include <string>
#include <map>
#include <vector>

#include "tl_types.hpp"
#include "tl_owp.hpp"
#include "tl_exception.hpp"
#include "tl_util.hpp"

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * Lyrics: That is "hardcoded" table TraceMain..Main, which was created
 * by Vladimir Alekseyev for submission field validation.
 * TL_TraceFields is collection of TL_TraceField class instances, and
 * each TL_TraceField class has some members, which describes actual
 * column in TraceMain..Trace table.
 *
 * During the Redesign, we are going to rid off dictionary fields,
 * so some members of TL_TraceField like table_name, or field_prime
 * are getting obsolete. Howvever, I'll leaving that information for
 * future generations ... commented
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

namespace _tl_ {

    /*  No setters for that class. Sorry guys :D
     */
class TL_TraceField {
public:
    TL_TraceField ();
    TL_TraceField ( const TL_TraceField & Field );
    ~TL_TraceField ();

    TL_TraceField & operator = ( const TL_TraceField & Field );
    TL_TraceField & operator = ( const TL_Owp & Desc );

    TL_Owp Compile () const;
    bool IsEmpty () const { return _owp . IsEmpty (); };

protected:
    void _reset ();
        /*  Use that method to initialize field object
         */
    void _set ( const TL_Owp & Owp );

private:
    void _validate () const;

    TL_Owp _owp;

public:
        /*  I left names of fields almost the same as Vladimir
         *  Alekseyev called them. I put here his comments also
         */

        /*  column name in the RFC
         */
    inline const std::string &Name () const { return _name; };

        /*  mandatory field or not  
            ( former presence )
         */
    inline bool IsMandatory () const { return _mandatory; };

        /*  if to log, which position to put into
         */
    inline int LogPos () const { return _log_pos; };

        /*  meaning of the field
         *  ( That one was deep )
         */
    inline const std::string & Description () const { return _description; };

        /*  can be specified under common section
         */
    inline bool IsCommon () const { return _common; };

        /*  can be ignored during parsing
         */
    inline bool CanBeIgnored () const { return _can_be_ignored; };

        /*  That does not have comment
         */
    inline bool CanBeSynonym () const { return _can_be_synonym; };

        /*  ( My comment ) Filed listed in RFC document
         */
    inline bool IsRfc () const { return _rfc; };

protected:
    std::string _name;
    bool _mandatory;

    int _log_pos;
    std::string _description;

    bool _common;
    bool _can_be_ignored;
    bool _can_be_synonym;
    bool _rfc;

    /* These are obsolete fields
     */
public:

        /*  table wich normalizes it
         */
    inline const std::string & TableName () const { return _table_name; };

        /*  primary filed name in Trace table
         */
    inline const std::string & FieldPrim () const { return  _field_prim; };

        /*  normailzed field (filed with the actual value)
         */
    inline const std::string & FieldNorm () const { return  _field_norm; };

        /*  what vocabulary to use ? static
         */
    inline bool IsStatic () const { return _static; };

        /*  allow to search by this field
         *  ( My note ) since VDB, no searches allowed
         */
    inline bool IsSearchable () const { return _searchable; };

        /*  There is no comment on that account, so ... obsolete
         */
    inline bool IsForDump () const { return _for_dump; };

    inline int FieldPos () const { return _field_pos; };

protected:
    std::string _table_name;
    std::string _field_prim;
    std::string _field_norm;
    bool _static;
    bool _searchable;
    bool _for_dump;
    int _field_pos;


    /* These are VDB related fields
     */
public:
        /*  name of field in VDB
         */
    inline const std::string & VdbName () const { return _vdb_name; };
        /*  size of field
         */
    inline size_t VdbSize () const { return _vdb_size; };

        /*  type of field ... std::string ... hmmm
         */
    inline const std::string & VdbType () const { return _vdb_type; };
        /*  elem_bits for GeneralWriter :: addColumn ()
         */

        /*  type to cast field
         */
    inline const std::string & VdbTypeCast () const { return _vdb_type_cast; };
        /*  column compression type
         */
    inline const std::string & VdbCompression () const { return _vdb_compression; };
        /*  column default value to set for writer
         */
    inline const std::string & VdbDefaultValue () const { return _vdb_default_value; };
        /*  comments used in schema file
         */
    inline const std::string & Comment () const { return _comment; };

protected:
    std::string _vdb_name;
    size_t _vdb_size;
    std::string _vdb_type;
    std::string _vdb_type_cast;
    std::string _vdb_compression;
    std::string _vdb_default_value;
    std::string _comment;

    /* These are accession and other flag fields
     */
public:
        /*  if field was deprecated
         */
    inline bool IsDeprecated () const { return _deprecated; };

protected:
    bool _deprecated;

};

typedef std::vector < TL_TraceField > TL_TFVec;
typedef TL_TFVec :: const_iterator TL_TFVecI;
typedef std::map < std::string, TL_TraceField, TL_SSCompare > TL_TFMap;
typedef TL_TFMap :: const_iterator TL_TFMapI;

class TL_TraceFields {

private:
    TL_TraceFields ();

public:
    static TL_TFMapI Begin ();
    static TL_TFMapI End ();

    static bool Has ( const std::string & FieldName );
    static bool Has ( const TL_TraceField & Field );

    static const TL_TraceField & Find ( const std::string & FieldName );

    static void Dump ();
};

}   /* namespace _tl_ */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *  Something for free.
 *  We are introducing new VDB TYPE for local time. But, we do avoid
 *  changing schema's headers. So, I put it here. BTW, full list of
 *  VDB TYPEs You ma find here : ncbi-vdb/interfaces/sra/types.h
 *
 *  Why: LOAD_DATE for submission is recorded as a localtime, which is
 *  sometime is day saving. So, we just convert localt time to GMT, and
 *  store it as uint64_t ( almost time_t ), just for sure
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

#define vdb_loctm_t     "loctm"     /* apparently uint64_t - 5HR  */

#endif /* _tl_tracefields_ */

