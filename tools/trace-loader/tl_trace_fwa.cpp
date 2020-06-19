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

#include <iostream>
#include <sstream>
#include <set>

#include <sra/types.h>
#include <insdc/sra.h>
#include <klib/log.h>

#include "tl_util.hpp"
#include "tl_names.hpp"

#include "tl_tracefields.hpp"
#include "tl_tracedata.hpp"
#include "tl_traceinfo.hpp"
#include "tl_trace_fwa.hpp"

#include <unistd.h>

using namespace std;
using namespace _tl_;
using namespace ncbi;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_StreamId - or field which comes with schema, and ...
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
namespace _tl_ {

class TL_StreamId {
public:
    TL_StreamId ();
    virtual ~TL_StreamId ();

    inline bool IsEmpty () const { return _name . empty (); };

    inline const string & Name () const { return _name; };
    inline int32_t StreamId () const { return _stream_id; };
    inline size_t ElemBits () const { return _elem_bits; };
    inline bool DefaultSet () const { return _default_set; };

    virtual void Add (
                    int TableId,
                    GeneralWriter * Writer,
                    TL_TraceFieldAdapter :: SMap & Map
                    ) = 0;
    virtual void SetDefault ( GeneralWriter * Writer ) = 0;
    virtual void Submit (
                    GeneralWriter * Writer,
                    const string & Data
                    ) = 0;
    virtual void Submit ( 
                    GeneralWriter * Writer,
                    const void * Data,
                    uint32_t ElemCount
                    ) = 0;
        /*  Do not know who need that method :LOL:
         */
    virtual bool IsAlias () const;

    virtual bool IsExternal () const;

protected:
    inline void _setName ( const string & Name )
                    { _name = TL_StringU :: ToUpper ( Name ); };
    inline void _setStreamId ( int32_t StreamId )
                    { _stream_id = StreamId; };
    inline void _setElemBits ( int32_t ElemBits )
                    { _elem_bits = ElemBits; };
    inline void _setDefaultSet ( bool DefaultSet )
                    { _default_set = DefaultSet; };
    inline void _addAlias ( const string & Alias )
                    { _aliases . insert ( TL_StringU :: ToUpper ( Alias ) ); }
    inline TL_SSet :: const_iterator _aliasesBegin ()
                    { return _aliases . begin (); };
    inline TL_SSet :: const_iterator _aliasesEnd ()
                    { return _aliases . end (); };
private:
    string _name;
    int32_t _stream_id;
    size_t _elem_bits;
    bool _default_set;
    TL_SSet _aliases;
};

}; /* namespace _tl_ */

TL_StreamId :: TL_StreamId ()
:   _name ( "" )
,   _stream_id ( 0 )
,   _elem_bits ( 8 )
,   _default_set ( false )
,   _aliases ()
{
}   /* TL_StreamId :: TL_StreamId () */

TL_StreamId :: ~TL_StreamId ()
{
    _name = "";
    _stream_id = 0;
    _elem_bits = 8;
    _default_set = false;
    _aliases . clear ();
}   /* TL_StreamId :: ~TL_StreamId () */

void
TL_StreamId :: Add (
                    int TableId,
                    GeneralWriter * Writer,
                    TL_TraceFieldAdapter :: SMap & Map
)
{
    throw TL_Exception ( "Calling unimplemented TL_StreamId :: Add ()" );
}   /* TL_StreamId :: Add () */

void
TL_StreamId :: SetDefault ( GeneralWriter * Writer )
{
    throw TL_Exception ( "Calling unimplemented TL_StreamId :: SetDefault ()" );
}   /* TL_StreamId :: SetDefault () */

void
TL_StreamId :: Submit ( GeneralWriter * Writer, const string & Data )
{
    throw TL_Exception ( "Calling unimplemented TL_StreamId :: Submit ()" );
}   /* TL_StreamId :: Submit () */

void
TL_StreamId :: Submit (
                    GeneralWriter * Writer,
                    const void * Data,
                    uint32_t ElemCount
)
{
    throw TL_Exception ( "Calling unimplemented TL_StreamId :: Submit ()" );
}   /* TL_StreamId :: Submit () */

bool
TL_StreamId :: IsAlias () const
{
    return true;
}   /* TL_StreamId :: IsAlias () */

bool
TL_StreamId :: IsExternal () const
{
    return false;
}   /* TL_StreamId :: IsExternal () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _AliasStreamId - stupid class which should be replaced with
 *_                  refcounts
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _AliasStreamId : public TL_StreamId {
public :
    _AliasStreamId ( TL_StreamId * Aliased );
    ~_AliasStreamId ();

    void Add (
            int TableId,
            GeneralWriter * Writer,
            TL_TraceFieldAdapter :: SMap & Map
            );
    void SetDefault ( GeneralWriter * Writer );
    void Submit (
            GeneralWriter * Writer,
            const string & Data
            );
    void Submit ( 
            GeneralWriter * Writer,
            const void * Data,
            uint32_t ElemCount
            );
    bool IsAlias () const;

private:
    TL_StreamId * _aliased;
};

_AliasStreamId :: _AliasStreamId ( TL_StreamId * Aliased )
:   TL_StreamId ()
,   _aliased ( Aliased )
{
}   /* _AliasStreamId :: _AliasStreamId () */

_AliasStreamId :: ~_AliasStreamId ()
{
    _aliased = NULL;
}   /* _AliasStreamId :: ~_AliasStreamId () */

void
_AliasStreamId :: Add (
                        int TableId,
                        GeneralWriter * Writer,
                        TL_TraceFieldAdapter :: SMap & Map
)
{
    if ( _aliased != NULL ) {
        _aliased -> Add ( TableId, Writer, Map );
    }
}   /* _AliasStreamId :: _AliasStreamId () */

void
_AliasStreamId :: SetDefault ( GeneralWriter * Writer )
{
    if ( _aliased != NULL ) {
        _aliased -> SetDefault ( Writer );
    }
}   /* _AliasStreamId :: _AliasStreamId () */

void
_AliasStreamId :: Submit ( GeneralWriter * Writer, const string & Data )
{
    if ( _aliased != NULL ) {
        _aliased -> Submit ( Writer, Data );
    }
}   /* _AliasStreamId :: _AliasStreamId () */

void
_AliasStreamId :: Submit (
                        GeneralWriter * Writer,
                        const void * Data,
                        uint32_t ElemCount
)
{
    if ( _aliased != NULL ) {
        _aliased -> Submit ( Writer, Data, ElemCount );
    }
}   /* _AliasStreamId :: Submit () */

bool
_AliasStreamId :: IsAlias () const
{
    return true;
}   /* _AliasStreamId :: IsAlias () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _DummyStreamId - dummy class for fields, which are internal but
 *_                  there should not be any action on them, like
 *_                  TRACE_FILE, BASES_FILE, etc
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _DummyStreamId : public TL_StreamId {
public:
    _DummyStreamId ( const string & Name );

    void Add (
            int TableId,
            GeneralWriter * Writer,
            TL_TraceFieldAdapter :: SMap & Map
            );
    void SetDefault ( GeneralWriter * Writer );
    void Submit (
            GeneralWriter * Writer,
            const string & Data
            );
    void Submit ( 
            GeneralWriter * Writer,
            const void * Data,
            uint32_t ElemCount
            );
};

_DummyStreamId :: _DummyStreamId ( const string & Name )
:   TL_StreamId ()
{
    _setName ( Name );
}   /* _DummyStreamId :: _DummyStreamId () */

void
_DummyStreamId :: Add (
                    int TableId,
                    GeneralWriter * Writer,
                    TL_TraceFieldAdapter :: SMap & Map
)
{
    string TheName = Name ();
    if ( TheName . empty () ) {
            /*  Not sure if it is error or something else
             */
        throw TL_Exception ( "Attempting to add column with empty name" );
    }

    if ( StreamId () != 0 ) {
        delete this;    /* It was allocated */

        return;     /* Already set */
    }

    if ( Map . find ( TheName ) != Map . end () ) {
        delete this;    /* It was allocated */

        return;     /* something wrong, but it is not an error.
                     * Prolly I should report about that situation
                     * JOJOBA
                     */
    }

    int SId = 333 + 222 + 111;

    _setStreamId ( SId );
    if ( SId != 0 ) {
        Map [ TheName ] = this;

        for ( TL_SSet :: const_iterator It = _aliasesBegin ();
                It != _aliasesEnd ();
                It ++
        ) {
            string DosEquis = TL_StringU :: ToUpper ( * It );
            if ( Map . find ( DosEquis ) == Map . end () ) {
                Map [ DosEquis ] = new _AliasStreamId ( this );
            }
        }
    }
    /* Nothing to do here ... it is dummy class */
}   /* _DummyStreamId :: Add () */

void
_DummyStreamId :: SetDefault ( GeneralWriter * W )
{
    /* Hoorray! Holiday! Nothing to do here :D */
}   /* _DummyStreamId :: Submit () */

void
_DummyStreamId :: Submit ( GeneralWriter * W, const string & D )
{
    /* Hoorray! Holiday! Nothing to do here :D */
}   /* _DummyStreamId :: Submit () */

void
_DummyStreamId :: Submit ( GeneralWriter * W, const void * D, uint32_t E)
{
    /* Hoorray! Holiday! Nothing to do here :D */
}   /* _DummyStreamId :: Submit () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _VdbStreamId - simple field with defined VDB type
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _VdbStreamId : public TL_StreamId {
public:
    static bool CalculateIsInteger ( const string & VdbType );
    static size_t CalculateElemBits ( const string & VdbType );

public:
    _VdbStreamId ( const string & VdbType );
    ~_VdbStreamId ();

    void Add (
            int TableId,
            GeneralWriter * Writer,
            TL_TraceFieldAdapter :: SMap & Map
            );
    void SetDefault ( GeneralWriter * Writer );
    void Submit (
            GeneralWriter * Writer,
            const string & Data
            );
    void Submit (
            GeneralWriter * Writer,
            const void * Data,
            uint32_t ElemCount
            );

private:
    void _set ( const string & VdbType );

private:
    string _vdb_type;

    bool _is_integer;
};

_VdbStreamId :: _VdbStreamId ( const string & VdbType )
:   TL_StreamId ()
,   _vdb_type ( "" )
,   _is_integer ( false )
{
    _set ( VdbType );
}   /* _VdbStreamId :: _VdbStreamId () */

_VdbStreamId :: ~_VdbStreamId ()
{
    _vdb_type = "";
    _is_integer = false;
}   /* _VdbStreamId :: ~_VdbStreamId () */

void
_VdbStreamId :: Add (
                    int TableId,
                    GeneralWriter * Writer,
                    TL_TraceFieldAdapter :: SMap & Map
)
{
    string TheName = Name ();
    if ( TheName . empty () ) {
            /*  Not sure if it is error or something else
             */
        throw TL_Exception ( "Attempting to add column with empty name" );
    }

    if ( StreamId () != 0 ) {
        delete this;    /* to avoid leaks */

        return;     /* Already set */
    }

    if ( Map . find ( TheName ) != Map . end () ) {
        delete this;    /* to avoid leaks */

        return;     /* something wrong, but it is not an error.
                     * Prolly I should report about that situation
                     * JOJOBA
                     */
    }

    int32_t SId = _is_integer 
                            ? Writer -> addIntegerColumn (
                                        TableId, TheName, ElemBits ()
                                        )
                            : Writer -> addColumn (
                                        TableId, TheName, ElemBits ()
                                        )
                            ;
    _setStreamId ( SId );
    if ( SId != 0 ) {
        Map [ TheName ] = this;

        for ( TL_SSet :: const_iterator It = _aliasesBegin ();
                It != _aliasesEnd ();
                It ++
        ) {
            string DosEquis = TL_StringU :: ToUpper ( * It );
            if ( Map . find ( DosEquis ) == Map . end () ) {
                Map [ DosEquis ] = new _AliasStreamId ( this );
            }
        }
    }
}   /* _VdbStreamId :: Add () */

void
_VdbStreamId :: SetDefault ( GeneralWriter * Writer )
{
    if ( ! DefaultSet () ) {
        if ( 0 < StreamId () ) {
            Writer -> columnDefault ( StreamId (), ElemBits (), "", 0 );
            _setDefaultSet ( true );
        }
    }
}   /* _VdbStreamId :: SetDefault () */

/*  Originally LOAD_DATE is coming as LocalTime, and we need to convert
 *  it to GMT.
 */
void _submit_localtime (
                        GeneralWriter * Writer,
                        int StreamId,
                        const string & Data
)
{
    uint64_t LocalTime = 0;
    uint64_t GMTime = 0;

    istringstream Out ( Data );
    Out >> LocalTime;
    if ( ! Out ) {
        throw TL_Exception ( string ( "Can not convert '" ) + Data + "' to LocalTime, smh" );
    }

    GMTime = TL_DateU :: LocalToGMT ( LocalTime );

    Writer -> write ( StreamId, sizeof ( uint64_t ) * 8, & GMTime, 1 );
}   /* _submit_localtime () */

template < class Type >
void _submit_data (
                GeneralWriter * Writer,
                int StreamId,
                const string & Data
)
{
    Type Var;

    istringstream Out ( Data );
    Out >> Var;
    if ( ! Out ) {
        throw TL_Exception ( string ( "Can not convert '" ) + Data + "' to integer, smh" );
    }

    Writer -> write ( StreamId, sizeof ( Type ) * 8, & Var, 1 );
}   /* template _submit_data () */

void
_VdbStreamId :: Submit ( GeneralWriter * Writer, const string & Data )
{
    if ( Name () . empty () ) {
        return;     /* Apparently is not an error */
    }

    int32_t SId = StreamId ();
    if ( SId <= 0 ) {
        return;     /* Apparently is not an error */
    }

    if ( _vdb_type . empty () ) {
        return;     /* Apparently is not an error */
    }

    if ( _vdb_type == vdb_uint64_t ) {
        _submit_data < uint64_t > ( Writer, SId, Data );
        return;
    }

    if ( _vdb_type == vdb_int64_t ) {
        _submit_data < int64_t > ( Writer, SId, Data );
        return;
    }

    if ( _vdb_type == vdb_uint32_t ) {
        _submit_data < uint32_t > ( Writer, SId, Data );
        return;
    }

    if ( _vdb_type == vdb_int32_t ) {
        _submit_data < int32_t > ( Writer, SId, Data );
        return;
    }

    if ( _vdb_type == vdb_uint16_t ) {
        _submit_data < uint16_t > ( Writer, SId, Data );
        return;
    }

    if ( _vdb_type == vdb_int16_t ) {
        _submit_data < int16_t > ( Writer, SId, Data );
        return;
    }

    if ( _vdb_type == vdb_uint8_t ) {
        _submit_data < uint8_t > ( Writer, SId, Data );
        return;
    }

    if ( _vdb_type == vdb_int8_t ) {
        _submit_data < int8_t > ( Writer, SId, Data );
        return;
    }

    if ( _vdb_type == vdb_float32_t ) {
        float Val;
        TL_StringU :: ToBin ( Data, & Val, sizeof ( Val ) );
        Writer -> write ( SId, sizeof ( Val ) * 8, & Val, 1 );
        return;
    }

    if ( _vdb_type == vdb_float64_t ) {
        double Val;
        TL_StringU :: ToBin ( Data, & Val, sizeof ( Val ) );
        Writer -> write ( SId, sizeof ( Val ) * 8, & Val, 1 );
        return;
    }

    if ( _vdb_type == vdb_loctm_t ) {
        _submit_localtime ( Writer, SId, Data );
        return;
    }

    if ( _vdb_type == vdb_ascii_t ) {
        Writer -> write (
                        SId,
                        ElemBits (),
                        Data . c_str (),
                        Data . size ()
                        );
        return;
    }

        /*  JOJOBA don't know what to do here ... may be that is an error
         */
    throw TL_Exception ( string ( "Submit () : Unknown format for filed '" ) + Name () + "'" );
}   /* _VdbStreamId :: Submit () */

void
_VdbStreamId :: Submit (
                    GeneralWriter * Writer,
                    const void * Data,
                    uint32_t ElemCount
)
{
    if ( Name () . empty () ) {
        return;     /* Apparently is not an error */
    }

    int32_t SId = StreamId ();
    if ( SId <= 0 ) {
        return;     /* Apparently is not an error */
    }

    if ( _vdb_type . empty () ) {
        return;     /* Apparently is not an error */
    }

    Writer -> write ( SId, ElemBits (), Data, ElemCount );

}   /* _VdbStreamId :: Submit () */

void
_VdbStreamId :: _set ( const string & VdbType )
{
    _vdb_type = VdbType;
    _is_integer = CalculateIsInteger ( VdbType );
    _setElemBits ( CalculateElemBits ( VdbType ) );
}   /* _VdbStreamId :: _set () */

bool
_VdbStreamId :: CalculateIsInteger ( const string & Type )
{
    if (    Type == vdb_uint64_t || Type == vdb_int64_t
        ||  Type == vdb_uint32_t || Type == vdb_int32_t
        ||  Type == vdb_uint16_t || Type == vdb_int16_t
        ||  Type == vdb_uint8_t || Type == vdb_int8_t
        ||  Type == vdb_loctm_t
        ) {
        return true;
    }

    return false;
}   /* _VdbStreamId :: CalculateIsInteger () */

size_t
_VdbStreamId :: CalculateElemBits ( const string & Type )
{
        /* All these constants taken from sra/types.h file
         */

        /*  Integers
         */
    if ( Type == vdb_uint64_t ) return sizeof ( uint64_t ) * 8;
    if ( Type == vdb_int64_t ) return sizeof ( int64_t ) * 8;

    if ( Type == vdb_uint32_t ) return sizeof ( uint32_t ) * 8;
    if ( Type == vdb_int32_t ) return sizeof ( int32_t ) * 8;

    if ( Type == vdb_uint16_t ) return sizeof ( uint16_t ) * 8;
    if ( Type == vdb_int16_t ) return sizeof ( int16_t ) * 8;

    if ( Type == vdb_uint8_t ) return sizeof ( uint8_t ) * 8;
    if ( Type == vdb_int8_t ) return sizeof ( int8_t ) * 8;

        /*  Floats
         */
    if ( Type == vdb_float32_t ) return sizeof ( float ) * 8;
    if ( Type == vdb_float64_t ) return sizeof ( double ) * 8;

        /*  LocalTime
         */
    if ( Type == vdb_loctm_t ) return sizeof ( uint64_t ) * 8;

        /*  Chars and other
         */
    return 8;
}   /* _VdbStreamId :: CalculateElemBits () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _IrregularStreamId - non-mandatory field which data is not external
 *_                      and comes from anywhere, or derivative
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _IrregularStreamId : public _VdbStreamId {
public: 
    _IrregularStreamId ( const string & Name, const string & VdbType );
    ~_IrregularStreamId ();

    bool IsExternal () const;
};

_IrregularStreamId :: _IrregularStreamId (
                                        const string & Name,
                                        const string & VdbType
)
:   _VdbStreamId ( VdbType )
{
        /*  Not sure if empty name is an error
         */
    _setName ( Name );

    _setStreamId ( 0 );

    _setDefaultSet ( false );
}   /* _IrregularStreamId :: _IrregularStreamId () */

_IrregularStreamId :: ~_IrregularStreamId ()
{
}   /* _IrregularStreamId :: ~_IrregularStreamId () */

bool
_IrregularStreamId :: IsExternal () const
{
    return false;
}   /* _IrregularStreamId :: IsExternal () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _RegularStreamId - non-mandatory field which data comes from
 *_                    external files, and TRACEINFO.tbl
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _RegularStreamId : public _VdbStreamId {
public :
    _RegularStreamId ( const TL_TraceField & Field );
    ~_RegularStreamId ();

    bool IsExternal () const;

private :
    void _set ( const TL_TraceField & Field );

};

_RegularStreamId :: _RegularStreamId ( const TL_TraceField & Field )
:   _VdbStreamId ( Field . VdbType () )
{
    _set ( Field );
}   /* _RegularStreamId :: _RegularStreamId () */

_RegularStreamId :: ~_RegularStreamId ()
{
}   /* _RegularStreamId :: ~_RegularStreamId () */

void
_RegularStreamId :: _set ( const TL_TraceField & Field )
{
    string Name = Field . Name ();

    if ( Name . empty () ) {
        /*  JOJOBA Prolly that is error 
         */
        return;
    }

    _setName ( Name );
    _setStreamId ( 0 );

    _setDefaultSet ( false );

}   /* _RegularStreamId :: _set () */

bool
_RegularStreamId :: IsExternal () const
{
    return true;
}   /* _RegularStreamId :: IsExternal () */

/******************************************************************
 *
 * Here we are ... we should add set of mandatory columns, and
 * set for them correct default valuse for sure. Thare list:
 *
 * CLIP_QUALITY_LEFT & CLIP_QUALITY_RIGHT - come with ZTR and SFF
 *                                          files. U32 type
 *                                          Added only if exists
 * CLIP_ADAPTER_LEFT & CLIP_ADAPTER_RIGHT - come with ZTR and SFF
 *                                          files. U32 type
 *                                          Added only if exists
 * CLIP_VECTOR_LEFT & CLIP_VECTOR_RIGHT - alias for 'CLIP_QUALITY'
 *
 * FLOW_CHARS & KEY_SEQUENCE - come with SFF file. ASCII type ( CHAR )
 *                             added only if exists 
 * NAME - mandatory field, and should be written ASCII type
 * 
 * PLATFORM - should be set as default value SRA_PLATFORM_CAPILLARY
 *            U8 type.
 * POSITION - peak_index is coming here ... U32
 * QUALITY - qualities PHRED. U8
 * READ - basecall. ASCII ( CHAR ) type
 * READ_FILTER - should be set default value SRA_READ_FILTER_PASS U8
 * READ_LEN - default #, #, #, # ... Usually it is "1, NUM_BASES, 0, 0"
 * READ_TYPE - should be set default value SRA_READ_TYPE_BIOLOGICAL
 * SIGNAL - here we are writing 4 columns of samples
 * SIGNAL_LEN - length of signal column
 *
 ******************************************************************/

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _ClipStreamId - mandatory field which data comes from both:
 *_                 external files, and trace file. If it is not
 *_                 set it should be set into '0' ( default )
 *_                 it is int32_t and defined in 'ncbi/clip.vschema'
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _ClipStreamId : public _VdbStreamId {
public:
    _ClipStreamId ( bool Left, bool Qual );

    void SetDefault ( GeneralWriter * Writer );

    bool IsExternal () const;

private:
    void _set ( bool Left, bool Qual );
};

_ClipStreamId :: _ClipStreamId ( bool Left, bool Qual )
:   _VdbStreamId ( vdb_int32_t )
{
    _set ( Left, Qual );
}   /* _ClipStreamId :: _ClipStreamId () */


void
_ClipStreamId :: SetDefault ( GeneralWriter * Writer )
{
    if ( ! DefaultSet () ) {
        if ( 0 < StreamId () ) {
            int32_t Var = 0;
            Writer -> columnDefault (
                                    StreamId (),
                                    ElemBits (),
                                    & Var,
                                    1
                                    );
            _setDefaultSet ( true );
        }
    }
}   /* _ClipStreamId :: SetDefault () */

void
_ClipStreamId :: _set ( bool Left, bool Qual )
{
    if ( Left ) {
        if ( Qual ) {
            _setName ( _TL_CLIP_QUALITY_LEFT );
            // JOJOBA _addAlias ( _TL_CLIP_VECTOR_LEFT );
        }
        else {
            _setName ( _TL_CLIP_ADAPTER_LEFT );
        }
    }
    else {
        if ( Qual ) {
            _setName ( _TL_CLIP_QUALITY_RIGHT );
            // JOJOBA _addAlias ( _TL_CLIP_VECTOR_RIGHT );
        }
        else {
            _setName ( _TL_CLIP_ADAPTER_RIGHT );
        }
    }

    _setStreamId ( 0 );
    _setElemBits ( 32 );
    _setDefaultSet ( false );
}   /* _ClipStreamId :: _set () */

bool
_ClipStreamId :: IsExternal () const
{
    return true;
}   /* _ClipStreamId :: IsExternal () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _CharStreamId - character oriented data which comes with TRACE
 *_                 files. Those are:
 *_                         FLOW_CHARS 
 *_                         KEY_SEQUENCE
 *_                         NAME ( TRACE_NAME )
 *_                         NAME ( PROGRAM_ID )
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _CharStreamId : public _VdbStreamId {
public:
    _CharStreamId ( const string & Name, const string & Alias = "" );

private:
    void _set ( const string & Name, const string & Alias );
};

_CharStreamId :: _CharStreamId (
                            const string & Name,
                            const string & Alias
)
:   _VdbStreamId ( vdb_ascii_t )
{
    _set ( Name, Alias );
}   /* _CharStreamId :: _CharStreamId () */

void
_CharStreamId :: _set ( const string & Name, const string & Alias )
{
    _setName ( Name );
    _addAlias ( Alias );
}   /* _CharStreamId :: _set () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _DateStreamId - There are two types of Date : string and time_t
 *_                 it doing both for fields :
 *_                         RUN_DATE 
 *_                         COLLECTION_DATE
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _DateStreamId : public _VdbStreamId {
public:
    _DateStreamId ( const string & Name );

    void Submit ( GeneralWriter * Writer, const string & Data );

    bool IsExternal () const;
};

_DateStreamId :: _DateStreamId ( const string & Name )
:   _VdbStreamId ( vdb_uint64_t )
{
    _setName ( Name );
}   /* _DateStreamId :: _DateStreamId () */

void
_DateStreamId :: Submit (
                GeneralWriter * Writer,
                const string & Data
                )
{
    uint64_t tTime = TL_DateU :: Date ( Data );
    if ( tTime == 0 ) {
        tTime = TL_StringU :: FromStr < uint64_t > ( Data );
    }

    _VdbStreamId :: Submit ( Writer, & tTime, 1 );
}   /* _DateStreamId :: Submit () */

bool
_DateStreamId :: IsExternal () const
{
    return true;
}   /* _DateStreamId :: IsExternal () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _ConstU8StreamId - That is fill field defined by name with some
 *_                    some default const ( U8 )
 *_ Work for:
 *_     PLATFORM -> SRA_PLATFORM_CAPILLARY (insdc/sra.h) (U8)
 *_     READ_FILTER -> SRA_READ_FILTER_PASS (insdc/insdc.vschema) (U8)
 *_     READ_TYPE -> SRA_READ_TYPE_BIOLOGICAL (insdc/insdc.vschema) (U8)
 *_     NREADS -> something strange asked by Kurt (U8)
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _ConstU8StreamId : public _VdbStreamId {
public:
    _ConstU8StreamId ( const string & Name, uint8_t Default );

    void Submit (
                GeneralWriter * Writer,
                const string & Data
                );
    void Submit ( 
                GeneralWriter * Writer,
                const void * Data,
                uint32_t ElemCount
                );
    void SetDefault ( GeneralWriter * Writer );
private:
    uint8_t _default_value;
};

_ConstU8StreamId :: _ConstU8StreamId (
                                    const string & Name,
                                    uint8_t Default
)
:   _VdbStreamId ( vdb_uint8_t )
{
    _setName ( Name );
    _default_value = Default;
}   /* _ConstU8StreamId :: _ConstU8StreamId () */

void
_ConstU8StreamId :: Submit (
                GeneralWriter * Writer,
                const string & Data
                )
{
    /* There is nothing to do ... Proceed with default value */
}   /* _ConstU8StreamId :: Submit () */

void
_ConstU8StreamId :: Submit ( 
                GeneralWriter * Writer,
                const void * Data,
                uint32_t ElemCount
                )
{
    /* There is nothing to do ... Proceed with default value */
}   /* _ConstU8StreamId :: Submit () */

void
_ConstU8StreamId :: SetDefault ( GeneralWriter * Writer )
{
    if ( ! DefaultSet () ) {
        if ( 0 < StreamId () ) {
            Writer -> columnDefault (
                                    StreamId (),
                                    ElemBits (),
                                    & _default_value,
                                    1
                                    );
            _setDefaultSet ( true );
        }
    }
}   /* _ConstU8StreamId :: SetDefault () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _ConstU32StreamId - That is fill field defined by name with some
 *_                    some default const ( U32 )
 *_ Work for:
 *_     READ_START -> something strange asked by Kurt (U8)
 *_
 *_ It is too many Const...StreamId's ... prolly I should parametrize it
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _ConstU32StreamId : public _VdbStreamId {
public:
    _ConstU32StreamId ( const string & Name, uint32_t Default );

    void Submit (
                GeneralWriter * Writer,
                const string & Data
                );
    void Submit ( 
                GeneralWriter * Writer,
                const void * Data,
                uint32_t ElemCount
                );
    void SetDefault ( GeneralWriter * Writer );
private:
    uint32_t _default_value;
};

_ConstU32StreamId :: _ConstU32StreamId (
                                    const string & Name,
                                    uint32_t Default
)
:   _VdbStreamId ( vdb_uint32_t )
{
    _setName ( Name );
    _default_value = Default;
}   /* _ConstU32StreamId :: _ConstU32StreamId () */

void
_ConstU32StreamId :: Submit (
                GeneralWriter * Writer,
                const string & Data
                )
{
    /* There is nothing to do ... Proceed with default value */
}   /* _ConstU32StreamId :: Submit () */

void
_ConstU32StreamId :: Submit ( 
                GeneralWriter * Writer,
                const void * Data,
                uint32_t ElemCount
                )
{
    /* There is nothing to do ... Proceed with default value */
}   /* _ConstU32StreamId :: Submit () */

void
_ConstU32StreamId :: SetDefault ( GeneralWriter * Writer )
{
    if ( ! DefaultSet () ) {
        if ( 0 < StreamId () ) {
            Writer -> columnDefault (
                                    StreamId (),
                                    ElemBits (),
                                    & _default_value,
                                    1
                                    );
            _setDefaultSet ( true );
        }
    }
}   /* _ConstU32StreamId :: SetDefault () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _PositionStreamId - POSITION is the same as "peak_index"
 *_                     Array of U32
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _PositionStreamId : public _VdbStreamId {
public:
    _PositionStreamId ();

    void Submit (
                GeneralWriter * Writer,
                const string & Data
                );
};

_PositionStreamId :: _PositionStreamId ()
:   _VdbStreamId ( vdb_uint32_t )
{
    _setName ( _TL_POSITION );
}   /* _PositionStreamId :: _PositionStreamId () */

void
_PositionStreamId :: Submit (
                            GeneralWriter * Writer,
                            const string & Data
)
{
    throw TL_Exception ( "That method could not be used for that stream 'POSITION')" );
}   /* _PositionStreamId :: Submit () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _QualityStreamId - QUALITY Array of U8
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _QualityStreamId : public _VdbStreamId {
public:
    _QualityStreamId ();

    void Submit (
                GeneralWriter * Writer,
                const string & Data
                );
};

_QualityStreamId :: _QualityStreamId ()
:   _VdbStreamId ( vdb_uint8_t )
{
    _setName ( _TL_QUALITY );
}   /* _QualityStreamId :: _QualityStreamId () */

void
_QualityStreamId :: Submit (
                            GeneralWriter * Writer,
                            const string & Data
)
{
    throw TL_Exception ( "That method could not be used for that stream 'QUALITY'" );
}   /* _QualityStreamId :: Submit () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _ReadStreamId - READ ( basecall ) Array of CHAR
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _ReadStreamId : public _VdbStreamId {
public:
    _ReadStreamId ();

    void Submit (
                GeneralWriter * Writer,
                const string & Data
                );
};

_ReadStreamId :: _ReadStreamId ()
:   _VdbStreamId ( vdb_ascii_t )
{
    _setName ( _TL_READ );
}   /* _ReadStreamId :: _ReadStreamId () */

void
_ReadStreamId :: Submit (
                            GeneralWriter * Writer,
                            const string & Data
)
{
    throw TL_Exception ( "That method could not be used for that stream 'READ'" );
}   /* _ReadStreamId :: Submit () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _ReadLenStreamId - READ_LEN contains 1, ##, 0, 0 ... or one 
 *_                    INSDC_coord_len (U32)
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _ReadLenStreamId : public _VdbStreamId {
public:
    _ReadLenStreamId ();
};

_ReadLenStreamId :: _ReadLenStreamId ()
:   _VdbStreamId ( vdb_uint32_t )
{
    _setName ( _TL_READ_LEN );
}   /* _ReadLenStreamId :: _ReadLenStreamId () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _SignalStreamId - SIGNAL ( traces ) Array of U16
 *_                   Right now we store these as is, but later ...
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _SignalStreamId : public _VdbStreamId {
public:
    _SignalStreamId ();

    void Submit (
                GeneralWriter * Writer,
                const string & Data
                );
};

_SignalStreamId :: _SignalStreamId ()
:   _VdbStreamId ( vdb_uint64_t )
{
    _setName ( _TL_SIGNAL );
}   /* _SignalStreamId :: _SignalStreamId () */

void
_SignalStreamId :: Submit (
                            GeneralWriter * Writer,
                            const string & Data
)
{
    throw TL_Exception ( "That method could not be used for that stream 'READ'" );
}   /* _SignalStreamId :: Submit () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _CommExtStreamId - COMMENTS and EXTENDED data are going here
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _CommExtStreamId : public _VdbStreamId {
public:
    _CommExtStreamId ( const string & Name );
};

_CommExtStreamId :: _CommExtStreamId ( const string & Name )
:   _VdbStreamId ( vdb_ascii_t )
{
    _setName ( Name );
}   /* _CommExtStreamId :: _CommExtStreamId () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _ConstU16StreamId - StreamId for stroring U16 for such fields as
 *_                     SIGNAL_LEN, TRACE_MAX_VALUE, BASES_(20, 40, 60)
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _ConstU16StreamId : public _VdbStreamId {
public:
    _ConstU16StreamId ( const string & Name );
};

_ConstU16StreamId :: _ConstU16StreamId ( const string & Name )
:   _VdbStreamId ( vdb_uint16_t )
{
    _setName ( Name );
}   /* _ConstU16StreamId :: _ConstU16StreamId () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceFieldAdapter
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_TraceFieldAdapter :: TL_TraceFieldAdapter ( const string & Output )
:   _output ( Output )
,   _writer ( NULL )
,   _table_id ( - 1 )
,   _stream_ids ()
,   _filler_capacity ( 0 )
,   _filler ( NULL )
{
}   /* TL_TraceFieldAdapter :: TL_TraceFieldAdapter () */

TL_TraceFieldAdapter :: ~TL_TraceFieldAdapter ()
{
    _Reset ();

    _filler_capacity = 0;
    if ( _filler != NULL ) {
        delete [] _filler;
    }
    _filler = NULL;
}   /* TL_TraceFieldAdapter :: ~TL_TraceFieldAdapter () */

void
TL_TraceFieldAdapter :: _Reset ()
{
    for ( SMap :: iterator It = _stream_ids . begin ();
            It != _stream_ids . end ();
            It ++
    ) {
        if ( It -> second != NULL ) {
            delete It -> second;
            It -> second = NULL;
        }
    }
    _stream_ids . clear ();

    _table_id = - 1;

    if ( _writer != NULL ) {
        delete _writer;
        _writer = NULL;
    }
}   /* TL_TraceFieldAdapter :: _Reset () */

void
TL_TraceFieldAdapter :: Init ( const string & Schema )
{
    _Reset ();

    _writer = _output . empty ()
                    ? new GeneralWriter ( STDOUT_FILENO, 128 * 1024 )
                    : new GeneralWriter ( _output )
                    ;

    _writer -> setSoftwareName ( "trace-load", "01.01.01" );
    _writer -> setRemotePath ( "TEST" );
    _writer -> useSchema ( Schema, "NCBI:sra:db:trace #1" );

    _table_id = _writer -> addTable ( "SEQUENCE" );
}   /* TL_TraceFieldAdapter :: Init () */

TL_StreamId *
TL_TraceFieldAdapter :: GetStreamId ( const string & Name, bool Throw )
{
    SMap :: iterator It = _stream_ids . find (
                                        TL_StringU :: ToUpper ( Name )
                                        );
    if ( It != _stream_ids . end () ) {
        return It -> second ;
    }

    if ( Throw ) {
        throw TL_Exception ( string ( "Can not find StreamId for '" ) + Name + "'" );
    }

    return NULL;
}   /* TL_TraceFieldAdapter :: GetStreamId () */

void
TL_TraceFieldAdapter :: AddColumn ( const TL_TraceField & Field )
{
    if ( _table_id < 0 ) {
        throw TL_Exception ( "TraceFieldAdapter :: AddColumn: table is not set" );
    }

    string ColumnName = Field . Name (); 
        /* JOJOBA !!! duplicate fields are going here, so temporary
         */
    if ( ColumnName . empty () ) {
        return;
    }

    if ( Field . IsDeprecated () ) {
        PLOGMSG (
            klogDebug,
            ( klogDebug,
            "[DEPRECATED FIELD] [$(field)]",
            "severity=debug,field=%s",
            ColumnName . c_str ()
            )
            );
        return;
    }

    if ( ! HasColumn ( ColumnName ) ) {
        ( new _RegularStreamId ( Field ) )
                            -> Add ( _table_id, _writer, _stream_ids );
    }

}   /* TL_TraceFieldAdapter :: AddColumn () */

bool
TL_TraceFieldAdapter :: HasColumn ( const string & ColumnName ) const
{
    if ( _stream_ids . find ( ColumnName ) != _stream_ids . end () ) {
        return true;
    }

    if ( _stream_ids . find ( TL_StringU :: ToUpper ( ColumnName ) ) != _stream_ids . end () ) {
        return true;
    }

    return false; 
}   /* TL_TraceFieldAdapter :: HasColumn () */

void
TL_TraceFieldAdapter :: _AddMandatoryColumns ()
{
    if ( _table_id < 0 ) {
        throw TL_Exception ( "TraceFieldAdapter :: AddColumn: table is not set" );
    }

        /*  Mandatory CLIP_XXXX_XXX columns
         */
    ( new _ClipStreamId ( true, true ) )
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _ClipStreamId ( false, true ) )
                            -> Add ( _table_id, _writer, _stream_ids );

/* We retain CLIP_VECTOR_XXXX information
    ( new _DummyStreamId ( _TL_CLIP_VECTOR_LEFT ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _DummyStreamId ( _TL_CLIP_VECTOR_RIGHT ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
*/

/* These two depricated accorting to SRA-5693
    ( new _CharStreamId ( _TL_FLOW_CHARS ) )
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _CharStreamId ( _TL_KEY_SEQUENCE ) )
                            -> Add ( _table_id, _writer, _stream_ids );
*/

// JOJOBA
    ( new _CharStreamId ( _TL_TRACE_NAME ) )
                            -> Add ( _table_id, _writer, _stream_ids );
/*
    ( new _DummyStreamId ( _TL_TRACE_NAME ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _CharStreamId ( "NAME" ) ) 
    // ( new _CharStreamId ( "name" ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
*/
// JOJOBA

    ( new _CharStreamId ( _TL_PROGRAM_ID ) )
                            -> Add ( _table_id, _writer, _stream_ids );

    ( new _ConstU8StreamId ( _TL_PLATFORM, SRA_PLATFORM_CAPILLARY ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _PositionStreamId () ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _QualityStreamId () ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _ReadStreamId () ) 
                            -> Add ( _table_id, _writer, _stream_ids );

#ifdef JOJOBA
/* by Ken request */
/* deprecated now */
    ( new _ConstU8StreamId ( _TL_READ_FILTER, SRA_READ_FILTER_PASS ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
/* Do not need that */
    ( new _ConstU8StreamId ( _TL_TRACE_FILTER, 0 ) )
                            -> Add ( _table_id, _writer, _stream_ids );
#endif /* JOJOBA */


    ( new _ReadLenStreamId () ) 
                            -> Add ( _table_id, _writer, _stream_ids );

/*  Since we using TRACE_END to synthesize name
    ( new _ConstU8StreamId ( _TL_READ_TYPE, SRA_READ_TYPE_BIOLOGICAL ) )
                            -> Add ( _table_id, _writer, _stream_ids );
 */
    ( new _IrregularStreamId ( _TL_READ_TYPE, vdb_uint8_t ) )
                            -> Add ( _table_id, _writer, _stream_ids );

#ifdef __NEW_CONTROL_FLAGS__
        /*  Since ROBERT told me that CONTROL_FLAGS supposed to be 
            corrected by me
         */
    ( new _IrregularStreamId ( _TL_CONTROL_FLAGS, vdb_uint16_t ) )
                            -> Add ( _table_id, _writer, _stream_ids );
#endif /* __NEW_CONTROL_FLAGS__ */

    ( new _SignalStreamId () ) 
                            -> Add ( _table_id, _writer, _stream_ids );

    ( new _ConstU16StreamId ( _TL_TRACE_MAX_VALUE ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _ConstU16StreamId ( _TL_BASES_20 ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _ConstU16StreamId ( _TL_BASES_40 ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _ConstU16StreamId ( _TL_BASES_60 ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );

    ( new _ConstU8StreamId ( _TL_NREADS, 1 ) )
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _ConstU32StreamId ( _TL_READ_START, 0 /* 1 */ ) )
                            -> Add ( _table_id, _writer, _stream_ids );

        /*  Two special cases for RUN_DATE and COLLECTION_DATE
         */
    ( new _DateStreamId ( _TL_RUN_DATE ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _DateStreamId ( _TL_COLLECTION_DATE ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );

    ( new _CommExtStreamId ( _TL_COMMENTS ) )
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _CommExtStreamId ( _TL_EXTENDED_DATA ) )
                            -> Add ( _table_id, _writer, _stream_ids );

        /*  Here are Dummies ... 
         */
/* I depricated these fields on the schema level
    ( new _DummyStreamId ( _TL_TRACE_FILE ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _DummyStreamId ( _TL_BASE_FILE ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _DummyStreamId ( _TL_QUAL_FILE ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
    ( new _DummyStreamId ( _TL_PEAK_FILE ) ) 
                            -> Add ( _table_id, _writer, _stream_ids );
*/
}   /* TL_TraceFieldAdaptre :: _AddMandatoryColumns () */

void
TL_TraceFieldAdapter :: _AddOtherColumns ( const TL_TraceInfo & Info )
{
    for ( size_t llp = 0; llp < Info . FieldQty (); llp ++ ) {
        AddColumn ( Info . Field ( llp ) );
    }
}   /* TL_TraceFieldAdapter :: _AddOtherColumns () */

void
TL_TraceFieldAdapter :: OpenWriter ( const TL_TraceInfo & Info )
{
    _AddMandatoryColumns ();
    _AddOtherColumns ( Info );

    _writer -> open ();

    _SetColumnDefaults ();
}   /* TL_TraceFieldAdapter :: OpenWriter () */

void
TL_TraceFieldAdapter :: _SetColumnDefaults ()
{
    for (
        SMap :: iterator It = _stream_ids . begin ();
        It != _stream_ids . end ();
        It ++
    ) {
        It -> second -> SetDefault ( _writer );
    }
}   /* TL_TraceFieldAdapter :: _SetColumnDefaults () */

void
TL_TraceFieldAdapter :: _FillStream (
                                    const string & StreamName,
                                    size_t Qty
)
{
    TL_StreamId * SId = GetStreamId ( StreamName );
    if ( SId == NULL ) {
            /* JOJOBA ... prolly that is error
             */
        return;
    }

    if ( 0 < Qty ) {
        size_t WriteSize = Qty * ( SId -> ElemBits () / 8 );

        if ( _filler_capacity < WriteSize ) {
            size_t NewCap = ( ( WriteSize / 1024 ) + 1 ) * 1024;

            char * NewFil = new char [ NewCap ];
            memset ( NewFil, 0, NewCap );

            if ( _filler != NULL ) {
                delete [] _filler;
            }

            _filler = NewFil;
            _filler_capacity = NewCap;
        }

        GetStreamId ( StreamName ) -> Submit ( _writer, _filler, Qty );
    }
}   /* TL_TraceFieldAdapter :: _FillStream () */

void
TL_TraceFieldAdapter :: WriteTraceData (
                                const TL_TraceDataProvider & Provider,
                                const TL_TraceInfo & TraceInfo,
                                size_t RowIdx
)
{
        /* TRACE_NAME */
    string Str = Provider . Name ();
    if ( Str . empty () ) {
        throw TL_Exception ( "WriteTraceData (): Empty TRACE_NAME" );
    }
// JOJOBA
    GetStreamId ( _TL_TRACE_NAME ) ->
                    Submit ( _writer, Str . c_str (), Str . size () );
/*
    GetStreamId ( "NAME" ) ->
    // GetStreamId ( "name" ) ->
                    Submit ( _writer, Str . c_str (), Str . size () );
*/

// JOJOBA

        /* PROGRAM_ID */
    Str = Provider . ProgramID ();
#ifdef JOJOBA
    if ( Str . empty () ) {
        throw TL_Exception ( "WriteTraceData (): Empty PROGRAM_ID" );
    }
    GetStreamId ( _TL_PROGRAM_ID ) ->
                    Submit ( _writer, Str . c_str (), Str . size () );
#else /* JOJOBA */
    if ( ! Str . empty () ) {
        GetStreamId ( _TL_PROGRAM_ID ) ->
                    Submit ( _writer, Str . c_str (), Str . size () );
    }
#endif /* JOJOBA */

        /* Bases : READ and READ_LEN */
    size_t ReadQty = Provider . Bases () . size ();
    if ( ReadQty != 0 ) {
        GetStreamId ( _TL_READ ) -> Submit (
                                            _writer,
                                            * Provider . Bases (),
                                            ReadQty
                                            );
        PLOGMSG (
            klogDebug,
            ( klogDebug,
            "[BASES WRR] [$(read_qty)] chars",
            "severity=debug,read_qty=%d",
            ReadQty
            )
            );

        GetStreamId ( _TL_READ_LEN ) -> Submit ( _writer, & ReadQty, 1 );
    }

        /* PeakIndex : POSITION */
    size_t Qty = Provider . PeakIndex () . size ();
    if ( Qty != 0 ) {
        GetStreamId ( _TL_POSITION ) -> Submit (
                                            _writer,
                                            * Provider . PeakIndex (),
                                            Qty
                                            );
        PLOGMSG (
            klogDebug,
            ( klogDebug,
            "[PEAK WRR] [$(read_qty)] 32 uint",
            "severity=debug,read_qty=%d",
            Qty
            )
            );
    }
    else {
        if ( ReadQty != 0 ) {
            _FillStream ( _TL_POSITION, ReadQty );
        }
    }

        /* Quals : QUALITY and BASES_(20/40/60) */
    Qty = Provider . ProbScores () . size ();
    if ( Qty != 0 ) {
        GetStreamId ( _TL_QUALITY ) -> Submit (
                                            _writer,
                                            * Provider . ProbScores (),
                                            Qty
                                            );
        PLOGMSG (
            klogDebug,
            ( klogDebug,
            "[QUAL WRR] [$(read_qty)] uchars",
            "severity=debug,read_qty=%d",
            Qty
            )
            );

        uint16_t Val = Provider . Bases20 ();
        GetStreamId ( _TL_BASES_20 ) -> Submit ( _writer, & Val, 1 );

        Val = Provider . Bases40 ();
        GetStreamId ( _TL_BASES_40 ) -> Submit ( _writer, & Val, 1 );

        Val = Provider . Bases60 ();
        GetStreamId ( _TL_BASES_60 ) -> Submit ( _writer, & Val, 1 );
    }
    else {
        if ( ReadQty != 0 ) {
            _FillStream ( _TL_QUALITY, ReadQty );
            _FillStream ( _TL_BASES_20, 1 );
            _FillStream ( _TL_BASES_40, 1 );
            _FillStream ( _TL_BASES_60, 1 );
        }
    }

        /* Traces : SIGNAL */
    Qty = Provider . SampleCombined () . size ();
    if ( Qty != 0 ) {
        GetStreamId ( _TL_SIGNAL ) -> Submit (
                                        _writer,
                                        * Provider . SampleCombined (),
                                        Qty
                                        );
        PLOGMSG (
            klogDebug,
            ( klogDebug,
            "[SIGN WRR] [$(read_qty)] 64 uint",
            "severity=debug,read_qty=%d",
            Qty
            )
            );

        Qty = Provider . MaxTraceVal ();
        GetStreamId ( _TL_TRACE_MAX_VALUE ) -> Submit ( _writer, & Qty, 1 );
    }

#ifdef _NO_FLOW_KEY_
// These two depricated accorting to SRA-5693
        /* FlowChars : FLOW_CHARS */
    Qty = Provider . FlowChars () . size ();
    if ( Qty != 0 ) {
        GetStreamId ( _TL_FLOW_CHARS ) -> Submit (
                                            _writer,
                                            * Provider . FlowChars (),
                                            Qty
                                            );
    }

        /* KeySequence : _TL_KEY_SEQUENCE */
    Qty = Provider . KeySequence () . size ();
    if ( Qty != 0 ) {
        GetStreamId ( _TL_KEY_SEQUENCE ) -> Submit (
                                            _writer,
                                            * Provider . KeySequence (),
                                            Qty
                                            );
    }
#endif /* _NO_FLOW_KEY_ */

        /* Special order : _TL_NREADS and _TL_READ_START */
    uint8_t V8 = 1;
    GetStreamId ( _TL_NREADS ) -> Submit ( _writer, & V8, 1 );

        /*  We should add Comments and ExtendedData if they are exists
         */
    Str = Provider . Comments ();
    if ( ! Str . empty () ) {
        GetStreamId ( _TL_COMMENTS ) ->
                    Submit ( _writer, Str );
    }

    Str = Provider . ExtendedData ();
    if ( Str . empty () ) {
        Str = TraceInfo . Value ( RowIdx, _TL_EXTENDED_DATA );
    }
    if ( ! Str . empty () ) {
        GetStreamId ( _TL_EXTENDED_DATA ) ->
                    Submit ( _writer, Str );
    }

    /* Apparently next two lines are not needed, and could be
     * misleading, but let it be here :LOL:
     */
    uint32_t V32 = 0;
    GetStreamId ( _TL_READ_START ) -> Submit ( _writer, & V32, 1 );

        /*  Initialize and write regular fields
         */
    for ( size_t llp = 0; llp < TraceInfo . FieldQty (); llp ++ ) {
        if ( TraceInfo . Field ( llp ) . IsDeprecated () ) {
                /* Party time */
            continue;
        }

        string Name = TraceInfo . Field ( llp ) . Name ();

        string Data = TraceInfo . Value ( RowIdx, Name );
        if ( ! Data . empty () ) {
            TL_StreamId * Id = GetStreamId ( Name );
            if ( Id -> IsExternal () ) {
                Id -> Submit ( _writer, Data );
            }
        }
    }

    _WrtieDerivatives ( TraceInfo, RowIdx );

        /*  Flueshing
         */
    _writer -> nextRow ( _table_id );
}   /* TL_TraceFieldAdapter :: WriteTraceData () */

    /* Those are for some values, which are derived from the other
     * values, like READ_TYPE, or CONTROL_FLAGS
     */
void
TL_TraceFieldAdapter :: _WrtieDerivatives (
                                        const TL_TraceInfo & TraceInfo,
                                        size_t RowIdx
)
{
        /*  READ_TYPE : as for SRA-5693 
         */
    uint8_t ReadType = SRA_READ_TYPE_BIOLOGICAL;
    string Value = TraceInfo . Value ( RowIdx, _TL_TRACE_END );

    if ( ! Value . empty () ) {
        if ( Value == "F" || Value == "FORWARD" ) {
            ReadType |= SRA_READ_TYPE_FORWARD;
        }
        if ( Value == "R" || Value == "REVERSE" ) {
            ReadType |= SRA_READ_TYPE_REVERSE;
        }
    }
    GetStreamId ( _TL_READ_TYPE )
                            -> Submit ( _writer, & ReadType, 1 );

#ifdef __NEW_CONTROL_FLAGS__
        /*  Next two will be about WITHDRAWAL
         */
    bool IsWithdraw  = TraceInfo . Value ( RowIdx, _TL_TRACE_END ) == "WITHDRAW";
        /*  CONTROL_FLAGS : as for SRA-5693 
         */
    uint16_t ControlFlags = IsWithdraw ? 2 : 0;
    GetStreamId ( _TL_CONTROL_FLAGS )
                            -> Submit ( _writer, & ControlFlags, 1 );
#endif /* __NEW_CONTROL_FLAGS__ */

}   /* TL_TraceFieldAdapter :: _WrtieDerivatives () */
