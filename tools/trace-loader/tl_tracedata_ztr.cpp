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
#include <list>
#include <sstream>

#include <arpa/inet.h>
#include <zlib.h>

#include <klib/log.h>

#include "tl_util.hpp"
#include "tl_tracedata.hpp"
#include "tl_tracedata_reader.hpp"

using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * This file contains reader of 'ZTR' file
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

#define ZTR_MAGIC           "\256ZTR\r\n\032\n"
#define ZTR_VERSION_MAJOR   1
#define ZTR_VERSION_MINOR   2

struct TL_ZTRHeader {
    unsigned char _magic [8];       /* 0xae5a54520d0a1a0a (be) */
    unsigned char _major;           /* ZTR_VERSION_MAJOR */
    unsigned char _minor;           /* ZTR_VERSION_MINOR */
};

namespace _tl_ {

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * ZTR Text
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

class TL_ZTRText {
public :
    TL_ZTRText ();
    TL_ZTRText ( const string & Identifier, const string & Value );
    TL_ZTRText ( const TL_ZTRText & Text );
    ~TL_ZTRText ();

    TL_ZTRText & operator = ( const TL_ZTRText & Text );

    void Set ( const string & Identifier, const string & Value );
    void Reset ();

    inline bool IsEmpty () const { return _identifier . empty (); };

    inline const string & Identifier () const { return _identifier; };
    inline const string & Value () const { return _value; };

private :
    string _identifier;
    string _value;
};

TL_ZTRText :: TL_ZTRText ()
:   _identifier ()
,   _value ()
{
}   /* TL_ZTRText :: TL_ZTRText () */

TL_ZTRText :: TL_ZTRText ( const string & Id, const string & Val )
:   _identifier ( Id )
,   _value ( Val )
{
}   /* TL_ZTRText :: TL_ZTRText () */

TL_ZTRText :: TL_ZTRText ( const TL_ZTRText & Text )
:   _identifier ( Text . _identifier )
,   _value ( Text . _value )
{
}   /* TL_ZTRText :: TL_ZTRText () */

TL_ZTRText :: ~TL_ZTRText ()
{
    Reset ();
}   /* TL_ZTRText :: ~TL_ZTRText () */

TL_ZTRText &
TL_ZTRText :: operator = ( const TL_ZTRText & Text )
{
    if ( this != & Text ) {
        _identifier = Text . _identifier;
        _value = Text . _value;
    }

    return * this;
}   /* TL_ZTRText :: operator = () */

void
TL_ZTRText :: Set ( const string & Identifier, const string & Value )
{
    _identifier = Identifier;
    _value = Value;
}   /* TL_ZTRText :: Set () */

void
TL_ZTRText :: Reset ()
{
    _identifier . clear ();
    _value . clear ();
}   /* TL_ZTRText :: Reset () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * ZTR Chunk
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

class TL_ZTRTraceDataReader;

class TL_ZTRChunk {
public :
    TL_ZTRChunk ();
    TL_ZTRChunk ( const TL_ZTRChunk & Chunk );
    TL_ZTRChunk (
                uint32_t Type,
                const char_a_t & MetaData,
                const char_a_t & Data
                );
    ~TL_ZTRChunk ();

    TL_ZTRChunk & operator = ( const TL_ZTRChunk & Chunk );

    void Reset ();

    void Uncompress ();

    bool ProcessText ( TL_ZTRText & Text );

    inline uint32_t Type () const { return _type; };

public :    /* Volodya's methods */
    void decode_samples4 ( TL_ZTRTraceDataReader & Reader ) const;
    void decode_samples ( TL_ZTRTraceDataReader & Reader ) const;
    void decode_bases ( TL_ZTRTraceDataReader & Reader ) const;
    void decode_peaks ( TL_ZTRTraceDataReader & Reader ) const;
    void decode_confidence4 ( TL_ZTRTraceDataReader & Reader ) const;
    void decode_confidence1 ( TL_ZTRTraceDataReader & Reader ) const;
    void decode_clips ( TL_ZTRTraceDataReader & Reader ) const;

private :
    uint32_t _type;
    char_a_t _meta_data;
    char_a_t _data;
};

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * ZTR Data Reader
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class TL_ZTRTraceDataReader : public TL_TraceDataReaderWithExport {
    friend class TL_ZTRChunk;
public:
    typedef list < TL_ZTRChunk > Lunk;
    typedef list < TL_ZTRText > Lext;

public:
    TL_ZTRTraceDataReader ();
    ~TL_ZTRTraceDataReader ();

    void Read ( const TL_TraceFile & File );
    void Export ( TL_Traces & Trace );

private:
    void _ReadCheckHeader ( const TL_TraceFile & File );
    void _ReadChunks ( const TL_TraceFile & File );
    void _UncompressChunks ();
    void _ProcessTextChunks ();
    void _DecodeText ();
    void _DecodeData ();

private:
    Lunk _chunks;
    Lext _texts;

    uint32_t _fields_presents;
};

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * ZTR Chunk Impelemene
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

#define ZTR_TYPE_HEADER 0xae5a5452 /* M-. Z T R */

#define ZTR_TYPE_SAMP   0x53414d50
#define ZTR_TYPE_SMP4   0x534d5034
#define ZTR_TYPE_BASE   0x42415345
#define ZTR_TYPE_BPOS   0x42504f53
#define ZTR_TYPE_CNF4   0x434e4634
#define ZTR_TYPE_CNF1   0x434e4631
#define ZTR_TYPE_CSID   0x43534944
#define ZTR_TYPE_TEXT   0x54455854
#define ZTR_TYPE_CLIP   0x434c4950
#define ZTR_TYPE_COMM   0x434f4d4d
#define ZTR_TYPE_CR32   0x43523332
#define ZTR_TYPE_UNKN   0xffffffff

// Format types
//
#define ZTR_FORM_RAW             0
#define ZTR_FORM_RLE             1
#define ZTR_FORM_ZLIB            2
#define ZTR_FORM_DELTA1     64
#define ZTR_FORM_DELTA2     65
#define ZTR_FORM_DELTA4     66
#define ZTR_FORM_DDELTA1    67
#define ZTR_FORM_DDELTA2    68
#define ZTR_FORM_DDELTA4    69
#define ZTR_FORM_16TO8      70
#define ZTR_FORM_32TO8      71
#define ZTR_FORM_FOLLOW1    72
#define ZTR_FORM_CHEB445    73
#define ZTR_FORM_ICHEB      74

TL_ZTRChunk :: TL_ZTRChunk ()
:   _type ( 0 )
,   _meta_data ()
,   _data ()
{
}   /* TL_ZTRChunk :: TL_ZTRChunk () */

TL_ZTRChunk :: TL_ZTRChunk ( const TL_ZTRChunk & Chunk )
:   _type ( Chunk . _type )
,   _meta_data ( Chunk . _meta_data )
,   _data ( Chunk . _data )
{
}   /* TL_ZTRChunk :: TL_ZTRChunk () */

TL_ZTRChunk :: TL_ZTRChunk (
                            uint32_t Type,
                            const char_a_t & MetaData,
                            const char_a_t & Data
)
:   _type ( Type )
,   _meta_data ( MetaData )
,   _data ( Data )
{
}   /* TL_ZTRChunk :: TL_ZTRChunk () */

TL_ZTRChunk :: ~TL_ZTRChunk ()
{
    Reset ();
}   /* TL_ZTRChunk :: ~TL_ZTRChunk () */

TL_ZTRChunk &
TL_ZTRChunk :: operator = ( const TL_ZTRChunk & Chunk )
{
    if ( this != & Chunk ) {
        _type = Chunk . _type;
        _meta_data = Chunk . _meta_data;
        _data = Chunk . _data;
    }

    return * this;
}   /* TL_ZTRChunk :: operator = () */

void
TL_ZTRChunk :: Reset ()
{
    _type = 0;
    _meta_data . reset ();
    _data . reset ();
}   /* TL_ZTRChunk :: Reset () */


char * s_ztr_unrle         (char*, int, int*);
char * s_ztr_inflate       (char*, int, int*);
char * s_ztr_recorrelate1  (char*, int, int*);
char * s_ztr_recorrelate2  (char*, int, int*);
char * s_ztr_recorrelate4  (char*, int, int*);
char * s_ztr_expand_8to16  (char*, int, int*);
char * s_ztr_expand_8to32  (char*, int, int*);
char * s_ztr_unfollow      (char*, int, int*);
char * s_ztr_ichebuncomp   (char*, int, int*);

void
TL_ZTRChunk :: Uncompress ()
{
        /*  Empty Chunk
         */
    if ( ! _data . size () ) {
        return;
    }

    char * Ptr = ( char * ) * _data;
    int Len = ( int ) _data . size ();
    char * NewPtr = NULL;
    int NewLen = 0;

    do {
        switch ( Ptr [ 0 ]) {
            case ZTR_FORM_RLE:
                PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> RLE]", "severity=debug" ) );
                NewPtr = s_ztr_unrle ( Ptr, Len, & NewLen );
                break;

            case ZTR_FORM_ZLIB:
                PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> ZLIB]", "severity=debug" ) );
                NewPtr = s_ztr_inflate ( Ptr, Len, & NewLen );
                break;

            case ZTR_FORM_DELTA1:
                PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> DELTA1]", "severity=debug" ) );
                NewPtr = s_ztr_recorrelate1 ( Ptr, Len, & NewLen );
                break;

            case ZTR_FORM_DELTA2:
                PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> DELTA2]", "severity=debug" ) );
                NewPtr = s_ztr_recorrelate2 ( Ptr, Len, & NewLen );
                break;

            case ZTR_FORM_DELTA4:
                PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> DELTA4]", "severity=debug" ) );
                NewPtr = s_ztr_recorrelate4 ( Ptr, Len, & NewLen );
                break;

            case ZTR_FORM_16TO8:
                PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> 16TO8]", "severity=debug" ) );
                NewPtr = s_ztr_expand_8to16 ( Ptr, Len, & NewLen );
                break;

            case ZTR_FORM_32TO8:
                PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> 32TO8]", "severity=debug" ) );
                NewPtr = s_ztr_expand_8to32 ( Ptr, Len, & NewLen );
                break;

            case ZTR_FORM_FOLLOW1:
                PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> FOLLOW1]", "severity=debug" ) );
                NewPtr = s_ztr_unfollow ( Ptr, Len, & NewLen );
                break;

            case ZTR_FORM_ICHEB:
                PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> ICHEB]", "severity=debug" ) );
                NewPtr = s_ztr_ichebuncomp ( Ptr, Len, & NewLen );
                break;

            case ZTR_FORM_RAW:
                PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> RAW]", "severity=debug" ) );
                default:
                break;
        }

        if ( NewPtr ) {
            _data . set ( NewLen, NewPtr );

            delete [] NewPtr;

            Ptr = ( char * ) * _data;
            Len = ( int ) _data . size ();
        }
    } while ( Ptr [ 0 ] != ZTR_FORM_RAW );
}   /* TL_ZTRChunk :: Uncompress () */

bool
TL_ZTRChunk :: ProcessText ( TL_ZTRText & Text )
{
    Text . Reset ();

    if ( _type != ZTR_TYPE_TEXT ) {
        return false;
    }

    if ( _data . size () == 0 ) {
        return false;
    }

    char * Ptr = ( char * ) * _data;
    char * End = Ptr + _data . size ();

    Ptr ++;

    while ( Ptr < End ) {
        char * Ident = Ptr;
        size_t Dat = strlen ( Ident ) + 1;
        if ( End <= Ptr + Dat ) {
            break;
        }
        Ptr += Dat;

        char * Value = Ident + strlen ( Ident ) + 1;

        if ( strlen ( Ident ) != 0 ) {
            Text . Set ( Ident, Value );
        }

        Dat = strlen ( Value ) + 1;
        if ( End < Ptr + Dat ) {
            break;
        }

        Ptr += Dat;
    }

    return ! Text . IsEmpty ();
}   /* TL_ZTRChunk :: ProcessText () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * Decoding methods by Alekseyev, slightly modified
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
void
TL_ZTRChunk :: decode_samples4 ( TL_ZTRTraceDataReader & Reader ) const
{
    if ( _data . size () == 0 ) {
        throw TL_Exception ( "decode_samples4 () : no data provided" );
    }

    unsigned char * Data = ( ( unsigned char * ) * _data ) + 2;
    size_t Size = ( _data . size () - 2 ) / ( 4 * sizeof ( uint16_t ) );
    if ( Size == 0 ) {
        throw TL_Exception ( "decode_samples4 () : invalid data provided" );
    }

    size_t Offset = Size * sizeof ( uint16_t );

    Reader . _samples_A . set ( Size, 0 );
    memmove ( * Reader . _samples_A, Data + ( 0 * Offset ), Offset );
    TL_McArch :: FromNet ( Reader . _samples_A );

    Reader . _samples_C . set ( Size, 0 );
    memmove ( * Reader . _samples_C, Data + ( 1 * Offset ), Offset );
    TL_McArch :: FromNet ( Reader . _samples_C );

    Reader . _samples_G . set ( Size, 0 );
    memmove ( * Reader . _samples_G, Data + ( 2 * Offset ), Offset );
    TL_McArch :: FromNet ( Reader . _samples_G );

    Reader . _samples_T . set ( Size, 0 );
    memmove ( * Reader . _samples_T, Data + ( 3 * Offset ), Offset );
    TL_McArch :: FromNet ( Reader . _samples_T );

}   /* TL_ZTRChunk :: decode_samples4 () */

void
TL_ZTRChunk :: decode_samples ( TL_ZTRTraceDataReader & Reader ) const
{
    if ( _data . size () == 0 ) {
        throw TL_Exception ( "decode_samples () : no data provided" );
    }

    if ( _meta_data . size () == 0 ) {
        throw TL_Exception ( "decode_samples () : no meta data provided" );
    }

    unsigned char * Data = ( ( unsigned char * ) * _data ) + 2;
    size_t Size = ( _data . size () - 2 ) / ( 4 * sizeof ( uint16_t ) );
    if ( Size == 0 ) {
        throw TL_Exception ( "decode_samples () : invalid data provided" );
    }

    Reader . _samples_A . set ( Size, 0 );
    Reader . _samples_C . set ( Size, 0 );
    Reader . _samples_G . set ( Size, 0 );
    Reader . _samples_T . set ( Size, 0 );

    switch ( _meta_data [ 0 ] ) {
        case 'A': case 'a':
            memmove (
                    * Reader . _samples_A,
                    Data,
                    Size * sizeof ( uint16_t )
                    );
            TL_McArch :: FromNet ( Reader . _samples_A );
            break;

        case 'C': case 'c':
            memmove (
                    * Reader . _samples_C,
                    Data,
                    Size * sizeof ( uint16_t )
                    );
            TL_McArch :: FromNet ( Reader . _samples_C );
            break;

        case 'T': case 't':
            memmove (
                    * Reader . _samples_T,
                    Data,
                    Size * sizeof ( uint16_t )
                    );
            TL_McArch :: FromNet ( Reader . _samples_T );
            break;

        case 'G': case 'g':
            memmove (
                    * Reader . _samples_G,
                    Data,
                    Size * sizeof ( uint16_t )
                    );
            TL_McArch :: FromNet ( Reader . _samples_G );
            break;

    }
}   /* TL_ZTRChunk :: decode_samples () */

void
TL_ZTRChunk :: decode_bases ( TL_ZTRTraceDataReader & Reader ) const
{
    if ( _data . size () == 0 ) {
        throw TL_Exception ( "decode_bases () : no data provided" );
    }

    size_t Size = _data . size () - 1;
    if ( Size == 0 ) {
        throw TL_Exception ( "decode_bases () : invalid data provided" );
    }

    const char * Data = * _data;
    Data ++;

    Reader . _bases . set ( Size, 0 );
    memmove ( * Reader . _bases, Data, Size);

       /*  After conversation with Kurt we decided to keep all
        *  bases character, and not substitute them for 'n'
        */
#ifdef _N_SUBSTITUTE_
    for ( size_t i = 0; i < Size; i ++ ) {
        switch ( Reader . _bases [ i ] ) {
            case 'A': case 'a':
            case 'C': case 'c':
            case 'G': case 'g':
            case 'T': case 't':
                Reader . _bases [ i ] =
                            (char) tolower ( Reader . _bases [ i ] );
                break;
            default:
                Reader . _bases [ i ] = 'n';
                break;
        }
    }
#else /* _N_SUBSTITUTE_ */
    for ( size_t i = 0; i < Size; i ++ ) {
        if ( isalpha ( Reader . _bases [ i ] ) ) {
            Reader . _bases [ i ] =
                            (char) tolower ( Reader . _bases [ i ] );
        }
        else {
            Reader . _bases [ i ] = 'n';
        }
    }
#endif /* _N_SUBSTITUTE_ */
}   /* TL_ZTRChunk :: decode_bases () */

void
TL_ZTRChunk :: decode_peaks ( TL_ZTRTraceDataReader & Reader ) const
{
        /* That is Volodya's comment */
        // no data in here or total screw up
    if ( _data . size () == 0 ) {
        throw TL_Exception ( "decode_peaks () : no data provided" );
    }

    if ( _data . size () == 4 ) {
        throw TL_Exception ( "decode_peaks () : no data provided" );
    }


        /* Volodya */ // ignore first 4 bytes
    size_t Size = ( _data . size () - 4 ) / sizeof ( uint32_t );

    if ( Reader . _bases . size () != Size ) {
        throw TL_Exception ( "decode_peaks () : invalid data provided" );
    }

    Reader . _peak_index . set ( Size, 0 );
    memmove (
            * Reader . _peak_index,
            ( * _data ) + 4,
            Size * sizeof ( uint32_t )
            );
    TL_McArch :: FromNet ( Reader . _peak_index );
}   /* TL_ZTRChunk :: decode_peaks () */

void
TL_ZTRChunk :: decode_confidence4 ( TL_ZTRTraceDataReader & Reader ) const
{
    if ( _data . size () == 0 ) {
        throw TL_Exception ( "decode_confidence4 () : no data provided" );
    }

    size_t Size = Reader . _bases . size ();

        /* Volodya */ // unable to decode this block before basecalls
    if ( Size == 0 ) {
        throw TL_Exception ( "decode_confidence4 () : invalid data provided" );
    }

    Reader . _prob_A . set ( Size, 0 );
    Reader . _prob_C . set ( Size, 0 );
    Reader . _prob_G . set ( Size, 0 );
    Reader . _prob_T . set ( Size, 0 );

    unsigned char * Data = ( unsigned char * ) ( * _data ) + 1;

    for ( size_t i = 0, j = Size; i < Size; i ++ ) {
        switch ( Reader . _bases [ i ] ) {
            case 'A': case 'a':
                Reader . _prob_A [i] = Data [i];
                Reader . _prob_C [i] = Data [j++];
                Reader . _prob_G [i] = Data [j++];
                Reader . _prob_T [i] = Data [j++];
                break;
            case 'C': case 'c':
                Reader . _prob_A [i] = Data [j++];
                Reader . _prob_C [i] = Data [i];
                Reader . _prob_G [i] = Data [j++];
                Reader . _prob_T [i] = Data [j++];
                break;
            case 'G': case 'g':
                Reader . _prob_A [i] = Data [j++];
                Reader . _prob_C [i] = Data [j++];
                Reader . _prob_G [i] = Data [i];
                Reader . _prob_T [i] = Data [j++];
                break;
            default:
                Reader . _prob_A [i] = Data [j++];
                Reader . _prob_C [i] = Data [j++];
                Reader . _prob_G [i] = Data [j++];
                Reader . _prob_T [i] = Data [i];
                break;
        }
    }

    Reader . _valid_scores = true;
}   /* TL_ZTRChunk :: decode_confidence4 () */

void
TL_ZTRChunk :: decode_confidence1 ( TL_ZTRTraceDataReader & Reader ) const
{
    if ( _data . size () == 0 ) {
        throw TL_Exception ( "decode_confidence1 () : no data provided" );
    }

    size_t Size = Reader . _bases . size ();

        /* Volodya */ // unable to decode this block before basecalls
    if ( Size == 0 ) {
        throw TL_Exception ( "decode_confidence1 () : invalid data provided" );
    }

    Reader . _prob_A . set ( Size, 0 );
    Reader . _prob_C . set ( Size, 0 );
    Reader . _prob_G . set ( Size, 0 );
    Reader . _prob_T . set ( Size, 0 );

    unsigned char * Data = ( unsigned char * ) ( * _data ) + 1;

    for ( size_t i = 0; i < Size; i ++ ) {
        switch ( Reader . _bases [ i ] ) {
            case 'A': case 'a':
                Reader . _prob_A [i] = Data [i];
                break;
            case 'C': case 'c':
                Reader . _prob_C [i] = Data [i];
                break;
            case 'G': case 'g':
                Reader . _prob_G [i] = Data [i];
                break;
            case 'T': case 't':
                Reader . _prob_T [i] = Data [i];
                break;
            default:
                Reader . _prob_A [i] = Data [i];
                Reader . _prob_C [i] = Data [i];
                Reader . _prob_G [i] = Data [i];
                Reader . _prob_T [i] = Data [i];
                break;
        }
    }

    Reader . _valid_scores = true;
}   /* TL_ZTRChunk :: decode_confidence1 () */

void
TL_ZTRChunk :: decode_clips ( TL_ZTRTraceDataReader & Reader ) const
{
    if ( _data . size () == 0 ) {
        throw TL_Exception ( "decode_clips () : no data provided" );
    }

    unsigned char * Data = ( unsigned char * ) * _data;

    Reader . _clip_quality_left = TL_a2ui ( Data );
    Reader . _clip_quality_right = TL_a2ui ( Data + 5 );
}   /* TL_ZTRChunk :: decode_clips () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * Uncompress methods by Volodya
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

//
//
char* s_ztr_unrle (char* comp, int comp_len, int* uncomp_len)
{
  int      val, count, out_len;
  char*    uncomp = 0;
  unsigned char * in     = (unsigned char *)comp+6;
  int      guard  = (unsigned char ) comp [5];

  assert (comp && uncomp_len);

  /* Allocate
   */
  out_len =
    ((unsigned char )comp[1] <<  0) +
    ((unsigned char )comp[2] <<  8) +
    ((unsigned char )comp[3] << 16) +
    ((unsigned char )comp[4] << 24);

  uncomp = new char [out_len];
  assert (uncomp);

  for (register int in_i=0, out_i=0; out_i < out_len; in_i++)
  {
    assert (in_i+6 <= comp_len);
    if (in[in_i] != guard)
    {
      /* When not 'guard' it's easy - just output this token */
      assert (out_i < out_len);
      uncomp [out_i++] = in [in_i];
    }
    else
    {
      /*
       * Found an 'guard' token. If next token is zero, then
       * we were simply escaping a real 'guard' token in the input
       * data, otherwise output a string of bytes.
       */
      count = in [++in_i];
      if (count != 0)
      {
        val = in [++in_i];
        assert (out_i+count <= out_len);
        memset (uncomp+out_i, val, count*sizeof(char));
        out_i+= count;
      }
      else
      {
        assert (out_i >= 0 && out_i < out_len);
        uncomp [out_i++] = guard;
      }
    }
  }
  *uncomp_len = out_len;

  return uncomp;
}   /* s_ztr_unrle () */
//
//
char* s_ztr_inflate (char* comp, int comp_len, int* uncomp_len)
{
  z_stream zstr;
  char*    uncomp;
  int      ulen;

  assert (comp && uncomp_len);

  /* Allocate
   */
  ulen =
    ((unsigned char )comp[1] <<  0) +
    ((unsigned char )comp[2] <<  8) +
    ((unsigned char )comp[3] << 16) +
    ((unsigned char )comp[4] << 24);

  uncomp = new char [ulen];
  assert (uncomp);

  /* Initialise zlib */
  zstr.zalloc    = (alloc_func) 0;
  zstr.zfree     = (free_func)  0;
  zstr.opaque    = (voidpf)     0;
  zstr.next_in   = (unsigned char *) comp+5;
  zstr.avail_in  = comp_len-5;
  zstr.next_out  = (unsigned char *) uncomp;
  zstr.avail_out = ulen;

  if (inflateInit(&zstr)       != Z_OK)         { delete [] uncomp; inflateEnd (&zstr); return 0; }
  if (inflate(&zstr, Z_FINISH) != Z_STREAM_END) { delete [] uncomp; inflateEnd (&zstr); return 0; }
  if (inflateEnd(&zstr)        != Z_OK)         { delete [] uncomp; return 0; }

  *uncomp_len = ulen;

  return uncomp;
}   /* s_ztr_inflate () */
//
//
char* s_ztr_recorrelate1 (char* comp, int comp_len, int* uncomp_len)
{
  int   z;
  int   u1 = 0, u2 = 0, u3 = 0;
  int   level = comp [1];
  char* uncomp;

  assert (comp && uncomp_len);

  uncomp = new char [comp_len-2];
  assert (uncomp);

  comp       += 2;
  comp_len   -= 2;
  *uncomp_len = comp_len;

  switch (level)
  {
    case 1:
      for (register int i=0;i<comp_len; i++)
      {
        z         = u1;
        uncomp[i] = (unsigned char )comp[i] + z;
        u1        = uncomp[i];
      }
      break;

    case 2:
      for (register int i=0; i<comp_len; i++)
      {
        z         = u1+u1-u2;
        uncomp[i] = (unsigned char )comp[i] + z;
        u2        = u1;
        u1        = uncomp[i];
      }
      break;

    case 3:
      for (register int i=0; i<comp_len; i++)
      {
        z         = u1+u1+u1-u2-u2-u2+u3;
        uncomp[i] = (unsigned char )comp[i] + z;
        u3        = u2;
        u2        = u1;
        u1        = uncomp[i];
      }
      break;
  }

  return uncomp;
}   /* s_ztr_recorrelate1 () */
//
//
char* s_ztr_recorrelate2 (char* comp,int comp_len, int* uncomp_len)
{
  int   z;
  int   u1 = 0, u2 = 0, u3 = 0;
  int   level = comp[1];
  char* uncomp;

  assert (comp && uncomp_len);

  uncomp = new char [comp_len-2];
  assert (uncomp);

  comp       += 2;
  comp_len   -= 2;
  *uncomp_len = comp_len;

  switch (level)
  {
    case 1:
      for (register int i=0; i<comp_len; i+=2)
      {
        z  = u1;
        u1 = (((unsigned char )comp[i] << 8) | (unsigned char )comp[i+1]) + z;
        uncomp[i  ] = (u1 >> 8) & 0xff;
        uncomp[i+1] = (u1 >> 0) & 0xff;
      }
      break;

    case 2:
      for (register int i=0; i<comp_len; i+=2)
      {
        z  = u1+u1-u2;
        u2 = u1;
        u1 = (((unsigned char )comp[i] << 8) | (unsigned char )comp[i+1]) + z;
        uncomp[i  ] = (u1 >> 8) & 0xff;
        uncomp[i+1] = (u1 >> 0) & 0xff;
      }
      break;

    case 3:
      for (register int i=0; i<comp_len; i+=2)
      {
        z  = u1+u1+u1-u2-u2-u2+u3;
        u3 = u2;
        u2 = u1;
        u1 = (((unsigned char )comp[i] << 8) | (unsigned char )comp[i+1]) + z;
        uncomp[i  ] = (u1 >> 8) & 0xff;
        uncomp[i+1] = (u1 >> 0) & 0xff;
      }
      break;
  }

  return uncomp;
}   /* s_ztr_recorrelate2 () */
//
//
char* s_ztr_recorrelate4 (char* comp, int comp_len, int* uncomp_len)
{
  int   z;
  int   u1 = 0, u2 = 0, u3 = 0;
  int   level = comp[1];
  char* uncomp;

  assert (comp && uncomp_len);

  uncomp = new char [comp_len-4];
  assert (uncomp);

  comp       += 4;
  comp_len   -= 4;
  *uncomp_len = comp_len;

  switch (level)
  {
    case 1:
      for (register int i=0; i<comp_len; i+=4)
      {
        z  = u1;
        u1 = z +(((unsigned char )comp[i  ] << 24) | ((unsigned char )comp[i+1] << 16) |
                 ((unsigned char )comp[i+2] <<  8) |  (unsigned char )comp[i+3]);
        uncomp[i  ] = (u1 >> 24) & 0xff;
        uncomp[i+1] = (u1 >> 16) & 0xff;
        uncomp[i+2] = (u1 >>  8) & 0xff;
        uncomp[i+3] = (u1 >>  0) & 0xff;
      }
      break;

    case 2:
      for (register int i=0; i<comp_len; i+=4)
      {
        z  = u1+u1 - u2;
        u2 = u1;
        u1 = z + (((unsigned char )comp[  i] << 24) | ((unsigned char )comp[i+1] << 16) |
                  ((unsigned char )comp[i+2] <<  8) |  (unsigned char )comp[i+3]);
        uncomp[i  ] = (u1 >> 24) & 0xff;
        uncomp[i+1] = (u1 >> 16) & 0xff;
        uncomp[i+2] = (u1 >>  8) & 0xff;
        uncomp[i+3] = (u1 >>  0) & 0xff;
      }
      break;
	
    case 3:
      for (register int i=0; i<comp_len; i+=4)
      {
        z  = u1+u1+u1-u2-u2-u2+u3;
        u3 = u2;
        u2 = u1;
        u1 = z + (((unsigned char )comp[i  ] << 24) | ((unsigned char )comp[i+1] << 16) |
                  ((unsigned char )comp[i+2] <<  8) |  (unsigned char )comp[i+3]);
        uncomp[i  ] = (u1 >> 24) & 0xff;
        uncomp[i+1] = (u1 >> 16) & 0xff;
        uncomp[i+2] = (u1 >>  8) & 0xff;
        uncomp[i+3] = (u1 >>  0) & 0xff;
      }
      break;
  }

  return uncomp;
}   /* s_ztr_recorrelate4 () */
//
//
char* s_ztr_expand_8to16 (char* comp, int comp_len, int* uncomp_len)
{
  register int i, j;
  char*        uncomp;

  assert (comp && uncomp_len);

  uncomp = new char [comp_len*2];
  assert (uncomp);

  for (i=0, j=1; j<comp_len; i+=2)
  {
    if (comp[j] >= 0)
    {
      uncomp[i  ] = 0;
      uncomp[i+1] = comp[j++];
    }
    else
    if (comp[j] != -128)
    {
      uncomp[i+1] = comp[j++];
      uncomp[i  ] = -1;
    }
    else
    {
      j++;
      uncomp[i  ] = comp[j++];
      uncomp[i+1] = comp[j++];
    }
  }

  *uncomp_len = i;

  return uncomp;
}   /* s_ztr_expand_8to16 () */
//
//
char* s_ztr_expand_8to32 (char* comp, int comp_len, int* uncomp_len)
{
  register int i, j;
  char*        uncomp;

  assert (comp && uncomp_len);

  uncomp = new char [comp_len*4];
  assert (uncomp);

  for (i=0, j=1; j<comp_len; i+=4)
  {
    if (comp[j] != -128)
    {
      uncomp[i  ] = comp[j] < 0 ? -1 : 0;
      uncomp[i+1] = comp[j] < 0 ? -1 : 0;
      uncomp[i+2] = comp[j] < 0 ? -1 : 0;
      uncomp[i+3] = comp[j++];
    }
    else
    {
      j++;
      uncomp[i  ] = comp[j++];
      uncomp[i+1] = comp[j++];
      uncomp[i+2] = comp[j++];
      uncomp[i+3] = comp[j++];
    }
  }

  *uncomp_len = i;

  return uncomp;
}   /* s_ztr_expand_8to32 () */
//
//
char* s_ztr_unfollow (char* comp, int comp_len, int* uncomp_len)
{
  char            next[256];
  unsigned  char* uncomp;

  assert (comp && uncomp_len);

  uncomp = new unsigned char [comp_len-257];
  assert (uncomp);

  /* Load next[] array
   */
  register int i;
  for (i=0; i<256; i++) next[i] = comp[i+1];

  /* Replace comp[x] with next[comp[x-1]] - comp[x]*/
  uncomp[0] = comp[i+1];

  comp_len -= 257;
  comp     += 257;
  for (i=1; i<comp_len; i++) uncomp[i] = next[uncomp[i-1]] - comp[i];

  *uncomp_len = i;

  return (char*) uncomp;
}   /* s_ztr_unfollow () */
//
//
#define CH1 150
#define CH2 105
#define _abs(n) ((n)>0?(n):-(n))
char* s_ztr_ichebuncomp (char* comp, int comp_len, int* uncomp_len)
{
  short* d16    = (short*)comp;
  int    nwords = comp_len/2-1;
  short* uncomp;

  assert (comp && uncomp_len);

  uncomp = new short [comp_len];
  assert (uncomp);

  // clock_started = clock();

/* Have big doubts that that code was tested on any platform different
 * than Linux: such constructions like: ntohs ( ntohs() ) are not good
 */
  /* Check for boundary cases
   */
  if (nwords <= 4)
  {
    switch (nwords)
    {
      case 4:
        uncomp[0] = ntohs (d16[1]);
        uncomp[1] = ntohs (ntohs(d16[2])+ntohs(uncomp[0]));
        uncomp[2] = ntohs (ntohs(d16[3])+ntohs(uncomp[1]));
        uncomp[3] = ntohs (ntohs(d16[4])+ntohs(uncomp[2]));
        break;
      case 3:
        uncomp[0] = ntohs (d16[1]);
        uncomp[1] = ntohs (ntohs(d16[2])+ntohs(uncomp[0]));
        uncomp[2] = ntohs (ntohs(d16[3])+ntohs(uncomp[1]));
        break;
      case 2:
        uncomp[0] = ntohs (d16[1]);
        uncomp[1] = ntohs (ntohs(d16[2])+ntohs(uncomp[0]));
        break;
      case 1:
        uncomp[0] = ntohs (d16[1]);
        break;
    }
    *uncomp_len = nwords*2;
    return (char*) uncomp;
  }

  // int frac[5] = {139,57,75,93,11};
#define _frac0 139
#define _frac1  57
#define _frac2  75
#define _frac3  93
#define _frac4  11
  // int fz [20] = {42,42,42,42,42,39,24,0,-24,-39,33,-12,-42,-12,33,24,-39,0,39,-24};
#define _fz0    42
#define _fz1    42
#define _fz2    42
#define _fz3    42
#define _fz4    42
#define _fz5    39
#define _fz6    24
#define _fz7    33
#define _fz8    12

  int dfac;

  /* First 4 values are just direct deltas */
  uncomp[0] = ntohs (d16[1]);
  uncomp[1] = ntohs (ntohs(d16[2])+ntohs(uncomp[0]));
  uncomp[2] = ntohs (ntohs(d16[3])+ntohs(uncomp[1]));
  uncomp[3] = ntohs (ntohs(d16[4])+ntohs(uncomp[2]));

  short* dptr  = uncomp;
  short* dptr2 = d16+5;

  int   p, dd, f0, f1, f2, f3, f4, coef[4];
  short diff;
  unsigned short d0, d1, d2, d3;

  /* Loop */

  for (register int i=4; i<nwords; i++)
  {
    d0 = ntohs(dptr[0]);
    d1 = ntohs(dptr[1]);
    d2 = ntohs(dptr[2]);
    d3 = ntohs(dptr[3]);

    f0 = d2 * _frac4 + d3 * _frac0;
    f1 = d2 * _frac3 + d3 * _frac1;
    f2 = d1 * _frac2 + d2 * _frac2;
    f3 = d0 * _frac1 + d1 * _frac3;
    f4 = d0 * _frac0 + d1 * _frac4;

    coef[0] = f0*_fz0 + f1*_fz1 + f2*_fz0  + f3*_fz0 + f4*_fz0;
    coef[1] = f0*_fz5 + f1*_fz6            - f3*_fz6 - f4*_fz5;
    coef[2] = f0*_fz7 - f1*_fz8 - f2*_fz0  - f3*_fz8 + f4*_fz7;
    coef[3] = f0*_fz6 - f1*_fz5            + f3*_fz5 - f4*_fz6;

    /*
     * computing p requires at most a temporary variable of
     * 24.1 * coef, but coef may be a full 32-bit integer.
     * If coef is sufficiently close to cause an integer overflow then
     * we scale it down.
     */
    int max = 0;
    if (max < _abs(coef[0])) max = _abs(coef[0]);
    if (max < _abs(coef[1])) max = _abs(coef[1]);
    if (max < _abs(coef[2])) max = _abs(coef[2]);
    if (max < _abs(coef[3])) max = _abs(coef[3]);

    if (max > (1<<26))
    {
      dfac = max/(1<<26) + 1;
      double one_dfac = 1. / dfac;
      // for (register int l=0; l<4; l++) coef[l] /= dfac;
      coef [0] = (int)((double)coef[0] * one_dfac);
      coef [1] = (int)((double)coef[1] * one_dfac);
      coef [2] = (int)((double)coef[2] * one_dfac);
      coef [3] = (int)((double)coef[3] * one_dfac);
    }
    else dfac = 1;

    dd = (coef[3]/3)*10+coef[2];
    p  = ((((dd/3)*10-coef[3]+coef[1])/3)*5-dd+coef[0]/2)/(CH1*CH2);
    p *= dfac;

    if (p < 0) p = 0;
    diff    = ntohs (*dptr2) + (short)p;
    dptr[4] = ntohs (diff);

    dptr++;
    dptr2++;
  }

  *uncomp_len = nwords*2;

  return (char*) uncomp;
}   /* s_ztr_ichebuncomp () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * ZTR Data Reader Impelemene
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_ZTRTraceDataReader :: TL_ZTRTraceDataReader ()
:   TL_TraceDataReaderWithExport ()
,   _fields_presents ( 0 )
{
}   /* TL_ZTRTraceDataReader :: TL_ZTRTraceDataReader () */

TL_ZTRTraceDataReader :: ~TL_ZTRTraceDataReader ()
{
    _fields_presents = 0;

    _chunks . clear ();
    _texts . clear ();
}   /* TL_ZTRTraceDataReader :: ~TL_ZTRTraceDataReader () */

void 
TL_ZTRTraceDataReader :: Read ( const TL_TraceFile & File ) 
{
    _ReadCheckHeader ( File );
    _ReadChunks ( File );
    _UncompressChunks ();
    _ProcessTextChunks ();

        /* Here we are */
    _DecodeText ();
    _DecodeData ();
}   /* TL_ZTRTraceReader :: Read () */

void 
TL_ZTRTraceDataReader :: Export ( TL_Traces & Traces ) 
{
        /* From volodya */
        // these have to be at least undecoded: samples, bases, base posititions, and qualities
        //
    if ( Traces . IsOverrideQuals () ) _fields_presents |= 8;
    if ( Traces . IsOverridePeaks () ) _fields_presents |= 4;

    if ( ( _fields_presents & 15 ) != 15 ) {
        throw TL_Exception ( "Can not esport TRACE data because not all fields presents" );
    }

    TL_TraceDataReaderWithExport :: Export ( Traces );
}   /* TL_ZTRTraceReader :: Export () */

void
TL_ZTRTraceDataReader :: _ReadCheckHeader ( const TL_TraceFile & File )
{
    if ( File . ContentSize () < sizeof ( TL_ZTRHeader ) ) {
        throw TL_Exception ( "Invalid ZTR file format" );
    }

    struct TL_ZTRHeader Header;
    memmove ( & Header, File . Content (), sizeof ( TL_ZTRHeader ) );

    if ( memcmp ( Header . _magic, ZTR_MAGIC, 8) != 0 ) {
        throw TL_Exception ( "Invalid ZTR file header" );
    }

    if ( Header . _major != ZTR_VERSION_MAJOR ) {
        throw TL_Exception ( "Invalid ZTR file version" );
    }
}   /* TL_ZTRTraceDataReader :: _ReadCheckHeader () */

const char *
_readUInt32 ( const char * Start, const char * End, uint32_t & Ret )
{
    if ( ( size_t ) ( End - Start ) < sizeof ( uint32_t ) ) {
        throw TL_Exception ( "Can not read uint32_t from file" );
    }

    uint32_t Type = 0;
    memmove ( & Type, Start, sizeof ( uint32_t ) );
    Ret = TL_McArch :: FromNet ( Type );

    return Start + sizeof ( uint32_t );
}   /* _readUInt32 () */

const char *
_readArr (
        const char * Start,
        const char * End,
        char_a_t & Ret
)
{
    uint32_t Size;
    const char * Next = _readUInt32 ( Start, End, Size );

    if ( ( size_t ) ( End - Start ) < Size ) {
        throw TL_Exception ( "Can not read char array from file" );
    }

    Ret . set ( Size, ( char * ) Next );

    return Next + Size;
}   /* _readArr () */

const char *
_skipArr ( const char * Start, const char * End )
{
    uint32_t Size;
    const char * Next = _readUInt32 ( Start, End, Size );

    if ( ( size_t ) ( End - Start ) < Size ) {
        throw TL_Exception ( "Can not read char array from file" );
    }

    return Next + Size;
}   /* _skipArr () */

void
TL_ZTRTraceDataReader :: _ReadChunks ( const TL_TraceFile & File )
{
    size_t Offset = sizeof ( struct TL_ZTRHeader );

    size_t ContentSize = File . ContentSize ();
    const char * Content = ( char * ) File . Content ();
    const char * ContentEnd = Content + ContentSize;
    Content += Offset;

    TL_ZTRChunk chunk;  // working horse
    TL_ZTRChunk pcurr;  // pointer to current
    TL_ZTRChunk pprev;  // pointer to previous in case need
                        // to drop current

    PLOGMSG( klogDebug, (klogDebug, "     [ZTR ----> READ CHUNKS]", "severity=debug" ) );

    while ( Content < ContentEnd ) {

        if ( ContentSize <= Offset ) {
            break;
        }

        if ( ContentSize < ( Offset + 4 ) ) {
            throw TL_Exception ( "Invalid ZTR file format" );
        }

        uint32_t Type = 0;
        Content = _readUInt32 ( Content, ContentEnd, Type );

        if ( Type == ZTR_TYPE_HEADER ) {
            break;
        }

        char_a_t MetaData;
        char_a_t Data;
        if (   Type == ZTR_TYPE_SAMP || Type == ZTR_TYPE_SMP4
            || Type == ZTR_TYPE_BASE || Type == ZTR_TYPE_BPOS
            || Type == ZTR_TYPE_CNF4 || Type == ZTR_TYPE_CNF1
            || Type == ZTR_TYPE_CSID || Type == ZTR_TYPE_TEXT
            || Type == ZTR_TYPE_CLIP
        ) {
            Content = _readArr ( Content, ContentEnd, MetaData );
            Content = _readArr ( Content, ContentEnd, Data );
        }
        else {
                // unknown chunk type, skipping
                //
            Type = ZTR_TYPE_UNKN;
            Content = _skipArr ( Content, ContentEnd );
            Content = _skipArr ( Content, ContentEnd );
        }

        TL_ZTRChunk Chunk ( Type, MetaData, Data );
        _chunks . insert ( _chunks . end (), Chunk );
    }

    PLOGMSG( klogDebug, (klogDebug, "     [ZTR <---- $(size) CHUNKS READED]", "severity=debug,size=%d", _chunks . size () ) );
}   /* TL_ZTRTraceDataReader :: _ReadChunks () */

void
TL_ZTRTraceDataReader :: _UncompressChunks ()
{
    PLOGMSG( klogDebug, (klogDebug, "     [ZTR ----> UNCOMP CHUNKS]", "severity=debug" ) );
    for ( 
        Lunk :: iterator It = _chunks . begin ();
        It != _chunks . end ();
        It ++
    ) {
        It -> Uncompress ();
    }
    PLOGMSG( klogDebug, (klogDebug, "     [ZTR <---- $(size) CHUNKS UNCOMP]", "severity=debug,size=%d", _chunks . size () ) );
}   /* TL_ZTRTraceDataReader :: _UncompressChunks () */

void
TL_ZTRTraceDataReader :: _ProcessTextChunks ()
{
    PLOGMSG( klogDebug, (klogDebug, "     [ZTR ----> READ TEXT]", "severity=debug" ) );
    size_t Processed = 0;
    for ( 
        Lunk :: iterator It = _chunks . begin ();
        It != _chunks . end ();
        It ++
    ) {
        TL_ZTRText Text;
        if ( It -> ProcessText ( Text ) ) {
            Processed ++;
            PLOGMSG( klogDebug, (klogDebug, "     [ZTR ----> $(ident)=$(value)]", "severity=debug,ident=%s,value=%s", Text . Identifier () . c_str (), Text . Value () . c_str () ) );

            _texts . insert ( _texts . end (), Text );
        }
    }
    PLOGMSG( klogDebug, (klogDebug, "     [ZTR <---- $(size) TEXT READ]", "severity=debug,size=%d", Processed ) );
}   /* TL_ZTRTraceDataReader :: _ProcessTextChunks () */

void
TL_ZTRTraceDataReader :: _DecodeText ()
{
    if ( _texts . size () == 0 ) {
        PLOGMSG( klogDebug, (klogDebug, "     [ZTR <---> NO TEXT TO DECODE]", "severity=debug" ) );
        return;
    }

    PLOGMSG( klogDebug, (klogDebug, "     [ZTR ----> DECODE TEXT]", "severity=debug" ) );
    stringstream Str;
    for ( 
        Lext :: iterator It = _texts . begin ();
        It != _texts . end ();
        It ++
    ) {
        Str << It -> Identifier () << "=" << It -> Value () << "\n";
    }

    _comments = Str . str ();

    PLOGMSG( klogDebug, (klogDebug, "     [ZTR <---- TEXT DECODED]", "severity=debug" ) );
}   /* TL_ZTRTraceDataReader :: _DecodeText () */

void
TL_ZTRTraceDataReader :: _DecodeData ()
{
    PLOGMSG( klogDebug, (klogDebug, "     [ZTR ----> DECODE DATA]", "severity=debug" ) );
    _fields_presents = 0;

    for ( 
        Lunk :: iterator It = _chunks . begin ();
        It != _chunks . end ();
        It ++
    ) {
        switch ( It -> Type () ) {
            case ZTR_TYPE_SMP4:
PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> SMP4]", "severity=debug" ) );
                It -> decode_samples4 ( * this );
                _fields_presents |= 1;
                break;

            case ZTR_TYPE_SAMP:
PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> SAMP]", "severity=debug" ) );
                It -> decode_samples ( * this );
                _fields_presents |= 1;
                break;

            case ZTR_TYPE_BASE:
PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> BASE]", "severity=debug" ) );
                It -> decode_bases ( * this );
                _fields_presents |= 2;
                break;

            case ZTR_TYPE_BPOS:
                It -> decode_peaks ( * this );
                _fields_presents |= 4;
                break;

            case ZTR_TYPE_CNF4:
                It -> decode_confidence4 ( * this );
                _fields_presents |= 8;
                break;

            case ZTR_TYPE_CNF1:
PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> CNF4]", "severity=debug" ) );
                It -> decode_confidence1 ( * this );
                _fields_presents |= 8;
                break;

            case ZTR_TYPE_CLIP:
PLOGMSG( klogDebug, (klogDebug, "          [ZTR ----> CLIP]", "severity=debug" ) );
                It -> decode_clips ( * this );
                _fields_presents |= 16;
                break;

            case ZTR_TYPE_TEXT:
            default:
                break;
        }
    }
    PLOGMSG( klogDebug, (klogDebug, "     [ZTR <---- DATA DECODED]", "severity=debug" ) );
}   /* TL_ZTRTraceDataReader :: _DecodeData () */

};  /* _tl_ namespace */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * Point of noreturn
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

void
_TraceDataReadZTR ( const TL_TraceFile & File, TL_Traces & Trace )
{
    TL_ZTRTraceDataReader Reader;

    Reader . Read ( File );

    Reader . Export ( Trace );
}   /* _TraceDataReadZTR () */
