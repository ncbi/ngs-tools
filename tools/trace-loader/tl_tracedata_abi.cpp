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

#include "tl_util.hpp"
#include "tl_tracedata.hpp"
#include "tl_tracedata_reader.hpp"

#include <stdio.h>  /* using sprintf once */

using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * This file contains reader of 'ABI' file
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

struct TL_ABIHeader
{
    uint32_t header_edge; // basic offset: 0 or 128
    uint32_t index_offset;
    uint32_t num_samples;
    uint32_t num_bases;
    uint32_t fwo;         // filter wheel order
    uint32_t dataA_offset;
    uint32_t dataC_offset;
    uint32_t dataG_offset;
    uint32_t dataT_offset;
    uint32_t bases_offset;
    uint32_t basepos_offset;
    uint32_t prob_offset;
    uint32_t comments_size;
    float  fspacing;
};

namespace _tl_ {

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * ABI Data Reader
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class TL_ABITraceDataReader : public TL_TraceDataReaderWithExport {
public:
    TL_ABITraceDataReader ();
    ~TL_ABITraceDataReader ();

    void Read ( const TL_TraceFile & File );
    // void Export ( TL_Traces & Traces );

private:
    uint32_t _ABIIndexEntryLW (
                                unsigned char * Ptr,
                                uint32_t Size,
                                uint32_t Label,
                                uint32_t Count,
                                int32_t LW,
                                uint32_t & Value
                                );

    uint32_t _ABIIndexEntryLW (
                                unsigned char * Ptr,
                                uint32_t Size,
                                const char * Label,
                                uint32_t Count,
                                int32_t LW,
                                uint32_t & Value
                                );
    string _ABIString (
                    unsigned char * Ptr,
                    uint32_t Size,
                    const char * Label,
                    uint32_t Count
                    );
    int32_t _ABIInt (
                    unsigned char * Ptr,
                    uint32_t Size,
                    const char * Label,
                    uint32_t Count,
                    char * Buffer,
                    uint32_t BufferSize
                    );

    void _ReadCheckHeader ( const TL_TraceFile & File );
    void _ReadSamples ( const TL_TraceFile & File );
    void _ReadBases ( const TL_TraceFile & File );
    void _ReadComments ( const TL_TraceFile & File );

    struct TL_ABIHeader _header;
};

TL_ABITraceDataReader :: TL_ABITraceDataReader ()
:   TL_TraceDataReaderWithExport ()
{
    memset ( & _header, 0, sizeof ( struct TL_ABIHeader ) );
}   /* TL_ABITraceDataReader :: TL_ABITraceDataReader () */

TL_ABITraceDataReader :: ~TL_ABITraceDataReader ()
{
    memset ( & _header, 0, sizeof ( struct TL_ABIHeader ) );
}   /* TL_ABITraceDataReader :: ~TL_ABITraceDataReader () */

void
TL_ABITraceDataReader :: Read ( const TL_TraceFile & File )
{
    _ReadCheckHeader ( File );
    _ReadSamples ( File );
    _ReadBases ( File );
    _ReadComments ( File );
}   /* TL_ABITraceDataReader :: Read () */

uint32_t
TL_ABITraceDataReader :: _ABIIndexEntryLW (
                                            unsigned char * Ptr,
                                            uint32_t Size,
                                            const char * Label,
                                            uint32_t Count,
                                            int32_t LW,
                                            uint32_t & Val
)
{
    uint32_t NewLabel = TL_a2ui ( ( unsigned char * ) Label );

    return _ABIIndexEntryLW ( Ptr, Size, NewLabel, Count, LW, Val );
}   /* TL_ABITraceDataReader :: _ABIIndexEntryLW () */

#define INDEX_ENTRYLENGTH 28

uint32_t
TL_ABITraceDataReader :: _ABIIndexEntryLW (
                                            unsigned char * Ptr,
                                            uint32_t Size,
                                            uint32_t Label,
                                            uint32_t Count,
                                            int32_t LW,
                                            uint32_t & Val
)
{

    unsigned char * ThePtr = NULL;
    uint32_t EntryNum = 0;
    uint32_t EntryLabel = 0;
    uint32_t EntryLW = 0;

    Val = 0;

    do {
        uint32_t SomeOffset =     _header . index_offset
                                + _header . header_edge
                                + ( EntryNum * INDEX_ENTRYLENGTH )
                                ;
        if ( SomeOffset + ( 2 * sizeof ( uint32_t ) ) > Size ) {
            return 0;
        };

        ThePtr = Ptr + SomeOffset;
        EntryLabel = TL_a2ui( ThePtr );

        ThePtr += sizeof ( uint32_t );
        EntryLW = TL_a2ui ( ThePtr );

        EntryNum++;
    } while ( ! ( EntryLabel == Label && EntryLW == Count ) );

    if ( LW >= 1 ) {
        Val = TL_a2ui ( ThePtr + ( (LW - 1) * sizeof ( uint32_t ) ) );
    }

    uint32_t Ret = _header . index_offset
                        + ( ( EntryNum - 1 ) * INDEX_ENTRYLENGTH );

    return Ret > Size ? 0 : Ret;
}   /* TL_ABITraceDataReader :: _ABIIndexEntryLW () */

#undef INDEX_ENTRYLENGTH

/*))    Original method did return bool ( always true ) now it string
 ((*/
string
TL_ABITraceDataReader :: _ABIString (
                                unsigned char * Ptr,
                                uint32_t Size,
                                const char * Label,
                                uint32_t Count
)
{
    uint32_t Length = 0;
    uint32_t Offset = _ABIIndexEntryLW (
                                        Ptr,
                                        Size,
                                        Label,
                                        Count,
                                        4,
                                        Length
                                        );
    if ( Length == 0 ) {
        return "";
    }

    if ( Length < 5 ) {
        Offset += 20;
    }
    else {
        if ( ! _ABIIndexEntryLW ( Ptr, Size, Label, Count, 5, Offset ) ) {
            return "";
        }
    }

        /* Huh ... short strings */
    Length = ( uint32_t ) ( * ( Ptr + Offset ) ); 

    string Ret ( "" );
    if ( Length != 0 ) {
        Ret = string ( ( const char * ) ( Ptr + Offset + 1 ), Length );
        for ( size_t llp = 0; llp < Length; llp ++ ) {
            if ( Ret [ llp ] == 0 ) {
                Ret . resize ( llp );
                break;
            }
        }
    }

    return Ret;
}   /* TL_ABITraceDataReader :: _ABIString () */

int32_t
TL_ABITraceDataReader :: _ABIInt (
                                unsigned char * Ptr,
                                uint32_t Size,
                                const char * Label,
                                uint32_t Count,
                                char * Buffer,
                                uint32_t BufferSize
)
{
    if ( _header . index_offset == 0 ) {
        memmove ( Buffer, Ptr, BufferSize );
        return ( int32_t ) BufferSize;
    }

    uint32_t Length = 0;
    uint32_t Offset = _ABIIndexEntryLW (
                                    Ptr,
                                    Size,
                                    Label,
                                    Count,
                                    4,
                                    Length
                                    );
    if ( ! Offset ) {
        return 1 * - 1;
    }

    if ( Length == 0 ) {
        return 0;
    }

        /* Determine offset
         */
    if ( Length <= 4) {
        Offset += 20;
    }
    else {
        if ( ! _ABIIndexEntryLW ( Ptr, Size, Label, Count, 5, Offset ) ) {
            return 0;
        }
    }

    if ( BufferSize < Length ) {
        Length = BufferSize;
    }

    memmove ( Buffer, Ptr + Offset, Length );

    return ( int32_t ) Length;
}   /* TL_ABITraceDataReader :: _ABIInt () */

#define UMSCHWANG(Where,What)   \
        if ( Where == 'C' ) _header . dataC_offset = What; else \
        if ( Where == 'A' ) _header . dataA_offset = What; else \
        if ( Where == 'G' ) _header . dataG_offset = What; else \
                            _header . dataT_offset = What
void
TL_ABITraceDataReader :: _ReadCheckHeader ( const TL_TraceFile & File )
{
    if ( File . ContentSize () < sizeof ( TL_ABIHeader ) ) {
        throw TL_Exception ( "Invalid ABI file format" );
    }

    unsigned char * Ptr = ( unsigned char * ) File . Content ();
    size_t Size = File . ContentSize ();

    uint32_t Magic = TL_a2ui ( Ptr );
    uint32_t Label = TL_a2ui ( ( unsigned char * ) "ABIF" );

    _header . header_edge = Magic == Label ? 0 : 128;
    _header . index_offset = TL_a2ui ( Ptr + _header . header_edge + 26 );

        // Get the number of points
    if ( ! _ABIIndexEntryLW (
                        Ptr, Size, "DATA", 9, 3, _header . num_samples
                        )
    ) {
        throw TL_Exception ( "Ivalid ABI file : NUM_SAMPLES undefined " );
    }

    if ( ! _ABIIndexEntryLW (
                        Ptr, Size, "PBAS", 1, 3, _header . num_bases
                        )
    ) {
        throw TL_Exception ( "Ivalid ABI file : NUM_BASES undefined " );
    }

        /*
         * The order of the DATA fields is determined by the field FWO_
         * Juggle around with data pointers to get it right
         */
        // Get the Filter Wheel Order (FWO_) field ...
    if ( ! _ABIIndexEntryLW (
                        Ptr, Size, "FWO_", 1, 5, _header . fwo
                        )
    ) {
            // Guess :     C A G T
        _header . fwo = 0x43414754;
    }

    uint32_t Offset;
    if ( ! _ABIIndexEntryLW ( Ptr, Size, "DATA", 9, 5, Offset ) ) {
        throw TL_Exception ( "Ivalid ABI file : OFFSET 9 undefined " );
    }
    UMSCHWANG ( ( ( _header . fwo >> 0x18 ) & 0xff ), Offset );

    if ( ! _ABIIndexEntryLW ( Ptr, Size, "DATA", 10, 5, Offset ) ) {
        throw TL_Exception ( "Ivalid ABI file : OFFSET 10 undefined " );
    }
    UMSCHWANG ( ( ( _header . fwo >> 0x10 ) & 0xff ), Offset );

    if ( ! _ABIIndexEntryLW ( Ptr, Size, "DATA", 11, 5, Offset ) ) {
        throw TL_Exception ( "Ivalid ABI file : OFFSET 11 undefined " );
    }
    UMSCHWANG ( ( ( _header . fwo >> 0x8 ) & 0xff ), Offset );

    if ( ! _ABIIndexEntryLW ( Ptr, Size, "DATA", 12, 5, Offset ) ) {
        throw TL_Exception ( "Ivalid ABI file : OFFSET 12 undefined " );
    }
    UMSCHWANG ( ( ( _header . fwo >> 0x0 ) & 0xff ), Offset );

    if ( ! _ABIIndexEntryLW (
                        Ptr, Size, "PBAS", 1, 5, _header . bases_offset
                        )
    ) {
        throw TL_Exception ( "Ivalid ABI file : BASES_OFFSET undefined " );
    }

    if ( ! _ABIIndexEntryLW (
                    Ptr, Size, "PLOC", 1, 5, _header . basepos_offset
                    )
    ) {
        throw TL_Exception ( "Ivalid ABI file : BASEPOS_OFFSET undefined " );
    }

    _ABIIndexEntryLW ( Ptr, Size, "PCON", 1, 5, _header . prob_offset );
}   /* TL_ABITraceDataReader :: _ReadCheckHeader () */

#undef UMSCHWANG

void
TL_ABITraceDataReader :: _ReadSamples ( const TL_TraceFile & File )
{
    if ( ! _header . num_samples ) {
        throw TL_Exception ( "Ivalid ABI file : NUM_SAMPLES is ZERO " );
    }

    size_t SizeToCopy = _header . num_samples * sizeof ( uint16_t );

    unsigned char * Ptr = ( unsigned char * ) File . Content ();
    size_t Size = File . ContentSize ();
    size_t Capacity = Size - _header . header_edge - SizeToCopy;

    if (   Capacity < _header . dataA_offset
        || Capacity < _header . dataC_offset
        || Capacity < _header . dataG_offset
        || Capacity < _header . dataT_offset
    ) {
        throw TL_Exception ( "Ivalid ABI file : invalid size " );
    }

    unsigned char * Base = Ptr + _header . header_edge;

    _samples_A . set ( _header . num_samples, 0 );
    memmove ( * _samples_A, Base + _header . dataA_offset, SizeToCopy );
    TL_McArch :: FromNet ( _samples_A );

    _samples_C . set ( _header . num_samples, 0 );
    memmove ( * _samples_C, Base + _header . dataC_offset, SizeToCopy );
    TL_McArch :: FromNet ( _samples_C );

    _samples_G . set ( _header . num_samples, 0 );
    memmove ( * _samples_G, Base + _header . dataG_offset, SizeToCopy );
    TL_McArch :: FromNet ( _samples_G );

    _samples_T . set ( _header . num_samples, 0 );
    memmove ( * _samples_T, Base + _header . dataT_offset, SizeToCopy );
    TL_McArch :: FromNet ( _samples_T );

        // Compute highest trace peak
    uint16_t MaxVal = 0;
    for ( uint32_t i = 0; i < _header . num_samples; i ++ ) {
        if ( MaxVal < _samples_A [ i ] ) MaxVal = _samples_A [ i ];
        if ( MaxVal < _samples_C [ i ] ) MaxVal = _samples_C [ i ];
        if ( MaxVal < _samples_G [ i ] ) MaxVal = _samples_G [ i ];
        if ( MaxVal < _samples_T [ i ] ) MaxVal = _samples_T [ i ];
    }

    _max_sample_value = MaxVal;
}   /* TL_ABITraceDataReader :: _ReadSamples () */

void
TL_ABITraceDataReader :: _ReadBases ( const TL_TraceFile & File )
{
    if ( ! _header . num_bases ) {
        throw TL_Exception ( "Ivalid ABI file : NUM_BASES is ZERO " );
    }

    size_t SizeToCopy = _header . num_bases * sizeof ( char );

    unsigned char * Ptr = ( unsigned char * ) File . Content ();
    size_t Size = File . ContentSize ();

    size_t Capacity = Size - _header . header_edge - SizeToCopy;

    if ( Capacity < _header . bases_offset ) {
        throw TL_Exception ( "Ivalid ABI file : invalid size " );
    }


        // get qualities
        //
    if ( _header . prob_offset != 0 ) {
        uint32_t Offset = _header . header_edge + _header . prob_offset;
        if ( Size < Offset + _header . num_bases) {
            throw TL_Exception ( "Ivalid ABI file : invalid PROB_OFFSET " );
        }

        _prob_A . set ( _header . num_bases, 0 );
        _prob_C . set ( _header . num_bases, 0 );
        _prob_G . set ( _header . num_bases, 0 );
        _prob_T . set ( _header . num_bases, 0 );

            // use slot A to read all the probabilities
            //
        memmove ( * _prob_A, Ptr + Offset, _header . num_bases * sizeof ( unsigned char ) );

        _valid_scores = true;
    }

        // get bases
        //
    _bases . set ( _header . num_bases, 0 );
    memmove ( * _bases, Ptr + _header . header_edge + _header . bases_offset, _header . num_bases * sizeof ( unsigned char ) );

    _peak_index . set ( _header . num_bases, 0 );
    uint32_t Offset = _header . header_edge + _header . basepos_offset;
    for ( uint32_t i = 0; i < _header . num_bases; i ++ ) {
        _peak_index [ i ] = ( uint32_t ) TL_a2us (
                        Ptr + Offset + ( i * sizeof ( unsigned short ) )
                        );
            // filter them out
        _bases [ i ] = tolower ( _bases [ i ] );
        switch ( _bases [ i ] ) {
            case 'a':
                break;
            case 'c':
                if ( _valid_scores ) {
                    _prob_C [ i ] = _prob_A [ i ];
                    _prob_A [ i ] = 0;
                }
                break;
            case 'g':
                if ( _valid_scores ) {
                    _prob_G [ i ] = _prob_A [ i ];
                    _prob_A [ i ] = 0;
                }
                break;
            case 't':
                if ( _valid_scores ) {
                    _prob_T [ i ] = _prob_A [ i ];
                    _prob_A [ i ] = 0;
                }
                break;
            default:
                if ( _valid_scores ) {
                    _prob_A [ i ] = 0;
                }
/*  After conversation with Kurt we decided to keep all
 *  bases character, and not substitute them for 'n'
 */
#ifdef _N_SUBSTITUTE_
                _bases [ i ] = 'n';
#else /* _N_SUBSTITUTE_ */
                if ( ! isalpha ( _bases [ i ] ) ) {
                    _bases [ i ] = 'n';
                }
#endif /* _N_SUBSTITUTE_ */
                break;
        }
    }
}   /* TL_ABITraceDataReader :: _ReadBases () */

void
TL_ABITraceDataReader :: _ReadComments ( const TL_TraceFile & File )
{
    unsigned char * Ptr = ( unsigned char * ) File . Content ();
    size_t Size = File . ContentSize ();

    stringstream Str;
    Str << "COMM=" << _ABIString ( Ptr, Size, "CMNT", 1 ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
    }

    Str << "NAME=" << _ABIString ( Ptr, Size, "SMPL", 1 ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
    }

    uint16_t UnSo = 0;
    _ABIInt ( Ptr, Size, "LANE", 1, ( char * ) & UnSo, sizeof ( uint16_t ) ); 
    Str << "LANE=" << UnSo << endl;
    if ( ! Str ) {
        throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
    }

        /*  Get Signal Strength Offset
         */
    uint32_t Offset;
    if ( _ABIIndexEntryLW ( Ptr, Size, "S/N%", 1, 5, Offset ) != 0 ) {
        int16_t Bases [ 4 ]; /* C, A, G, T */

#define UMSCHWANG(Where,What)   \
        if ( Where == 'C' ) Bases [ 0 ] = What; else \
        if ( Where == 'A' ) Bases [ 1 ] = What; else \
        if ( Where == 'G' ) Bases [ 2 ] = What; else \
                            Bases [ 3 ] = What
        unsigned char * Data = Ptr + Offset;
        UMSCHWANG ( ( ( _header . fwo >> 0x18 ) & 0xff ), TL_a2us ( Data ) );
        UMSCHWANG ( ( ( _header . fwo >> 0x10 ) & 0xff ), TL_a2us ( Data + 2 ) );
        UMSCHWANG ( ( ( _header . fwo >> 0x8 ) & 0xff ), TL_a2us ( Data + 4 ) );
        UMSCHWANG ( ( ( _header . fwo >> 0x0 ) & 0xff ), TL_a2us ( Data + 6 ) );
#undef UMSCHWANG

        Str << "Sign=" << "A=" << Bases [ 1 ] << ",C=" << Bases [ 0 ] << ",G=" << Bases [ 2 ] << ",T=" << Bases [ 3 ] << endl;
        if ( ! Str ) {
            throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
        }
    }

        /*  Get the spacing.. it's a float but don't worry yet
         */
    uint32_t UI4 = 0;
    _ABIInt ( Ptr, Size, "SPAC", 1, ( char * ) & UI4, sizeof ( uint32_t ) );
        /*  Converting to float only once ... no need separate function
         */
    union { int32_t i; float f; } uif;
    uif . i = UI4;
    _header. fspacing = uif . f;
    {   /*  dont want to use stream formatters setw, setf and
         *  restore them after one time using :LOL
         */
        char BB [ 16 ];
        sprintf ( BB, "%-6.2f", _header. fspacing );

        Str << "SPAC=" << BB << endl;
        if ( ! Str ) {
            throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
        }
    }

        /*  Get primer position
         */
    if ( _ABIIndexEntryLW ( Ptr, Size, "PPOS", 1, 5, UI4 ) ) {
        Str << "PRIM=" << ( UI4 >> 0x10 ) << endl;
        if ( ! Str ) {
            throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
        }
    }

        /*  RUND/RUNT
         */
    uint32_t Offset1;
    if ( _ABIIndexEntryLW ( Ptr, Size, "RUND", 1, 5, Offset1 ) ) {
        uint32_t Offset2, Offset3, Offset4;
        _ABIIndexEntryLW ( Ptr, Size, "RUND", 2, 5, Offset2 );
        _ABIIndexEntryLW ( Ptr, Size, "RUNT", 1, 5, Offset3 );
        _ABIIndexEntryLW ( Ptr, Size, "RUNT", 2, 5, Offset4 );

        tm t;
        memset ( & t, 0, sizeof ( t ) );
        t.tm_mday  =   ( Offset1 >> 0    ) & 0xff;
        t.tm_mon   = ( ( Offset1 >> 0x8  ) & 0xff ) - 1;
        t.tm_year  =   ( Offset1 >> 0x10 ) - 1900;
        t.tm_hour  =     Offset3 >> 0x18;
        t.tm_min   =   ( Offset3 >> 0x10 ) & 0xff;
        t.tm_sec   =   ( Offset3 >> 0x8  ) & 0xff;
        t.tm_isdst = -1;
        mktime   (&t);

        char BB [ 1024 ];
        strftime ( BB, sizeof ( BB ), "%a %d %b %T %Y", & t);

        Str << "DATE=" << BB << " to ";

        t.tm_mday  =   ( Offset2 >> 0    ) & 0xff;
        t.tm_mon   = ( ( Offset2 >> 0x8  ) & 0xff ) - 1;
        t.tm_year  =   ( Offset2 >> 0x10 ) - 1900;
        t.tm_hour  =     Offset4 >> 0x18;
        t.tm_min   =   ( Offset4 >> 0x10 ) & 0xff;
        t.tm_sec   =   ( Offset4 >> 0x8  ) & 0xff;
        t.tm_isdst = -1;
        mktime   (&t);

        strftime ( BB, sizeof ( BB ), "%a %d %b %T %Y", & t);

        Str << BB << endl;
        if ( ! Str ) {
            throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
        }

        sprintf (
                BB,
                "%04d%02d%02d.%02d%02d%02d - %04d%02d%02d.%02d%02d%02d",
                ( Offset1 >> 0x10 ),
                ( Offset1 >> 0x8  ) & 0xff,
                ( Offset1 >> 0x0  ) & 0xff,
                ( Offset3 >> 0x18 ),
                ( Offset3 >> 0x10 ) & 0xff,
                ( Offset3 >> 0x8  ) & 0xff,
                ( Offset2 >> 0x10 ),
                ( Offset2 >>  0x8 ) & 0xff,
                ( Offset2 >> 0x0  ) & 0xff,
                ( Offset4 >> 0x18 ),
                ( Offset4 >> 0x10 ) & 0xff,
                ( Offset4 >> 0x8  ) & 0xff
                );
        Str << "RUND=" << BB << endl;
        if ( ! Str ) {
            throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
        }
    }

        /*  Get Dye Primer Offset
         */
    Str << "DYEP=" << _ABIString ( Ptr, Size, "PDMF", 1 ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
    }

        /*  Get Machine Name Offset
         */
    Str << "MACH=" << _ABIString ( Ptr, Size, "MCHN", 1 ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
    }

        /*  Machine model
         */
    Str << "MODL=" << _ABIString ( Ptr, Size, "MODL", 1 ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
    }

        /*  Matrix file
         */
    Str << "MTXF=" << _ABIString ( Ptr, Size, "MTXF", 1 ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
    }

        /*  Base calling version
         */
    Str << "BCAL=" << _ABIString ( Ptr, Size, "SPAC", 2 ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
    }

        /*  Software versions
         */
    Str << "VER1=" << _ABIString ( Ptr, Size, "SVER", 1 ) << endl;
    Str << "VER2=" << _ABIString ( Ptr, Size, "SVER", 2 ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
    }

        /*  Get Gel Name Offset
         */
    Str << "GELN=" << _ABIString ( Ptr, Size, "GELN", 1 ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "_ReadComments () : Can not write stringstream" );
    }

    _comments = Str . str ();
}   /* TL_ABITraceDataReader :: _ReadComments () */

}; /* namespace _tl_ */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * Point of noreturn
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

void
_TraceDataReadABI ( const TL_TraceFile & File, TL_Traces & Trace )
{
    TL_ABITraceDataReader Reader;

    Reader . Read ( File );

    Reader . Export ( Trace );
}   /* _TraceDataReadABI () */

