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

using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * This file contains reader of 'SCF' file
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
struct TL_SCFHeader {
    uint32_t magic_number;       /* SCF_MAGIC */
    uint32_t num_samples;        /* Number of elements in Samples
                                  * matrix
                                  */
    uint32_t samples_offset;     /* Byte offset from start of file */
    uint32_t num_bases;          /* Number of bases in Bases matrix */
    uint32_t bases_left_clip;    /* OBSOLETE: No. bases in left clip
                                  * (vector)
                                  */
    uint32_t bases_right_clip;   /* OBSOLETE: No. bases in right clip
                                  * (qual)
                                  */
    uint32_t bases_offset;       /* Byte offset from start of file */
    uint32_t comments_size;      /* Number of bytes in Comment section
                                  */
    uint32_t comments_offset;    /* Byte offset from start of file */
    char     version[4];         /* "version.revision" */
    uint32_t sample_size;        /* precision of samples (in bytes) */
    uint32_t code_set;           /* uncertainty codes used */
    uint32_t private_size;       /* size of private data, 0 if none */
    uint32_t private_offset;     /* Byte offset from start of file */
    uint32_t spare[18];          /* Unused */
};

struct TL_SCFBases {
    uint32_t peak_index;    /* Index into Samples matrix for
                             * base position
                             */
    unsigned char probC;    /* Probability of it being a  C */
    unsigned char probA;    /* Probability of it being an A */
    unsigned char probG;    /* Probability of it being a  G */
    unsigned char probT;    /* Probability of it being a  T */
    char base;              /* Base called */
    unsigned char prob_sub; /* Probability of this base call being
                             * a substitution for another base
                             */
    unsigned char prob_ins; /* Probability of it being an overcall */
    unsigned char prob_del; /* Probability of an undercall at this
                             * point (extra base between this base
                             * and the previous base)
                             */
};

namespace _tl_ {

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * SCF Reader
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class TL_SCFTraceDataReader : public TL_TraceDataReaderWithExport {
public:
    TL_SCFTraceDataReader ();
    ~TL_SCFTraceDataReader ();

    void Read ( const TL_TraceFile & File );

private:
    void _ReadCheckHeader ( const TL_TraceFile & File );
    void _ReadSamples ( const TL_TraceFile & File );
    void _DeltaSamples (
                    const uint16_a_t & In,
                    uint16_a_t & Out,
                    uint16_t & MaxVal
                    );
    void _DeltaSamples (
                    const uchar_a_t & In,
                    uint16_a_t & Out,
                    uint16_t & MaxVal
                    );
    void _ReadBases ( const TL_TraceFile & File );
    void _ReadComments ( const TL_TraceFile & File );
    void _ReadPrivateData ( const TL_TraceFile & File );

    struct TL_SCFHeader _header;
    float _version;
};

TL_SCFTraceDataReader :: TL_SCFTraceDataReader ()
:   TL_TraceDataReaderWithExport ()
{
    memset ( & _header, 0, sizeof ( struct TL_SCFHeader ) );
    _version = 0.0f;
}   /* TL_SCFTraceDataReader :: TL_SCFTraceDataReader () */

TL_SCFTraceDataReader :: ~TL_SCFTraceDataReader ()
{
    memset ( & _header, 0, sizeof ( struct TL_SCFHeader ) );
    _version = 0.0f;
}   /* TL_SCFTraceDataReader :: ~TL_SCFTraceDataReader () */

void
TL_SCFTraceDataReader :: Read ( const TL_TraceFile & File )
{
    _ReadCheckHeader ( File );
    _ReadSamples ( File );
    _ReadBases ( File );
    _ReadComments ( File );
    _ReadPrivateData ( File );
}   /* TL_SCFTraceDataReader :: Read () */

void
TL_SCFTraceDataReader :: _ReadCheckHeader ( const TL_TraceFile & File )
{
    if ( File . ContentSize () < sizeof ( TL_SCFHeader ) ) {
        throw TL_Exception ( "Invalid SCF file format" );
    }

    memmove ( & _header, File . Content (), sizeof ( TL_SCFHeader ) );

    _header . magic_number     =
                    TL_McArch :: FromNet ( _header . magic_number );
    _header . num_samples      =
                    TL_McArch :: FromNet ( _header . num_samples );
    _header . samples_offset   =
                    TL_McArch :: FromNet ( _header . samples_offset );
    _header . num_bases        =
                    TL_McArch :: FromNet ( _header . num_bases );
    _header . bases_left_clip  =
                    TL_McArch :: FromNet ( _header . bases_left_clip );
    _header . bases_right_clip =
                    TL_McArch :: FromNet ( _header . bases_right_clip );
    _header . bases_offset     =
                    TL_McArch :: FromNet ( _header . bases_offset );
    _header . comments_size    =
                    TL_McArch :: FromNet ( _header . comments_size );
    _header . comments_offset  =
                    TL_McArch :: FromNet ( _header . comments_offset );
    _header . sample_size      =
                    TL_McArch :: FromNet ( _header . sample_size );
    _header . code_set         =
                    TL_McArch :: FromNet ( _header . code_set );
        /* Volodya forgot about it? */
    _header . private_size     =
                    TL_McArch :: FromNet ( _header . private_size );
    _header . private_offset   =
                    TL_McArch :: FromNet ( _header . private_offset );

        /* Version */
    _version = atof ( _header . version );
}   /* TL_SCFTraceDataReader :: _ReadCheckHeader () */

void
TL_SCFTraceDataReader :: _ReadSamples ( const TL_TraceFile & File )
{
    if ( ! _header . num_samples ) {
        throw TL_Exception ( "Ivalid SCF file : NUM_SAMPLES is ZERO " );
    }

    unsigned char * Ptr = ( unsigned char * ) File . Content ();
    size_t Size = File . ContentSize ();

    _samples_A . set ( _header . num_samples, 0 );
    _samples_C . set ( _header . num_samples, 0 );
    _samples_G . set ( _header . num_samples, 0 );
    _samples_T . set ( _header . num_samples, 0 );

    if ( _version < 2.9 ) {
        if ( 1 < _header . sample_size ) {
            uint32_t Size2Copy = _header . num_samples * 4 * sizeof ( uint16_t );

            if ( Size - _header . samples_offset < Size2Copy ) { 
                throw TL_Exception ( "Ivalid SCF file : invalid size " );
            }

            uint16_a_t Samples ( _header . num_samples * 4 , 0 );
            memmove ( * Samples, Ptr + _header . samples_offset, Size2Copy );
            _max_sample_value = 0;
            for ( uint32_t i = 0; i < _header . num_samples; i ++ ) {
                _samples_A [ i ] =
                    TL_McArch :: FromNet ( Samples [ ( i * 4 ) + 0 ] );
                _samples_C [ i ] =
                    TL_McArch :: FromNet ( Samples [ ( i * 4 ) + 1 ] );
                _samples_G [ i ] =
                    TL_McArch :: FromNet ( Samples [ ( i * 4 ) + 2 ] );
                _samples_T [ i ] =
                    TL_McArch :: FromNet ( Samples [ ( i * 4 ) + 3 ] );

                if ( _max_sample_value < _samples_A [ i ] )
                    _max_sample_value = _samples_A [ i ];
                if ( _max_sample_value < _samples_C [ i ] )
                    _max_sample_value = _samples_C [ i ];
                if ( _max_sample_value < _samples_G [ i ] )
                    _max_sample_value = _samples_G [ i ];
                if ( _max_sample_value < _samples_T [ i ] )
                    _max_sample_value = _samples_T [ i ];
            }
        }
        else {
            uint32_t Size2Copy = _header . num_samples * 4 * sizeof ( unsigned char );

            if ( Size - _header . samples_offset < Size2Copy ) { 
                throw TL_Exception ( "Ivalid SCF file : invalid size " );
            }

            uchar_a_t Samples ( _header . num_samples * 4 , 0 );
            memmove ( * Samples, Ptr + _header . samples_offset, Size2Copy );
            _max_sample_value = 0;
            for ( uint32_t i = 0; i < _header . num_samples; i ++ ) {
                _samples_A [ i ] = Samples [ ( i * 4 ) + 0 ];
                _samples_C [ i ] = Samples [ ( i * 4 ) + 1 ];
                _samples_G [ i ] = Samples [ ( i * 4 ) + 2 ];
                _samples_T [ i ] = Samples [ ( i * 4 ) + 3 ];

                if ( _max_sample_value < _samples_A [ i ] )
                    _max_sample_value = _samples_A [ i ];
                if ( _max_sample_value < _samples_C [ i ] )
                    _max_sample_value = _samples_C [ i ];
                if ( _max_sample_value < _samples_G [ i ] )
                    _max_sample_value = _samples_G [ i ];
                if ( _max_sample_value < _samples_T [ i ] )
                    _max_sample_value = _samples_T [ i ];
            }
        }
    }
    else {
        if ( 1 < _header . sample_size ) {
            uint32_t Size2Copy = _header . num_samples * 4 * sizeof ( uint16_t );

            if ( Size - _header . samples_offset < Size2Copy ) { 
                throw TL_Exception ( "Ivalid SCF file : invalid size " );
            }

            uint16_a_t Samples ( _header . num_samples, 0 );
            size_t MvSize = _header . num_samples * sizeof ( uint16_t );

            _max_sample_value = 0;
            unsigned char * Data = Ptr;

            memmove ( * Samples, Data, MvSize );
            _DeltaSamples ( Samples, _samples_A, _max_sample_value );
            Data += MvSize;

            memmove ( * Samples, Data, MvSize );
            _DeltaSamples ( Samples, _samples_C, _max_sample_value );
            Data += MvSize;

            memmove ( * Samples, Data, MvSize );
            _DeltaSamples ( Samples, _samples_G, _max_sample_value );
            Data += MvSize;

            memmove ( * Samples, Data, MvSize );
            _DeltaSamples ( Samples, _samples_T, _max_sample_value );
        }
        else {
            uint32_t Size2Copy = _header . num_samples * 4 * sizeof ( unsigned char );

            if ( Size - _header . samples_offset < Size2Copy ) { 
                throw TL_Exception ( "Ivalid SCF file : invalid size " );
            }

            uchar_a_t Samples ( _header . num_samples, 0 );
            size_t MvSize = _header . num_samples * sizeof ( unsigned char );

            _max_sample_value = 0;
            unsigned char * Data = Ptr;

            memmove ( * Samples, Data, MvSize );
            _DeltaSamples ( Samples, _samples_A, _max_sample_value );
            Data += MvSize;

            memmove ( * Samples, Data, MvSize );
            _DeltaSamples ( Samples, _samples_C, _max_sample_value );
            Data += MvSize;

            memmove ( * Samples, Data, MvSize );
            _DeltaSamples ( Samples, _samples_G, _max_sample_value );
            Data += MvSize;

            memmove ( * Samples, Data, MvSize );
            _DeltaSamples ( Samples, _samples_T, _max_sample_value );
        }
    }
}   /* TL_SCFTraceDataReader :: _ReadSamples () */

/*))    Note, MaxVal is accumulative, so do not reset it 
 ((*/
void
TL_SCFTraceDataReader :: _DeltaSamples (
                                    const uint16_a_t & In,
                                    uint16_a_t & Out,
                                    uint16_t & MaxVal
)
{
    Out . set ( In . size (), 0 );  /* really we do not need that */

    uint16_t Sample1 = 0;
    uint16_t Sample2 = 0;

    for ( uint32_t i = 0; i < In . size (); i ++ ) {
        Sample1 += TL_McArch :: FromNet ( In [ i ] );

        Sample2 = Sample1 + Sample2;

        Out [ i ] = Sample2;
        if ( MaxVal < Sample2 ) MaxVal = Sample2;
    }
}   /* TL_SCFTraceDataReader :: _DeltaSamples () */

/*))    Note, MaxVal is accumulative, so do not reset it 
 ((*/
void
TL_SCFTraceDataReader :: _DeltaSamples (
                                const uchar_a_t & In,
                                uint16_a_t & Out,
                                uint16_t & MaxVal
)
{
    Out . set ( In . size (), 0 );  /* really we do not need that */

    uint16_t Sample1 = 0;
    uint16_t Sample2 = 0;

    for ( uint32_t i = 0; i < In . size (); i ++ ) {
        Sample1 += ( uint16_t ) In [ i ];

        Sample2 = Sample1 + Sample2;

        Out [ i ] = Sample2;
        if ( MaxVal < Sample2 ) MaxVal = Sample2;
    }
}   /* TL_SCFTraceDataReader :: _DeltaSamples () */

void
TL_SCFTraceDataReader :: _ReadBases ( const TL_TraceFile & File )
{
    if ( ! _header . num_bases ) {
        throw TL_Exception ( "Ivalid SCF file : NUM_BASES is ZERO " );
    }

    unsigned char * Ptr = ( unsigned char * ) File . Content ();
    size_t Size = File . ContentSize ();

    uint32_t Size2Copy = _header . num_bases * sizeof ( struct TL_SCFBases );

    if ( Size - _header . bases_offset < Size2Copy ) {
        throw TL_Exception ( "Ivalid SCF file : invalid size " );
    }

    _prob_A . set ( _header . num_bases, 0 );
    _prob_C . set ( _header . num_bases, 0 );
    _prob_G . set ( _header . num_bases, 0 );
    _prob_T . set ( _header . num_bases, 0 );
    _bases . set ( _header . num_bases, 0 );
    _peak_index . set ( _header . num_bases, 0 );

    if ( _version < 2.9 ) {
        tl_array < struct TL_SCFBases > Bases ( _header . num_bases, 0 );
        memmove ( * Bases, Ptr + _header . bases_offset, Size2Copy );

        for ( uint32_t i = 0; i < _header . num_bases; i ++ ) {
            char * PChar = ( char * ) & ( Bases [ i ] . peak_index );
            _peak_index [ i ] = TL_a2ui ( ( unsigned char * ) PChar );
            _prob_A [ i ] = Bases [ i ] . probA;
            _prob_C [ i ] = Bases [ i ] . probC;
            _prob_G [ i ] = Bases [ i ] . probG;
            _prob_T [ i ] = Bases [ i ] . probT;
            _bases [ i ] = Bases [ i ] . base;
        }
    }
    else {
        unsigned char * Data = Ptr + _header . bases_offset;

        memmove ( * _peak_index, Data, _header . num_bases * sizeof ( uint32_t ) );
        TL_McArch :: FromNet ( _peak_index );

        Data += _header . num_bases * sizeof ( uint32_t );
        memmove ( * _prob_A , Data, _header . num_bases );

        Data += _header . num_bases;
        memmove ( * _prob_C , Data, _header . num_bases );

        Data += _header . num_bases;
        memmove ( * _prob_G , Data, _header . num_bases );

        Data += _header . num_bases;
        memmove ( * _prob_T , Data, _header . num_bases );

        Data += _header . num_bases;
        memmove ( * _bases , Data, _header . num_bases );

        if ( 3.0f < _version ) {

            Data += _header . num_bases;
            memmove ( * _prob_sub , Data, _header . num_bases );

            Data += _header . num_bases;
            memmove ( * _prob_ins , Data, _header . num_bases );

            Data += _header . num_bases;
            memmove ( * _prob_del , Data, _header . num_bases );

            _extra_probs = true;
        }
    }

        /*  Common for every version
         *  After conversation with Kurt we decided to keep all
         *  bases character, and not substitute them for 'n'
         */
    for ( uint32_t i = 0; i < _header . num_bases; i ++ ) {
#ifdef _N_SUBSTITUTE_
        switch ( _bases [ i ] ) {
            case 'A': case 'a': case 'C': case 'c':
            case 'G': case 'g': case 'T': case 't':
                _bases [ i ] = tolower ( _bases [ i ] );
                break;
            default:
                _bases [ i ] = 'n';
                break;
        }
#else /* _N_SUBSTITUTE_ */
            /*  some characters are recected on schema level
             */
        if ( isalpha ( _bases [ i ] ) ) {
            _bases [ i ] = tolower ( _bases [ i ] );
        }
        else {
            _bases [ i ] = 'n';
        }
#endif /* _N_SUBSTITUTE_ */
    }

    _valid_scores = true;
}   /* TL_SCFTraceDataReader :: _ReadBases () */

void
TL_SCFTraceDataReader :: _ReadComments ( const TL_TraceFile & File )
{
    if ( ! _header . comments_size ) {
       /* That is normal situation
        */
        return;
    }

    unsigned char * Ptr = ( unsigned char * ) File . Content ();
    size_t Size = File . ContentSize ();

    if ( Size - _header . comments_offset < _header . comments_size ) {
        throw TL_Exception ( "Ivalid SCF file : invalid comments size " );
    }

    _comments = string (
                        ( char * ) ( Ptr + _header . comments_offset ),
                        _header . comments_size
                        );
}   /* TL_SCFTraceDataReader :: _ReadComments () */

void
TL_SCFTraceDataReader :: _ReadPrivateData ( const TL_TraceFile & File )
{
    if ( ! _header . private_size ) {
       /* That is normal situation
        */
        return;
    }

    unsigned char * Ptr = ( unsigned char * ) File . Content ();
    size_t Size = File . ContentSize ();

    if ( Size - _header . private_offset < _header . private_size ) {
        throw TL_Exception ( "Ivalid SCF file : invalid private size " );
    }

    _private_data . set ( _header . private_size, 0 );
    memmove (
            * _private_data,
            Ptr + _header . private_offset,
            _header . private_size
            );
}   /* TL_SCFTraceDataReader :: _ReadPrivateData () */

};  /* _tl_ namespace */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * Point of noreturn
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

void
_TraceDataReadSCF ( const TL_TraceFile & File, TL_Traces & Trace )
{
    TL_SCFTraceDataReader Reader;

    Reader . Read ( File );

    Reader . Export ( Trace );
}   /* _TraceDataReadSCF () */

