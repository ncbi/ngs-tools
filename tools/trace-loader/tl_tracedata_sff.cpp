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
 * This file contains reader of 'SFF' file
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
struct TL_SFFHeader {
    uint32_t magic;                 // magic_number 0
    int32_t version;                // version 4
    uint64_t index_offset;          // index_offset 8
    uint32_t index_size;            // index_length 16
    int32_t num_reads;              // number_of_reads 20
    int16_t header_length;          // header_length 24
    int16_t key_length;             // key_length 26
    int16_t num_flows_per_read;     // number_of_flows_per_read 28
    char flowgram_format_code;      // flowgram_format_code 30
};

namespace _tl_ {

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * SFF Reader
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class TL_SFFTraceDataReader : public TL_TraceDataReaderWithExport {
public:
    TL_SFFTraceDataReader (
                        const string & Name,
                        const string & ProgramID
                        );
    ~TL_SFFTraceDataReader ();

    void Read ( const TL_TraceFile & File );

private:
    void ReadCheckHeader ( const TL_TraceFile & File );
    void ReadSamples ( const TL_TraceFile & File );
    void ReadBases ( const TL_TraceFile & File );
    void ReadComments ( const TL_TraceFile & File );

    string _name;
    string _program_id;

    TL_SFFHeader _header;

    char_a_t _flow_chars;
    char_a_t _key_sequence;
    map < string, uint32_t > _read_map;

    uint32_t _read_offset;
    uint32_t _num_samples;
    uint32_t _num_bases;
};

TL_SFFTraceDataReader :: TL_SFFTraceDataReader (
                                            const string & Name,
                                            const string & ProgramID
)
:   _name ( Name )
,   _program_id ( ProgramID )
,   _flow_chars ()
,   _key_sequence ()
,   _read_map ()
,   _read_offset ( 0 )
,   _num_samples ( 0 )
,   _num_bases ( 0 )
{
}   /* TL_SFFTraceDataReader :: TL_SFFTraceDataReader () */

TL_SFFTraceDataReader :: ~TL_SFFTraceDataReader ()
{
}   /* TL_SFFTraceDataReader :: ~TL_SFFTraceDataReader () */

void
TL_SFFTraceDataReader :: Read ( const TL_TraceFile & File )
{
    ReadCheckHeader ( File );
}   /* TL_SFFTraceDataReader :: Read () */

uint16_t 
_get_short ( unsigned char * Ptr )
{
    uint16_t Var;
    memmove ( & Var, Ptr, sizeof ( uint16_t ) );

    return TL_McArch :: FromNet ( Var );
}   /* _get_short () */

uint32_t 
_get_long ( unsigned char * Ptr )
{
    uint32_t Var;
    memmove ( & Var, Ptr, sizeof ( uint32_t ) );

    return TL_McArch :: FromNet ( Var );
}   /* _get_short () */

static const uint16_t _scFlowgramBytes = 2;

void
TL_SFFTraceDataReader :: ReadCheckHeader ( const TL_TraceFile & File )
{
    if ( File . ContentSize () < sizeof ( TL_SFFHeader ) ) {
        throw TL_Exception ( "Invalid SFF file size" );
    }

    memmove ( & _header, File . Content (), sizeof ( TL_SFFHeader ) );

    _header . version = TL_McArch :: ToNet ( _header . version );
    if ( _header . version != 1 ) {
        throw TL_Exception ( "Invalid SFF file format version" );
    }

    if ( _header . flowgram_format_code != 1 ) {
        throw TL_Exception ( "Invalid SFF file format code" );
    }

    // index_offset
    // index_size

    _header . num_reads = TL_McArch :: FromNet ( _header . num_reads );
    _header . header_length = TL_McArch :: FromNet ( _header . header_length );
    _header . key_length = TL_McArch :: FromNet ( _header . key_length );
    _header . num_flows_per_read = TL_McArch :: FromNet ( _header . num_flows_per_read );

    size_t hdrSize = sizeof ( struct TL_SFFHeader );
    hdrSize += _header . num_flows_per_read + _header . key_length;

                                              // Volodya comment
    hdrSize = ( ( hdrSize + 7u ) / 8u ) * 8u; // bump to next 8-byte
                                              // boundary

    if ( File . ContentSize () < hdrSize ) {
        throw TL_Exception ( "Invalid SFF file size" );
    }

    if ( _header . header_length != ( int32_t ) hdrSize ) {
        throw TL_Exception ( "Invalid SFF file header lenght" );
    }

    unsigned char * Ptr = ( unsigned char * ) File . Content ();

        /*  Reading actual data
         */
    _flow_chars . set ( _header . num_flows_per_read, 0 );
    memmove (
            * _flow_chars,
            Ptr + sizeof ( struct TL_SFFHeader ),
            _header . num_flows_per_read
            );
    for ( uint32_t i = 0; i < _flow_chars . size (); i ++ ) {
        _flow_chars [ i ] = tolower ( _flow_chars [ i ] );
    }

    _key_sequence . set ( _header . key_length, 0 );
    memmove (
            * _key_sequence,
            Ptr + sizeof ( struct TL_SFFHeader )
                            + _header . num_flows_per_read,
            _header . key_length
            );
    for ( uint32_t i = 0; i < _key_sequence . size (); i ++ ) {
        _key_sequence [ i ] = tolower ( _key_sequence [ i ] );
    }

        // Volydya comments
        // iterate through reads
    uint32_t Size = File . ContentSize ();
    size_t Offset = hdrSize;
    size_t NumReads = _header . num_reads;
    while ( NumReads != 0 ) {
        NumReads --;

            // read header sanity checks
        if ( Size <= Offset + 16 ) {
            throw TL_Exception ( "Invalid SFF file size" );
        }

        uint16_t NameLength = _get_short ( Ptr + Offset + 2 );
        if ( Size <= Offset + 16 + NameLength ) {
            throw TL_Exception ( "Invalid SFF file size" );
        }

            // insert name and offset into the map
        string Name ( ( char * ) ( Ptr + Offset + 16 ), NameLength );
        _read_map [ Name ] = Offset;

            // Volodya comment ... and what to do then ???
            /// todo: check if not unique

        Offset += ( ( 16 + NameLength + 7u ) / 8u ) * 8u;
        Offset += _header . num_flows_per_read * _scFlowgramBytes;
        Offset += 3 * _get_long ( Ptr + Offset + 4 );

        Offset = ( ( Offset + 7u ) / 8u ) * 8u;
    }


    map < string, uint32_t > :: const_iterator It =
                                            _read_map . find ( _name );
    if ( It == _read_map . end () ) {
        throw TL_Exception ( string ( "Can not find '" ) + _name + "' entry in SFF stream" );
    }

    _read_offset = It -> second;

    _num_bases = _get_long ( Ptr + _read_offset + 4 );
    _num_samples = _header . num_flows_per_read;

    _clip_quality_left = _get_short ( Ptr + _read_offset + 8 );
    _clip_quality_right = _get_short ( Ptr + _read_offset + 10 );
}   /* TL_SFFTraceDataReader :: ReadCheckHeader () */

void
TL_SFFTraceDataReader :: ReadSamples ( const TL_TraceFile & File )
{
    unsigned char * Ptr = ( unsigned char * ) File . Content ();

    _samples_A . set ( _num_samples, 0 );
    _samples_C . set ( _num_samples, 0 );
    _samples_G . set ( _num_samples, 0 );
    _samples_T . set ( _num_samples, 0 );

    uint16_t Offset = _get_short ( Ptr + _read_offset );
    unsigned char * Data = Ptr + _read_offset + Offset;

    for ( uint32_t i = 0; i < _num_samples; i ++ ) {
        switch ( _flow_chars [ i ] ) {
            case 'a': _samples_A [ i ] = _get_short ( Data ); break;
            case 'c': _samples_C [ i ] = _get_short ( Data ); break;
            case 'g': _samples_G [ i ] = _get_short ( Data ); break;
            case 't': _samples_T [ i ] = _get_short ( Data ); break;
        }

        Data += sizeof ( uint16_t );
    }
}   /* TL_SFFTraceDataReader :: ReadSamples () */

void
TL_SFFTraceDataReader :: ReadBases ( const TL_TraceFile & File )
{
    unsigned char * Ptr = ( unsigned char * ) File . Content ();

    _prob_A . set ( _num_bases, 0 );
    _prob_C . set ( _num_bases, 0 );
    _prob_G . set ( _num_bases, 0 );
    _prob_T . set ( _num_bases, 0 );
    _bases . set ( _num_bases, 0 );
    _peak_index . set ( _num_bases, 0 );

    uint16_t Offset = _get_short ( Ptr + _read_offset );
    unsigned char * Data = Ptr + _read_offset + Offset;
    Data += _num_samples * sizeof ( uint16_t );

    int32_t Index = -1;
    for ( uint32_t i = 0; i < _num_bases; i ++ ) {
        Index += ( int32_t ) * Data;
        _peak_index [ i ] = ( uint32_t ) Index;
        Data ++;
    }

    for ( uint32_t i = 0; i < _num_bases; i ++ ) {
        _bases [ i ] = tolower ( * Data );
        Data ++;
    }

    _valid_scores = true;
    for ( uint32_t i = 0; i < _num_bases; i++ ) {
        switch ( _bases [ i ] ) {
            case 'a' : _prob_A [ i ] = * Data;
            case 'c' : _prob_C [ i ] = * Data;
            case 'g' : _prob_G [ i ] = * Data;
            case 't' : _prob_T [ i ] = * Data;
        }
        Data ++;
    }
}   /* TL_SFFTraceDataReader :: ReadBases () */

void
TL_SFFTraceDataReader :: ReadComments ( const TL_TraceFile & File )
{
    stringstream Str;

    Str << "NAME=" << _name << endl;
    if ( ! Str ) {
        throw TL_Exception ( "Can not read comment" );
    }

    Str << "CONV=" << _program_id << endl;
    if ( ! Str ) {
        throw TL_Exception ( "Can not read comment" );
    }

    Str << "SFFVersion=0001" << endl;
    if ( ! Str ) {
        throw TL_Exception ( "Can not read comment" );
    }

    Str << "SFFKey=" << string ( * _key_sequence, _key_sequence . size () ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "Can not read comment" );
    }

    Str << "SFFFlowChars=" << string ( * _flow_chars, _flow_chars . size () ) << endl;
    if ( ! Str ) {
        throw TL_Exception ( "Can not read comment" );
    }

    _comments = Str . str ();
}   /* TL_SFFTraceDataReader :: ReadComments () */

};  /* _tl_ namespace */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * Point of noreturn
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

void
_TraceDataReadSFF ( const TL_TraceFile & File, TL_Traces & Traces )
{
    TL_SFFTraceDataReader Reader (
                                Traces . Name (),
                                Traces . ProgramID ()
                                );

    Reader . Read ( File );

    Reader . Export ( Traces );
}   /* _TraceDataReadSFF () */

