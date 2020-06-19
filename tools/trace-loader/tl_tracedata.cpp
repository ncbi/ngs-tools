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

#include <sstream>
#include <iostream>

#include <kfs/directory.h>
#include <kfs/file.h>
#include <kfs/mmap.h>
#include <kfs/tar.h>
#include <kfs/bzip.h>
#include <klib/rc.h>
#include <klib/namelist.h>

#include "tl_tracedata.hpp"
#include "tl_traceinfo.hpp"
#include "tl_names.hpp"
#include "tl_log.hpp"


using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceData
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_TraceData :: TL_TraceData ()
:   _config ( NULL )
{
}   /* TL_TraceData :: TL_TraceData () */

TL_TraceData :: ~TL_TraceData ()
{
    _config = NULL;
}   /* TL_TraceData :: ~Tl_TraceData () */

void
TL_TraceData :: Init ( TL_TraceConfig * Config )
{
    _config = Config;
}   /* TL_TraceData :: Init () */

void
TL_TraceData :: Dispose ()
{
    _config = NULL;
}   /* TL_TraceData :: Dispose () */

void
TL_TraceData :: LoadFile (
                        TL_TraceFile & File,
                        const string & ThePath
) const
{
    if ( ! IsGood () ) {
        throw TL_Exception ( "TL_TraceData :: LoadFile (): object in invalid state" );
    }

    string Path = _config -> NormalizePath ( ThePath );

    if ( ! TL_Io :: IsFile ( Path ) ) {
        throw TL_Exception ( string ( "TL_TraceData :: FileSize (): can not stat file '" ) + Path + "'" );
    }

    struct KDirectory * NatDir = NULL;
    rc_t RCt = KDirectoryNativeDir ( & NatDir );
    if ( RCt == 0 ) {
        const struct KFile * TheFile = NULL;

        RCt = KDirectoryOpenFileRead (
                                ( const struct KDirectory * ) NatDir,
                                &TheFile,
                                Path . c_str ()
                                );
        if ( RCt == 0 ) {
            File . Init ( TheFile );

            KFileRelease ( TheFile );
        }

        KDirectoryRelease ( NatDir );
    }

    if ( RCt != 0 ) {
        throw TL_Exception ( RCt );
    }


}   /* TL_TraceData :: LoadFile () */

bool
TL_TraceData :: FileExists ( const string & ThePath ) const
{
    if ( ! IsGood () ) {
        throw TL_Exception ( "TL_TraceData :: FileExists (): object in invalid state" );
    }

    return TL_Io :: IsFile ( _config -> NormalizePath ( ThePath ) );
}   /* TL_TraceData :: FileExists () */

uint64_t
TL_TraceData :: FileSize ( const string & ThePath ) const
{
    if ( ! IsGood () ) {
        throw TL_Exception ( "TL_TraceData :: FileSize (): object in invalid state" );
    }

    string Path = _config -> NormalizePath ( ThePath );

    if ( ! TL_Io :: IsFile ( Path ) ) {
        throw TL_Exception ( string ( "TL_TraceData :: FileSize (): can not stat file '" ) + Path + "'" );
    }

    uint64_t RetVal = 0;

    struct KDirectory * NatDir = NULL;
    rc_t RCt = KDirectoryNativeDir ( & NatDir );

    if ( RCt == 0 ) {
        int PathType = KDirectoryPathType ( NatDir, Path . c_str () );
        if ( PathType == kptFile ) {
            RCt = KDirectoryFileSize ( NatDir, & RetVal, Path . c_str () );
        }
        else {
            RetVal = RC ( rcApp, rcPath, rcAccessing, rcDirEntry, rcInvalid );
        }

        KDirectoryRelease ( NatDir );
    }

    if ( RCt != 0 ) {
        throw TL_Exception ( RCt );
    }

    return RetVal;
}   /* TL_TraceData :: FileSize () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_TraceFile
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_TraceFile :: TL_TraceFile ()
:   _content ( NULL )
,   _content_size ( 0 )
,   _allocated_content ( false )
,   _map ( NULL )
{
}   /* TL_TraceFile :: TL_TraceFile () */

TL_TraceFile :: TL_TraceFile ( const struct KFile * File )
:   _content ( NULL )
,   _content_size ( 0 )
,   _allocated_content ( false )
,   _map ( NULL )
{
    Init ( File );
}   /* TL_TraceFile :: TL_TraceFile () */

TL_TraceFile :: ~TL_TraceFile ()
{
    TL_TRY {
        Dispose ();
    }
    TL_CATCH_R
}   /* TL_TraceFile :: ~TL_TraceFile () */

void
TL_TraceFile :: Init ( const struct KFile * File )
{
    Dispose ();

    if ( File == NULL ) {
        throw TL_Exception ( "TraceFile: NULL file passed for Init()" );
    }

    if ( ! TL_GZipReader () . Uncompress ( File, & _content, & _content_size ) ) {
        const struct KMMap * Map = NULL;
        rc_t RCt = KMMapMakeRead ( & Map, File );
        if ( RCt != 0 ) {
            throw TL_Exception ( "Can not map file", RCt );
        }

        RCt = KMMapSize ( Map, & _content_size );
        if ( RCt == 0 ) {
            RCt = KMMapAddrRead ( Map, ( const void ** ) & _content );
            if ( RCt == 0 ) {
                _map = Map;
            }
        }

        if ( RCt != 0 ) {
            _content_size = 0;
            _content = NULL;
            _map = NULL;
            KMMapRelease ( Map );

            throw TL_Exception ( "Can not init TraceFile" );
        }
    }
    else {
        // This is a trap
    }

    _detectFileType ();
}   /* TL_TraceFile :: Init () */

void
TL_TraceFile :: Dispose ()
{
    if ( _allocated_content ) {
        _content_size = 0;

        if ( _content != NULL ) {
            char * Ptr = ( char * ) _content;
            delete [] Ptr;
            _content = NULL;
        }
    }

    if ( _map != NULL ) {
        KMMapRelease ( _map );
        _map = NULL;
    }
    
}   /* TL_TraceFile :: Dispose () */

static const struct TL_FormatDesc {
                TL_TraceFile :: FType _type;
                uint32_t _offset;
                const char _magic [ 15 ];
                string _description;
} _scKnownFormats [] = {
      { TL_TraceFile :: UNK,     0, "", "UNKNOWN" }            /* It is always first!!! */
    , { TL_TraceFile :: SCF,     0, ".scf", "SCF" }
    , { TL_TraceFile :: ZTR,     0, "\256ZTR\r\n\032\n", "ZTR" }
    , { TL_TraceFile :: ABI,     0, "ABIF", "ABI" }
    , { TL_TraceFile :: ABI,   128, "ABIF", "ABI" }
    , { TL_TraceFile :: SFF,     0, ".sff", "SCF" }
    , { TL_TraceFile :: RCF,     0, ".rcf", "RCF" }
    , { TL_TraceFile :: CTF,     1, "\007\375\343\000", "CTF" } /* mieg */
    , { TL_TraceFile :: TTFF,    0, "\364TFF\n\r\032\n", "TTFF" }
    , { TL_TraceFile :: ALF,   518, "ALF ", "ALF" }
    , { TL_TraceFile :: SCF,     0, "\234\330\300\000", "SCF" } /* Amersham variant */
    , { TL_TraceFile :: EXP,     0, "ID   ", "EXP" }
    , { TL_TraceFile :: ALF,     0, "ALF ", "ALF" }             /* Added by newer alfsplit programs */
    , { TL_TraceFile :: ALF,     0, "\021G\021G", "ALF" }       /* Pharmacia's alfsplit equiv */
    , { TL_TraceFile :: ALF,  1546, "X-axis", "ALF" }           /* Good guestimation if all else fails */
};

void
TL_TraceFile :: _detectFileType ()
{
    _type = UNK;
    _description = "UNKNOWN";

    char * Ptr = ( char * ) _content;
    size_t Size = _content_size;

    for ( uint32_t i = 1; i < sizeof ( _scKnownFormats ) / sizeof ( struct TL_FormatDesc ); i ++ ) {
        size_t MagicLen = strlen ( _scKnownFormats [ i ] . _magic );
        size_t MagicOff = _scKnownFormats [ i ] . _offset;

        if ( Size < MagicLen + MagicOff ) {
            continue;
        }

        if ( memcmp (
                    Ptr + MagicOff,
                    _scKnownFormats [ i ] . _magic,
                    MagicLen
                    ) == 0
        ) {
            _type = _scKnownFormats [ i ] . _type;
            _description = _scKnownFormats [ i ] . _description;

            break;
        }
    }
}   /* TL_TraceFile :: _detectFileType () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_Bases
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_Bases :: TL_Bases ()
:   _bases ()
,   _peak_index ()
,   _prob_A ()
,   _prob_C ()
,   _prob_G ()
,   _prob_T ()
,   _valid_scores ( false )
,   _prob_sub ()
,   _prob_ins ()
,   _prob_del ()
,   _extra_probs ( false )
,   _prob ()
,   _bases_20 ( 0 )
,   _bases_40 ( 0 )
,   _bases_60 ( 0 )
{
}   /* TL_Bases :: TL_Bases () */

TL_Bases :: ~TL_Bases ()
{
    Reset ();
}   /* TL_Bases :: ~TL_Bases () */

#define BOUILLABAISSE(Name)     if ( Name != NULL ) { delete [] Name; Name = NULL; }

void 
TL_Bases :: Reset ()
{
    _bases . reset ();
    _peak_index . reset ();
    _prob_A . reset ();
    _prob_C . reset ();
    _prob_G . reset ();
    _prob_T . reset ();
    _valid_scores = false;
    _prob_sub . reset ();
    _prob_ins . reset ();
    _prob_del . reset ();
    _extra_probs = false;
    _prob . reset ();
    _bases_20 = 0;
    _bases_40 = 0;
    _bases_60 = 0;
}   /* TL_Bases :: Reset () */

void
TL_Bases :: Prepare ()
{
        /*  UnPreparing first
         */
    _bases_20 = 0;
    _bases_40 = 0;
    _bases_60 = 0;
    _prob . reset ();

        /*  Preparing second
         */
    size_t NumBases = _bases . size ();
    if ( NumBases == 0 ) {
        /*  It is not an error
        throw TL_Exception ( "TL_Bases :: Prepare (): Empty bases array" );
        */
        return;
    }

    if (    _prob_A . size () != NumBases
        ||  _prob_C . size () != NumBases
        ||  _prob_G . size () != NumBases
        ||  _prob_T . size () != NumBases
    ) {
        throw TL_Exception ( "TL_Bases :: MakeProbScores (): Ivalid prob score quantities" );
    }

    _prob . set ( NumBases, 0 );

    for ( size_t llp = 0; llp < NumBases; llp ++ ) {
        switch ( _bases [ llp ] ) {
            case 'A': case 'a': _prob [ llp ] = _prob_A [ llp ]; break;
            case 'C': case 'c': _prob [ llp ] = _prob_C [ llp ]; break;
            case 'G': case 'g': _prob [ llp ] = _prob_G [ llp ]; break;
            case 'T': case 't': _prob [ llp ] = _prob_T [ llp ]; break;
        }

            /*  As Volodya said  */
        if ( 127 < _prob [ llp ] ) _prob [ llp ] = 0;

            /* bases statistic */
        if ( ( unsigned char ) 20 <= _prob [ llp ] ) _bases_20 ++;
        if ( ( unsigned char ) 40 <= _prob [ llp ] ) _bases_40 ++;
        if ( ( unsigned char ) 60 <= _prob [ llp ] ) _bases_60 ++;
    }
}   /* TL_Bases :: Prepare () */

void
TL_Bases :: ReadBaseFile ( const TL_TraceFile & File )
{
    const char * Orig = ( const char * ) File . Content ();
    size_t Size = File . ContentSize ();

        /*  First line is always comment, and actual data starts after
         */
    const char * Start = Orig;
    const char * End = Orig + Size;
    const char * Pos = Start;

    while ( Pos < End ) {
        if ( * Pos == '\n' ) {
            break;
        }
        Pos ++;
    }

    if ( End <= Pos ) {
        throw TL_Exception ( "Invalid FASTA file" );
    }

        /*  Skippin' '\n'
         */
    Pos ++;

        /*  WARNING: somebody will think that it is cool to
         *           upack bases directly to _bases array.
         *           It is not right, cuz there could be symbols
         *           to skip ( '\n' for example ). 
         *           So let it be as is.
         */

    char * Bases = new char [ End - Pos ];
    char * BasesPtr = Bases;

        /*  Some majic transformations here
         *  After conversation with Kurt we decided to keep all 
         *  bases character, and not substitute them for 'n'
         */
#ifdef _N_SUBSTITUTE_
    while ( Pos < End ) {
        switch ( * Pos ) {
            case 'A': case 'a':
            case 'C': case 'c':
            case 'G': case 'g':
            case 'T': case 't':
                    * BasesPtr = (char) tolower ( * Pos );
                    BasesPtr ++;
                    break;
            case '\n':
                    break;
            default:
                    * BasesPtr = 'n';
                    BasesPtr ++;
                    break;
        }
        Pos ++;
    }
#else /* _N_SUBSTITUTE_ */
    while ( Pos < End ) {
        if ( * Pos != '\n' ) {
            if ( isalpha ( * Pos ) ) {
                * BasesPtr = (char) tolower ( * Pos );
            }
            else {
                * BasesPtr = 'n';
            }
            BasesPtr ++;
        }
        Pos ++;
    }
#endif /* _N_SUBSTITUTE_ */

    _bases . set ( BasesPtr - Bases, Bases );

    delete [] Bases;
}   /* TL_Bases :: ReadBaseFile () */

/*  Just read all ints and fill a vector 
 *  We need it to read Peak and Qual files.
 */
static
size_t
_SplitLineBy ( const char * Start, const char * End, TL_IVec & Vec )
{
    Vec . clear ();

    const char * Pos = Start;

    while ( Pos < End ) {
            /*  Moving to first non ' ' or '\n' character 
             */
        if ( * Pos != ' ' && * Pos != '\n' ) {
            Vec . insert ( Vec . end (), atoi ( Pos ) );

                /*  Moving to first ' ' of '\n' character 
                 */
            while ( Pos < End ) {
                if ( * Pos == ' ' || * Pos == '\n' ) {
                    break;
                }
                Pos ++;
            }
        }

        Pos ++;
    }

    return Vec . size ();
}   /* _SplitLineBy () */

void
TL_Bases :: ReadQualFile ( const TL_TraceFile & File )
{
    if ( ! _bases ) {
        throw TL_Exception ( "ReadQualFile (): Bases aren't set yet" );
    }

    const char * Orig = ( const char * ) File . Content ();
    size_t Size = File . ContentSize ();

        /*  First line is always comment, and actual data starts after
         */
    const char * Start = Orig;
    const char * End = Orig + Size;
    const char * Pos = Start;

    while ( Pos < End ) {
        if ( * Pos == '\n' ) {
            break;
        }
        Pos ++;
    }

    if ( End <= Pos ) {
        throw TL_Exception ( "Invalid QUAL file" );
    }

    TL_IVec Vec;
    if ( _SplitLineBy ( Pos, End, Vec ) != _bases . size () ) {
        throw TL_Exception ( "Invalid QUALITY qty" );
    }

    _prob_A . set ( _bases . size (), 0 );
    _prob_C . set ( _bases . size (), 0 );
    _prob_G . set ( _bases . size (), 0 );
    _prob_T . set ( _bases . size (), 0 );

    for ( size_t llp = 0; llp < _bases . size (); llp ++ ) {
        switch ( _bases [ llp ] ) {
            case 'A': case 'a':
                _prob_A [ llp ] = ( unsigned char ) Vec [ llp ];
                break;
            case 'C': case 'c':
                _prob_C [ llp ] = ( unsigned char ) Vec [ llp ];
                break;
            case 'G': case 'g':
                _prob_G [ llp ] = ( unsigned char ) Vec [ llp ];
                break;
            case 'T': case 't':
                _prob_T [ llp ] = ( unsigned char ) Vec [ llp ];
                break;
        }
    }

    _valid_scores = true;
}   /* TL_Base :: ReadQualFile () */

void
TL_Bases :: ReadPeakFile ( const TL_TraceFile & File )
{
    if ( ! _bases ) {
        throw TL_Exception ( "ReadQualFile (): Bases aren't set yet" );
    }

    const char * Orig = ( const char * ) File . Content ();
    size_t Size = File . ContentSize ();

        /*  First line is always comment, and actual data starts after
         */
    const char * Start = Orig;
    const char * End = Orig + Size;
    const char * Pos = Start;

    while ( Pos < End ) {
        if ( * Pos == '\n' ) {
            break;
        }
        Pos ++;
    }

    if ( End <= Pos ) {
        throw TL_Exception ( "Invalid PEAK file" );
    }

    TL_IVec Vec;
    if ( _SplitLineBy ( Pos, End, Vec ) != _bases . size () ) {
        throw TL_Exception ( "Invalid PEAK qty" );
    }

    _peak_index . set ( _bases . size (), 0 );
    for ( size_t llp = 0; llp < _bases . size (); llp ++ ) {
        _peak_index [ llp ] = ( uint32_t ) Vec [ llp ];
    }
}   /* TL_Base :: ReadPeakFile () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_Samples
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_Samples :: TL_Samples ()
:   _sample_A ()
,   _sample_C ()
,   _sample_G ()
,   _sample_T ()
,   _sample_combined ()
,   _max_trace_val ( 0 )
,   _flow_chars ()
,   _key_sequence ()
{
}   /* TL_Samples :: TL_Samples () */

TL_Samples :: TL_Samples ( const TL_Samples & Samples )
:   _sample_A ( Samples . _sample_A )
,   _sample_C ( Samples . _sample_C )
,   _sample_G ( Samples . _sample_G )
,   _sample_T ( Samples . _sample_T )
,   _sample_combined ( Samples . _sample_combined )
,   _max_trace_val ( Samples . _max_trace_val )
,   _flow_chars ( Samples . _flow_chars )
,   _key_sequence ( Samples . _key_sequence )
{
}   /* TL_Samples :: TL_Samples () */

TL_Samples :: ~TL_Samples ()
{
    TL_TRY {
        Reset ();
    }
    TL_CATCH_R
}   /* TL_Samples :: ~TL_Samples () */

TL_Samples &
TL_Samples :: operator = ( const TL_Samples & Samples )
{
    if ( this != & Samples ) {
        _sample_A = Samples . _sample_A;
        _sample_C = Samples . _sample_C;
        _sample_G = Samples . _sample_G;
        _sample_T = Samples . _sample_T;
        _sample_combined = Samples . _sample_combined;
        _max_trace_val = Samples . _max_trace_val;
        _flow_chars = Samples . _flow_chars;
        _key_sequence = Samples . _key_sequence;
    }

    return * this;
}   /* TL_Samples :: ooperator = () */

void
TL_Samples :: Reset ()
{
    _sample_A . reset ();
    _sample_C . reset ();
    _sample_G . reset ();
    _sample_T . reset ();
    _sample_combined . reset ();
    _max_trace_val = 0;
    _flow_chars . reset ();
    _key_sequence . reset ();
}   /* TL_Samples :: Reset () */

void
TL_Samples :: Prepare ()
{
        /*  UnPrepare first
         */
    _max_trace_val = 0;
    _sample_combined . reset ();

        /*  Prepare second
         */
    size_t Size = _sample_A . size ();
    if (   Size != _sample_C . size ()
        || Size != _sample_G . size ()
        || Size != _sample_T . size ()
    ) {
        throw TL_Exception ( "TL_Bases :: MakeMaxTraceVal (): Ivalid sample quantities" );
    }

    _max_trace_val = 0;
    _sample_combined . set ( Size, 0 );

    union {
        uint16_t in [ 4 ];
        uint64_t ou;
    } WooHoo;

    for ( size_t llp = 0; llp < Size; llp ++ ) {
        if ( _max_trace_val < _sample_A [ llp ] )
            _max_trace_val = _sample_A [ llp ];
        if ( _max_trace_val < _sample_C [ llp ] )
            _max_trace_val = _sample_C [ llp ];
        if ( _max_trace_val < _sample_G [ llp ] )
            _max_trace_val = _sample_G [ llp ];
        if ( _max_trace_val < _sample_T [ llp ] )
            _max_trace_val = _sample_T [ llp ];

        WooHoo . in [ 0 ] = _sample_A [ llp ];
        WooHoo . in [ 1 ] = _sample_C [ llp ];
        WooHoo . in [ 2 ] = _sample_G [ llp ];
        WooHoo . in [ 3 ] = _sample_T [ llp ];
        _sample_combined [ llp ] = WooHoo . ou;

    }
}   /* TL_Samples :: Prepare () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_Traces
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_Traces :: TL_Traces ()
:   _name ( "" )
,   _program_id ( "" )
,   _clip_quality_right ( 0 )
,   _clip_quality_left ( 0 )
,   _clip_adapter_right ( 0 )
,   _clip_adapter_left ( 0 )
,   _override_peaks ( false )
,   _override_quals ( false )
,   _format ( 0 )
,   _bases ()
,   _samples ()
,   _comments ( "" )
,   _private_data ()
{
}   /* TL_Traces :: TL_Traces () */

TL_Traces :: ~TL_Traces ()
{
    TL_TRY {
        Reset ();
    }
    TL_CATCH_R
}   /* TL_Traces :: ~TL_Traces () */

void
TL_Traces :: Reset ()
{
    _clip_quality_right = 0;
    _clip_quality_left = 0;
    _clip_adapter_right = 0;
    _clip_adapter_left = 0;

    _override_peaks = false;
    _override_quals = false;

    _format = 0;

    _bases . Reset ();
    _samples . Reset ();

    _comments . clear ();

    _private_data . reset ();

    _program_id . clear ();
    _name . clear ();
}   /* TL_Traces :: Reset () */

void
TL_Traces :: Prepare ()
{
    _bases . Prepare ();
    _samples . Prepare ();
}   /* TL_Traces :: Prepare () */

void _TraceDataReadZTR ( const TL_TraceFile & File, TL_Traces & Trace );
void _TraceDataReadABI ( const TL_TraceFile & File, TL_Traces & Trace );
void _TraceDataReadSCF ( const TL_TraceFile & File, TL_Traces & Trace );
void _TraceDataReadSFF ( const TL_TraceFile & File, TL_Traces & Trace );

void
TL_Traces :: ReadTraces ( 
                        const string & Name,
                        const string & ProgramID,
                        const TL_TraceFile & File
)
{
    Reset ();

    _name = Name;
    _program_id = ProgramID;

    PLOGMSG(
        klogDebug,
        (klogDebug,
        "[TRACE] NAME [$(name)] PROGRAM_ID [$(prog_id)] FORMAT [$(format)]",
        "severity=debug,name=%s,prog_id=%s,format=%s",
        _name . c_str(),
        _program_id . c_str (),
        File . Description () . c_str ()
        )
    );

    switch ( File . Type () ) {
        case TL_TraceFile :: ZTR:
            _TraceDataReadZTR ( File, * this );
            break;
        case TL_TraceFile :: ABI:
            _TraceDataReadABI ( File, * this );
            break;
        case TL_TraceFile :: SCF:
            _TraceDataReadSCF ( File, * this );
            break;
        case TL_TraceFile :: SFF:
            _TraceDataReadSFF ( File, * this );
            break;
        case TL_TraceFile :: UNK:
            throw TL_Exception ( "UNKNOWN trace format" );
            break;
        default:
            throw TL_Exception ( string ( "Trace format \"" ) + File . Description () + "\" is not supported" );
    }
}   /* TL_Traces :: ReadTraces () */


/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_Trace_DP
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_Trace_DP :: TL_Trace_DP ()
:   TL_TraceDataProvider ()
,   _bases ()
,   _traces ()
,   _name ( "" )
,   _program_id ( "" )
{
}   /* TL_Trace_DP :: TL_Trace_DP () */

TL_Trace_DP :: ~TL_Trace_DP ()
{
    Reset ();
}   /* TL_Trace_DP :: ~TL_Trace_DP () */

void
TL_Trace_DP :: Set (
                    const TL_TraceData & TraceData,
                    const TL_TraceInfo & TraceInfo,
                    size_t RowIdx
)
{
    _name = TraceInfo . Value ( RowIdx, _TL_TRACE_NAME );
    if ( _name . empty () ) {
        throw TL_Exception ( "TL_Trace_DP: Empty name for row" );
    }

    _program_id = TraceInfo . Value ( RowIdx, _TL_PROGRAM_ID );
    if ( _program_id . empty () ) {
        throw TL_Exception ( "TL_Trace_DP: Empty program_id for row" );
    }

    string Val = TraceInfo . Value ( RowIdx, _TL_BASE_FILE );
    if ( ! Val . empty () ) {
        PLOGMSG(
                klogDebug,
                (klogDebug,
                "[LOADING] [BASE] FILE [$(file)]",
                "severity=debug,file=%s",
                Val . c_str()
                )
        );
        TL_TraceFile File;
        TraceData . LoadFile ( File, Val );
        _bases . ReadBaseFile ( File );

        Val = TraceInfo . Value ( RowIdx, _TL_QUAL_FILE );
        if ( ! Val . empty () ) {
            PLOGMSG(
                    klogDebug,
                    (klogDebug,
                    "[LOADING] [QUAL] FILE [$(file)]",
                    "severity=debug,file=%s",
                    Val . c_str()
                    )
            );
            TraceData . LoadFile ( File, Val );
            _bases . ReadQualFile ( File );

            _traces . OverrideQuals ( true );
        }

        Val = TraceInfo . Value ( RowIdx, _TL_PEAK_FILE );
        if ( ! Val . empty () ) {
            PLOGMSG(
                    klogDebug,
                    (klogDebug,
                    "[LOADING] [PEAK] FILE [$(file)]",
                    "severity=debug,file=%s",
                    Val . c_str()
                    )
            );
            TraceData . LoadFile ( File, Val );
            _bases . ReadPeakFile ( File );

            _traces . OverridePeaks ( true );
        }
    }

        /*  Here we are loading traces NOTE TRACELESS ???
         *  Under ordinary circumstances it is impossible situation
         */
        /* JOJOBA */
    Val = TraceInfo . Value ( RowIdx, _TL_TRACE_FILE );
    if ( ! Val . empty () ) {
        PLOGMSG(
                klogDebug,
                (klogDebug,
                "[LOADING] [TRACE] FILE [$(file)]",
                "severity=debug,file=%s",
                Val . c_str()
                )
        );
        TL_TraceFile File;
        TraceData . LoadFile ( File, Val );
        _traces . ReadTraces ( _name, _program_id, File );
    }

    _bases . Prepare ();
    _traces . Prepare ();

}   /* TL_Trace_DP :: Set () */

void
TL_Trace_DP :: Reset ()
{
    _name . clear ();
    _program_id . clear ();
    _traces . Reset ();
    _bases . Reset ();
}   /* TL_Trace_DP :: Reset () */

const std::string &
TL_Trace_DP :: Name () const
{
    return _name;
}   /* TL_Trace_DP :: Name () */

const std::string &
TL_Trace_DP :: ProgramID () const
{
    return _program_id;
}   /* TL_Trace_DP :: ProgramID () */

const char_a_t &
TL_Trace_DP :: Bases () const
{
    return _bases . Bases () . size ()
                                    ? _bases . Bases ()
                                    : _traces . Bases () . Bases ()
                                    ;
}   /* TL_Trace_DP :: Bases () */

const uint32_a_t &
TL_Trace_DP :: PeakIndex () const
{
    return _bases . PeakIndex () . size ()
                                    ? _bases . PeakIndex ()
                                    : _traces . Bases () . PeakIndex ()
                                    ;
}   /* TL_Trace_DP :: PeakIndex () */

const uchar_a_t &
TL_Trace_DP :: ProbScores () const
{
    return _bases . ProbScores () . size ()
                                    ? _bases . ProbScores ()
                                    : _traces . Bases () . ProbScores ()
                                    ;
}   /* TL_Trace_DP :: ProbScores () */

uint16_t
TL_Trace_DP :: Bases20 () const
{
    return _bases . ProbScores () . size ()
                                    ? _bases . Bases20 ()
                                    : _traces . Bases () . Bases20 ()
                                    ;
}   /* TL_Trace_DP :: Bases20 () */

uint16_t
TL_Trace_DP :: Bases40 () const
{
    return _bases . ProbScores () . size ()
                                    ? _bases . Bases40 ()
                                    : _traces . Bases () . Bases40 ()
                                    ;
}   /* TL_Trace_DP :: Bases40 () */

uint16_t
TL_Trace_DP :: Bases60 () const
{
    return _bases . ProbScores () . size ()
                                    ? _bases . Bases60 ()
                                    : _traces . Bases () . Bases60 ()
                                    ;
}   /* TL_Trace_DP :: Bases60 () */

const uint64_a_t &
TL_Trace_DP :: SampleCombined () const
{
    return _traces . Samples() . SampleCombined ();
}   /* TL_Trace_DP :: SampleCombined () */

uint16_t
TL_Trace_DP :: MaxTraceVal () const
{
    return _traces . Samples () . MaxTraceVal ();
}   /* TL_Trace_DP :: MaxTraceVal () */

uint32_t
TL_Trace_DP :: ClipQualityRight () const
{
    return _traces . ClipQualityRight ();
}   /* TL_Trace_DP :: ClipQualityRight () */

uint32_t
TL_Trace_DP :: ClipQualityLeft () const
{
    return _traces . ClipQualityLeft ();
}   /* TL_Trace_DP :: ClipQualityLeft () */

const char_a_t &
TL_Trace_DP :: FlowChars () const
{
    return _traces . Samples () . FlowChars ();
}   /* TL_Trace_DP :: FlowChars () */

const char_a_t &
TL_Trace_DP :: KeySequence () const
{
    return _traces . Samples () . KeySequence ();
}   /* TL_Trace_DP :: KeySequence () */

const string &
TL_Trace_DP :: Comments () const
{
    return _traces . Comments ();
}   /* TL_Trace_DP :: Comments () */

const string &
TL_Trace_DP :: ExtendedData () const
{
    return _traces . ExtendedData ();
}   /* TL_Trace_DP :: ExtendedData () */
