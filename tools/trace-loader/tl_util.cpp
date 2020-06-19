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
#include <fstream>
#include <algorithm>

#include "tl_util.hpp"
#include "tl_log.hpp"

#include <kfs/mmap.h>
#include <klib/rc.h>
#include <kfs/file.h>
#include <kfs/gzip.h>
#include <kfs/directory.h>

#include <ctype.h>  /* isspace () */

#include <arpa/inet.h>

using namespace std;
using namespace _tl_;

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_KfsMap 
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

TL_KfsMap :: TL_KfsMap ()
:   _handle ( NULL )
,   _path ( "" )
,   _size ( 0 )
,   _addr ( NULL )
{
}   /* TL_KfsMap :: TL_KfsMap () */

TL_KfsMap :: TL_KfsMap ( const TL_KfsMap & Map )
:   _handle ( NULL )
,   _path ( "" )
,   _size ( 0 )
,   _addr ( NULL )
{
    _DoMap ( Map );
}   /* TL_KfsMap :: TL_KfsMap () */

TL_KfsMap :: TL_KfsMap ( const string & Path )
:   _handle ( NULL )
,   _path ( "" )
,   _size ( 0 )
,   _addr ( NULL )
{
    _DoMap ( Path );
}   /* TL_KfsMap :: TL_KfsMap () */

TL_KfsMap :: ~TL_KfsMap ()
{
    TL_TRY {
        _DoUnMap ();
    }
    TL_CATCH_R
}   /* TL_KfsMap :: ~Tl_KfsMap () */

TL_KfsMap &
TL_KfsMap :: operator = ( const TL_KfsMap & Map )
{
    if ( this != & Map ) {
        _DoMap ( Map );
    }

    return * this;
}   /* TL_KfsMap :: operator = () */

TL_KfsMap &
TL_KfsMap :: operator = ( const string & Path )
{
    _DoMap ( Path );

    return * this;
}   /* TL_KfsMap :: operator = () */

void
TL_KfsMap :: Map ( const string & Path )
{
    _DoMap ( Path );
}   /* TL_KfsMap :: Map () */

void
TL_KfsMap :: UnMap ()
{
    _DoUnMap ();
}   /* TL_KfsMap :: UnMap () */

void
TL_KfsMap :: _DoMap ( const TL_KfsMap & Map )
{
    if ( IsMapped () ) {
        throw TL_Exception ( "Attempt to mmap already mmapped mmap" );
    }

    if ( Map . IsMapped () ) {
        const struct KMMap * kMap =
                                ( const struct KMMap * ) Map . _handle;

        rc_t RCt = KMMapAddRef ( kMap );
        if ( RCt != 0 ) {
            throw TL_Exception ( RCt );
        }

        _handle = ( void * ) kMap;
        _path = Map . _path;
        _size = Map . _size;
        _addr = Map . _addr;
    }

}   /* TL_KfsMap :: _DoMap () */

void
TL_KfsMap :: _DoMap ( const string & Path )
{
    if ( IsMapped () ) {
        throw TL_Exception ( "Attempt to mmap already mmapped mmap" );
    }

    struct KDirectory * NatDir = NULL;
    rc_t RCt = KDirectoryNativeDir ( & NatDir );
    if ( RCt == 0 ) {
        const struct KFile * File = NULL;
        RCt = KDirectoryOpenFileRead ( NatDir, & File, Path . c_str () );
        if ( RCt == 0 ) {
            const struct KMMap * Map = NULL;
            RCt = KMMapMakeRead ( & Map, File );
            if ( RCt == 0 ) {
                const void * Addr = NULL;
                RCt = KMMapAddrRead ( Map, & Addr );
                if ( RCt == 0 ) {
                    size_t Size = 0;
                    RCt = KMMapSize ( Map, & Size );
                    if ( RCt == 0 ) {
                        _handle = ( void * ) Map;
                        _path = Path;
                        _size = Size;
                        _addr = ( void * ) Addr;
                    }
                }

                if ( RCt != 0 ) {
                    KMMapRelease ( Map );
                }
            }

            KFileRelease ( File );
        }

        KDirectoryRelease ( NatDir );
    }

    if ( RCt != 0 ) {
        _handle = NULL;
        _size = 0;
        _path . clear ();
        _addr = NULL;

        throw TL_Exception ( RCt );
    }
}   /* TL_KfsMap :: _DoMap () */

void
TL_KfsMap :: _DoUnMap ()
{
    if ( IsMapped () ) {
        const struct KMMap * Map = ( const struct KMMap * ) _handle;
        KMMapRelease ( Map );
        _addr = NULL;
        _size = 0;
        _handle = NULL;
        _path . clear ();
    }
}   /* TL_KfsMap :: _DoUnMap () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_LineReader 
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

TL_LineReader :: TL_LineReader ( const string & Path )
:   _map ()
,   _line ( "" )
,   _line_no ( 0 )
,   _pos ( 0 )
{
    if ( ! Path . empty () ) {
        _DoOpen ( Path );
    }
}   /* TL_LineReader :: TL_LineReader () */

TL_LineReader :: ~TL_LineReader ()
{
    TL_TRY {
        _DoClose ();
    }
    TL_CATCH_R
}   /* TL_LineReader :: ~TL_LineReader () */

void
TL_LineReader :: Open ( const string & Path )
{
    if ( IsOpen () ) {
        throw TL_Exception ( "Attemp to open openedn LineReader" );
    }

    _DoOpen ( Path );
}   /* TL_LineReader :: Open () */

void
TL_LineReader :: Close ()
{
    if ( IsOpen () ) {
        _DoClose ();
    }
}   /* TL_LineReader :: Close () */

void
TL_LineReader :: Reset ()
{
    if ( ! IsOpen () ) {
        throw TL_Exception ( "Attempt to reset unopened reader" );
    }

    _line_no = 0;
    _line . clear ();
    _pos = 0;
}   /* TL_LineReader :: Reset () */

bool
TL_LineReader :: NextLine ()
{
    if ( ! IsOpen () ) {
        throw TL_Exception ( "Attempt to get next line from unopened reader" );
    }

        /*  Here we are looking to new line
         */
    const char * Bg = ( ( const char * ) _map . Get () ) + _pos;
    const char * En = Bg + ( _map . Size () - _pos );
    const char * Ps = Bg;

    while ( Ps < En ) {
        if ( * Ps == '\n' ) {
            break;
        }

        Ps ++;
    }

    if ( Ps < En || Ps != Bg ) {
        _line = string ( Bg, Ps - Bg );
        _line_no ++;
        _pos += ( Ps - Bg ) + 1;

        return true;
    }

    return false;
}   /* TL_LineReader :: NextLine () */

void
TL_LineReader :: _DoOpen ( const string & Path )
{
    if ( IsOpen () ) {
        throw TL_Exception ( "Attempt to open already opened reader" );
    }

    _map = Path;
    _line . clear ();
    _line_no = 0;
    _pos = 0;

}   /* TL_LineReader :: _DoOpen () */

void
TL_LineReader :: _DoClose ()
{
    _map . UnMap ();
    _line . clear ();
    _line_no = 0;
    _pos = 0;
}   /* TL_LineReader :: _DoClose () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_DateU 
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

/*)) Awkward class: there are two typese of time ( really four )
/\  described in RFC: run_date ( 2000-10-28) and
\/   collection_date ( Mar 2 2006 12:00AM ). However, I see it
/\  in many different formats, like run_date ( 05-Jan-2013 )
\/   and collection_date ( Thu Feb  6 16:33:00 2014 )
/\
\/   Ok, I run some script on TRACEINFO.tbl files for last couple
/\  years, and there is statistic:
\/
/\  run_dates are all in format :
\/       Fri Apr XX XX:XX:XX XXXX
/\      XX-APR-XXXX
\/   collect_dates are in format :
/\      XX-Aug-XXXX
\/       XXXX
/\
\/   Which is slightry different from RFC, so, I believe, I need
/\   to parse 6 different formats : 2 from RFC and 4 from real life
\/   Also, we may nod distinct run from collect date by format :/
((*/
string
TL_DateU :: RunDate ( uint64_t Time )
{
        /*  Do it as RFC prescribed : 2000-10-28
         */
    struct tm Tm;
    time_t TimeT = ( time_t ) Time;

    struct tm * RetTm = localtime_r ( & TimeT, & Tm );

    if ( RetTm == NULL ) {
            /* Invalid time */
        return TL_StringU :: EmptyString ();
    }

    char BB [ 16 ];

    if ( strftime ( BB, sizeof ( BB ), "%F", & Tm ) == 0 ) {
            /* Conversion error */
        return TL_StringU :: EmptyString ();
    }

    return BB;
}   /* TL_DateU :: RunDate () */

string
TL_DateU :: CollectionDate ( uint64_t Time )
{
        /*  Do it as RFC prescribed : Mar 2 2006 12:00AM
         */
    struct tm Tm;
    time_t TimeT = ( time_t ) Time;

    struct tm * RetTm = localtime_r ( & TimeT, & Tm );

    if ( RetTm == NULL ) {
            /* Invalid time */
        return TL_StringU :: EmptyString ();
    }

    char BB [ 32 ];

    if ( strftime ( BB, sizeof ( BB ), "%b %d %Y %I:%M%p", & Tm ) == 0 ) {
            /* Conversion error */
        return TL_StringU :: EmptyString ();
    }

    return BB;
}   /* TL_DateU :: CollectionDate () */

uint64_t
TL_DateU :: Date ( const string & Time )
{
    /*  Date could be in these formats :
     *      XXXX
     *      XXXX-XX-XX
     *      XX-APR-XXXX
     *      Apr XX XXXX XX:XXAM
     *      Fri Apr XX XX:XX:XX XXXX
     */

    string Str = TL_StringU :: Prune ( Time );

    TL_SVec Vec;
    size_t Q = TL_StringU :: Tokenize ( Time, ' ', Vec );

    const char * Format = NULL;
    if ( Q == 1 ) {
        if ( Str . length () == 4 ) {
                /*  First "XXXX" */
            Format = "%Y";
        }
        else {
            if ( Str . length () == 10 ) {
                    /*  First "XXXX-XX-XX" */
                Format = "%F";
            }
            else {
                if ( Str . length () == 11 ) {
                        /*  Second "XX-APR-XXXX" */
                    Format = "%d-%b-%Y";
                } 
            }
        }
    }
    else {
        if ( Q == 5 ) {
                /*  Last "Fri Apr XX XX:XX:XX XXXX" */
            Format = "%a %b %d %H:%M:%S %Y";
        }
        else {
            if ( Q == 4 ) {
                    /*  Lastest "Apr XX XXXX XX:XXAM" */
                Format = "%b %d %Y %I:%M%p";
            }
        }
    }

    if ( Format == NULL ) {
            /* Invalid format */
        return 0;
    }

    struct tm Tm;
    memset ( & Tm, 0, sizeof ( Tm ) ); /* Don't remove that :D */

    const char * P = strptime ( Str . c_str (), Format, & Tm );
    if ( P != Str . c_str () + Str . length () ) {
            /* Dat dis error */
        return 0;
    }

    return mktime ( & Tm );
}   /* TL_DateU :: Date () */

uint64_t
TL_DateU :: LocalToGMT ( uint64_t LocalTime )
{
    struct tm sTM;
    memset ( & sTM, 0, sizeof ( struct tm ) );

    if ( gmtime_r ( ( const time_t * ) & LocalTime, & sTM ) == NULL ) {
            /* Das ist errror :| */
        return - 1;
    }

    return mktime ( & sTM );
}   /* TL_DateU :: LocalToGMT () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_StringU 
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

TL_StringU :: TL_StringU ()
{
    throw TL_Exception ( "You are not allowed to create instance of that class" );
}   /* TL_StringU :: TL_StringU () */

string
TL_StringU :: Prune ( const string & Line, bool Front, bool End )
{
    if ( ! Front && ! End ) {
        return Line;
    }

    const char * Bg = Line . c_str ();
    const char * En = Bg + Line . length ();

        /*  Trimming front
         */
    if ( Front ) {
        while ( Bg < En ) {
            if ( ! isspace ( * Bg ) ) { break; }
            Bg ++;
        }
    }

        /*  Trimming end
         */
    if ( End ) {
        while ( Bg < En ) {
            if ( ! isspace ( * ( En - 1 ) ) ) { break; }
            En --;
        }
    }

    return string ( Bg, En - Bg );
}   /* TL_StringU :: Prune () */

static
string
_CutEndLine ( const string & Line )
{
    string Str = Line;

    string :: iterator Pos = Str . begin ();
    while ( Pos != Str . end () ) {
        char Ch = * Pos;
        if ( Ch == ';' || Ch == '\n' || Ch == '\r' || Ch == '\n' ) {
            * Pos = 0;
            break;
        }

        Pos ++;
    }

    return Str;
}   /* _CutEndLine () */

bool
TL_StringU :: Split (
                    const string & Line,
                    char Delim,
                    string & Key,
                    string & Value
)
{
    Key . clear ();
    Value . clear ();

    const char * Bg = Line . c_str ();
    const char * Ps = Bg;
    const char * En = Bg + Line . length ();

    while ( Ps < En ) {
        if ( * Ps == Delim ) {
            string TheKey = Prune ( string ( Bg, Ps - Bg ) );
            if ( TheKey . empty () ) {
                break;
            }

            Key = TheKey;
            Ps ++;
            Value = Prune ( _CutEndLine ( string ( Ps, En - Ps ) ) );

            return true;
        }

        Ps ++;
    }

    return false;
}   /* TL_StringU :: Split () */

size_t
TL_StringU :: Tokenize (
                        const string & Line,
                        char Delim,
                        TL_SVec & Tokens,
                        bool AddEmptyTokens
)
{
        /*  First we are clearing vector of Tokens
         */
    Tokens . clear ();

    const char * Bg = Line . c_str ();
    const char * Ps = Bg;
    const char * En = Ps + Line . length ();

    while ( Ps < En ) {
        if ( * Ps == Delim ) {
            if ( Bg != Ps || AddEmptyTokens ) {
                Tokens . insert (
                                Tokens . end (),
                                string ( Bg, Ps - Bg )
                                );
            }

            Bg = Ps + 1;
        }

        Ps ++;
    }

    if ( Bg != Ps ) {
        Tokens . insert ( Tokens . end (), string ( Bg, Ps - Bg ) );
    }

    return Tokens . size ();
}   /* TL_StringU :: Tokenize () */

string
TL_StringU :: ToUpper ( const string & Line )
{
    string RetVal = Line;
    std :: transform ( RetVal . begin (), RetVal . end (), RetVal . begin (), :: toupper );
    return RetVal;
}   /* TL_StringU :: ToUpper () */

string
TL_StringU :: ToLower ( const string & Line )
{
    string RetVal = Line;
    std :: transform ( RetVal . begin (), RetVal . end (), RetVal . begin (), :: tolower );
    return RetVal;
}   /* TL_StringU :: ToLower () */

string
TL_StringU :: ToPrintable ( const string & Line )
{
    stringstream Str;

    for (
        string :: const_iterator Cit = Line . begin ();
        Cit != Line . end ();
        Cit ++
    ) {
        if ( isprint ( * Cit ) ) {
            Str << * Cit;
        }
    }

    return Str . str ();
}   /* TL_StringU :: ToPrintable () */

string
TL_StringU :: FromBin ( const void * Data, size_t Size )
{
    stringstream Str;

    for ( size_t llp = 0; llp < Size; llp ++ ) {
        unsigned int F = ( ( unsigned char * ) Data ) [ llp ];
        Str << setfill ( '0' ) << setw ( 2 ) << hex << F;
        if ( ! Str ) {
            throw TL_Exception ( "Can not convert binary to string" );
        }
    }

    return Str . str ();
}   /* TL_StringU :: FromBin () */

void
TL_StringU :: ToBin ( const string & Src, void * Data, size_t Size )
{
    if ( Src . length () != Size << 1 ) {
        throw TL_Exception ( "Failed to convert string to binary" );
    }

    for ( size_t llp = 0; llp < Size; llp ++ ) {
        stringstream Str;
        Str << Src . substr ( llp * 2, 2 );
        unsigned int CH = 0;
        Str >> setw ( 2 ) >> hex >> CH;
        ( ( unsigned char * ) Data ) [ llp ] = CH;
    }
}   /* TL_StringU :: ToBin () */

bool
TL_StringU :: StartsWith ( const string & Line, const string & With )
{
    if ( With . length () <= Line . length () ) {
        return With == Line . substr ( 0, With . length () );
    }
    return false;
}   /* TL_StringU :: StartsWith () */

bool
TL_StringU :: EndsWith ( const string & Line, const string & With )
{
    if ( With . length () <= Line . length () ) {
        return With == Line . substr ( Line . length () - With . length () );
    }
    return false;
}   /* TL_StringU :: EndsWith () */

static const string _sStringU_emptyString ("");

const string &
TL_StringU :: EmptyString ()
{
    return _sStringU_emptyString;
}   /* TL_StringU :: EmptyString () */

bool
TL_StringU :: IsPrint ( const string & Line )
{
    for ( string :: const_iterator It = Line . begin ();
        It != Line . end ();
        It ++ 
    ) {
        if ( ! isprint ( * It ) && ! isspace ( * It ) ) {
            return false;
        }
    }

    return true;
}

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_Io 
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_Io :: TL_Io ()
{
    throw TL_Exception ( "You are not allowed to create instance of that class" );
}   /* TL_Io :: TL_Io () */

static
int
_getPathType ( const string & Path )
{
    rc_t RCt;
    struct KDirectory * Dir;
    int RetVal;

    RCt = 0;
    Dir = NULL;
    RetVal = kptNotFound;

    RCt = KDirectoryNativeDir ( & Dir );
    if (RCt == 0 ) {
        RetVal = KDirectoryPathType ( Dir, Path . c_str () );

        KDirectoryRelease ( Dir );
    }

    return RetVal; 
}   /* _getPathType () */

bool
TL_Io :: Exists ( const string & Path )
{
    return _getPathType ( Path ) != kptNotFound;
}   /* TL_Io :: Exists () */

bool
TL_Io :: IsDir ( const string & Path )
{
    return _getPathType ( Path ) == kptDir;
}   /* TL_Io :: Exists () */

bool
TL_Io :: IsFile ( const string & Path )
{
    return _getPathType ( Path ) == kptFile;
}   /* TL_Io :: Exists () */

uint64_t
TL_Io :: FileSize ( const std::string & Path )
{
    struct KDirectory * Dir = NULL;
    uint64_t Size = 0;

    if ( KDirectoryNativeDir ( & Dir ) == 0 ) {
        if ( KDirectoryPathType ( Dir, Path . c_str () ) == kptFile ) {
            if ( KDirectoryFileSize ( Dir, & Size, Path . c_str () ) != 0 ) {
                Size = 0;
            }
        }

        KDirectoryRelease ( Dir );
    }

    return Size;
}   /* TL_Io :: FileSize () */

void
TL_Io :: Remove ( const string & Path )
{
    rc_t RCt;
    struct KDirectory * Dir;

    RCt = 0;
    Dir = NULL;

    RCt = KDirectoryNativeDir ( & Dir );
    if ( RCt == 0 ) {
        RCt = KDirectoryRemove ( Dir, true, Path . c_str () );
        KDirectoryRelease ( Dir );
    }

    if ( RCt != 0 ) {
        throw TL_Exception ( string ( "Can not remove directory \"" ) + Path + "\"", RCt );
    }
}   /* TL_Io :: Exists () */


void
TL_Io :: MkDir ( const string & Path )
{
    rc_t RCt;
    struct KDirectory * Dir;

    RCt = 0;
    Dir = NULL;

    RCt = KDirectoryNativeDir ( & Dir );
    if ( RCt == 0 ) {
        RCt = KDirectoryCreateDir ( Dir, 0775, kcmCreate, Path . c_str () );
        KDirectoryRelease ( Dir );
    }

    if ( RCt != 0 ) {
        throw TL_Exception ( string ( "Can not create directory \"" ) + Path + "\"", RCt );
    }
}   /* TL_Io :: Exists () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_GZipReader 
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

/*
    #define TL_DEF_INC_SIZE 65536
*/
#define TL_DEF_INC_SIZE 16384

TL_GZipReader :: TL_GZipReader ()
:   _buffer ( NULL )
,   _buffer_size ( 0 )
{
}   /* TL_GZipReader :: GZipReader () */

TL_GZipReader :: ~TL_GZipReader ()
{
    if ( _buffer != NULL ) {
        delete [] _buffer;
        _buffer = NULL;
    }
    _buffer_size = 0;
}   /* TL_GZipReader :: ~TL_GZipReader () */

/*  Problem is: Kfs GZip file does not support KFileSize method.
 *  That's why.
 *
 */
bool
TL_GZipReader :: Uncompress (
                        const struct KFile * File,
                        void ** Buffer,
                        size_t * BufferSize

)
{
    if ( File == NULL || Buffer == NULL || BufferSize == NULL ) {
        throw TL_Exception ( "TL_GZipReader :: ReadAll : Invalid parameters" );
    }

    * Buffer = NULL;
    * BufferSize = 0;

        /*  First we are trying to get actual file size
         */
    uint64_t Size = 0;
    rc_t RCt = KFileSize ( File, & Size );
    if ( RCt != 0 ) {
            /*  Sometimes size of file is unknown :LOL:
             */
        Size = TL_DEF_INC_SIZE;
    }
    else {
            /*  Some archivers could compress three times :LOL:
             */
        Size *= 3;
    }

    const struct KFile * CompressedFile = NULL;

    RCt = KFileMakeGzipForRead ( & CompressedFile, File );
    if ( RCt == 0 ) {
        _realloc ( Size );
        uint64_t Offset = 0;

        size_t Readed = 0;
        while ( true ) {
            _realloc ( Offset + TL_DEF_INC_SIZE );

            RCt = KFileRead (
                            CompressedFile,
                            Offset,
                            _buffer + Offset,
                            TL_DEF_INC_SIZE,
                            & Readed
                            );
            if ( RCt != 0 ) {
                KFileRelease ( CompressedFile );
                return false;
            }

            if ( Readed == 0 ) {
                char * RetBuf = new char [ Offset ];

                memmove ( RetBuf, _buffer, sizeof ( char ) * Offset );

                * Buffer = ( void * ) RetBuf;
                * BufferSize = Offset;
                break;
            }

            Offset += Readed;
        }

        KFileRelease ( CompressedFile );
    }
    else {
        return false;
    }

    return true;
}   /* TL_GZipReader :: Uncompress () */

void
TL_GZipReader :: _realloc ( size_t Size )
{
    const int _cIncSize = TL_DEF_INC_SIZE;

    if (  _buffer_size <= Size ) {
        size_t NBS = ( ( Size / _cIncSize ) + 1 ) * _cIncSize;

        char * NBF = new char [ NBS ];

        if ( 0 < _buffer_size ) {
            memmove ( NBF, _buffer, sizeof ( char ) * _buffer_size );

            delete [] _buffer;
            _buffer = NULL;
            _buffer_size = 0;
        }

        _buffer = NBF;
        _buffer_size = NBS;
    }
}   /* TL_GZipReader :: _realloc () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_McArch 
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

#define DEF_UND_ENDIAN    0
#define DEF_BIG_ENDIAN    4321
#define DEF_LITTLE_ENDIAN 1234
#define DEF_PDP_ENDIAN    3412

static int _sEndian = DEF_UND_ENDIAN;

void
_calcEndian ()
{
    union { unsigned char c[4]; uint32_t v; } ue;

    ue.c[3]=1; ue.c[2]=2; ue.c[1]=3; ue.c[0]=4;
    switch ( ue.v ) {
        case 0x04030201L: _sEndian = DEF_BIG_ENDIAN; break;
        case 0x01020304L: _sEndian = DEF_LITTLE_ENDIAN; break;
        case 0x03040102L: _sEndian = DEF_PDP_ENDIAN; break;
    }
}   /* _calcEndian () */

TL_McArch :: TL_McArch ()
{
    throw TL_Exception ( "Should not call that" );
}   /* TL_McArch :: TL_McArch () */

bool
TL_McArch :: IsBig ()
{
    if ( _sEndian == DEF_UND_ENDIAN ) {
        _calcEndian ();
    }
    return _sEndian == DEF_BIG_ENDIAN;
}   /* TL_McArch :: IsBig () */

bool
TL_McArch :: IsLittle ()
{
    if ( _sEndian == DEF_UND_ENDIAN ) {
        _calcEndian ();
    }
    return _sEndian == DEF_LITTLE_ENDIAN;
}   /* TL_McArch :: IsLittle () */

bool
TL_McArch :: IsPdp ()
{
    if ( _sEndian == DEF_UND_ENDIAN ) {
        _calcEndian ();
    }
    return _sEndian == DEF_PDP_ENDIAN;
}   /* TL_McArch :: IsPdp () */

uint32_t
TL_McArch :: ToNet ( uint32_t HostLong )
{
    return htonl ( HostLong );
}   /* TL_McArch :: ToNet () */

int32_t
TL_McArch :: ToNet ( int32_t HostLong )
{
    return htonl ( HostLong );
}   /* TL_McArch :: ToNet () */

uint16_t
TL_McArch :: ToNet ( uint16_t HostShort )
{
    return htons ( HostShort );
}   /* TL_McArch :: ToNet () */

int16_t
TL_McArch :: ToNet ( int16_t HostShort )
{
    return htons ( HostShort );
}   /* TL_McArch :: ToNet () */

uint32_t
TL_McArch :: FromNet ( uint32_t NetLong )
{
    return ntohl ( NetLong );
}   /* TL_McArch :: FromNet () */

int32_t
TL_McArch :: FromNet ( int32_t NetLong )
{
    return ntohl ( NetLong );
}   /* TL_McArch :: FromNet () */

uint16_t
TL_McArch :: FromNet ( uint16_t NetShort )
{
    return ntohs ( NetShort );
}   /* TL_McArch :: FromNet () */

int16_t
TL_McArch :: FromNet ( int16_t NetShort )
{
    return ntohs ( NetShort );
}   /* TL_McArch :: FromNet () */

void
TL_McArch :: FromNet ( uint32_a_t & NetLongArray )
{
    if ( ! IsBig () ) {
        for ( uint32_t i = 0; i < NetLongArray . size (); i ++ ) {
            NetLongArray [ i ] = FromNet ( NetLongArray [ i ] );
        }
    }
}   /* TL_McArch :: FromNet () */

void
TL_McArch :: FromNet ( uint16_a_t & NetLongArray )
{
    if ( ! IsBig () ) {
        for ( uint32_t i = 0; i < NetLongArray . size (); i ++ ) {
            NetLongArray [ i ] = FromNet ( NetLongArray [ i ] );
        }
    }
}   /* TL_McArch :: FromNet () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_a2u[s,i] somehow like that
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

uint32_t
TL_a2ui ( unsigned char * B )
{
    union { uint32_t v; unsigned char c [ 4 ]; } d;

    if ( TL_McArch :: IsBig () ) {
        d . c [ 0 ] = B [ 0 ];
        d . c [ 1 ] = B [ 1 ];
        d . c [ 2 ] = B [ 2 ];
        d . c [ 3 ] = B [ 3 ];
    }
    else {
        d . c [ 3 ] = B [ 0 ];
        d . c [ 2 ] = B [ 1 ];
        d . c [ 1 ] = B [ 2 ];
        d . c [ 0 ] = B [ 3 ];
    }

    return d . v;
}   /* _a2ui () */

uint16_t
TL_a2us ( unsigned char * B )
{
    union { uint16_t v; unsigned char c [2]; } d;

    if ( TL_McArch :: IsBig () ) {
        d . c [0] = B [0];
        d . c [1] = B [1];

    }
    else {
        d . c [1] = B [0];
        d . c [0] = B [1];
    }

    return d . v;
}

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_OwpVector
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

static const TL_Owp _sTL_OwpEmptyOwp ( "" );

TL_OwpVector :: TL_OwpVector ()
:   _path ( "" )
,   _vec ()
{
}   /* TL_OwpVector :: TL_OwpVector () */

TL_OwpVector :: ~TL_OwpVector ()
{
    Dispose ();
}   /* TL_OwpVector :: ~TL_OwpVector() */

static
bool
_nextOwp ( TL_LineReader & Reader, TL_OVec & Vec )
{
    bool RetVal = false;
    TL_Owp Owp;

    while ( Reader . NextLine () ) {
        string Line = TL_StringU :: Prune ( Reader . Line () );

            /*  Commentary
             */
        if ( Line [ 0 ] == '#' ) {
            continue;
        }

        if ( Line . empty () ) {
            continue;
        }

        if ( Line == TL_OWP_DELIM ) {
            RetVal = true;
            break;
        }

        string fKey, fValue;

        if ( ! TL_StringU :: Split ( Line, '=', fKey, fValue ) ) {
            throw TL_Exception ( string ( "invalid line format '" ) + Line + "'" );
        }

        if ( fKey == TL_NAME_TAG ) {
            Owp . SetName ( fValue );
        }
        else {
            Owp . SetValue ( fKey, fValue );
        }
    }

    if ( ! Owp . IsEmpty () ) {
        Vec . insert ( Vec . end (), Owp );
        RetVal = true;
    }

    return RetVal;
}   /* _nextOwp () */

void
TL_OwpVector :: Load ( const string & Path, TL_OVec & Vec )
{
    TL_LineReader Reader ( Path );

    while ( _nextOwp ( Reader, Vec ) );

    Reader . Close ();
}   /* TL_OwpVector :: Load () */

void
TL_OwpVector :: Load ( const string & Path )
{
    Dispose ();

    _path = Path;

    Load ( _path, _vec );
}   /* TL_OwpVector :: Load () */

void
TL_OwpVector :: Store ( const string & Path, const TL_OVec & Vec )
{
    ofstream Out ( Path . c_str () );
    if ( ! Out ) {
        throw TL_Exception ( string ( "Can not open file for write \"" ) + Path + "\"" );
    }

    for (
        TL_OVec :: const_iterator It = Vec . begin ();
        It != Vec . end ();
        It ++
    ) {
        Out << It -> ToString () << endl << endl;
        if ( ! Out ) {
            throw TL_Exception ( string ( "Can not write file \"" ) + Path + "\"" );
        }
    }

    Out . close ();
    if ( ! Out ) {
        throw TL_Exception ( string ( "Can not close file \"" ) + Path + "\"" );
    }
}   /* TL_OwpVector :: Store () */


void
TL_OwpVector :: Store ( const string & Path )
{
    string Fsol = Path . empty () ? _path : Path;
    if ( Fsol . empty () ) {
        throw TL_Exception ( "TL_OwpVector :: Store() : empty path passed" );
    }

    Store ( Fsol, _vec );
}   /* TL_OwpVector :: Store () */

void
TL_OwpVector :: Dispose ()
{
    _vec . clear ();
    _path . clear ();
}   /* TL_OwpVector :: Dispose () */

const TL_Owp &
TL_OwpVector :: Get ( size_t Idx ) const
{
    if ( 0 <= Idx && Idx < _vec . size () ) {
        return _vec [ Idx ];
    }
    return _sTL_OwpEmptyOwp;
}   /* TL_OwpVector :: Get () */

const TL_Owp &
TL_OwpVector :: operator [] ( size_t Idx ) const
{
    return Get ( Idx );
}   /* TL_OwpVector :: operator [] () */

void
TL_OwpVector :: Add ( const TL_Owp & Owp )
{
    _vec . insert ( _vec . end (), Owp );
}   /* TL_OwpVector :: Add () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_OwpVector
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_OwpLR :: TL_OwpLR ()
:   _reader ()
,   _owp ( _sTL_OwpEmptyOwp )
,   _owp_no ( ( size_t ) ~0 )
{
}   /* TL_OwpLR :: TL_OwpLR () */

TL_OwpLR :: ~TL_OwpLR ()
{
}   /* TL_OwpLR :: ~TL_OwpLR () */

void
TL_OwpLR :: Open ( const std::string & Path )
{
    Close ();

    _owp_no =  ( size_t ) ~0;
    _owp = _sTL_OwpEmptyOwp;

    _reader . Open ( Path );
}   /* TL_OwpLR :: Open () */

void
TL_OwpLR :: Close ()
{
    _reader . Close ();

    _owp = _sTL_OwpEmptyOwp;
    _owp_no = ( size_t ) ~0;
}   /* TL_OwpLR :: Close () */

void
TL_OwpLR :: Reset ()
{
    _owp_no = ( size_t ) ~0;
    _owp = _sTL_OwpEmptyOwp;

    _reader . Reset ();
}   /* TL_OwpLR :: Reset () */

bool
TL_OwpLR :: NextOwp ()
{
        /* First we are checking if reader was never called
         * in that case we are skipping to delimiter "<<--OWP-->>"
         */
    if ( _owp_no == ( size_t ) ~0 ) {
        bool HasDelimiter = false;
        while ( _reader . NextLine () ) {
            string Line = TL_StringU :: Prune ( _reader . Line () );
            if ( Line == TL_OWP_DELIM ) {
                HasDelimiter = true;
                break;
            }
        }
            /*  If we did not find delimiter that means file is 
             *  formatted badly
             */
        if ( ! HasDelimiter ) {
            throw TL_Exception ( "invalid OWP file format : missed delimiter" );
        }
    }

        /*  Here we are packing an OWP, reading till the next delimiter
         */
    bool RetVal = false;
    TL_Owp Owp;
    while ( _reader . NextLine () ) {
        string Line = TL_StringU :: Prune ( _reader . Line () );

            /*  Commentary
             */
        if ( Line [ 0 ] == '#' ) {
            continue;
        }

        if ( Line . empty () ) {
            continue;
        }

        if ( Line == TL_OWP_DELIM ) {
            break;
        }

        string fKey, fValue;

        if ( ! TL_StringU :: Split ( Line, '=', fKey, fValue ) ) {
            throw TL_Exception ( string ( "invalid line format '" ) + Line + "'" );
        }

        if ( fKey == TL_NAME_TAG ) {
            Owp . SetName ( fValue );
        }
        else {
            Owp . SetValue ( fKey, fValue );
        }
    }

    if ( ! Owp . IsEmpty () ) {
        _owp = Owp;
        _owp_no ++;

        RetVal = true;
    }

    return RetVal;
}   /* TL_OwpLR :: NextOwp () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *
 * TL_MemB
 *
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_MemB :: TL_MemB ()
:   _size ( 0 )
,   _capacity ( 0 )
,   _data ( NULL )
{
}   /* TL_MemB :: TL_MemB () */

TL_MemB :: ~TL_MemB ()
{
    if ( _data != NULL ) {
        delete [] _data;
    }
    _data = NULL;

    _size = 0;
    _capacity = 0;
}   /* TL_MemB :: ~TL_MemB () */

void
TL_MemB :: Reset ()
{
    _size = 0;
}   /* TL_MemB :: Reset () */

void
TL_MemB :: Add ( void * Data, size_t Size )
{
    if ( Size == 0 || Data == NULL ) {
        return;
    }

        /* _reallocing */
    _realloc ( _size + Size );

    memmove ( _data + _size, Data, Size );
    _size += Size;
}   /* TL_MemB :: Add () */

void
TL_MemB :: _realloc ( size_t Size )
{
    if ( _capacity < Size ) {
        const size_t JoJoBa = 2048;
        size_t NewCap = ( ( Size / JoJoBa ) + 1 ) * JoJoBa;
        NewCap *= 2;

        char * NewDat = new char [ NewCap ];
        memset ( NewDat, 0, NewCap );

        if ( _data != NULL ) {
            if ( _size != 0 ) {
                memmove ( NewDat, _data, _size );
            }
 
            delete [] _data;
        }

        _data = NewDat;
        _capacity = NewCap;
    }
}   /* TL_MemB :: _realloc () */

void *
TL_MemB :: Reserve ( size_t Size )
{
    if ( Size != 0 ) {
        _realloc ( Size );

        _size = Size;
        memset ( _data, 0, sizeof ( char ) * Size );

    }

    return _data;
}   /* TL_MemB :: Reserve () */

