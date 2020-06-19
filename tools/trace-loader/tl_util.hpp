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

#ifndef _tl_util_
#define _tl_util_

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "tl_exception.hpp"
#include "tl_types.hpp"
#include "tl_owp.hpp"

#include <kfs/file.h>

namespace _tl_ {

/*))    Simple wrapper for KMMap ... just basic methods.
  ((*/
class TL_KfsMap {
public:
    TL_KfsMap ();
    TL_KfsMap ( const TL_KfsMap & Map );
    TL_KfsMap ( const std::string & Path );
    ~TL_KfsMap ();

    TL_KfsMap & operator = ( const TL_KfsMap & Map );
    TL_KfsMap & operator = ( const std::string & Path );

    void Map ( const std::string & Path );
    void UnMap ();

    inline bool IsMapped () const { return _handle != NULL; };

    inline const std::string & Path () const { return _path; };
    inline size_t Size () const { return _size; };

    inline const void * Get () const { return _addr; };

private:
    void _DoMap ( const TL_KfsMap & Map );
    void _DoMap ( const std::string & Path );
    void _DoUnMap ();

private:
    void * _handle;

    std::string _path;
    size_t _size;
    const void * _addr;
};  /* TL_KfsMap */

/*))    Simple line reader ... just basic
  ||    Just open file as map and "while ( NextLine () ) GetLine ();"
  ||    Warning: it cut off '\n'
  ((*/
class TL_LineReader {
public:
    TL_LineReader ( const std::string & Path = "" );
    ~TL_LineReader ();

    void Open ( const std::string & Path );
    void Close ();
    inline bool IsOpen () const { return _map . IsMapped (); };

    void Reset ();
    bool NextLine ();
    inline const std::string & Line () const { return _line; };
    inline size_t LineNo () const { return _line_no; };

private:
    void _DoOpen ( const std::string & Path );
    void _DoClose ();

private:
    TL_KfsMap _map;

    std::string _line;
    size_t _line_no;

    size_t _pos;
};

    /*)) Awkward class: there are two typese of time ( really four )
     //  described in RFC: run_date ( 2000-10-28) and
    //   collection_date ( Mar 2 2006 12:00AM ). However, I see it 
     //  in many different formats, like run_date ( 05-Jan-2013 )
    //   and collection_date ( Thu Feb  6 16:33:00 2014 )
     //
    //   Ok, I run some script on TRACEINFO.tbl files for last couple
     //  years, and there is statistic:
    //
     //  run_dates are all in format :
    //       Fri Apr XX XX:XX:XX XXXX
     //      XX-APR-XXXX
    //   collect_dates are in format :
     //      X-Aug-XXXX
    //       XXXX
     //
    //   Which is slightry different from RFC, so, I believe, I need
     //  to parse 6 different formats : 2 from RFC and 4 from real life
    //   Also, we may nod distinct run from collect date by format :/
   ((*/
class TL_DateU {
public:
        /*  These two returns RFC format or empty string in error
         */
    static std::string RunDate ( uint64_t TimeT );
    static std::string CollectionDate ( uint64_t Time );

        /*  This returns time_t or ZERO in error
         */
    static uint64_t Date ( const std::string & Time );

        /*  This converts local time_t to GMT time_t
         */
    static uint64_t LocalToGMT ( uint64_t LocalTime );
};

class TL_StringU {
public:
        /*  Removes 'space' characters from begin and end of line
         */
    static std::string Prune (
                        const std::string & Line,
                        bool Front = true,
                        bool End = true
                        );

        /*  splits line into two parts by given delimiter
         */
    static bool Split (
                        const std::string & Line,
                        char Delim,
                        std::string & Key,
                        std::string & Value
                        );

        /*  splits line of text by given delimiter character
            WARNING: it will clear vector from previous tokens
            RETURNS: num elements in array
         */
    static size_t Tokenize (
                        const std::string & Line,
                        char Delim,
                        TL_SVec & Tokens,
                        bool AddEmptyTokens = false
                        );

        /*  Kinda not cool, will reimplement later
         */
    static std::string ToUpper ( const std::string & Line );
    static std::string ToLower ( const std::string & Line );
    static std::string ToPrintable ( const std::string & Line );

        /*  Kinda not cool, will reimplement later
         */
    static std::string FromBin ( const void * Src, size_t SrcSize );
    static void ToBin (
                        const std::string & Src,
                        void * Dst,
                        size_t DstSize
                        );

        /*  Kinda stupid, but sometimes needed
         */
    static bool StartsWith (
                        const std::string & Line,
                        const std::string & With
                        );
    static bool EndsWith (
                        const std::string & Line,
                        const std::string & With
                        );

    template < class Type > static std::string ToStr (
                                                    Type Data,
                                                    int Prec = 0
                                                    );
    template < class Type > static Type FromStr ( std::string Data );

        /*  Kinda cool thing to have it up and ready
         */
    static const std::string & EmptyString ();

        /*  In some cases it could be printable :LOL:
         *  There will be dups from 'ctype.h' effu
         */
    static bool IsPrint ( const std::string & Line );

private:
    TL_StringU ();
};

class TL_Io {
public:
    static bool Exists ( const std::string & Path );
    static bool IsDir ( const std::string & Path );
    static bool IsFile ( const std::string & Path );

        /*  For everything what is not a file that returns 0
         */
    static uint64_t FileSize ( const std::string & Path );

        /*  Always forse removing
         */
    static void Remove ( const std::string & Path );

    static void MkDir ( const std::string & Path );

private:
    TL_Io ();
};

/*  For (u)int8_t use cast to (u)int32_t, or it will be formatted as 
 *  a char :LOL:
 */
template < class Type >
std::string
TL_StringU :: ToStr ( Type Data, int Prec )
{
    std::ostringstream Out;

    if ( Prec != 0 ) {
        Out << std :: setprecision ( Prec );
        if ( ! Out . good () ) {
            throw TL_Exception ( "Can not convert data to string" );
        }
    }

    Out << ( Type ) Data;
    if ( ! Out . good () ) {
        throw TL_Exception ( "Can not convert data to string" );
    }

    return Out . str ();
}   /* TL_StringU :: ToStr () */

template < class Type >
Type
TL_StringU :: FromStr ( std::string Data )
{
    Type Var;

    std::istringstream Out ( Data );

    Out >> Var;
        /* eof() and good() are diamedral */
    if ( ! Out . eof () && ! Out . good () ) {
        throw TL_Exception ( std::string ( "Can not convert string to data [" ) + Data + "]" );
    }

    return Var;
}   /* TL_StringU :: FromStr () */

    /*  That class reads zip file, uncompress it and stores in a memory
     *  NOTE: Small files, please! 
     */
class TL_GZipReader {
public :
    TL_GZipReader ();
    ~TL_GZipReader ();

        /*  Will return false if file is not compressed
         */
    bool Uncompress (
                    const struct KFile * CompressedFile,
                    void ** Buffer,
                    size_t * BufferSize
                    );

private :
    void _realloc ( size_t Size );

    char * _buffer;
    size_t _buffer_size;
};

}   /* namespace _tl_ */

/*))    Some kind of integer gettin' from ASCII :LOL:
 ((     Ha-ha, these functions do not make check for length of buffer,
  ))    so, plese call a2ui with 4 char lenght, and a2us with 2 char.
 ((*/
uint32_t TL_a2ui ( unsigned char * Buf );
uint16_t TL_a2us ( unsigned char * Buf );

namespace _tl_ {

/*  Just a template :D
 *  NOTE: it makes array copy, it is Your task to delete initial data
 */
template < class Type >
class tl_array {
public:
    tl_array ();
    tl_array ( size_t Size, Type * Arr );
    tl_array ( size_t Size, int Fill );
    tl_array ( const tl_array & Arr );
    ~tl_array ();

    void reset ();
    void set ( size_t Size, Type * Arr );
    void set ( size_t Size, int Fill );
    void set ( const tl_array & Arr );

    tl_array & operator = ( const tl_array & Arr );

    inline bool empty () const { return _size == 0; };
    inline bool operator ! () const { return empty (); };

    inline size_t size () const { return _size; };
    inline const Type * array () const { return _array; };

    inline const Type * operator -> () const { return _array; };
    inline Type * operator -> () { return _array; };
    inline const Type * operator * () const { return _array; };
    inline Type * operator * () { return _array; };

        /*  Do I need to check Index here ?
         */
    inline const Type & operator [] ( int Index ) const
                                    { return _array [ Index ]; };
    inline Type & operator [] ( int Index )
                                    { return _array [ Index ]; };

private:
    void _reset ();
    void _set ( size_t Size, Type * Arr );
    void _set ( size_t Size, int Fill );
    void _set ( const tl_array & Arr );

    size_t _size;
    Type * _array;
};

template < class Type >
tl_array < Type > :: tl_array ()
:   _size ( 0 )
,   _array ( NULL )
{
}   /* tl_array < Type > :: tl_array () */

template < class Type >
tl_array < Type > :: tl_array ( size_t Size, Type * Arr )
:   _size ( 0 )
,   _array ( NULL )
{
    _set ( Size, Arr );
}   /* tl_array < Type > :: tl_array () */

template < class Type >
tl_array < Type > :: tl_array ( size_t Size, int Fill )
:   _size ( 0 )
,   _array ( NULL )
{
    _set ( Size, Fill );
}   /* tl_array < Type > :: tl_array () */

template < class Type >
tl_array < Type > :: tl_array ( const tl_array < Type > & Arr )
:   _size ( 0 )
,   _array ( NULL )
{
    _set ( Arr );
}   /* tl_array < Type > :: tl_array () */

template < class Type >
tl_array < Type > :: ~tl_array ()
{
    _reset ();
}   /* tl_array < Type > :: ~tl_array () */

template < class Type >
void
tl_array < Type > :: reset ()
{
    _reset ();
}   /* tl_array < Type > :: set () */

template < class Type >
void
tl_array < Type > :: set ( size_t Size, Type * Arr )
{
    _set ( Size, Arr );
}   /* tl_array < Type > :: set () */

template < class Type >
void
tl_array < Type > :: set ( size_t Size, int Fill )
{
    _set ( Size, Fill );
}   /* tl_array < Type > :: set () */

template < class Type >
void
tl_array < Type > :: set ( const tl_array < Type > & Arr )
{
    _set ( Arr );
}   /* tl_array < Type > :: set () */

template < class Type >
tl_array < Type > &
tl_array < Type > :: operator = ( const tl_array < Type > & Arr )
{
    if ( this != & Arr ) {
        _set ( Arr );
    }

    return * this;
}   /* tl_array < Type > :: operator = () */

template < class Type >
void
tl_array < Type > :: _reset ()
{
    if ( _array != NULL ) {
        delete [] _array;
        _array = NULL;
    }

    _size = 0;
}   /* tl_array < Type > :: _reset () */

template < class Type >
void
tl_array < Type > :: _set ( size_t Size, Type * Arr )
{
    _reset ();

    if ( Size != 0 ) {
        _array = new Type [ Size ];
        memmove ( _array, Arr, sizeof ( Type ) * Size );
        _size = Size;
    }
}   /* tl_array < Type > :: _set () */

template < class Type >
void
tl_array < Type > :: _set ( size_t Size, int Fill )
{
    _reset ();

    if ( Size != 0 ) {
        _array = new Type [ Size ];
        memset ( _array, Fill, sizeof ( Type ) * Size );
        _size = Size;
    }
}   /* tl_array < Type > :: _set () */

template < class Type >
void
tl_array < Type > :: _set ( const tl_array < Type > & Arr )
{
    _set ( Arr . _size, Arr . _array );
}   /* tl_array < Type > :: _set () */

/*  Something used often
 */
typedef tl_array < unsigned char > uchar_a_t;
typedef tl_array < char > char_a_t;
typedef tl_array < uint16_t > uint16_a_t;
typedef tl_array < int16_t > int16_a_t;
typedef tl_array < uint32_t > uint32_a_t;
typedef tl_array < int32_t > int32_a_t;
typedef tl_array < uint64_t > uint64_a_t;
typedef tl_array < int64_t > int64_a_t;

}   /* namespace _tl_ */

    /*  Sorry the order of that stuff is quite stochastic
     */
namespace _tl_ {

/*))    McArch .. machine architect stuff
  ((*/
class TL_McArch {
public :
    static bool IsBig ();
    static bool IsLittle ();
    static bool IsPdp ();

        /*  I think that it is possible to make template,
         *  but ... do not like that idea :D
         */
    static uint32_t ToNet ( uint32_t HostLong );
    static int32_t ToNet ( int32_t HostLong );
    static uint16_t ToNet ( uint16_t HostShort );
    static int16_t ToNet ( int16_t HostShort );
    static uint32_t FromNet ( uint32_t NetLong );
    static int32_t FromNet ( int32_t NetLong );
    static uint16_t FromNet ( uint16_t NetShort );
    static int16_t FromNet ( int16_t NetShort );

        /*  Warning!!! that wil modify input array
         */
    static void FromNet ( uint32_a_t & NetLong );
    static void FromNet ( uint16_a_t & NetShort );

private :
    TL_McArch ();
};


    /*  Simple archive of objects with properties
     */
class TL_OwpVector {
public:
    static void Load ( const std::string & Path, TL_OVec & Vec );
    static void Store ( const std::string & Path, const TL_OVec & Vec );

public:
    TL_OwpVector ();
    ~TL_OwpVector ();

    void Load ( const std::string & Path );
    void Store ( const std::string & Path = "" );
    void Dispose ();

    size_t Size () const { return _vec . size (); };
    const TL_Owp & Get ( size_t Idx ) const;
    const TL_Owp & operator [] ( size_t Idx ) const;

    void Add ( const TL_Owp & Owp );

private:
    std::string _path;

    TL_OVec _vec;
};

/*  Linear Reader for OWP : maps file and reads TL_Owp one by one
 */
class TL_OwpLR {
public :
    TL_OwpLR ();
    ~TL_OwpLR ();

    void Open ( const std::string & Path );
    void Close ();
    inline bool IsOpen () const { return _reader . IsOpen (); };

    void Reset ();
    bool NextOwp ();
    inline const TL_Owp & Owp () const { return _owp; };
    inline size_t OwpNo () const { return _owp_no; };

private :
    TL_LineReader _reader;

    TL_Owp _owp;
    size_t _owp_no;
};


    /* Some simple memory buffer
     */
class TL_MemB {
public:
    TL_MemB ();
    ~TL_MemB ();

    void Reset ();
    void Add ( void * Data, size_t Size );

        /*  that method is to initialize and modify data
         *  if Size == 0, then there will be no any reallocation
         *  but if array is empty - you will get NULL :D
         */
    void * Reserve ( size_t Size = 0 );

    inline const void * Data () const { return _data; };
    inline size_t Size () const { return _size; };

private:
    void _realloc ( size_t Size );

    size_t _size;
    size_t _capacity;
    char * _data;
};

}   /* namespace _tl_ */

#endif /* _tl_util_ */

