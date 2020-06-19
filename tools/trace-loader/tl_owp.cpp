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

#include "tl_owp.hpp"
#include "tl_util.hpp"

using namespace std;
using namespace _tl_;

static const string _sTL_OwpEmptyString ( "" );

TL_Owp :: TL_Owp ( const string & Name )
:   _name ( Name )
,   _map ()
{
}   /* TL_Owp :: TL_Owp () */

TL_Owp :: TL_Owp ( const TL_Owp & Owp )
:   _name ( Owp . _name )
,   _map ( Owp . _map )
{
}   /* TL_Owp :: TL_Owp () */

TL_Owp :: ~TL_Owp ()
{
    Clear (); 
}   /* TL_Owp :: ~TL_Owp () */

TL_Owp &
TL_Owp :: operator = ( const string & Name )
{
    Clear ();

    _name = Name;

    return * this;
}   /* TL_Owp :: operator = () */

TL_Owp &
TL_Owp :: operator = ( const TL_Owp & Owp )
{
    if ( this != & Owp ) {

        Clear ();

        _name = Owp . _name;
        _map = Owp . _map;
    }
    return * this;
}   /* TL_Owp :: operator = () */

void
TL_Owp :: Clear ()
{
    _name . clear ();
    _map . clear ();
}   /* TL_Owp :: Clear () */

bool
TL_Owp :: operator == ( const TL_Owp & Owp ) const
{
    if ( Name () != Owp . Name () ) {
        return false;
    }

    if ( _map . size () != Owp . _map . size () ) {
        return false;
    }

    for ( 
        TL_SSMap :: const_iterator It = _map . begin ();
        It != _map . end ();
        It ++
    ) {
        if ( Value ( It -> first ) != Owp . Value ( It -> first ) ) {
            return false;
        }
    }

    return true;
}   /* TL_Owp :: operator == () */

bool
TL_Owp :: Has ( const string & Key ) const
{
    return _map . find ( Key ) != _map . end ();
}   /* TL_Owp :: Has () */

const string &
TL_Owp :: Value ( const string & Key ) const
{
    TL_SSMap :: const_iterator It = _map . find ( Key );
    if ( It != _map . end () ) {
        return It -> second;
    }

    return _sTL_OwpEmptyString;
}   /* TL_Owp :: Value () */

void
TL_Owp :: SetValue ( const string & Key, const string & Value )
{
    _map [ Key ] = Value;
}   /* TL_Owp :: SetValue () */

void
TL_Owp :: EraseValue ( const string & Key )
{
    _map . erase ( Key );
}   /* TL_Owp :: EraseValue () */

bool
TL_Owp :: ListKeys ( TL_SVec & Keys ) const
{
    Keys . clear ();

    for (
        TL_SSMap :: const_iterator It = _map . begin ();
        It != _map . end ();
        It ++
    ) {
        Keys . insert ( Keys . end (), It -> first );
    }

    return Keys . size () != 0;
}   /* TL_Owp :: ListKeys () */

TL_Owp &
TL_Owp :: operator -= ( const TL_Owp & Owp )
{
    if ( this == & Owp ) {
        Clear ();
    }
    else {
        TL_SVec Keys;
        if ( ListKeys ( Keys ) ) {
            for (
                TL_SVec :: const_iterator It = Keys . begin ();
                It != Keys . end ();
                It ++
            ) {
                if ( ! Owp . Has ( * It ) ) {
                    _map . erase ( * It );
                }
            }
        }
    }

    return * this;
}   /* TL_Owp :: operator - () */

TL_Owp &
TL_Owp :: operator += ( const TL_Owp & Owp )
{
    if ( this != & Owp ) {
        if ( ! Owp . _name . empty () ) {
            _name = Owp . _name;
        }

        for ( 
            TL_SSMap :: const_iterator It = Owp . _map . begin ();
            It != Owp . _map . end ();
            It ++
        ) {
            _map [ It -> first ] = It -> second;
        }
    }

    return * this;
}   /* TL_Owp :: operator + () */

string
TL_Owp :: ToString () const
{
    stringstream Str;

    Str << TL_OWP_DELIM << endl;
    Str << TL_NAME_TAG << " = " << _name << endl;

    for (
        TL_SSMap :: const_iterator It = _map . begin ();
        It != _map . end ();
        It ++
    ) {
        Str << It -> first << " = " << It -> second << endl;
    }

    Str << endl;

    return Str . str ();
}   /* TL_Owp :: ToString () */

