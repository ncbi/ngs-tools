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
#include <stack>

#include <math.h>

#include <klib/log.h>

#include "tl_vrules.hpp"
#include "tl_tracefields.hpp"

using namespace std;
using namespace _tl_;

/*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
 *  The VRule and syntax.
 *  Very simple: VRule consists from two parts : condition and action
 *
 *  Condition is a expression written in some language.
 *  Action is predefined set of checks.
 *
 *  Expression could be described formally as ( sorry I am not good
 *  in BNF, and it is very schematic approach ) :
 *  
 *  <rule> ::= <expr> | "ANY"
 *  <expr> ::=   "(" <expr> ")"
 *              | ! < expr >
 *              | <expr> <log_op> <expr>
 *              | <field>
 *              | <field> <oper> <value>
 *  <log_op> ::= <AND> | <OR>
 *  <oper> ::= "<" | "<=" | "==" | ">=" | ">" | "!="
 *  <field> ::= <name of field>
 *  <value> ::= <name of field> | "\"" <al-num> "\""
 *  "ANY" - reserved word for any condition
 *
 *  Action is defined by it's name, and some additional parameters
 *  are added
 *
 *  Rules are stored as vector of TL_Owp objects:
 *     Each rule has unique name
 *     Each rule has section 'condition' with expression
 *     Each rule has section 'action' with name of action
 *     All necessary parameters are stored as properties of TL_Owp
 *     Two rules with same conditions applies by order of appearance
 *     The expression is parsed and executed from left to right
 *
 *  Priorit by descend
 *
 *     !
 *     <= < == > >= != 
 *     AND OR
 *
 *  List of available actions:
 *     CheckFieldExists
 *     CheckFieldValue
 *     CheckFieldNotExists
 *
 *  Will extend it with time
 *
 *=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*
 *  Implementation:
 *  Because rules are very simple and, prolly, short, will do three
 *  four scanning of each line :LOL:
 *=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ Weird stuff
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
#define _VR_NAME              "name"
#define _VR_CONDITION         "condition"
#define _VR_ACTION            "action"
#define _VR_ACTION_FTYPE      "check_field_type"
#define _VR_ACTION_FPRESENCE  "check_field_presence"
#define _VR_YES_FIELDS        "yes_fields"
#define _VR_NO_FIELDS         "no_fields"
#define _VR_YES_VALUE          1
#define _VR_NO_VALUE           0
#define _VR_TYPE              "type"
#define _VR_TYPE_DATE         "date"
#define _VR_TYPE_FLOAT        "float"
#define _VR_TYPE_INT          "int"
#define _VR_TYPE_STRING       "string"
#define _VR_MAX               "max"
#define _VR_MIN               "min"
#define _VR_FIELD             "field"
#define _VR_ALPHABET          "alphabet"
#define _VR_OTHER             "other"

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_  Simple parser ... very dubious, but works.
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_  token - generally we do not have any other type of token
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _VR_t {
public:
    typedef enum {
        eUN = -1,   /* type invalid */
        eANY,       /* type ANY */
        eWR,        /* type word: or literal */
        eSTR,       /* type string: or literal in "" */
        eLB,        /* type left bracket */
        eRB,        /* type right bracket */
        eOP,        /* type operation */
        eLOP,       /* type logical operation */
        eNG         /* type O-negative */
        } eType;
public:
    _VR_t ()
            : _type ( eUN ), _token ( "" ) {};
    _VR_t ( eType Tp, const string & Tkn = "" )
            : _type ( Tp ), _token ( Tkn ) {};
    _VR_t ( const _VR_t & T )
            : _type ( T . _type ), _token ( T . _token ) {};
    ~_VR_t () { _type = eUN; _token . clear (); };

    _VR_t & operator = ( const _VR_t & T )
            { if ( this != & T ) {
                  _type = T . _type; _token = T . _token;
              }
              return * this;
            };

    inline void Invalidate () { _type = eUN; _token . clear (); };

    inline eType Type () const { return _type; };
    inline const string & Token () const { return _token; };

    inline bool gOOD () const { return _type != eUN; };

    string ToString () const;

    static string T2S ( eType T );

private:
    eType _type;
    string _token;
};

string
_VR_t :: T2S ( eType T )
{
    switch ( T ) {
        case eUN: return "eUN";
        case eANY: return "eANY";
        case eWR: return "eWR";
        case eSTR: return "eSTR";
        case eLB: return "eLB";
        case eRB: return "eRB";
        case eOP: return "eOP";
        case eLOP: return "eLOP";
        case eNG: return "eNG";
        default:
            return "Unknown";
    }
}   /* _VR_t :: T2S () */

string
_VR_t :: ToString () const
{
    stringstream Str;

    Str << "[" << _token << "](" << T2S ( _type ) << ")";

    return Str . str ();
}   /* _VR_t :: ToString () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_  token puller pulls token from string and returns it's type
 *_  and may be token.
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _VR_tp {
public:
    _VR_tp () : _line ( "" ), _curr ( NULL ) {};
    ~_VR_tp () { Dispose (); };

    void Init ( const string & Line )
                { _line = Line; _curr = Line . c_str (); };
    void Rewind () { _curr =  _line . c_str (); };
    void Dispose () { _line . clear (); _curr = NULL; };

    bool Pull ( _VR_t & Token );

    void PullAll ( vector < _VR_t > & Tokens );

private:
    string _line;
    const char * _curr;
};

/*  Not really smart, but simple :LOL:
 */
bool
_VR_tp :: Pull ( _VR_t & Token )
{
    Token . Invalidate ();

    const char * end = _line . c_str () + _line . length ();

        /*  skipping trailing spaces
         */
    while ( _curr < end ) {
        if ( ! isspace ( * _curr ) ) {
            break;
        }
        _curr ++;
    }

    if ( end <= _curr ) {
        return false;
    }

        /*  First are brakets, negatives and comparision operations
         */
    if ( * _curr == '(' ) {
        Token = _VR_t ( _VR_t :: eLB, "(" );
        _curr ++;
        return true;
    }

    if ( * _curr == ')' ) {
        Token = _VR_t ( _VR_t :: eRB, ")" );
        _curr ++;
        return true;
    }

    if ( * _curr == '!' ) {
        _curr ++;
        if ( _curr < end ) {
            if ( * _curr == '=' ) {
                Token = _VR_t ( _VR_t :: eOP, "!=" );
                _curr ++;
                return true;
            }
        }

        Token = _VR_t ( _VR_t :: eNG, "!" );
        return true;
    }

        /*  Operations : "<=" "<"
         */
    if ( * _curr == '<' ) {
        _curr ++;
        if ( _curr < end ) {
            if ( * _curr == '=' ) {
                Token = _VR_t ( _VR_t :: eOP, "<=" );
                _curr ++;
                return true;
            }
        }

        Token = _VR_t ( _VR_t :: eOP, "<" );
        return true;
    }

        /*  Operations : ">=" ">"
         */
    if ( * _curr == '>' ) {
        _curr ++;
        if ( _curr < end ) {
            if ( * _curr == '=' ) {
                Token = _VR_t ( _VR_t :: eOP, ">=" );
                _curr ++;
                return true;
            }
        }

        Token = _VR_t ( _VR_t :: eOP, ">" );
        return true;
    }

        /*  Operations : "=="
         */
    if ( * _curr == '=' ) {
        _curr ++;
        if ( _curr < end ) {
            if ( * _curr == '=' ) {
                Token = _VR_t ( _VR_t :: eOP, "==" );
                _curr ++;
                return true;
            }
        }
    }

        /*  String
         */
    if ( * _curr == '"' ) {
        _curr ++;
        const char * word = _curr;

            /* skipping to end of string */
        while ( _curr < end ) {
            if ( * _curr == '"' ) {
                break;
            }
            _curr ++;
        }

        if ( end <= _curr ) {
            throw TL_Exception ( string ( "_VR_tp::Pull(): ERROR: Unbalanced '\"\' character in line [" ) + _line + "]" );
        }

        Token = _VR_t ( _VR_t :: eSTR, string ( word, _curr - word ) );
        _curr ++;

        return true;
    }

        /*  Word, Field Value, OR and AND
         */
        /* Just skipping till first space character */
    const char * word = _curr;
    while ( _curr < end ) {
        char Ch = * _curr;
        if ( isspace ( Ch ) ) {
            break;
        }

            /*  Comparision operations ... */
        if ( Ch == '<' || Ch == '>'
            || Ch == '=' || Ch == '!'
            || Ch == '"'
            || Ch == '(' || Ch == ')' ) {
            break;
        }
        _curr ++;
    }

    size_t Len = _curr - word;
    string Tok ( word, Len );
    if ( Tok == "AND" ) {
        Token = _VR_t ( _VR_t :: eLOP, "AND" );
    }
    else {
        if ( Tok == "OR" ) {
            Token = _VR_t ( _VR_t :: eLOP, "OR" );
        }
        else {
            if ( Tok == "ANY" ) {
                Token = _VR_t ( _VR_t :: eANY, "ANY" );
            }
            else {
                Token = _VR_t ( _VR_t :: eWR, string ( word, Len ) );
            }
        }
    }

    if ( isspace ( * _curr ) ) {
        _curr ++;
    }

    return true;
}   /* _VR_tp :: Pull () */

void
_VR_tp :: PullAll ( vector < _VR_t > & Tokens )
{
    Tokens . clear ();

    Rewind ();

    bool AnySet = false;
    _VR_t T;
    while ( Pull ( T ) ) {
        if ( T . Type () == _VR_t :: eANY ) {
            AnySet = true;
        }

        Tokens . insert ( Tokens . end (), T );
    }

    if ( Tokens . size () != 1 && AnySet ) {
        throw TL_Exception ( "_VR_tp::PullAll(): Condition \"ANY\" could be set only stand alone" );
    }

}   /* _VR_tp :: PullAll () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_  Parsing tree. It will be one pass and recursion
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _VR_ex {
public :
    _VR_ex () {};
    virtual ~_VR_ex () {};

    virtual bool Validate ( const TL_Owp & Row ) const = 0;

    virtual string ToString () const = 0;
};

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
    /*  Just a field, technically field name
     */
class _VR_ex_ANY : public _VR_ex {
public :
    _VR_ex_ANY ( const _VR_t & Token );
    ~_VR_ex_ANY ();

    bool Validate ( const TL_Owp & Row ) const;
    string ToString () const;
};

_VR_ex_ANY :: _VR_ex_ANY ( const _VR_t & Token )
{
    if ( Token . Type () != _VR_t :: eANY ) {
        throw TL_Exception ( string ( "_VR_ex_ANY: Invalid token type [" ) + Token . Token () + "]" );
    }
}   /* _VR_ex_ANY :: _VR_ex_ANY () */

_VR_ex_ANY :: ~_VR_ex_ANY ()
{
}   /* _VR_ex_ANY :: ~_VR_ex_ANY () */

bool
_VR_ex_ANY :: Validate ( const TL_Owp & Row ) const
{
    return true;
}   /* _VR_ex_ANY :: Validate () */

string
_VR_ex_ANY :: ToString () const
{
    return "ANY";
}   /* _VR_ex_ANY :: ToString () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
    /*  Just a field, technically field name
     */
class _VR_ex_WR : public _VR_ex {
public :
    _VR_ex_WR ( const _VR_t & Token );
    ~_VR_ex_WR ();

    bool Validate ( const TL_Owp & Row ) const;
    string ToString () const;
private :
    string _field;
};

_VR_ex_WR :: _VR_ex_WR ( const _VR_t & Token )
:   _field ( "" )
{
        /* I know it is not good to throw here
         */
    if ( Token . Type () != _VR_t :: eWR ) {
        throw TL_Exception ( string ( "_VR_ex_WR: Invalid token type [" ) + Token . Token () + "]" );
    }

        /* Here, before we will create that rule, we should
         * be sure that that field does exists
         */
    _field = Token . Token ();

    if ( ! TL_TraceFields :: Has ( _field ) ) {
        throw TL_Exception ( string ( "_VR_ex_WR: Can not find field with name \"" ) + _field + "\"" );
    }

    PLOGMSG(
            klogDebug,
            (klogDebug,
            "_VR_ex_WR ($(field))",
            "severity=debug,field=%s",
            _field . c_str ()
            )
            );
}   /* _VR_ex_WR :: _VR_ex_WR () */

_VR_ex_WR :: ~_VR_ex_WR ()
{
    _field . clear ();
}   /* _VR_ex_WR :: ~_VR_ex_WR () */

bool
_VR_ex_WR :: Validate ( const TL_Owp & Row ) const
{
    return Row . Has ( _field );
}   /* _VR_ex_WR :: Validate () */

string
_VR_ex_WR :: ToString () const
{
    return _field;
}   /* _VR_ex_WR :: ToString () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
    /*  Simple comparision operation: always <field> <op> <value>
     */
class _VR_ex_OP : public _VR_ex {
private:
    typedef enum {
            eUND = -1,
            eL,     // "<"
            eLE,    // "<="
            eE,     // "=="
            eGE,    // ">=" 
            eG,     // ">"
            eNE     // "!="
            } eOp;
public:
    _VR_ex_OP (
                    const _VR_t & Left,
                    const _VR_t & Op,
                    const _VR_t & Right
                    );
    ~_VR_ex_OP ();

    bool Validate ( const TL_Owp & Row ) const;
    string ToString () const;

private :
    string _left;  /* Always name */
    eOp _op;
    string _r_value;  /* Right operand as value */
    string _r_field;  /* Right operand as field name */
};

_VR_ex_OP :: _VR_ex_OP (
                    const _VR_t & Left,
                    const _VR_t & Op,
                    const _VR_t & Right
)
:   _left ( "" )
,   _op ( eUND )
,   _r_value ( "" )
,   _r_field ( "" )
{
        /* I know it is not good to throw here
         */
    if ( Left . Type () != _VR_t :: eWR ) {
        throw TL_Exception ( string ( "_VR_ex_OP: Invalid left token type [" ) + Left . Token () + "]" );
    }
    _left = Left . Token ();

    if ( ! TL_TraceFields :: Has ( _left ) ) {
        throw TL_Exception ( string ( "_VR_ex_WR: Can not find field with name \"" ) + _left + "\"" );
    }

    _op = eUND;

    if ( Op . Type () != _VR_t :: eOP ) {
        throw TL_Exception ( string ( "_VR_ex_OP: Invalid operation token type [" ) + Op . Token () + "]" );
    }
    if ( Op . Token () == "<" ) { _op = eL; }
    else {
        if ( Op . Token () == "<=" ) { _op = eLE; }
        else {
            if ( Op . Token () == "==" ) { _op = eE; }
            else {
                if ( Op . Token () == ">=" ) { _op = eGE; }
                else {
                    if ( Op . Token () == ">" ) { _op = eG; }
                    else {
                        if ( Op . Token () == "!=" ) { _op = eNE; }
                        else {
                            throw TL_Exception ( string ( "_VR_ex_OP: Invalid operation token type [" ) + Op . Token () + "]" );
                        }
                    }
                }
            }
        }
    }

    if ( Right . Type () == _VR_t :: eWR ) {
        _r_field = Right . Token ();

        if ( ! TL_TraceFields :: Has ( _r_field ) ) {
            throw TL_Exception ( string ( "_VR_ex_WR: Can not find field with name \"" ) + _r_field + "\"" );
        }
    }
    else {
        _r_value = Right . Token ();
    }

    PLOGMSG(
            klogDebug,
            (klogDebug,
            "_VR_ex_OP ($(field) $(token) $(right))",
            "severity=debug,field=%s,token=%s,right=%s",
            Left . Token () . c_str (),
            Op . Token () . c_str (),
            Right . Token () . c_str ()
            )
            );
}   /* _VR_ex_OP :: _VR_ex_OP () */

_VR_ex_OP :: ~_VR_ex_OP ()
{
    _op = eUND;
}   /* _VR_ex_OP :: ~_VR_ex_OP () */

/*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
 * NOTE ... right now it compares only text values - no INTs or whatever
 * IMPORTANT : it is case insensitive comparision
 *#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*/
bool
_VR_ex_OP :: Validate ( const TL_Owp & Row ) const
{
        /*  First we should get value of field on the left
         */
    string LValue = TL_StringU :: ToUpper ( Row . Value ( _left ) );

    string RValue = "";

    if ( ! _r_value . empty () ) {
        RValue = _r_value;
    }
    else {
        if ( ! _r_field . empty () ) {
            RValue = Row . Value ( _r_field );
        }
        else {
            return false;
        }
    }

    RValue = TL_StringU :: ToUpper ( RValue );

    bool RetVal = false;

    switch ( _op ) {
        case eL :   RetVal = LValue < RValue;   break;
        case eLE :  RetVal = LValue <= RValue;  break;
        case eE :   RetVal = LValue == RValue;  break;
        case eGE :  RetVal = LValue >= RValue;  break;
        case eG :   RetVal = LValue > RValue;   break;
        case eNE :  RetVal = LValue != RValue;  break;
        default:    RetVal = false;             break;
    }

    return RetVal;
}   /* _VR_ex_OP :: Validate () */

string
_VR_ex_OP :: ToString () const 
{
    string Ret = "(";
    Ret += _left;

    switch ( _op ) {
        case eL : Ret += " < "; break;
        case eLE : Ret += " <= "; break;
        case eE : Ret += " == "; break;
        case eGE : Ret += " >= "; break;
        case eG : Ret += " > "; break;
        case eNE : Ret += " != "; break;
        default: return "INVALID COMPARISION";
    }

    if ( ! _r_field . empty () ) {
        Ret += _r_field;
    }
    else {
        if ( ! _r_value . empty () ) {
            Ret += "\"" + _r_value + "\"";
        }
        else {
            return "INVALID COMPARISION";
        }
    }

    Ret += ")";

    return Ret;
}   /* _VR_ex_OP :: ToString () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
    /*  Negation
     */
class _VR_ex_NG : public _VR_ex {
public:
    _VR_ex_NG ( _VR_ex * Ex );
    ~_VR_ex_NG ();

    bool Validate ( const TL_Owp & Row ) const;
    string ToString () const;

private :
    _VR_ex * _ex;
};

_VR_ex_NG :: _VR_ex_NG ( _VR_ex * Ex )
:   _ex ( Ex )
{
    if ( _ex == NULL ) {
        throw TL_Exception ( "_VR_ex_NG : NULL expression passed" );
    }
    PLOGMSG(
            klogDebug,
            (klogDebug,
            "_VR_ex_NG",
            "severity=debug"
            )
            );
}   /* _VR_ex_NG :: _VR_ex_NG () */

_VR_ex_NG :: ~_VR_ex_NG ()
{
    if ( _ex != NULL ) {
        delete _ex;
        _ex = NULL;
    }
}   /* _VR_ex_NG :: ~_VR_ex_NG () */

bool
_VR_ex_NG :: Validate ( const TL_Owp & Row ) const
{
    return ! _ex -> Validate ( Row );
}   /* _VR_ex_NG :: Validate () */

string
_VR_ex_NG :: ToString () const
{
    return "!" + _ex -> ToString ();
}   /* _VR_ex_NG :: ToString () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
    /*  Simple LOGICAL operation: always AND or OR
     */
class _VR_ex_LOP : public _VR_ex {
private:
    typedef enum { eUND = -1, eAND, eOR } eLOp;
public:
    _VR_ex_LOP ( _VR_ex * Left, _VR_ex * Right, const _VR_t & Op );
    ~_VR_ex_LOP ();

    bool Validate ( const TL_Owp & Row ) const;
    string ToString () const;

private :
    _VR_ex * _left;
    _VR_ex * _right;
    eLOp _op;
};

_VR_ex_LOP :: _VR_ex_LOP (
                        _VR_ex * Left,
                        _VR_ex * Right,
                        const _VR_t & Op
)
:   _left ( Left )
,   _right ( Right )
,   _op ( eUND )
{
    if ( _left == NULL ) {
        throw TL_Exception ( "_VR_ex_LOP : NULL left expression passed" );
    }
    if ( _right == NULL ) {
        throw TL_Exception ( "_VR_ex_LOP : NULL right expression passed" );
    }

    if ( Op . Token () == "AND" ) {
        _op = eAND;
    }
    else {
        if ( Op . Token () == "OR" ) {
            _op = eOR;
        }
        else {
            throw TL_Exception ( string ( "_VR_ex_LOP : unknown logical operation [" ) + Op . Token () + "]" );
        }
    }
    PLOGMSG(
            klogDebug,
            (klogDebug,
            "_VR_ex_LOP ($(op))",
            "severity=debug,op=%s",
            Op . Token () . c_str ()
            )
            );
}   /* _VR_ex_LOP :: _VR_ex_LOP () */

_VR_ex_LOP :: ~_VR_ex_LOP ()
{
    _op = eUND;

    if ( _left != NULL ) {
        delete _left;
        _left = NULL;
    }
    if ( _right != NULL ) {
        delete _right;
        _right = NULL;
    }
}   /* _VR_ex_LOP :: ~_VR_ex_LOP () */

bool
_VR_ex_LOP :: Validate ( const TL_Owp & Row ) const
{
    if ( _op == eAND ) {
        return _right -> Validate ( Row )
                && _left -> Validate ( Row );
    }

    if ( _op == eOR ) {
        return _right -> Validate ( Row )
                || _left -> Validate ( Row );
    }

    // throw TL_Exception ( "_VR_ex_LOP : unknown logical operation" );
    return false;
}   /* _VR_ex_LOP :: Validate () */

string
_VR_ex_LOP :: ToString () const
{
    string Ret = "(";

    Ret += _right -> ToString ();

    switch ( _op ) {
        case eOR : Ret += " OR "; break;
        case eAND : Ret += " AND "; break;
        default : return "InVaLiD logical operation";  

    }

    Ret += _left -> ToString ();

    Ret += ")";

    return Ret;
}   /* _VR_ex_LOP :: ToString () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _VR_p
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _VR_p {
public :
    _VR_p ();
    ~_VR_p ();

    _VR_ex * Parse ( const string & Expression ) const;

private :
    _VR_ex * _compile ( vector < _VR_t > Tokens ) const;
    _VR_ex * _readEx (
                    vector < _VR_t > Tokens,
                    size_t & Pos,
                    bool Expand = true
                    ) const;
    _VR_ex * _readActualEx (
                    vector < _VR_t > Tokens,
                    size_t & Pos
                    ) const;

    _VR_ex * _ex;
};

_VR_p :: _VR_p ()
{
}   /* _VR_p :: _VR_p () */

_VR_p :: ~_VR_p ()
{
}   /* _VR_p :: ~_VR_p () */

_VR_ex *
_VR_p :: Parse ( const string & Expression ) const
{
        /*  Read all tokens
         */
    _VR_tp Tp;
    Tp . Init ( Expression );
    vector < _VR_t > Tokens;
    Tp . PullAll ( Tokens );

        /*  Compile Expression
         */
    _VR_ex * Ex = _compile ( Tokens );
    if ( Ex == NULL ) {
        throw TL_Exception ( "ERROR can not parse rule :LOL:" );
    }

    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[ IN] [$(expr)]",
            "severity=debug,expr=%s",
            Expression . c_str ()
            )
            );
    PLOGMSG(
            klogDebug,
            (klogDebug,
            "[OUT] [$(expr)]",
            "severity=debug,expr=%s",
            Ex -> ToString () . c_str ()
            )
            );

    return Ex;
}   /* _VR_p :: Parse () */

_VR_ex *
_VR_p :: _compile ( vector < _VR_t > Tokens ) const
{
    size_t TokQty = Tokens . size ();
    size_t Pos = 0;

    if ( TokQty == 1 ) {
        if ( Tokens [ 0 ] . Type () == _VR_t :: eANY ) {
            return new _VR_ex_ANY ( Tokens [ 0 ] );
        }
    }

    _VR_ex * Left = NULL;

    while ( Pos < TokQty ) {
        if ( Left == NULL ) {
            Left = _readEx ( Tokens, Pos );
        }
        else {
            _VR_t T = Tokens [ Pos ];

            Pos ++;

            if ( TokQty <= Pos ) {
                throw TL_Exception ( string ( "operation \"" ) + T . Token () + "\" requires expression" );
            }

            _VR_ex * Right = _readEx ( Tokens, Pos, false );
            if ( Right == NULL ) {
                throw TL_Exception ( string ( "NULL expression passed for operation \"" ) + T . Token () + "\"" );
            }

            Left = new _VR_ex_LOP ( Left, Right, T );
        }
    }

    return Left;
}   /* _VR_p :: _compile () */

_VR_ex *
_VR_p :: _readEx (
                vector < _VR_t > Tokens,
                size_t & Pos,
                bool Expand
) const
{
    size_t TokQty = Tokens . size ();

    if ( TokQty <= Pos ) {
            /* Kinda weird, but not error
             */
        return NULL;
    }

    _VR_ex * RetVal = NULL;

    switch ( Tokens [ Pos ] . Type () ) {
        case _VR_t :: eLB :
            Pos ++;
            RetVal = _readEx ( Tokens, Pos );
            if ( Tokens [ Pos ] . Type () != _VR_t :: eRB ) {
                throw TL_Exception ( "ERROR unbalansed brackets :LOL:" );
            }
            Pos ++;
            break;

        case _VR_t :: eNG :
            Pos ++;
            RetVal = new _VR_ex_NG ( _readEx ( Tokens, Pos, false ) );
            break;

        case _VR_t :: eWR :
            RetVal = _readActualEx ( Tokens, Pos );
            break;

        default:
            throw TL_Exception ( string ( "ERROR: unexpected token [" ) + Tokens [ Pos ] . Token () + "]" );
    }


    while ( Expand && Pos < TokQty ) {
        if ( Tokens [ Pos ] . Type () != _VR_t :: eLOP ) {
            break;
        }

        _VR_t T = Tokens [ Pos ];
        Pos ++;
        RetVal = new _VR_ex_LOP ( RetVal, _readEx ( Tokens, Pos, false ), T );
    }

    return RetVal;
}   /* _VR_p :: _readEx () */

_VR_ex *
_VR_p :: _readActualEx ( vector < _VR_t > Tokens, size_t & Pos ) const
{
    size_t TokQty = Tokens . size ();

    if ( TokQty <= Pos ) {
        throw TL_Exception ( "ERROR: unexpected position" );
    }

        /* Here we are */
    if ( Tokens [ Pos ] . Type () != _VR_t :: eWR ) {
        throw TL_Exception ( "ERROR: unexpected token type" );
    }

    size_t TokLeft = TokQty - Pos;

    if ( 3 <= TokLeft ) {
            /*  Simple comparision : three tokens with possible bracket
             */
        if ( Tokens [ Pos + 1 ] . Type () == _VR_t :: eOP ) {
            _VR_ex * Ex = new _VR_ex_OP (
                                        Tokens [ Pos ],
                                        Tokens [ Pos + 1 ],
                                        Tokens [ Pos + 2 ]
                                        );
            Pos += 3;
            return Ex;
        }
    }

        /*  Just a field
         */
    Pos ++;
    return new _VR_ex_WR ( Tokens [ Pos - 1 ] );
}   /* _VR_p :: _readActualEx () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _VR_ch - actually class will validate content of TraceInfo row
 *_          object
 *_
 *_ There will be couple classes
 *_     _VR_ch_FieldType - do basics checks on field type and value
 *_     _VR_ch_FieldPresence - check that field present or not
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _VR_ch {
public :
    static _VR_ch * Create ( const TL_Owp & Owp );

    virtual void Init ( const TL_Owp & Owp );

            /*  if it returns false, Message could contain reason why
             */
    virtual bool Validate (
                        const TL_Owp & Row,
                        string & Message
                        ) const = 0;

    inline const string & Name () const { return _name; };

protected :
    _VR_ch ();

public :
    virtual ~_VR_ch ();

private :
    string _name;
};

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
    /*  Simple TYPE checks
     */
    /*  There are four types recognized : date, float, int, string.
     */
class _VR_ch_FieldType : public _VR_ch {
public:
    typedef _VR_ch PAPA;
    typedef enum {
        eInvalid = 0,
        eInt,
        eDate,
        eFloat,
        eString
    }   eType;

public:
    _VR_ch_FieldType ();
    ~_VR_ch_FieldType ();

    virtual void Init ( const TL_Owp & Owp );
    void Dispose ();

    bool Validate ( const TL_Owp & Row, string & Message ) const;

private:
    bool _validateInt ( const TL_Owp & Row, string & Message ) const;
    bool _validateDate ( const TL_Owp & Row, string & Message ) const;
    bool _validateFloat ( const TL_Owp & Row, string & Message ) const;
    bool _validateString ( const TL_Owp & Row, string & Message ) const;

private:
    string _field;
    eType _type;
    int32_t _min;
    bool _min_set;
    int32_t _max;
    bool _max_set;
    string _alphabet;
    string _other;

};

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
    /*  Simple TYPE checks
     */
class _VR_ch_FieldPresence : public _VR_ch {
public:
    typedef _VR_ch PAPA;

public:
    _VR_ch_FieldPresence ();
    ~_VR_ch_FieldPresence ();

    virtual void Init ( const TL_Owp & Owp );

    bool Validate ( const TL_Owp & Row, string & Message ) const;

private:
    TL_SIMap _yes_no;
};

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
    /*  Basic class implementation
     */
_VR_ch :: _VR_ch ()
{
}   /* _VR_ch :: _VR_ch () */

_VR_ch :: ~_VR_ch ()
{
    _name . clear ();
}   /* _VR_ch :: ~_VR_ch () */

_VR_ch *
_VR_ch :: Create ( const TL_Owp & Owp )
{
    string Action = Owp . Value ( _VR_ACTION );
    if ( Action . empty () ) {
        throw TL_Exception ( string ( "ERROR: Action is not defined for validation rule \"" ) + Owp . Name () + "\"" );
    }

    _VR_ch * RetVal = NULL;

    if ( Action == _VR_ACTION_FTYPE ) {
        RetVal = new _VR_ch_FieldType ();
    }

    if ( Action == _VR_ACTION_FPRESENCE ) {
        RetVal = new _VR_ch_FieldPresence ();
    }

    if ( RetVal == NULL ) {
        throw TL_Exception ( string ( "ERROR: Unknown action \"" ) + Action + "\" for validation rule \"" + Owp . Name () + "\"" );
    }

    RetVal -> Init ( Owp );

    return RetVal;
}   /* _VR_ch :: Create () */

void
_VR_ch :: Init ( const TL_Owp & Owp )
{
    _name = Owp . Name ();
}   /* _VR_ch :: Init () */

bool
_VR_ch :: Validate ( const TL_Owp & Row, string & Message ) const
{
    Message . clear ();

    throw TL_Exception ( "ERROR: calling unimplemented 'Validate' method" );
}   /* _VR_ch :: Validate () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
    /*  FieldType check implementation
     */
_VR_ch_FieldType :: _VR_ch_FieldType ()
:   _VR_ch ()
,   _field ( "" )
,   _type ( eInvalid )
,   _min ( 0 )
,   _min_set ( false )
,   _max ( 0 )
,   _max_set ( false )
,   _alphabet ( "" )
,   _other ( "" )
{
}   /* _VR_ch_FieldType :: _VR_ch_FieldType () */

_VR_ch_FieldType :: ~_VR_ch_FieldType ()
{
    Dispose ();
}   /* _VR_ch_FieldType :: ~_VR_ch_FieldType () */

void
_VR_ch_FieldType :: Init ( const TL_Owp & Owp )
{
        /* Some scolding here */
    string Action = Owp . Value ( _VR_ACTION );
    if ( Action != _VR_ACTION_FTYPE ) {
        throw TL_Exception ( "_VR_ch_FieldType: Inproper action" );
    }

    PAPA :: Init ( Owp );

    Dispose ();

    _field = Owp . Value ( _VR_FIELD );
    if ( _field . empty () ) {
        throw TL_Exception ( string ( "_VR_ch_FieldType: field name is not defined for rule \"" ) + Name () + "\"" );
    }

    if ( ! TL_TraceFields :: Has ( _field ) ) {
        throw TL_Exception ( string ( "_VR_ch_FieldType: Can not find field with name \"" ) + _field + "\" for rule \"" + Name () + "\"" );
    }

        /* Determine type of */
    string Str = Owp . Value ( _VR_TYPE );
    if ( Str == _VR_TYPE_DATE ) { _type = eDate; }
    else {
        if ( Str == _VR_TYPE_FLOAT ) { _type = eFloat; }
        else {
            if ( Str == _VR_TYPE_INT ) { _type = eInt; }
            else {
                if ( Str == _VR_TYPE_STRING ) { _type = eString; }
                else {
                    _type = eInvalid;
                    throw TL_Exception ( string ( "_VR_ch_FieldType: Invalid field type \"" ) + Str + "\"" );
                }
            }
        }
    }

        /* Reading Min-Max */
    Str = Owp . Value ( _VR_MAX );
    if ( ! Str . empty () ) {
        _max = TL_StringU :: FromStr < int32_t > ( Str );
        _max_set = true;
    }
    Str = Owp . Value ( _VR_MIN );
    if ( ! Str . empty () ) {
        _min = TL_StringU :: FromStr < int32_t > ( Str );
        _min_set = true;
    }

    _alphabet = Owp . Value ( _VR_ALPHABET );
    _other = Owp . Value ( _VR_OTHER );
}   /* _VR_ch_FieldType :: Init () */

void
_VR_ch_FieldType :: Dispose ()
{
    _field . clear ();
    _type = eInvalid;
    _min = 0;
    _min_set = false;
    _max = 0;
    _max_set = false;
    _alphabet . clear ();
    _other . clear ();
}   /* _VR_ch_FieldType :: Dispose () */

bool
_VR_ch_FieldType :: Validate (
                            const TL_Owp & Row,
                            string & Message
) const
{
        /*  If field is not present, that is not an error
         */
    if ( Row . Has ( _field ) ) {
            /*  Here we are :D
             */
        switch ( _type ) {
            case eInt : return _validateInt ( Row, Message );
            case eDate : return _validateDate ( Row, Message );
            case eFloat : return _validateFloat ( Row, Message );
            case eString : return _validateString ( Row, Message );
            default : break;
        }
    }

    return true;
}   /* _VR_ch_FieldType :: Validate () */

bool
_VR_ch_FieldType :: _validateInt (
                                    const TL_Owp & Row,
                                    string & Message
) const
{
    Message . clear ();

        /*  While validation of string suppose empty string, that
         *  is not true for int values
         */
    string Str = Row . Value ( _field );
    if ( Str . empty () ) {
        stringstream Out;
        Out << "[_validateInt] field [" << _field << "]";
        Out << " contains empty value ( suppose to be INT )";
        Message = Out . str ();

        return false;
    }

    for ( size_t llp = 0; llp < Str . length (); llp ++ ) {
        if ( ! isdigit ( Str . at ( llp ) ) ) {
            stringstream Out;
            Out << "[_validateInt] field [" << _field << "]";
            Out << " value contains non digit characters ( suppose to be INT )";
            Message = Out . str ();

            return false;
        }
    }

    int32_t Val = 0;
    try {
        Val = TL_StringU :: FromStr < int32_t > ( Str );
    }
    catch ( ... ) {
        stringstream Out;
        Out << "[_validateInt] field [" << _field << "]";
        Out << " value can not be converted to int ( suppose to be INT )";
        Message = Out . str ();

        return false;
    }

        /* Min/Max val */
    if ( _min_set && Val < _min ) {
        stringstream Out;
        Out << "[_validateInt] field [" << _field << "]";
        Out << " value suppose to be greater or equal to [";
        Out << TL_StringU :: ToStr < int32_t > ( _min ) << "]";
        Message = Out . str ();

        return false;
    }

    if ( _max_set && _max < Val ) {
        stringstream Out;
        Out << "[_validateInt] field [" << _field << "]";
        Out << " value suppose to be less or equal to [";
        Out << TL_StringU :: ToStr < int32_t > ( _max ) << "]";
        Message = Out . str ();

        return false;
    }

    return true;
}   /* _VR_ch_FieldType :: _validateInt () */

bool
_VR_ch_FieldType :: _validateDate (
                                    const TL_Owp & Row,
                                    string & Message
) const
{
    Message . clear ();

    /* JOJOBA */
    /* Rules file contains 2 different types of date : run_date and
     * collection_date.
     * We know two more dates : load_date and update_date
     * run_date is in format 05-Jan-2013; 11-Jun-2008
     * collection_date is like "Thu Feb  6 16:33:00 2014"
     *                         "Mon Feb 24 13:48:00 2014"
     */
    /* Somebody here did not have any check */
        /*  While validation of string suppose empty string, that
         *  is not true for date values
         */
    string Str = Row . Value ( _field );
    if ( Str . empty () ) {
        stringstream Out;
        Out << "[_validateDate] field [" << _field << "]";
        Out << " contains empty value ( suppose to be DATE )";
        Message = Out . str ();

        return false;
    }

    if ( TL_DateU :: Date ( Str ) == 0 ) {
        stringstream Out;
        Out << "[_validateDate] field [" << _field << "]";
        Out << " contans invalid date [" << Str << "]";
        Message = Out . str ();

        return false;
    }

    return true;
}   /* _VR_ch_FieldType :: _validateDate () */

bool
_VR_ch_FieldType :: _validateFloat (
                                    const TL_Owp & Row,
                                    string & Message
) const
{
    Message . clear ();

        /*  While validation of string suppose empty string, that
         *  is not true for float values
         */
    string Str = Row . Value ( _field );
    if ( Str . empty () ) {
        stringstream Out;
        Out << "[_validateFloat] field [" << _field << "]";
        Out << " contains empty value ( suppose to be FLOAT )";
        Message = Out . str ();

        return false;
    }

    float Val = 0;
    try {
        Val = TL_StringU :: FromStr < float > ( Str );
    }
    catch ( ... ) {
        stringstream Out;
        Out << "[_validateFloat] field [" << _field << "]";
        Out << " value can not be converted to float [" << Str << "]";
        Message = Out . str ();

        return false;
    }

        /* Min/Max val */
    if ( _min_set && Val < ( float ) _min ) {
        stringstream Out;
        Out << "[_validateFloat] field [" << _field << "]";
        Out << " value suppose to be greater or equal to ";
        Out << "[" << TL_StringU :: ToStr < int32_t > ( _min ) << "]";
        Message = Out . str ();

        return false;
    }

    if ( _max_set && ( float ) _max < Val ) {
        stringstream Out;
        Out << "[_validateFloat] field [" << _field << "]";
        Out << " value suppose to be less or equal to ";
        Out << "[" << TL_StringU :: ToStr < int32_t > ( _max ) << "]";
        Message = Out . str ();

        return false;
    }

    return true;
}   /* _VR_ch_FieldType :: _validateFloat () */

bool
_VR_ch_FieldType :: _validateString (
                                    const TL_Owp & Row,
                                    string & Message
) const
{
    Message . clear ();

    string Str = Row . Value ( _field );

/* Really, another study of Michael codes shows 
   that the string could be empty

    if ( Str . empty () ) {
        stringstream Out;
        Out << "[_validateString] field [" << _field << "]";
        Out << " contains empty value ( suppose to be STRING )";
        Message = Out . str ();

        return false;
    }
*/

        /* Min/Max len */
    if ( _min_set && ( int32_t ) Str . length () < _min ) {
        stringstream Out;
        Out << "[_validateString] field [" << _field << "]";
        Out << " value length suppose to be greater or equal to ";
        Out << "[" << TL_StringU :: ToStr < int32_t > ( _min ) << "]";
        Message = Out . str ();

        return false;
    }

    if ( _max_set && _max < ( int32_t ) Str . length () ) {
        stringstream Out;
        Out << "[_validateString] field [" << _field << "]";
        Out << " value length suppose to be less or equal to ";
        Out << "[" << TL_StringU :: ToStr < int32_t > ( _max ) << "]";
        Message = Out . str ();

        return false;
    }

    if ( ! Str . empty () ) {
            /* Alphabet */
        if ( ! _alphabet . empty () ) {
            bool AB [ 256 ];
            memset ( AB, 0, sizeof ( AB ) );
            for ( size_t llp = 0; llp < _alphabet . length (); llp ++ ) {
                AB [ ( unsigned char ) _alphabet . at ( llp ) ] = true;
            }

            for ( size_t llp = 0; llp < Str . length (); llp ++ ) {
                if ( ! AB [ ( unsigned char ) Str . at ( llp ) ] ) {
                    stringstream Out;
                    Out << "[_validateString] field [" << _field << "]";
                    Out << " contans non-alphabet characters";
                    Out << "[" << Str << "] [" << ( char ) Str . at ( llp ) << "] [" << _alphabet << "]";
                    Message = Out . str ();

                    return false;
                }
            }
        }
    }

    return true;
}   /* _VR_ch_FieldType :: _validateString () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
    /*  FieldPresence implementation
     */
_VR_ch_FieldPresence :: _VR_ch_FieldPresence ()
:   _VR_ch ()
,   _yes_no ()
{
}   /* _VR_ch_FieldPresence :: _VR_ch_FieldPresence () */

_VR_ch_FieldPresence :: ~_VR_ch_FieldPresence ()
{
    _yes_no . clear ();
}   /* _VR_ch_FieldPresence :: ~_VR_ch_FieldPresence () */

void
_VR_ch_FieldPresence :: Init ( const TL_Owp & Owp )
{
        /* Some scolding here */
    string Action = Owp . Value ( _VR_ACTION );
    if ( Action != _VR_ACTION_FPRESENCE ) {
        throw TL_Exception ( "_VR_ch_FieldPresence: Inproper action" );
    }

    PAPA :: Init ( Owp );

    _yes_no . clear ();

    string Fields = Owp . Value ( _VR_YES_FIELDS );
    TL_SVec Tokens;
    size_t Qty = TL_StringU :: Tokenize ( Fields, ' ', Tokens );
    if ( Qty != 0 ) {
        for ( size_t llp = 0; llp < Qty; llp ++ ) {
            string Field = Tokens [ llp ];
            if ( ! TL_TraceFields :: Has ( Field ) ) {
                throw TL_Exception ( string ( "_VR_ch_FieldPresence: Can not find field with name \"" ) + Field + "\"" );
            }

            _yes_no [ Field ] = 1;
        }
    }

    Fields = Owp . Value ( _VR_NO_FIELDS );
    Tokens . clear ();
    Qty = TL_StringU :: Tokenize ( Fields, ' ', Tokens );
    if ( Qty != 0 ) {
        for ( size_t llp = 0; llp < Qty; llp ++ ) {
            string Field = Tokens [ llp ];

            if ( _yes_no . find ( Field ) != _yes_no . end () ) {
                throw TL_Exception ( string ( "_VR_ch_FieldPresence: Ambiguous rule: field \"" ) + Field + "\" belongs to both Yes and No" );
            }

            if ( ! TL_TraceFields :: Has ( Field ) ) {
                throw TL_Exception ( string ( "_VR_ch_FieldPresence: Can not find field with name \"" ) + Field + "\"" );
            }

            _yes_no [ Field ] = 0;
        }
    }
}   /* _VR_ch_FieldPresence :: Init () */

bool
_VR_ch_FieldPresence :: Validate (
                                const TL_Owp & Row,
                                string & Message
) const
{
    Message . clear ();

    for ( 
        TL_SIMap :: const_iterator It = _yes_no . begin ();
        It != _yes_no . end ();
        It ++
    ) {
        string Field = It -> first;

        string Val = Row . Value ( Field );
        if ( It -> second == 0 ) {
            /* No field should be */
            if ( ! Val . empty () ) {
                stringstream Out;
                Out << "[_VR_ch_FieldPresence] field [" << Field << "]";
                Out << " should not present according to rule ";
                Out << "[" << Name () << "]";
                Message = Out . str ();

                return false;
            }
        }
        else {
            if ( Val . empty () ) {
                stringstream Out;
                Out << "[_VR_ch_FieldPresence] field [" << Field << "]";
                Out << " should be present according to rule ";
                Out << "[" << Name () << "]";
                Message = Out . str ();

                return false;
            }
        }
    }

    return true;
}   /* _VR_ch_FieldPresence :: Validate () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _VRule
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

 /*))
  //    Because rules will be short, we are doing recursion ... :LOL:
 ((*/
class _VRule {
public:
    _VRule ();
    ~_VRule ();

        /*  Here the parsing is going on
         */
    void Init ( const TL_Owp & Owp );

        /*  Second You may validate that rule
         *  if validation failed, message may contain reason why
         */
    bool Validate ( const TL_Owp & Row, string & Message );

private:
    _VR_ex * _condition;
    _VR_ch * _validator;
};

_VRule :: _VRule ()
:   _condition ( NULL )
,   _validator ( NULL )
{
}   /* _VRule :: _VRule () */

_VRule :: ~_VRule ()
{
    if ( _condition != NULL ) {
        delete _condition;
        _condition = NULL;
    }

    if ( _validator != NULL ) {
        delete _validator;
        _validator = NULL;
    }
}   /* _VRule :: ~_VRule () */

void
_VRule :: Init ( const TL_Owp & Owp )
{
    string Cond = Owp . Value ( _VR_CONDITION );
    if ( Cond . empty () ) {
        throw TL_Exception ( "ERROR baad rule :LOL:" );
    }

    _condition = _VR_p () . Parse ( Cond );
    if ( _condition == NULL ) {
        throw TL_Exception ( "_VRule :: NULL condition set" );
    }

    _validator = _VR_ch :: Create ( Owp );
    if ( _validator == NULL ) {
        throw TL_Exception ( "_VRule :: NULL validator set" );
    }

}   /* _VRule :: _Init () */

bool
_VRule :: Validate ( const TL_Owp & Row, string & Message )
{
    Message . clear ();

    return _condition -> Validate  ( Row )
                            ? _validator -> Validate ( Row, Message )
                            : true
                            ;
}   /* _VRule :: _Validate () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ _VRulesRep
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
class _VRulesRep {
private:
    typedef vector < _VRule * > rVec;
public:
    _VRulesRep ();
    ~_VRulesRep ();

    void Init ( const std::string & Path );
    void Dispose ();

        /*  If validation failed, Message could contain reason why
         */
    bool Validate ( const TL_Owp & Row, string & Message );

private:

    rVec _rules;
};

static _VRulesRep _sVRuleRepInstance;

_VRulesRep :: _VRulesRep ()
{
}   /* _VRulesRep :: _VRulesRep () */

_VRulesRep :: ~_VRulesRep ()
{
    Dispose ();
}   /* _VRulesRep :: ~_VRulesRep () */

string
_resolvePaath ( const string & GenPath, const string & Path )
{
    if ( Path . at ( 0 ) == '/' ) {
        return Path;
    }

    size_t Pos = GenPath . rfind ( '/' );
    if ( Pos == string :: npos ) {
        return Path;
    }

    return GenPath . substr ( 0, Pos ) + '/' + Path;
}   /* _resolvePaath () */

void
_VRulesRep :: Init ( const std::string & Path )
{
    stack < string > toLoad;
    TL_SSet  loadNames;

    toLoad . push ( Path );
    loadNames . insert ( Path );

    while ( ! toLoad . empty () ) {
        string theP = toLoad . top ();
        toLoad . pop ();

        PLOGMSG(
                klogDebug,
                (klogDebug,
                "[Loading V-Rules]",
                "severity=debug,file=%s",
                theP . c_str ()
                )
                );

        TL_OwpVector OVec;
        OVec . Load ( theP );

        for ( size_t llp = 0; llp < OVec . Size (); llp ++ ) {
            string Name = OVec . Get ( llp ) . Name ();
            if ( Name == "INIT" ) {
                TL_SVec Keys;
                OVec . Get ( llp ) . ListKeys ( Keys );
                for ( size_t ppl = 0; ppl < Keys . size (); ppl ++ ) {
                    string Key = Keys [ ppl ];
                    string Val = OVec . Get ( llp ) . Value ( Key );

                    if ( loadNames . find ( Val ) == loadNames . end () ) {
                        loadNames . insert ( Val );
                        toLoad . push ( _resolvePaath ( Path, Val ) );
                    }
                }
            }
            else {
                _VRule * Rule = new _VRule;
                Rule -> Init ( OVec . Get ( llp ) );
                _rules . insert ( _rules . end (), Rule );
            }
        }
    }
}   /* _VRulesRep :: Init () */

void
_VRulesRep :: Dispose ()
{
    for ( size_t llp = 0; llp < _rules . size (); llp ++ ) {
        if ( _rules [ llp ] != NULL ) {
            delete _rules [ llp ];
            _rules [ llp ] = NULL;
        }
    }

    _rules . clear ();
}   /* _VRulesRep :: Dispose () */

bool
_VRulesRep :: Validate ( const TL_Owp & Row, string & Message )
{
    Message . clear ();

    for ( size_t llp = 0; llp < _rules . size (); llp ++ ) {
        if ( _rules [ llp ] != NULL ) {
            if ( ! _rules [ llp ] -> Validate ( Row, Message ) ) {
                return false;
            }
        }
    }
    return true;
}   /* _VRulesRep :: Validate () */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
 *_
 *_ TL_VRules
 *_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
TL_VRules :: TL_VRules ()
{
}   /* TL_VRules :: TL_VRules () */

TL_VRules :: ~TL_VRules ()
{
    Dispose ();
}   /* TL_VRules :: ~TL_VRules () */

void
TL_VRules :: Init ( const std::string & Path )
{
    _sVRuleRepInstance . Init ( Path );
}   /* TL_VRules :: Init () */

void
TL_VRules :: Dispose ()
{
    _sVRuleRepInstance . Dispose ();
}   /* TL_VRules :: Dispose () */

bool
TL_VRules :: Validate ( const TL_Owp & Row, string & Message )
{
    return _sVRuleRepInstance . Validate ( Row, Message );
}   /* TL_VRules :: Validate () */

