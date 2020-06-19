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

#ifndef _tl_ows_
#define _tl_ows_

namespace _tl_ {

/*  Object With State: Ready, Good, Bad ...
 */
class TL_Ows {
public:
    typedef enum {
        kBad = - 1,
        kReady,
        kGood
    } _St;

public:
    TL_Ows ( _St state = kReady );
    ~TL_Ows ();

    inline _St State () const { return _state; };
    inline bool IsStateGood () const { return _state == kGood; };
    inline bool IsStateBad () const { return _state == kBad; };
    inline bool IsStateReady () const { return _state == kReady; };

protected:
    inline void SetState ( _St state ) const  { _state = state; };
    inline void SetStateGood () const { _state = kGood; };
    inline void SetStateBad () const { _state = kBad; };
    inline void SetStateReady () const { _state = kReady; };

private:
    mutable _St _state;
};

}   /* namespace _tl_ */

#endif /* _tl_ows_ */

