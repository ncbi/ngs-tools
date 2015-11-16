/*===========================================================================
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
* ===========================================================================
*
*/

#ifndef _hpp_ngs_filter_
#define _hpp_ngs_filter_

#include "ngs/StringRef.hpp"
#include <vector>

namespace ngs
{
    class Read;
    class Fragment;
    class ReferenceSequence;
}

namespace fastrq
{
    /*======================================================================
     * Filter
     */
    struct FastRQFilter
    {
        bool pass ( const ngs :: Read & read, const ngs :: String & bases ) const
        { return pass ( read, bases . data (), bases . size () ); }
        bool pass ( const ngs :: Read & read, const ngs :: StringRef & bases ) const
        { return pass ( read, bases . data (), bases . size () ); }

        bool pass ( const ngs :: Fragment & frag, const ngs :: String & bases ) const
        { return pass ( frag, bases . data (), bases . size () ); }
        bool pass ( const ngs :: Fragment & frag, const ngs :: StringRef & bases ) const
        { return pass ( frag, bases . data (), bases . size () ); }

        bool pass ( const ngs :: ReferenceSequence & ref, const ngs :: String & bases ) const
        { return pass ( ref, bases . data (), bases . size () ); }
        bool pass ( const ngs :: ReferenceSequence & ref, const ngs :: StringRef & bases ) const
        { return pass ( ref, bases . data (), bases . size () ); }

        FastRQFilter (){}
        virtual ~FastRQFilter (){}

        virtual bool pass ( const ngs :: Read & read, const char * bases, size_t blen ) const = 0;
        virtual bool pass ( const ngs :: Fragment & frag, const char * bases, size_t blen ) const = 0;
        virtual bool pass ( const ngs :: ReferenceSequence & ref, const char * bases, size_t blen ) const = 0;
    };

    FastRQFilter * makeLengthFilter ( size_t min_length );
    FastRQFilter * makeNCountFilter ( uint32_t max_count );

    /*======================================================================
     * FilterPurse
     */
    struct FilterPurse : FastRQFilter
    {
        virtual bool pass ( const ngs :: Read & read, const char * bases, size_t blen ) const;
        virtual bool pass ( const ngs :: Fragment & frag, const char * bases, size_t blen ) const;
        virtual bool pass ( const ngs :: ReferenceSequence & ref, const char * bases, size_t blen ) const;

        void addFilter ( FastRQFilter * filter );

        FilterPurse () {}
        ~FilterPurse ();

        std :: vector < FastRQFilter * > filters;
    };

    FilterPurse * makeFilterPurse ();
}

#endif //  _hpp_ngs_filter_
