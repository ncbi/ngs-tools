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

#include "filter.hpp"

using namespace ngs;

namespace fastrq
{
    /*======================================================================
     * FilterPurse
     */

    FilterPurse :: ~FilterPurse ()
    {
        int count = filters . size ();
        
        for ( int i = 0; i < count; ++ i )
            delete filters [ i ];
        
        filters . clear ();
    }

    FilterPurse * makeFilterPurse ()
    {
        return new FilterPurse ();
    }

    void FilterPurse :: addFilter ( FastRQFilter * filter )
    {
        if ( filter != 0 )
            filters . push_back ( filter );
    }

    bool FilterPurse :: pass ( const Read & read, const char * bases, size_t blen ) const
    {
        if ( bases == 0 )
            return false;

        if ( blen == 0 )
            return false;

        size_t f_count = filters . size ();

        for ( size_t i = 0; i < f_count; ++ i )
        {
            FastRQFilter *filter = filters [ i ];
            if ( filter -> pass ( read, bases, blen ) == false )
                return false;
        }
        return true;
    }

    bool FilterPurse :: pass ( const Fragment & frag, const char * bases, size_t blen ) const
    {
        if ( bases == 0 )
            return false;

        if ( blen == 0 )
            return false;

        size_t f_count = filters . size ();

        for ( size_t i = 0; i < f_count; ++ i )
        {
            FastRQFilter *filter = filters [ i ];
            if ( filter -> pass ( frag, bases, blen ) == false )
                return false;
        }
        return true;
    }

    bool FilterPurse :: pass ( const ReferenceSequence & ref, const char * bases, size_t blen ) const
    {
        if ( bases == 0 )
            return false;

        if ( blen == 0 )
            return false;

        size_t f_count = filters . size ();

        for ( size_t i = 0; i < f_count; ++ i )
        {
            FastRQFilter *filter = filters [ i ];
            if ( filter -> pass ( ref, bases, blen ) == false )
                return false;
        }
        return true;
    }

    /*======================================================================
     * Length Filter
     */
    struct LengthFilter : FastRQFilter
    {
        virtual bool pass ( const Read & read, const char * bases, size_t blen ) const;
        virtual bool pass ( const Fragment & frag, const char * bases, size_t blen ) const;
        virtual bool pass ( const ReferenceSequence & ref, const char * bases, size_t blen ) const;

        LengthFilter ( size_t _min_length )
        : min_length ( _min_length )
        {
        }
           
        size_t min_length;
    };

    FastRQFilter * makeLengthFilter ( size_t min_length )
    {
        return new LengthFilter ( min_length );
    }

    bool LengthFilter :: pass ( const Read & read, const char * bases, size_t blen ) const
    {
        return ( min_length > blen );
    }

    bool LengthFilter :: pass ( const Fragment & frag, const char * bases, size_t blen ) const
    {
        return ( min_length > blen );
    }

    bool LengthFilter :: pass ( const ReferenceSequence & ref, const char * bases, size_t blen ) const
    {
        return ( min_length > blen );
    }

    /*======================================================================
     * N-Count Filter
     */
    struct NCountFilter : FastRQFilter
    {
        virtual bool pass ( const Read & read, const char * bases, size_t blen ) const;
        virtual bool pass ( const Fragment & frag, const char * bases, size_t blen ) const;
        virtual bool pass ( const ReferenceSequence & ref, const char * bases, size_t blen ) const;

        bool int_pass ( const char *bases, size_t blen ) const;
        
        NCountFilter ( uint32_t _max_count )
        : max_count ( _max_count )
        {
        }
           
        size_t max_count;
    };

    FastRQFilter * makeNCountFilter ( uint32_t max_count )
    {
        return new NCountFilter ( max_count );
    }

    bool NCountFilter :: int_pass ( const char *bases, size_t blen ) const
    {
        if ( bases == 0 )
            return false;

        uint32_t count = 0;
        for ( size_t i = 0; i < blen; ++ i )
        {
            if ( bases [ i ] == 'N' )
                ++ count;

            if ( count == max_count )
                return false;
        }

        return true;
    }

    bool NCountFilter :: pass ( const Read & read, const char * bases, size_t blen ) const
    {
        return int_pass ( bases, blen );
    }

    bool NCountFilter :: pass ( const Fragment & frag, const char * bases, size_t blen ) const
    {
        return int_pass ( bases, blen );
    }

    bool NCountFilter :: pass ( const ReferenceSequence & ref, const char * bases, size_t blen ) const
    {
        return int_pass ( bases, blen );
    }

}
