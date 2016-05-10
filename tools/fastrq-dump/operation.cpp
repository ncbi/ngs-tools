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

#include "operation.hpp"

#include <ngs/StringRef.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>

#include "formatter.hpp"
#include "filter.hpp"

using namespace ngs;

namespace fastrq
{
     /*======================================================================
     * FastRQOperation
     */
    FastRQOperation :: FastRQOperation ()
    {}
    FastRQOperation :: ~FastRQOperation ()
    {}

    void FastRQOperation :: run ( ngs :: ReadCollection & run, const FastRQFilter * filter, const FastRQFormatter * fmt ) const
    {
        throw __FUNCTION__;
    }

    void FastRQOperation :: run ( ngs :: ReadCollection & run, uint64_t start_row, uint64_t num_rows,
                                  const FastRQFilter * filter, const FastRQFormatter * fmt ) const
    {
        throw  __FUNCTION__;
    }

    void FastRQOperation :: run ( ngs :: ReferenceSequence & ref, const FastRQFilter * filter, const FastRQFormatter * fmt ) const
    {
        throw  __FUNCTION__;
    }
    void FastRQOperation :: run ( ngs :: ReferenceSequence & ref, uint64_t start_row, uint64_t num_rows,
                                  const FastRQFilter * filter, const FastRQFormatter * fmt ) const
    {
        throw  __FUNCTION__;
    }

     /*======================================================================
     * ReadOperation
     */

    struct ReadOperation : FastRQOperation
    {
        // filter ok NULL
        virtual void run ( ReadCollection & run, const FastRQFilter * filter, const FastRQFormatter * fmt ) const;
        virtual void run ( ReadCollection & run, uint64_t start_row, uint64_t num_rows,
                           const FastRQFilter * filter, const FastRQFormatter * fmt ) const {}

        ReadOperation ()
        : start_row ( 0 ), num_rows ( 0 )
        {
        }
        
        ReadOperation ( uint64_t _start_row, uint64_t _num_rows )
        : start_row ( _start_row ), num_rows ( _num_rows )
        {
        }

        uint64_t start_row, num_rows;
    };

    FastRQOperation * makeReadOperation ()
    {
        return new ReadOperation ();
    }

    FastRQOperation * makeReadOperation ( uint64_t start_row, uint64_t num_rows )
    {
        return new ReadOperation ( start_row, num_rows );
    }

    void ReadOperation :: run ( ReadCollection & run, const FastRQFilter * filter, const FastRQFormatter * fmt ) const
    {
        uint64_t spotID = 1;
        String run_name = run . getName ();
        
        ReadIterator reads = run . getReads ( Read :: all );
        while ( reads . nextRead () )
        {
            StringRef bases = reads . getReadBases ();
            if ( filter == 0 )
                fmt -> dump ( run_name, spotID, reads, bases );
            else
            {
                if ( filter -> pass ( reads, bases ) )
                    fmt -> dump ( run_name, spotID, reads, bases );
            }
            
            ++ spotID;
        }
    }

    /*======================================================================
     * FragmentOperation
     */
    struct RefSeqOperation : FastRQOperation
    {
        // filter ok NULL
        virtual void run ( ReferenceSequence & ref, const FastRQFilter * filter, const FastRQFormatter * fmt ) const;
        virtual void run ( ReferenceSequence & ref, uint64_t start_row, uint64_t num_rows,
                           const FastRQFilter * filter, const FastRQFormatter * fmt ) const {}

        
        RefSeqOperation ()
        : start_row ( 0 ), num_rows ( 0 )
        {
        }
        
        RefSeqOperation ( uint64_t _start_row, uint64_t _num_rows )
        : start_row ( _start_row ), num_rows ( _num_rows )
        {
        }

        uint64_t start_row, num_rows;
    };

    FastRQOperation * makeRefSeqOperation ()
    {
        return new RefSeqOperation ();
    }

    void RefSeqOperation :: run ( ReferenceSequence & ref, const FastRQFilter * filter, const FastRQFormatter * fmt ) const
    {
    }

}
