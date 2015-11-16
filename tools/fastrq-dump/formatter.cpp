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

#include "formatter.hpp"

#include <ngs/ReadIterator.hpp>

using namespace ngs;

namespace fastrq
{
    /*======================================================================
     * Fastq Formatter
     */
    struct FastQFormatter : FastRQFormatter
    {
        virtual void dump ( String run_name, uint64_t spotID, const ReadIterator &reads, const StringRef &bases ) const;

        FastQFormatter ()
        {
        }
    };

    FastRQFormatter * makeFastQFormatter ()
    {
        return new FastQFormatter ();
    }

    void FastQFormatter :: dump ( String run_name, uint64_t spotID, const ReadIterator &reads, const StringRef &bases ) const
    {
        StringRef read_name = reads . getReadName ();
        //StringRef read_id = reads . getReadId ();
        StringRef qualities = reads . getReadQualities ();
        
        std :: cout 
            << '@'
            << run_name
            << '.'
            << spotID
            //<< read_id
            << ' '
            << read_name
            << " length="
            << bases . size ()
            << '\n'
            << bases
            << '\n'
            << '+'
            << run_name
            << '.'
            << spotID
            //<< read_id
            << ' '
            << read_name
            << " length="
            << qualities . size ()
            << '\n'
            << qualities
            << '\n'
            ;
    }

    /*======================================================================
     * Fasta Formatter
     */
    struct FastAFormatter : FastRQFormatter 
    {
        virtual void dump ( String run_name, uint64_t spotID, const ReadIterator &reads, const StringRef &bases ) const;

        FastAFormatter ()
        {
        }
           
    };

    FastRQFormatter * makeFastAFormatter ()
    {
        return new FastAFormatter ();
    }

    void FastAFormatter :: dump ( String run_name, uint64_t spotID, const ReadIterator &reads, const StringRef &bases ) const
    {
        const size_t limit = 70;

        StringRef read_name = reads . getReadName ();

        if ( bases . size () > 0 )
        {
            StringRef base_str = bases;

            std :: cout 
                << '>'
                << run_name
                << '.'
                << spotID
                << ' '
                << read_name
                << " length="
                << base_str . size ()
                << '\n'
                ;
            
            while ( base_str.size () > limit )
            {
                StringRef sub = base_str.substr(0,limit);

                base_str = base_str.substr ( limit );


                std :: cout << sub << '\n';
            }
            std :: cout << base_str << '\n';
        }        
    }

}
