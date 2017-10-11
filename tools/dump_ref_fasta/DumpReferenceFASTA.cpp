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

#include <ngs/ncbi/NGS.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>


#include <math.h>
#include <iostream>
#include <sstream>
#include <limits>

using namespace ngs;
using namespace std;

const String COLON = std::string( ":" );
const String DASH  = std::string( "-" );
const String DOT   = std::string( "." );

class Range
{
    public :
        uint64_t offset;
        uint64_t len;

        Range() : offset( 0 ), len( 0 ) {}
        
        Range( const String & spec )
        {
            const string::size_type dash_pos = spec.find( DASH );
            if ( dash_pos == string::npos )
            {
                const string::size_type dot_pos = spec.find( DOT );
                if ( dot_pos == string::npos )
                {
                    offset = to_uint64_t( spec );
                    len = 0;
                }
                else
                {
                    const String offset_string = spec.substr( 0, dot_pos );
                    offset = to_uint64_t( offset_string );
                    const String len_string = spec.substr( dot_pos + 1 );
                    len = len_string.empty() ? std::numeric_limits< uint64_t >::max() : to_uint64_t( len_string );
                }
            }
            else
            {
                const String offset_string = spec.substr( 0, dash_pos );
                offset = to_uint64_t( offset_string );
                const String end_string = spec.substr( dash_pos + 1 );
                uint64_t end = end_string.empty() ? std::numeric_limits< uint64_t >::max() : to_uint64_t( end_string );
                len = offset <= end ? end - offset : offset - end;
                if ( offset > end ) { offset = end; }
            }
            
            /* we expect the coordinates 1-based, but we are using internally 0-based offsets */
            if ( offset > 0 ) { offset--; }
        }
        
        static uint64_t to_uint64_t( const String & num )
        {
            uint64_t res = 0;
            std::istringstream ss( num );
            ss >> res;
            return res;
        }

        void adjust_by_reflen( const uint64_t reflen )
        {
            if ( offset == 0 )
            {
                if ( len > reflen || len == 0 ) { len = reflen; }
            }
            else
            {
                if ( offset >= reflen ) { offset = reflen - 1; }
                
                if ( len == 0 )
                {
                    len = reflen - offset;
                }
                else
                {
                    if ( offset + len > reflen ) { len = reflen - offset; }
                }
            }
        }
};

class Slice
{
    public :
        String name;
        Range range;
    
        Slice() {}
        
        // spec = 'name[:from[-to]]' or 'name[:from[.count]]'
        Slice( const String & spec )
        {
            const string::size_type colon_pos = spec.find( COLON );
            if ( colon_pos == string::npos )
            {
                name = spec;
            }
            else
            {
                name = spec.substr( 0, colon_pos );
                const String range_string = spec.substr( colon_pos + 1 );
                range = Range( range_string );
            }
        }
};

const uint64_t CHUNK_SIZE = 5000;
const uint64_t LINE_LEN = 70;

class DumpReferenceFASTA
{
    public :
        static void process ( const Reference & ref, const Range & range )
        {
            cout << '>' << ref.getCanonicalName () << ':' << range.offset + 1 << '.' << range.len << '\n';
            try
            {
                size_t line = 0;
                uint64_t count = 0;
                uint64_t offset = range.offset;
                uint64_t left = range.len;
                
                while ( left > 0 )
                {
                    StringRef chunk = ref.getReferenceChunk( offset, left > CHUNK_SIZE ? CHUNK_SIZE : left );
                    size_t chunk_len = chunk.size ();
                    for ( size_t chunk_idx = 0; chunk_idx < chunk_len; )
                    {
                        StringRef chunk_line = chunk.substr( chunk_idx, LINE_LEN - line );
                        line += chunk_line.size ();
                        chunk_idx += chunk_line.size ();

                        cout << chunk_line;
                        if ( line >= LINE_LEN )
                        {
                            cout << '\n';
                            line = 0;
                        }
                    }
                    offset += chunk_len;                    
                    count += chunk_len;
                    left -= chunk_len;
                }
                if ( line != 0 )
                    cout << '\n';
            }
            catch ( ErrorMsg x )
            {
                cerr <<  x.toString() << '\n';            
            }
        }
    
    static void run( const String & acc, Slice * slices, int n )
    {
        // open requested accession using SRA implementation of the API
        ReadCollection run = ncbi::NGS::openReadCollection( acc );
        Reference ref = run.getReference( slices[ 0 ].name );
        for ( int i = 0; i < n; ++i )
        {
            if ( i > 0 && !slices[ i ].name.empty() )
            {
                ref = run.getReference( slices[ i ].name );
            }
            slices[ i ].range.adjust_by_reflen( ref.getLength() );
            process( ref, slices[ i ].range );
        }
    }

    static void run( const String & acc )
    {

        // open requested accession using SRA implementation of the API
        ReadCollection run = ncbi::NGS::openReadCollection( acc );
        ReferenceIterator refs = run.getReferences();
        while ( refs.nextReference () )
        {
            Slice slice( refs.getCanonicalName() );
            slice.range.adjust_by_reflen( refs.getLength() );
            process( refs, slice.range );
        }
    }
};

static void print_help ( void )
{
    cout << "Usage: dump-ref-fasta accession [ reference[ slice ] ] ... " << endl;
}

static void print_version ( void )
{
    cout << "dump-ref-fata : 2.9.0" << endl;
}


int main ( int argc, char const *argv[] )
{
    if ( argc < 2 )
    {
        print_help ();
    }
    else try
    {
        ncbi::NGS::setAppVersionString ( "DumpReferenceFASTA.1.0.0" );
        const String arg1( argv[ 1 ] );
        if ( arg1 == "-h" || arg1 == "--help" )
        {
            print_help ();
        }
        else if ( arg1 == "-v" || arg1 == "--version" )
        {
            print_version ();
        }
        else
        {
            if ( argc > 2 )
            {
                int n = ( argc - 2 );
                Slice * slices = new Slice[ n ];
                for ( int i = 0; i < n; ++i )
                    slices[ i ] = Slice( argv[ 2 + i ] );
                DumpReferenceFASTA::run( arg1, slices, n );
                delete [] slices;
            }
            else
            {
                DumpReferenceFASTA::run( arg1 );
            }
        }
        return 0;
    }
    catch ( ErrorMsg &x )
    {
        cerr <<  x.toString() << '\n';
    }
    catch ( exception &x )
    {
        cerr <<  x.what() << '\n';
    }
    catch ( ... )
    {
        cerr <<  "unknown exception\n";
    }

    return 10;
}
