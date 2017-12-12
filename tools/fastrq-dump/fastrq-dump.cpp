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
#include <ngs/StringRef.hpp>

#include "settings.hpp"
#include "filter.hpp"
#include "formatter.hpp"
#include "operation.hpp"

#include <kfc/rsrc.h>
#include <kapp/main.h>

#include <iostream>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

using namespace ngs;

namespace fastrq
{
    static
    FastRQFormatter * make_formatter ( const FastRQSettings &settings )
    {
        if ( settings . fmt == fmt_fasta )
            return makeFastAFormatter ();

        return makeFastQFormatter ();
    }

    static
    FastRQFilter * make_filters ( const FastRQSettings &settings )
    {
        FastRQFilter * filter;
        FilterPurse * purse = makeFilterPurse ();

        if ( settings . min_length > 0 )
        {
            filter = makeLengthFilter ( settings . min_length );
            purse -> addFilter ( filter );

            filter = 0;
        }

        if ( settings . n_count > 0 )
        {
            filter = makeNCountFilter ( settings . n_count );
            purse -> addFilter ( filter );

            filter = 0;
        }

        return purse;
    }

    static
    FastRQOperation * make_operation ( const FastRQSettings &settings )
    {
        if ( settings . op == op_refseq )
            return makeRefSeqOperation ();

        return makeReadOperation ();
    }

    rc_t CC Usage ( struct Args const * args )
    {
    }

    static
    void handle_help ()
    {
        std :: cout
            << '\n'
            << "Usage:\n"
            << "  " << UsageDefaultName << " [options] <accession>"
            << "\n\n"
            << "Options:\n"
            << "  -l|-L|--spot_id_length <value>   Minimum spot_id length.\n"
            << "  -n|-N|--max_n_count <value>      Filter out reads with more than <value> 'N' bases.\n"
            << "  --fasta                          Output in FASTA format.\n"
            << "  -h|-?|--help                     Output brief explanation for the program. \n"
            << "  -V|--version                     Display the version of the program then\n"
            << "                                   quit.\n"
            ;
        HelpVersion ( UsageDefaultName, KAppVersion () );
    }

    static
    const char * nextArg ( int & i, int argc, const char * argv [] )
    {
        if ( ++ i >= argc )
            throw "aaah!";

        return argv [ i ];
    }

    static
    const char * nextArg ( const char * & arg_ref, int & i, int argc, const char * argv [] )
    {
        const char * arg = arg_ref;
        arg_ref = "\0";

        if ( arg [ 0 ] != 0 )
            return arg;

        return nextArg ( i, argc, argv );
    }

    static
    long smart_atoi ( const char * text )
    {
        char * end;
        long val = strtol ( text, & end, 0 );
        if ( end [ 0 ] != 0 )
            throw "bad integer specification";
        return val;
    }

    static
    bool parse_cmdline ( int argc, const char *argv [], FastRQSettings &settings )
    {
        for ( int i = 1; i < argc; ++ i )
        {
            const char * arg = argv [ i ];
            if ( arg [ 0 ] != '-' )
            {
                // have an input run
                argv [ ++ settings . num_accessions ] = arg;
            }
            else do switch ( ( ++ arg ) [ 0 ] )
            {
            case 'h':
            case '?':
                handle_help ();
                return true;
            case 'l':
            case 'L': // Get the minimum spot_id length
                ++ arg;
                settings . min_length = ( size_t ) smart_atoi ( nextArg ( arg, i, argc, argv ) );
                break;
            case 'n':
            case 'N':
                ++ arg;
                settings . n_count = ( uint32_t ) smart_atoi ( nextArg ( arg, i, argc, argv ) );
                break;
            case 'V':
                HelpVersion ( UsageDefaultName, KAppVersion () );
                return true;
            case '-':
                ++ arg;
                if ( strcmp ( arg, "fasta"  ) == 0 )
                    settings . fmt = fmt_fasta;
                else if ( strcmp ( arg, "spot_id_length"  ) == 0 )
                    settings . min_length = ( size_t ) atoi ( nextArg ( i, argc, argv ) );
                else if ( strcmp ( arg, "max_n_count" ) == 0 )
                    settings . n_count = ( uint32_t ) atoi ( nextArg ( i, argc, argv ) );
                else if ( strcmp ( arg, "help"  ) == 0 )
                {
                    handle_help ();
                    return true;
                }
                else if ( strcmp ( arg, "version"  ) == 0 )
                {
                    HelpVersion ( UsageDefaultName, KAppVersion () );
                    return true;
                }
                else
                {
                    throw "Invalid Argument";
                }

                arg = "\0";

                break;
            default:
                throw "Invalid argument";
            }
            while ( arg [ 1 ] != 0 );
        }

        return false;
    }

    static
    int run ( int argc, const char *argv[] )
    {
        FastRQSettings settings;

        if ( parse_cmdline ( argc, argv, settings ) )
            return 0;

        settings.validate ();

        FastRQOperation *operation = make_operation ( settings );
        if ( operation == NULL )
            throw "Failed to make operation";

        // NULL OK
        FastRQFilter *filters = make_filters ( settings );

        // NULL OK
        FastRQFormatter *formatter = make_formatter ( settings );

        for ( uint32_t i = 1; i <= settings . num_accessions; ++ i )
        {
            const char *accession = argv [ i ];
            ReadCollection rCollection = ncbi :: NGS :: openReadCollection ( accession );

            operation -> run ( rCollection, filters, formatter );
        }

        delete ( operation );
        delete ( filters );
        delete ( formatter );

        return 0;
    }
}

extern "C"
{
    const char UsageDefaultName[] = "fastrq-dump";

    rc_t CC UsageSummary (const char * progname)
    {   // this is not used at this point, see handle_help()
        return 0;
    }

    rc_t CC Usage ( struct Args const * args )
    {   // this is not used at this point, see handle_help()
        return 0;
    }

    rc_t CC KMain ( int argc, char *argv [] )
    {
        try
        {
            return fastrq :: run ( argc, (const char**)argv );
        }
        catch ( ErrorMsg & x )
        {
            std :: cerr <<  x.toString () << '\n';
            return -1;
        }
        catch ( std :: exception & x )
        {
            std :: cerr <<  x.what () << '\n';
            return -1;
        }
        catch ( const char x [] )
        {
            std :: cerr <<  x << '\n';
            return -1;
        }
        catch ( ... )
        {
            std :: cerr <<  "unknown exception\n";
            return -1;
        }

        return 0;
    }
}

