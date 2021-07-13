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

#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <stdexcept>
#include "log.h"
#include "missing_cpp_features.h"

struct Args
{
    using FileList = std::list <std::string>;
    FileList files;
    std::string spot_filter_file;

    bool unaligned_only = false;
    bool hide_counts = false, compact = false;

    int optimization_ultrafast_skip_reader = 0;

    Args(int argc, char const *argv[])
    {
        int arg = 1;
        auto have_arg = [&]() -> bool { return arg < argc; };
        auto pop_arg = [&]() -> std::string {
            if (have_arg()) {
                auto const result = argv[arg];
                ++arg;
                return std::string(result);
            }
            fail("need more args");
        };
        std::string file;
        auto have_file = false;

        while (have_arg()) {
            auto const &arg = pop_arg();

            if (arg.empty()) {
USAGE_ERROR:
                std::string reason = "unexpected argument: \"" + arg + "\"";
                fail(reason.c_str());
            }
            if (arg.front() == '-') {
                if (arg == "-unaligned_only")
                    unaligned_only = true;
                else if (arg == "-spot_filter")
                    spot_filter_file = pop_arg();
                else if (arg == "-optimization_ultrafast_skip_reader")
                    optimization_ultrafast_skip_reader = std::stoi(pop_arg());
                else
                    goto USAGE_ERROR;
            }
            else if (have_file)
                goto USAGE_ERROR;
            else {
                file = arg;
                have_file = true;
            }
        }

        // exactly one should exist
        if (!have_file)
            fail("nothing to process");

        if (ends_with(file, ".list"))
            files = load_list(file);
        else
            files.push_back(file);

        if (files.empty())
            fail("loaded empty list of files to process");
    }

    static std::list<std::string> load_list(const std::string &filename)
    {
        std::ifstream f(filename);
        if (f.fail())
            throw std::runtime_error(std::string("cannot open list file ") + filename);

        std::list<std::string> items;

        while (!f.eof())
        {
            std::string s;
            f >> s;            
            if (f.fail())
                break;

            items.push_back(s);
        }

        return items;
    }

    static void fail [[noreturn]] (const char* reason = "invalid arguments")
    {
        LOG(reason);
        print_usage();
        exit(1);
    }

    static void print_usage()
    {
        std::cerr <<
            "need fasta/accession or file containing list of fasta/accession\n"
            "examples:\n"
            "    reads_from inputs.list # NB. extension is .list\n"
            "    reads_from chr1.fasta # NB. extension is .fasta, .fa, or .fna\n"
            "    reads_from SRR000001 # not one of the above\n"
            << std::endl;
    }
};
