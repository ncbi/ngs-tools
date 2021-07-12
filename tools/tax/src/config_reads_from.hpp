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

struct Config
{
    std::list <std::string> contig_files;
    std::string spot_filter_file;

    bool unaligned_only = false;
    bool hide_counts = false, compact = false;

    int optimization_ultrafast_skip_reader = 0;

    Config(int argc, char const *argv[])
    {
        std::string contig_file;

        std::list<std::string> args;
        for (int i = 1; i < argc; ++i)
            args.push_back(std::string(argv[i]));

        while (!args.empty()) {
            auto arg = pop_arg(args);

            if (arg == "-unaligned_only")
                unaligned_only = true;
            else if (arg == "-spot_filter")
                spot_filter_file = pop_arg(args);
            else if (arg == "-optimization_ultrafast_skip_reader")
                optimization_ultrafast_skip_reader = std::stoi(pop_arg(args));
            else if (arg.empty() || arg[0] == '-' || !contig_file.empty())
            {
                std::string reason = "unexpected argument: " + arg;
                fail(reason.c_str());
            }
            else
                contig_file = arg;
        }

        // exactly one should exist
        if (contig_file.empty()) // == contig_files.empty())
            fail("please provide either contig file or list");

        if (ends_with(contig_file, ".list"))
            contig_files = load_list(contig_file);
        else
            contig_files.push_back(contig_file);

        if (contig_files.empty())
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

    static void fail(const char* reason = "invalid arguments")
    {
        print_usage();
        LOG(reason);
        exit(1);
    }

    static void print_usage()
    {
        std::cerr << "need <contig fasta, accession or .list file of fasta/accessions>" << std::endl;
    }

private:

    static std::string pop_arg(std::list<std::string>& args)
    {
        if (args.empty())
            fail("need more args");

        std::string arg = args.front();
        args.pop_front();
        return arg;
    }
};
