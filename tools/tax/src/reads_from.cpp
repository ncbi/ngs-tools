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

#include "config_reads_from.hpp"

#include <iostream>
#include <chrono>
#include <thread>
#include <list>

const std::string VERSION = "0.1";

typedef uint64_t hash_t;

#include "reader.h"
#include "missing_cpp_features.h"

using namespace std;
using namespace std::chrono;

static int process(Args::FileList const &files, Reader::Params const &params)
{
    auto const before = high_resolution_clock::now();
    auto const multifile = files.size() > 1;

    for (auto &file : files)
    {
        LOG(file);

        try {
            auto reader = Reader::create(file, params);
            Reader::Fragment frag;

            while (reader->read(&frag)) {
                cout << '>' << frag.spotid << '\n' <<
                               frag.bases << '\n';
            }
        }
        catch (std::exception const &e)
        {
            LOG(e.what());
            if (!multifile)
                throw e;
        }
    }

    LOG("total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count());

    return 0;
}

static int process(Args const &args)
{
    Reader::Params params;

    params.filter_file = args.spot_filter_file;
    params.ultrafast_skip_reader = args.optimization_ultrafast_skip_reader;
    params.unaligned_only = args.unaligned_only;

    return process(args.files, params);
}

int main(int argc, char const *argv[])
{
#ifdef __GLIBCXX__
    std::set_terminate(__gnu_cxx::__verbose_terminate_handler);
#endif

    LOG("reads_from version " << VERSION);

    return process(Args(argc, argv));
}

