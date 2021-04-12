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

#include "aligns_to_job.h"
#include "hash.h"
#include "seq_transform.h"
#include <map>
#include "omp_adapter.h"
#include "missing_cpp_features.h"

struct ManyJobs : public Job
{
    static std::unique_ptr<Job> create_job(const std::string &db)
    {
        if (ends_with(db, ".db"))
            return make_unique<DBJob>(db);
        if (ends_with(db, ".dbs"))
            return make_unique<DBSBasicJob>(db);
        if (ends_with(db, ".dbss"))
            throw std::runtime_error("dbss is not supported for <many> option at this moment");
            //return make_unique<DBSSJob>(db);

        throw std::runtime_error(std::string("cannot identify database type of ") + db);
    }

    std::list< std::unique_ptr<Job> > jobs;

	ManyJobs(const std::string &many)
	{
        auto files = split(many, ',');
        if (files.empty())
            throw std::runtime_error("no databases to load");

        for (auto &file : files)
            jobs.push_back(create_job(file));
	}

	virtual void run(const std::string &filename, IO::Writer &writer, const Config &config) override
	{
        if (writer.out_f.size() != jobs.size())
            throw std::runtime_error("you need as many output files as jobs");

		Job::run_for_matcher(filename, config.spot_filter_file, config.unaligned_only, [&](const std::vector<Reader::Fragment> &chunk){ match_and_print_chunk(chunk, writer); } );
	}

    virtual void match_and_print_chunk(const std::vector<Reader::Fragment> &chunk, IO::Writer &_writer)
    {
        int stream_id = 0;
        for (auto &job : jobs)
        {
            IO::Writer writer(_writer, stream_id);
            job->match_and_print_chunk(chunk, writer);
            stream_id++;
        }
    }
};

