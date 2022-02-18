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

#ifndef ALIGNS_TO_JOB_H_INCLUDED
#define ALIGNS_TO_JOB_H_INCLUDED

#include "config_align_to.h"
#include <time.h>
#include <thread>
#include "log.h"
#include "reader.h"
#include "fasta_reader.h"
#include "io.h"

struct Job
{
	virtual void run(const std::string &contig_filename, IO::Writer &writer, const Config &config) = 0;
//    virtual void match_and_print_chunk(const std::vector<Reader::Fragment> &chunk, IO::Writer &writer);

	template <class Matcher, class Printer, class MatchId>
    static void match_and_print(const std::vector<Reader::Fragment> &chunk, Printer &print, Matcher &matcher)
    {
        std::vector<MatchId> matched_ids;
        matched_ids.reserve(chunk.size()); // todo: tune

        for (size_t seq_id = 0; seq_id < chunk.size(); ++seq_id) 
        {
            // auto const &spotid = chunk[seq_id].spotid;
            auto const &bases = chunk[seq_id].bases;
            if (auto const m = matcher(bases)) {
                matched_ids.emplace_back((int)seq_id, m);
            }
        }

        #pragma omp critical (output)
        {
            print(chunk, matched_ids);
        }
    }     

    template <class MatchAndPrint>
	static void run_for_matcher(const std::string &contig_filename, const std::string &spot_filter_file, bool unaligned_only, int ultrafast_skip_reader, MatchAndPrint &&match_and_print)
	{
		Progress progress;
        Reader::Params params;
        params.filter_file = spot_filter_file;
        params.ultrafast_skip_reader = ultrafast_skip_reader;
        params.unaligned_only = unaligned_only;
        auto reader = Reader::create(contig_filename, params);

        #pragma omp parallel
        {
            std::vector<Reader::Fragment> chunk;
            bool done = false;
            while (!done) {
                #pragma omp critical (read)
                {
                    done = !reader->read_many(chunk);
                    progress.report(reader->progress());
                }

                match_and_print(chunk);
            }

        }

        progress.report(1, true); // always report 100%, needed by pipeline for proper progress report

        Reader::SourceStats total_stats;
        if (unaligned_only) {
            auto unaligned_stats = reader->stats();
            LOG("unaligned spot count: " << unaligned_stats.spot_count);
            LOG("unaligned read count: " << unaligned_stats.read_count);
        
            Reader::Params total_params;
            total_params.thread_count = 0;
            if (FastaReader::is_fasta(contig_filename)) {
                total_stats = unaligned_stats;
            } else {
                total_stats = Reader::create(contig_filename, total_params)->stats();
            }
        } else {
            total_stats = reader->stats();
        }
        
        LOG("total spot count: " << total_stats.spot_count);
        LOG("total read count: " << total_stats.read_count);
	}

	virtual size_t db_kmers() const { return 0; }
    virtual ~Job() {}

private:
	struct Progress
	{
        time_t last_timestamp;
		int last_reported;
		Progress() : last_timestamp(0), last_reported(-1) {}

		void report(float progress, bool force = false)
		{
            assert(progress >= 0 && progress <= 1);
            int percent = 100 * progress;
            if (percent != last_reported) {
                auto timestamp = time(nullptr);
                if (force || (timestamp > last_timestamp + 5)) {
                    LOG(percent << "% processed");
                    last_reported = percent;
                    last_timestamp = timestamp;
                }
            }
		}
	};

};

#endif
