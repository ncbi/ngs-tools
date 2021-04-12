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

#include <time.h>
#include <thread>
#include "log.h"
#include "reader.h"
#include "fasta_reader.h"

struct BasicMatchId
{
	int seq_id;

	BasicMatchId(int seq_id, int matches) : seq_id(seq_id){}
	bool operator < (const BasicMatchId &b) const { return seq_id < b.seq_id; }
};

struct BasicPrinter
{
	std::ostream &out_f;
	BasicPrinter(std::ostream &out_f) : out_f(out_f){}

	void operator() (const std::vector<Reader::Fragment> &processing_sequences, const std::vector<BasicMatchId> &ids)
	{
		for (auto seq_id : ids)
			out_f << processing_sequences[seq_id.seq_id].spotid << std::endl;
	}
};

struct Job
{
	virtual void run(const std::string &contig_filename, std::ostream &out_f) = 0;

	template <class Matcher, class Printer, class MatchId = BasicMatchId>
	static void run(const std::string &contig_filename, Printer &print, Matcher &matcher, size_t min_sequence_len, const std::string &spot_filter_file, bool unaligned_only)
	{
		Progress progress;
        Reader::Params params;
        params.filter_file = spot_filter_file;
        params.split_non_atgc = true;
        params.unaligned_only = unaligned_only;
        auto reader = Reader::create(contig_filename, params);

        #pragma omp parallel
        {
            std::vector<MatchId> matched_ids;
            std::vector<Reader::Fragment> chunk;
            bool done = false;
            while (!done) {
                #pragma omp critical (read)
                {
                    done = !reader->read_many(chunk);
                    progress.report(reader->progress());
                }

                matched_ids.clear();
                for (size_t seq_id = 0; seq_id < chunk.size(); ++seq_id) {
                    auto& spotid = chunk[seq_id].spotid;
                    auto& bases = chunk[seq_id].bases;
                    //processed_spots.insert(spotid);
                    if (bases.size() >= min_sequence_len) {
                        if (auto m = matcher(bases)) {
                            //identified_spots.insert(spotid);
                            matched_ids.push_back(MatchId(seq_id, m));
                        }
                    }
                }

                #pragma omp critical (output)
                {
                    print(chunk, matched_ids);
                }
            }

        }

        progress.report(1, true); // always report 100%, needed by pipeline for proper progress report

        Reader::SourceStats total_stats;
        if (unaligned_only) {
            auto unaligned_stats = reader->stats();
            LOG("unaligned spot count: " << unaligned_stats.spot_count);
            LOG("unaligned read count: " << unaligned_stats.frag_count());
        
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
        LOG("total read count: " << total_stats.frag_count());
	}

	virtual size_t db_kmers() const { return 0;}

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
