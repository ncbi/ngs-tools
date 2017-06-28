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

#include <iostream>
#include <chrono>
#include <thread>
#include <array>
#include "omp_adapter.h"
#include "kmers.h"
#include "kmer_io.h"
#include "kmer_hash.h"
#include "ready_seq.h"
#include "filename_meta.h"
#include "config_build_index.h"
#include "file_list_loader.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.35";

size_t weight(size_t kmers_count)
{
	return kmers_count*16;
}

int calculate_window_size_(size_t filesize, bool eukaryota, bool virus)
{
	if (virus)
		return 200;

	if (!eukaryota)
		return 2000;

	return 8000;
}

int calculate_window_size(size_t filesize, bool eukaryota, bool virus, int window_divider, int min_window_size)
{
	return std::max(min_window_size, calculate_window_size_(filesize, eukaryota, virus) / window_divider);
}

void process_window(Kmers &kmers, const char *s, int len, tax_id_t tax_id, int kmer_len)
{
	if (len < kmer_len)
		return;

	const int THREADS = 16;
	struct ThreadFinding
	{
		KmerHash::hash_of_hash_t min_hash;
		int has_kmer_pos, min_hash_pos;
		ThreadFinding() : min_hash(std::numeric_limits<size_t>::max()), has_kmer_pos(-1), min_hash_pos(-1){}
	};

	array<ThreadFinding, THREADS> thread_findings;
	bool has_kmer_found = false;

	#pragma omp parallel num_threads(THREADS)
	for (int i = omp_get_thread_num(); !has_kmer_found && i <= len - kmer_len; i += omp_get_num_threads())
	{
		hash_t kmer = KmerIO::kmer_from(s, i, kmer_len);

		auto thread_id = omp_get_thread_num();

		auto h = KmerHash::hash_of(kmer); // todo: can be optimized

		if (h < thread_findings[thread_id].min_hash)
		{
			thread_findings[thread_id].min_hash = h;
			thread_findings[thread_id].min_hash_pos = i;
		}
	}

	int chosen_kmer_pos = -1;

	{
		auto min_hash = thread_findings[0].min_hash;
		int min_hash_pos = thread_findings[0].min_hash_pos;

		for (int i=1; i<THREADS; i++)
			if (thread_findings[i].min_hash < min_hash)
			{
				min_hash = thread_findings[i].min_hash;
				min_hash_pos = thread_findings[i].min_hash_pos;
			}

		chosen_kmer_pos = min_hash_pos;

		if (chosen_kmer_pos < 0)
			throw std::runtime_error("cannot find min hash");
	}

	kmers.add_kmer(KmerIO::kmer_from(s, chosen_kmer_pos, kmer_len), tax_id);
}

void process_clean_string(Kmers &kmers, p_string p_str, int window_size, tax_id_t tax_id, int kmer_len)
{
	for (int start = 0; start <= p_str.len - window_size; start += window_size)
    {
        int from = std::max(0, start - (kmer_len - 1));
        int to = std::min(start + window_size, p_str.len);
		process_window(kmers, p_str.s + from, to - from, tax_id, kmer_len);
    }
}

size_t add_kmers(Kmers &kmers, const string &filename, tax_id_t tax_id, int window_size, int kmer_len)
{
	Fasta fasta(filename);

	size_t seq_index = 0;
	size_t total_size = 0;
	const int DOT_INTERVAL = 128;

	ReadySeq loading_seq, processing_seq;
	load_sequence(&fasta, &processing_seq);

	while (!processing_seq.seq.empty())
	{
		std::thread loading_thread(load_sequence, &fasta, &loading_seq);

		total_size += processing_seq.seq.size();

		for (auto &clean_string : processing_seq.clean_strings)
			process_clean_string(kmers, clean_string, window_size, tax_id, kmer_len);

		seq_index++;
		if (seq_index % DOT_INTERVAL == 0)
			cerr << ".";

		loading_thread.join();
		swap(processing_seq, loading_seq); // todo: move?
	}

	if (seq_index >= DOT_INTERVAL)
		cerr << endl;

	return total_size;
}


int main(int argc, char const *argv[])
{
	LOG("build_index version " << VERSION);
	ConfigBuildIndex config(argc, argv);
	LOG("window divider: " << config.window_divider);
	LOG("kmer len: " << config.kmer_len);

	auto before = high_resolution_clock::now();

	FileListLoader file_list(config.file_list);

	TaxIdTree tax_id_tree;
	TaxIdTreeLoader::load_tax_id_tree(tax_id_tree, config.tax_parents_file);

	Kmers kmers(tax_id_tree);
	size_t total_size = 0;
	for (auto &file_list_element : file_list.files)
	{
		auto window_size = calculate_window_size(file_list_element.filesize, FilenameMeta::is_eukaryota(file_list_element.filename), FilenameMeta::is_virus(file_list_element.filename), config.window_divider, config.min_window_size);
		auto tax_id = FilenameMeta::tax_id_from(file_list_element.filename);
		LOG(file_list_element.filesize << "\t" << window_size << "\t" << tax_id << "\t" << file_list_element.filename);
		total_size += add_kmers(kmers, file_list_element.filename, tax_id, window_size, config.kmer_len);
		{
			auto seconds_past = std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count();
			if (seconds_past < 1)
				seconds_past = 1;

			size_t megs = total_size/1000000;
			LOG("processed size " << megs << "M = " << (total_size/1000)/seconds_past << "K/sec, kmers: " << kmers.storage.size()/1000 << "K, compression rate " << total_size/std::max(size_t(1), weight(kmers.storage.size())));
		}
	}

	KmerIO::print_kmers(kmers, config.kmer_len);

	LOG("total time (min) " << std::chrono::duration_cast<std::chrono::minutes>( high_resolution_clock::now() - before ).count());
}
