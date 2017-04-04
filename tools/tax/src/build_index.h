#include <iostream>
#include <chrono>
#include <thread>
#include <list>
#include <iomanip>
#include <algorithm>
#include <map>
#include "config_build_index.h"
#include <omp.h>
#include <math.h>
#include <set>
#include "p_string.h"
#include "seq_cleaner.h"
#include "file_list_loader.h"
#include "kmers.h"
#include "fasta.h"
#include "hash.h"
#include <functional>
#include "seq_transform.h"

using namespace std;
using namespace std::chrono;

//#include "time.h"
void print_current_time() // todo: move out
{
	auto t = std::time(nullptr);
	auto timeinfo = std::localtime(&t);
	const int BUFFER_SIZE = 256;
	char buffer[BUFFER_SIZE];
	strftime(buffer, BUFFER_SIZE, "%m/%d/%Y %H:%M:%S", timeinfo);
	cerr << "time is " << buffer << endl;
}


//const int KMER_LEN = 32;

hash_t kmer_from(const char *s, int from, int len)
{
	hash_t kmer = Hash<hash_t>::hash_of(s + from, len);
	return seq_transform<hash_t>::min_hash_variant(kmer, len);
}

uint64_t fnv1_hash (void *key, int n_bytes)
{
    unsigned char *p = (unsigned char *)key;
    uint64_t h = 14695981039346656037UL;
    
    for (int i = 0; i < n_bytes; i++)
        h = (h * 1099511628211) ^ p[i];

    return h;
}

hash_t hash_of(hash_t hash)
{
	return fnv1_hash(&hash, sizeof(hash));
}
#define STOP_WHEN_FOUND_KMER 0

void process_window(Kmers &kmers, const char *s, int len, tax_id_t tax_id, int kmer_len)
{
//	int kmer_len = KMER_LEN; //kmers.kmer_len;
	if (len < kmer_len)
		return;

	const int THREADS = 16;
	struct ThreadFinding
	{
		size_t min_hash;
		int has_kmer_pos, min_hash_pos;
		ThreadFinding() : min_hash(std::numeric_limits<size_t>::max()), has_kmer_pos(-1), min_hash_pos(-1){}
	};

	array<ThreadFinding, THREADS> thread_findings;
	bool has_kmer_found = false;

	#pragma omp parallel num_threads(THREADS)
	for (int i = omp_get_thread_num(); !has_kmer_found && i <= len - kmer_len; i+=THREADS)
	{
		hash_t kmer = kmer_from(s, i, kmer_len);

		auto thread_id = omp_get_thread_num();

#if STOP_WHEN_FOUND_KMER
		if (kmers.has_kmer(kmer))
		{
			thread_findings[thread_id].has_kmer_pos = i;
			has_kmer_found = true;
			break;
		}
#endif

		std::size_t h = hash_of(kmer); // todo: can be optimized
		if (h < thread_findings[thread_id].min_hash)
		{
			thread_findings[thread_id].min_hash = h;
			thread_findings[thread_id].min_hash_pos = i;
		}
	}

	int chosen_kmer_pos = -1;

#if STOP_WHEN_FOUND_KMER
	if (has_kmer_found)
	{
		for (int i=0; i<THREADS; i++)
			if (thread_findings[i].has_kmer_pos >= 0)
			{
				chosen_kmer_pos = thread_findings[i].has_kmer_pos;
				break;
			}

		if (chosen_kmer_pos < 0)
			throw "has_kmer_found but cannot find for which thread";
	}
	else
#endif

	{
		size_t min_hash = thread_findings[0].min_hash;
		int min_hash_pos = thread_findings[0].min_hash_pos;

		for (int i=1; i<THREADS; i++)
			if (thread_findings[i].min_hash < min_hash)
			{
				min_hash = thread_findings[i].min_hash;
				min_hash_pos = thread_findings[i].min_hash_pos;
			}

		chosen_kmer_pos = min_hash_pos;

		if (chosen_kmer_pos < 0)
			throw "cannot find min hash";
	}

	kmers.add_kmer(kmer_from(s, chosen_kmer_pos, kmer_len), tax_id);
}

#define OVERLAPPING_WINDOWS 1

void process_clean_string(Kmers &kmers, p_string p_str, int window_size, tax_id_t tax_id, int kmer_len)
{
	for (int start = 0; start <= p_str.len - window_size; start += window_size)
    {
#if OVERLAPPING_WINDOWS
        int from = std::max(0, start - (kmer_len - 1));
        int to = std::min(start + window_size, p_str.len);
		process_window(kmers, p_str.s + from, to - from, tax_id, kmer_len);
#else
		process_window(kmers, p_str.s + start, window_size, tax_id, kmer_len);
#endif
    }
}

int calculate_window_size_(size_t filesize, bool eukariota, bool virus)
{
	size_t K = 1000;
	size_t M = 1000*K;
	size_t G = 1000*M;

	if (filesize < 2087819 || virus)
		return 200;

	if (filesize < 2932668)
		return 1000;

	if (filesize < 4774843 || !eukariota)
		return 2000;

	if (filesize < 500*M)
		return 4000;

	return 8000;
//	return (filesize/(4*G)) * 8000;
}

int calculate_window_size(size_t filesize, bool eukariota, bool virus, int window_divider)
{
	const int MIN_WINDOW_SIZE = 64;
	return std::max(MIN_WINDOW_SIZE, calculate_window_size_(filesize, eukariota, virus) / window_divider);
}
// todo: switch to seq_loader
struct ReadySeq
{
	string seq;
	SeqCleaner::p_strings clean_strings;
};

void swap(ReadySeq &a, ReadySeq &b)
{
	swap(a.seq, b.seq);
	swap(a.clean_strings, b.clean_strings);
}

void load_sequence(Fasta *_fasta, ReadySeq *_seq)
{
	Fasta &fasta = *_fasta;
	ReadySeq &seq = *_seq;

	seq.seq.clear(); // for better performance clear instad of constructor
	seq.clean_strings.clear();

	if (!fasta.get_next_sequence(seq.seq))
		return;

	SeqCleaner cleaner(seq.seq);
	seq.clean_strings = move(cleaner.clean_strings);
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

string str_kmer(hash_t kmer, int kmer_len)
{
	return Hash<hash_t>::str_from_hash(kmer, kmer_len);
}

#if 1
void print_kmers(const Kmers &kmers, int kmer_len)
{
	for (auto &kmer : kmers.storage)
	{
		auto tax_id = kmer.second;
		if (tax_id != TaxIdTree::ROOT)
			cout << str_kmer(kmer.first, kmer_len) << '\t' << tax_id << endl;
	}
}

#else
void print_kmers(const Kmers &kmers)
{
	for (auto &kmer : kmers.storage)
	{
		cout << str_kmer(kmer.first);
		for (auto &id : kmer.second)
			cout << '\t' << id;
		cout << endl;
	}
}
#endif

tax_id_t tax_id_from(const string &filename)
{
	auto to = filename.find_last_of('.');
	auto from = filename.find_last_of('/') + 1;
	return stoi(filename.substr(from, to - from));
}
