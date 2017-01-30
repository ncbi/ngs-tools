#include <iostream>
#include <chrono>
#include <thread>
#include <list>
#include <iomanip>
#include <algorithm>
#include <map>
#include "config_compress.h"
#include "text_loader_mt.h"
#include "kmers_loading.h"
#include "translation.h"
#include <omp.h>
#include <math.h>
#include "stringn.h"
#include "kmers.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.1";

#define TRANSLATE_SEQ 0

#if TRANSLATE_SEQ
const int NUCLS_PER_ACID = 3;
#else
const int NUCLS_PER_ACID = 1;
#endif

void print_current_time()
{
	auto t = std::time(nullptr);
	auto timeinfo = std::localtime(&t);
	const int BUFFER_SIZE = 256;
	char buffer[BUFFER_SIZE];
	strftime(buffer, BUFFER_SIZE, "%m/%d/%Y %H:%M:%S", timeinfo);
	cerr << "time is " << buffer << endl;
}

struct Sequence
{
	string nucleotides;
	string frame0, frame1, frame2;

	bool empty() const { return nucleotides.empty(); }
	size_t size() const { return nucleotides.size(); }
};

#if TRANSLATE_SEQ
Kmer get_translated_kmer(const Kmers &kmers, const string &seq, int kmer_len, int from)
{
	if (kmer_len == 1)
		return NON_ACID_1_MER;

	if (kmer_len == 2)
		return NON_ACID_2_MER;

	if (kmer_len %3 != 0) // todo: remove
	{
		throw "kmer_len %3 != 0";
		return NON_ACID_X_MER;
	}

	const char *translated;
	get_translation(seq, from, &translated);
	return kmers.get_kmer(translated, kmer_len/3);
}

#else

void print_kmer(const Kmer &kmer);

Kmer get_translated_kmer(const Kmers &kmers, const Sequence &seq, int kmer_len, int from)
{
	//cout << "|||";
	//print_kmer(kmers.get_kmer("-", 1));
	//cout << "|||" << endl;
	auto kmer = kmers.get_kmer(&seq.nucleotides[from], kmer_len);
	if (kmer_len == 1 && kmer.empty())
	{
		cout << "-----------" << endl;
		cout << seq.size() << " " << from << " " << seq.nucleotides << endl;
		cerr << "not found seq " << seq.nucleotides[from];
		//print_string(&seq.nucleotides[from], kmer_len);
		cerr << "| for kmer len " << kmer_len << endl;
	}

	return kmer;
}

#endif


struct Compression
{
	Compression *left;
	Kmer kmer;
	double weight;

	Compression(Compression *left = nullptr, const Kmer &kmer = MISSING_KMER) : left(left), kmer(kmer)
	{
		weight = left ? left->weight : 0;
		weight += kmer.bits;
	}

	bool empty() const
	{
		return kmer.empty();
	}
};

//Compression combine(Compression *left, Kmer &kmer)
//{
//	return Compression(left, kmer);
//}
//

template <class It, class End>
void safe_increment(It &it, const End &end)
{
	if (it!=end)
		it++;
}

typedef list<Compression> CompressedResult;
double print_compression(const Compression &compression);

void compress_seq(const Sequence &seq, const Kmers &kmers, CompressedResult &compressions)
{
	compressions.clear();

	if (seq.empty())
		return;
		//throw "compression string is empty";

//	cout << "seq size " << seq.size() << endl;
	for (int to_inclusive = 0; to_inclusive < int(seq.size()); to_inclusive++)
	{
		Compression best_compression;
		auto left_part = compressions.rbegin();

		for (int kmer_len = 1; kmer_len <= kmers.max_kmer_len()*NUCLS_PER_ACID; kmer_len++, safe_increment(left_part, compressions.rend()))
		{
			int left_part_to_inclusive = to_inclusive - kmer_len;
			if (left_part_to_inclusive < -1)
				break;

#if TRANSLATE_SEQ
			if (kmer_len > NUCLS_PER_ACID && ((kmer_len % NUCLS_PER_ACID) != 0) ) // todo: optimize - increment by 3 - - do kmer_len += 2
				continue;
#endif

			auto kmer = get_translated_kmer(kmers, seq, kmer_len, left_part_to_inclusive + 1);
			if (kmer.empty())
				continue;

//			auto left_part = left_part_to_inclusive >=0 ? &compression[left_part_to_inclusive] : nullptr;
			auto compressed = Compression(left_part == compressions.rend() ? nullptr : &*left_part, kmer);
//			print_compression(compressed);

			if (!compressed.empty() && (best_compression.empty() || compressed.weight < best_compression.weight))
				best_compression = compressed;
		}

//		print_compression(best_compression);
//		print_kmer(best_compression.kmer);
//		cout << endl;
		if (best_compression.empty())
		{
			cout << "seq " << seq.nucleotides << endl;
			cerr << "bad letter: " << to_inclusive << " " << seq.nucleotides[to_inclusive] << endl;
			throw "best_compression is empty";
		}

		compressions.push_back(best_compression);
	}

//	print_compression(*compressions.rbegin());
}

struct CompressedResults
{
	std::mutex result_mutex;

	typedef std::pair<CompressedResult, int> ResultWithId;

	list<ResultWithId> results;

	CompressedResult& reserve(int seq_id) // todo: keep order
	{
//		return;
		std::lock_guard<std::mutex> lock(result_mutex);

//		if (results.size() <= seq_id)
//			results.resize(seq_id + 1);
		results.push_back(std::make_pair(CompressedResult(), seq_id));
		return results.rbegin()->first;

//		results[seq_id] = std::move(r);
	}
};

bool operator < (const CompressedResults::ResultWithId &a, const CompressedResults::ResultWithId &b)
{
	return a.second < b.second;
}

string stringn(const char *s, int len)
{
	if (!s)
		return "null";

	return string(s, s + len);
}

void print_kmer(const Kmer &kmer)
{
	if (kmer.empty())
		cout << "[]";
	else
		cout << stringn(kmer.seq, kmer.len) << " "; // << std::setprecision(2) << kmer.bits; // << ")";
}

double print_compression(const Compression &compression)
{
	const Compression *c = &compression;
	list<Kmer> kmers;
//	for (auto &c : compression)
	for(; c; c = c->left)
		kmers.push_back(c->kmer);

	kmers.reverse();

	double weight = 0;
	for (auto &kmer : kmers)
	{
		print_kmer(kmer);
		weight += kmer.bits;
	}

	cout << endl;
	return weight;
//	cerr << " = " << std::setprecision(10) << weight << endl;
}

void print_compressions(const CompressedResults &compressions)
{
	double weight = 0;
	for (auto &result : compressions.results)
		weight += print_compression(*result.first.rbegin());

	cerr << "total weight " << std::setprecision(10) << weight << endl;
}

void compress(TextLoaderMTNoStore &seq_loader, CompressedResults &compressions, const Kmers &kmers)
{
	Sequence seq;
	size_t seq_id = 0;
	while (seq_loader.load_next_sequence(seq.nucleotides, seq_id))
	{
//		cout << "seq " << seq.size() << " " << seq.nucleotides << endl;
#if TRANSLATE_SEQ
		translation::translate(seq.nucleotides, 0, seq.frame0);
		translation::translate(seq.nucleotides, 1, seq.frame1);
		translation::translate(seq.nucleotides, 2, seq.frame2);
#endif
		CompressedResult &r = compressions.reserve(seq_id);
		compress_seq(seq, kmers, r);
		cerr << ".";

//		compressions.save_result(, seq_id);
	}
}

struct PStringCounted
{
	Kmers::p_string p_str;
	size_t count;

	PStringCounted(Kmers::p_string p_str, size_t count) : p_str(p_str), count(count) {}
	bool operator < (const PStringCounted &ps) const { return count > ps.count; }
};

void dump_frequency(size_t kmer_len, const vector<PStringCounted> &len_freq)
{
	string filename = string("stat") + std::to_string(kmer_len) + string(".txt");
	ofstream f(filename);
	f.flush();
	for (auto &ps : len_freq)
		f << string(ps.p_str.s, ps.p_str.s + ps.p_str.len) << '|' << ps.count << endl;
}

void dump_frequency(const CompressedResults &compressions, const Kmers &kmers, int min_coverage)
{
	vector< map<Kmers::p_string, size_t> > freq(kmers.max_kmer_len() + 1);

	for (auto &result : compressions.results)
	{
		const Compression *c = &*result.first.rbegin();
		for(; c; c = c->left)
			freq[c->kmer.len][Kmers::p_string(c->kmer.seq, c->kmer.len)] ++;
	}

	for (size_t kmer_len = 1; kmer_len <= kmers.max_kmer_len(); kmer_len++)
	{
		vector<PStringCounted> len_freq;

		auto &storage = kmers.kmers[kmer_len];
		for (auto &it : storage)
		{
			auto p_str = it.first;
			size_t count = freq[kmer_len][p_str]; // todo: optimize - do not create key if do not exists
			if (kmer_len == 1 || count >= min_coverage)
				len_freq.push_back(PStringCounted(p_str, count));
		}

		sort(len_freq.begin(), len_freq.end());
		dump_frequency(kmer_len, len_freq);
	}
}


const int THREADS = 16;

int main(int argc, char const *argv[])
{
    try
    {
		ConfigCompress config(argc, argv);
		cerr << "compress version " << VERSION << endl;

		print_current_time();
		auto before = high_resolution_clock::now();

		TextLoaderMTNoStore seq_loader(config.reference);

		Kmers kmers;
		KmerLoading kmer_loading(kmers, config.kmer_path_mask);
		cerr << "kmers loaded at (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;

		CompressedResults compressions;

		#pragma omp parallel num_threads(THREADS)
		{
			compress(seq_loader, compressions, kmers);
		}

		cerr << "compression time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;

		compressions.results.sort();

		cerr << "results sorted at (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;

		print_compressions(compressions);
		cerr << "dumping frequencies" << endl;
		dump_frequency(compressions, kmers, config.min_dump_coverage);

		cerr << "total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;

		exit(0); // dont want to wait for KmerMaps destructors
        return 0;
    }
    catch ( exception & x )
    {
        cerr << x.what() << endl;
//		cerr << "exit 3" << endl;
		return 3;
    }
    catch ( string & x )
    {
        cerr << x << endl;
//		cerr << "exit 4" << endl;
		return 4;
    }
    catch ( const char * x )
    {
        cerr << x << endl;
//		cerr << "exit 5" << endl;
		return 5;
    }
    catch ( ... )
    {
        cerr << "unknown exception" << endl;
//		cerr << "exit 6" << endl;
		return 6;
    }
}
