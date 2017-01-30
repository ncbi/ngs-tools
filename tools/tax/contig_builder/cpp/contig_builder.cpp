#include "kmer_map.h"
#include "kmer_loader.h"
#include "seq_transform.h"
#include "contig_builder.h"
#include "begins.h"
#include "coverage.h"
#include "contig_cutter.h"
#include "reduce.h"
#include <iostream>
#include <chrono>
#include "contig.h"
#include "config.h"
#include <iomanip>
#include <ctime>

using namespace std;
using namespace std::chrono;

#include "print_time.h"

const string VERSION = "0.19";

typedef std::list<string> Strings;

void print_coverage(const std::vector<int> &cov)
{
	for (auto c : cov)
		cout << c << ", ";
	cout << endl;
}

template <class KmerMap>
std::list<string> cut_contigs(const KmerMap &kmers, const string &contig)
{
    const int COVERAGE_MEDIAN_FILTER_RADIUS = 5;
	auto coverage = reduce::median_filter(Coverage<KmerMap>(contig, kmers), COVERAGE_MEDIAN_FILTER_RADIUS);
//	print_coverage(coverage);
	return ContigCutter::cut(contig, coverage, kmers.kmer_len);
}


template <class KmerMap>
typename KmerMap::hash_t restore_orientation(typename KmerMap::hash_t hash, KmerMap &kmers)
{
	bool orig_complement = false, orig_reverse = false;
	kmers.get_original_compl_rev(hash, &orig_complement, &orig_reverse);
	return seq_transform<typename KmerMap::hash_t>::apply_transformation(hash, kmers.kmer_len, orig_reverse, orig_complement);
}

template <class MainKmerMap>
Strings build_contigs(MainKmerMap &kmers, KmerMap16 &kmers16, int MIN_SEQUENCE_LEN)
{
	auto before = high_resolution_clock::now();
	Begins<MainKmerMap> begins(kmers);
	cerr << "building begins time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;
	std::cerr << "mem usage before build_contigs (G) " << mem_usage()/1000000000 << std::endl;

	typename MainKmerMap::hash_t hash = 0;

	before = high_resolution_clock::now();

	std::chrono::duration<double> perf_begin, perf_contig;

	Strings contigs;

	auto before_loop = high_resolution_clock::now();
	while (true)
	{
		before = high_resolution_clock::now();

		if (!begins.next(&hash))
			break;

		hash = restore_orientation(hash, kmers);

		perf_begin += high_resolution_clock::now() - before;

		before = high_resolution_clock::now();
		string contig = ContigBuilder::get_next_contig(kmers, hash);

		if (contig.length() >= MIN_SEQUENCE_LEN)
		{
#if 0
			auto cutted_contigs = cut_contigs(kmers16, contig);
			for (auto &c : cutted_contigs)
			{
				if (c.size() < MIN_SEQUENCE_LEN)
					continue;

				contigs.push_back(c);
			}
#else
			contigs.push_back(contig);
#endif
		}
		perf_contig += high_resolution_clock::now() - before;
	}

	cerr << "building contigs time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before_loop ).count() << endl;
	cerr << "perf begin   (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( perf_begin ).count() << endl;
	cerr << "perf contig  (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( perf_contig ).count() << endl;

	return contigs;
}

template <class KmerMap>
double percent_of_run(double coverage_sum, const KmerMap &kmers)
{
	auto w = kmers.total_weight();
	return ( w > 0 ? coverage_sum/w : 0 ) * 100.0;
}

template <class KmerMap>
double print_seqs(const list<string> &seqs, const KmerMap &kmers, const string &desc)
{
    // todo: parallel for ?
	Contigs contigs;
	for (auto &seq : seqs)
	{
		auto coverage = Coverage<KmerMap>(seq, kmers);
		double coverage_sum = 0; //std::accumulate(coverage.begin(), coverage.end(), 0);
		for (auto c : coverage)
			coverage_sum += c;

		contigs.push_back(Contig(seq, percent_of_run(coverage_sum, kmers), coverage.empty() ? 0 : coverage_sum/coverage.size() ));
	}

	contigs.sort();

	int index = 0;
	double sum_percent = 0;
	for (auto &c :contigs)
	{
		cout << ">" << index << desc << c.data_percent << "%_cov_" << int(0.5 + c.average_coverage) << "_len_" << c.seq.size() << endl;
		cout << c.seq << endl;
		sum_percent += c.data_percent;
		index ++;
	}

	return sum_percent;
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);
	ngs::String acc = config.accession;
	cerr << "assembler version " << VERSION << endl;
	cerr << "accession: " << acc << endl;
	cerr << "-unaligned_only: " << config.unaligned_only << endl;
	cerr << "-min_contig_len: " << config.min_contig_len << endl;
	cerr << "-filter_file: " << config.filter_file << endl;
    cerr << "-exclude_filter: " << config.exclude_filter << endl;

	print_current_time();
	auto before = high_resolution_clock::now();

	KmerMap64 kmers64;
	KmerMap16 kmers16;

	KmerLoader loader(kmers64, kmers16, config.unaligned_only, config.filter_file, config.exclude_filter);
	loader.load(acc);

    Strings contigs_seqs = build_contigs(kmers64, kmers16, config.min_contig_len);
    double contig_percent = print_seqs(contigs_seqs, kmers64, "_");

	cerr << "reported contigs % " << contig_percent << endl;
	cerr << "reported contigs count " << contigs_seqs.size() << endl;
//	cerr << "reported sum % " << contig_percent + cont_percent << endl;
	cerr << "total time (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count() << endl;

//		exit(0); // dont want to wait for KmerMaps destructors
    return 0;
}
