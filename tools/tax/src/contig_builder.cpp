#include "kmer_map.h"
#include "kmer_loader.h"
#include "seq_transform.h"
#include "contig_builder.h"
#include "begins.h"
#include "coverage.h"
#include <iostream>
#include <chrono>
#include "contig.h"
#include "config_contig_builder.h"
#include <iomanip>
#include <ctime>

using namespace std;
using namespace std::chrono;

const string VERSION = "0.19";

typedef std::list<string> Strings;

void print_coverage(const std::vector<int> &cov)
{
	for (auto c : cov)
		cout << c << ", ";
	cout << endl;
}

template <class KmerMap>
typename KmerMap::hash_t restore_orientation(typename KmerMap::hash_t hash, KmerMap &kmers) // todo: remove ?
{
	bool orig_complement = false, orig_reverse = false;
	kmers.get_original_compl_rev(hash, &orig_complement, &orig_reverse);
	return seq_transform<typename KmerMap::hash_t>::apply_transformation(hash, kmers.kmer_len, orig_reverse, orig_complement);
}

template <class MainKmerMap>
Strings build_contigs(MainKmerMap &kmers, int MIN_SEQUENCE_LEN)
{
	auto before = high_resolution_clock::now();
	Begins<MainKmerMap> begins(kmers);
	LOG("building begins time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before ).count());
	LOG("mem usage before build_contigs (G) " << mem_usage()/1000000000);

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
			contigs.push_back(contig);

		perf_contig += high_resolution_clock::now() - before;
	}

	LOG("building contigs time is (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( high_resolution_clock::now() - before_loop ).count());
	LOG("perf begin   (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( perf_begin ).count());
	LOG("perf contig  (ms) " << std::chrono::duration_cast<std::chrono::milliseconds>( perf_contig ).count());

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
	LOG("assembler version " << VERSION);
	LOG("accession: " << acc);
	LOG("-unaligned_only: " << config.unaligned_only);
	LOG("-min_contig_len: " << config.min_contig_len);
	LOG("-filter_file: " << config.filter_file);
    LOG("-exclude_filter: " << config.exclude_filter);

	auto before = high_resolution_clock::now();

	KmerMap32 kmers;

	KmerLoader loader(kmers, config.unaligned_only, config.filter_file, config.exclude_filter);
	loader.load(acc);

    Strings contigs_seqs = build_contigs(kmers, config.min_contig_len);
    double contig_percent = print_seqs(contigs_seqs, kmers, "_");

	LOG("reported contigs % " << contig_percent);
	LOG("reported contigs count " << contigs_seqs.size());
//	LOG("reported sum % " << contig_percent + cont_percent);
	LOG("total time (s) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count());

    return 0;
}
