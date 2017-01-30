#include <iostream>
#include <chrono>
#include <thread>
#include <list>
#include <iomanip>
#include <algorithm>
#include <map>
#include "config_min_coverage.h"
#include "text_loader_mt.h"
#include <omp.h>
#include <math.h>
#include <set>
#include "p_string.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.1";
#include "time.h"

typedef p_string Kmer;

struct Sequence
{
	string seq;
	vector<bool> used;

	bool used_mark(size_t pos) const
	{
		if (pos >= used.size()) // todo: removes
			throw "used_mark pos >= used.size()";

		return used[pos];
	}

	void mark_used(size_t pos)
	{
		if (pos >= used.size()) // todo: removes
			throw "used_mark pos >= used.size()";

		used[pos] = true;
	}

	size_t length() const 
	{
		return seq.length();
	}
};

struct KmerPositions
{
	struct Position
	{
		Sequence *seq;
		size_t pos; // can use int, actually
		Position(Sequence *seq, size_t pos) : seq(seq), pos(pos){}

		bool operator < (const Position &p) const
		{
			return pos == p.pos ? seq < p.seq : pos < p.pos;
		}
	};
	
	typedef set<Position> Positions;

	std::map<p_string, Positions> storage;

	void add(Sequence *seq, size_t pos, int kmer_len)
	{
		auto key = p_string(&seq->seq[pos], kmer_len);
//		print_string(key);
//		cout << endl;
		auto &at = storage[key];
		at.insert(Position(seq, pos));
//		cout << "added seq " << seq << endl;
//		cout << seq->seq << endl;
	}

	size_t coverage_of(const Kmer &kmer) const
	{
		auto it = storage.find(kmer); // p_string(kmer.s, kmer.len));
		if (it == storage.end())
			return 0;

		return it->second.size();
	}

	void remove(const Kmer &kmer)
	{
		auto it = storage.find(kmer); // p_string(kmer.s, kmer.len));
		if (it == storage.end())
			throw "removing non existing element from set";

		storage.erase(it);
	}

	Positions &get_positions(const Kmer &kmer)
	{
		auto it = storage.find(kmer); //p_string(kmer.s, kmer.len));
		if (it == storage.end())
			throw "cannot find kmer in storage";

		return it->second;
	}
};

size_t find_first_used_kmer_mark(size_t from, int len, Sequence *seq)
{
	size_t to = min(seq->length(), from + len);
//	cout << "looking for kmer mark from " << from << " to " << to << " ";
	for (size_t pos = from; pos < to; pos++)
		if (seq->used_mark(pos))
		{
//			cout << "found at " << pos << endl;
			return pos;
		}

//	cout << "not found" << endl;
	return 0;
}

size_t find_last_used_kmer_mark(size_t from, int len, Sequence *seq)
{
	auto pos = find_first_used_kmer_mark(from, len, seq);
	if (!pos)
		return 0;

	pos++;
	for (; pos < seq->length(); pos++)
		if (!seq->used_mark(pos))
			return pos - 1;

	return pos - 1;
}

size_t find_next_free_pos(size_t from, int kmer_len, int lookup_len, Sequence *seq, const KmerPositions &kmers)
{
//	const int LOOKUP_LEN = kmer_len;
	while (true)
	{
		size_t to = std::min(from + lookup_len, seq->length() - kmer_len);
		if (from >= to)
			return 0;

//		cout << "find free pos from " << from << " to " << to << endl;

		auto used_kmer_mark = find_last_used_kmer_mark(from, kmer_len + lookup_len, seq);
		if (used_kmer_mark != 0)
		{
			from = used_kmer_mark + 1;
			continue;
		}

		int max_cov = 0;
		size_t best_pos = 0;
		for (size_t pos = from; pos < to; pos++)
		{
			auto current_cov = kmers.coverage_of(Kmer(&seq->seq[pos], kmer_len));
			if (current_cov >= max_cov)
			{
				max_cov = current_cov;
				best_pos = pos;
			}
		}

		return best_pos;
	}
}

void mark_as_used(Sequence *seq, size_t from, int kmer_len)
{
	size_t to = min(seq->length(), from + kmer_len);
	for (size_t pos = from; pos < to; pos++)
		seq->mark_used(pos);
}

void mark_kmer(const Kmer &kmer, size_t pos, int kmer_len, Sequence *seq, KmerPositions &kmers)
{
	if (pos + kmer_len > seq->length())
		throw "choose kmer:: pos + kmer_len > seq->length()";

	if (find_first_used_kmer_mark(pos, kmer_len, seq)) // todo: remove
		throw "reused sequence";

//	auto kmer = Kmer(seq->seq, pos, kmer_len);
	auto &poss = kmers.get_positions(kmer);
	for (auto &pos : poss)
	{
//		cout << "seq " << seq << " pos seq " << pos.seq << endl;
//		mark_as_used(seq, pos.pos, kmer_len);
		mark_as_used(pos.seq, pos.pos, kmer_len);
//		cout << "seq is " << seq->seq << endl;
//		cout << "pos seq is " << pos.seq->seq << endl;
	}

	kmers.remove(kmer);
}

struct KmerBuilder
{
	std::list<Sequence> sequences;
	KmerBuilder(TextLoaderSTNoStore &seq_loader, KmerPositions &kmers, int kmer_len)
	{
		sequences.push_back(Sequence());
		while (seq_loader.load_next_sequence(sequences.rbegin()->seq))
		{
			Sequence *seq = &*sequences.rbegin();
//			cout << seq->seq << endl;
			seq->used.resize(seq->seq.size());
			fill_kmers(kmers, seq, kmer_len);
			sequences.push_back(Sequence());
			cerr << ".";
		}

		sequences.pop_back(); // todo: think
	}

	void fill_kmers(KmerPositions &kmers, Sequence *seq, int kmer_len)
	{
		for (size_t pos = 0; pos <= seq->length() - kmer_len; pos++)
			kmers.add(seq, pos, kmer_len);
	}
};

template <class OnKmerFound, class OnSeqCompressed>
void compress(KmerBuilder &kmer_builder, KmerPositions &kmers, int kmer_len, int lookup_len, OnKmerFound &&kmer_found, OnSeqCompressed &&seq_compressed)
{
//	while (Sequence *seq = kmer_builder.next_sequence())
	cerr << "sequences: " << kmer_builder.sequences.size() << endl;
	//for (auto &it : kmer_builder.sequences)
	for (auto it = kmer_builder.sequences.begin(); it!=kmer_builder.sequences.end(); it++)
	{
		Sequence *seq = &*it;
//		cout << "next seq " << seq->length() << endl;
//		cout << seq->seq << endl;

		size_t pos = 1; // 0 is invalid pos for everything
		for (; pos = find_next_free_pos(pos, kmer_len, lookup_len, seq, kmers); pos += kmer_len)
		{
			auto kmer = Kmer(&seq->seq[pos], kmer_len);
			kmer_found(seq, pos, kmer, kmers);
			mark_kmer(kmer, pos, kmer_len, seq, kmers);
		}

		seq_compressed(seq);
		//cerr << ".";
	}
}

void print_used(const vector<bool> &used)
{
	for (int i=0; i<int(used.size()); i++)
		if (used[i])
			cout << 'x';
		else
			cout << '-';

	cout << endl;
}

