#include <iostream>
#include <chrono>
#include <thread>
#include <omp.h>
#include <algorithm>
#include "freq_db.h"
#include "freq_db_io.h"
#include "hash.h"
#include "p_string.h"
#include "seq_transform.h"
#include "seq_loader.h"
#include "file_list_loader.h"

using namespace std;
using namespace std::chrono;

hash_t get_kmer_from(const char *s, int from, int len)
{
//	return Hash<hash_t>::hash_of(s + from, len);
	return Hash<hash_t>::hash_of_amino(s + from, len);
//	return seq_transform<hash_t>::min_hash_variant2(kmer, len);
}

void get_window_from_to(int &from, int &to, int kmer_from, int line_len, int kmer_len, int window_r)
{
	from = std::max(0, kmer_from - window_r);
	to = std::min(line_len - kmer_len + 1, kmer_from + window_r + 1);
}

template <class Lambda>
void do_for_left_and_right_sides(p_string line, int kmer_from, int kmer_len, int window_r, Lambda &&lambda)
{
	int from, to;
	get_window_from_to(from, to, kmer_from, line.len, kmer_len, window_r);
	for (int pos = from; pos <= kmer_from - kmer_len; pos++)
        lambda(pos);
    
	for (int pos = kmer_from + kmer_len; pos < to; pos++)
        lambda(pos);
}

struct WindowBits
{
    struct CompressedAmino
    {
        int frame;
        double bits;
        bool intra;

        CompressedAmino(int frame = -1, double bits = 1000000000, bool intra = false) : frame(frame), bits(bits), intra(intra){}
    };

	map<int, CompressedAmino> windows;

	double sum_bits_with_frame_penalties() 
    {
        const double FRAME_SHIFT_PENALTY = 0;
        return sum_bits() + frame_shifts() * FRAME_SHIFT_PENALTY;
    }

    int sum_intra_windows() const
    {
        int sum = 0;
		for (auto &it : windows)
            if (it.second.intra)
                sum++;

        return sum;
    }

private:
	double sum_bits() const
	{
		double sum = 0;
		for (auto &it : windows)
			sum += it.second.bits;

		return sum;
	}

	int frame_shifts() 
	{
		int shifts = 0;
		for (auto &it : windows)
        {
            int window_number = it.first;
            if (window_number == 0)
                continue;

            if (windows.find(window_number - 1) == windows.end())
                throw "windows.find(window_number - 1) == windows.end()";

            int frame = it.second.frame;
            if (frame != windows[window_number - 1].frame)
                shifts++;
        }

		return shifts;
	}
};

double bits_for_freq(double freq)
{
	if (freq <= 0)
		throw "bits for freq: freq <=0 ";
	
	return log2(1.0/freq);
}

double intra_bits_for(p_string text, int pos, hash_t kmer, Frequences &freqs)
{
    const int INTRA_PENALTY = 0; // todo: tune
	if (freqs.freqs.find(kmer) == freqs.freqs.end())
		return 13; //log2(AMINOACIDS^3) // todo: make it right bits_for_freq(1.0/freqs.total) + INTRA_PENALTY; // todo + penalty ?

	return bits_for_freq(freqs.freqs[kmer]) + INTRA_PENALTY;
}

double combine_freqs_max(const vector<double> &found_freqs)
{
	if (found_freqs.empty())
		throw "found_freqs.empty()";

	double max_freq = 0;
	for (auto f : found_freqs)
	{
		max_freq = std::max(max_freq, f);
//		cout << "freq " << f << endl;
	}
		
//	cout << "combined freq is " << max_freq << endl;
	return max_freq;
}

double combine_freqs_mul(const vector<double> &found_freqs)
{
	if (found_freqs.empty())
		throw "found_freqs.empty()";

	double err = 1.0;
	for (auto f : found_freqs)
	{
		err *= (1.0 - f);
//		max_freq = std::max(max_freq, f);
//		cout << "freq " << f << endl;
	}
		
//	cout << "combined freq is " << 1.0 - err << endl;
	return 1.0 - err;
}

double combine_freqs(const vector<double> &found_freqs)
{
    return combine_freqs_mul(found_freqs);
//    return combine_freqs_max(found_freqs);
}

bool same_frame(int pos1, int pos2)
{
    const int FRAME_NUCLEOTIDES = 3;
    return ( (pos1 - pos2) % FRAME_NUCLEOTIDES ) == 0;
}

const double INF_BITS = 1000000; // todo: make it right

double inter_bits_for(p_string text, int kmer_from, hash_t kmer, int kmer_len, PredictedFrequences &pred_freqs, int window_r)
{
//    return 1000000; // todo: remove
	//int from = 0, to = 0;
	//int kmer_len = kmer.size();
	//get_window_from_to(from, to, kmer_from, text.size(), kmer.size(), window_r);
	vector<double> found_freqs;

    do_for_left_and_right_sides(text, kmer_from, kmer_len, window_r, [&](int pos)
		{
            if (!same_frame(pos, kmer_from)) // todo: remove or optomoze
                return;

			if (pos + kmer_len > text.len) // todo: remove
				throw "pos + kmer_len > int(text.size()) 2";

			auto with_kmer = get_kmer_from(text.s, pos, kmer_len);
			if (pred_freqs.find(with_kmer) == pred_freqs.end())
				return;

			auto &one_kmer_freq = pred_freqs[with_kmer];
			if (one_kmer_freq.find(kmer) == one_kmer_freq.end())
				return;

			double freq = one_kmer_freq[kmer];
			found_freqs.push_back(freq);
//			cout << with_kmer << " - " << freq << endl;
		});

	if (found_freqs.empty())
		return INF_BITS; // todo: make it right

	return bits_for_freq(combine_freqs(found_freqs));
}

struct BitsForResult
{
    double bits;
    bool intra;
    BitsForResult(double bits, bool intra) : bits(bits), intra(intra){}
};

BitsForResult bits_for(p_string text, int pos, hash_t kmer, int kmer_len, Frequences &freqs, PredictedFrequences &pred_freqs, int inter_window_r)
{
	if (pos + kmer_len > text.len)
		throw "pos + kmer_len > text.size()";

	auto intra = intra_bits_for(text, pos, kmer, freqs);
	auto inter = inter_bits_for(text, pos, kmer, kmer_len, pred_freqs, inter_window_r);

    if (intra < inter)
        return BitsForResult(intra, true);

    return BitsForResult(inter, false);
//	cout << Hash<hash_t>::str_from_hash(kmer, kmer_len) << " intra " << intra << " inter " << inter << endl;
//	return std::min(intra, inter);
}

void compress(p_string text, int pos, WindowBits &window_bits, int kmer_len, Frequences &freqs, PredictedFrequences &pred_freqs, int inter_window_r, bool rev_compl)
{
	const int COMPRESSION_WINDOW = 3; //inter_window_r; // todo: think ? 8 ? kmer_len * 2 ? parameter?

	auto kmer = get_kmer_from(text.s, pos, kmer_len);
	auto bits_result = bits_for(text, pos, kmer, kmer_len, freqs, pred_freqs, inter_window_r);
    auto bits = bits_result.bits;

	int window = pos / COMPRESSION_WINDOW;
    int frame  = pos - window * COMPRESSION_WINDOW;

    WindowBits::CompressedAmino compressed(frame, bits, bits_result.intra);

	if (window_bits.windows.find(window) == window_bits.windows.end())
		window_bits.windows[window] = compressed;
	else
        if (window_bits.windows[window].bits > bits)
		    window_bits.windows[window] = compressed;

    //if (!rev_compl) // && frame == 2)
    //{
    //    cout << "at pos " << pos << " frame " << frame << " " << bits << endl;
    //    cout << "window_bits.windows[window]: " << window_bits.windows[window] << endl;
    //}
}

struct WindowBitsBothDirections
{
    WindowBits direct, rev_complement;
};

bool file_exists(const string &filename)
{
    ifstream f(filename);
    return f.good();
}

struct CompressResult
{
    double bits, intra_part;
    CompressResult(double bits, double intra_part) : bits(bits), intra_part(intra_part){}
};

CompressResult compress(const string &freq_filename, const string &fasta_filename)
{
    if (!file_exists(freq_filename)) // todo: remove
        return CompressResult(-1, 0);

    int kmer_len = 0, window_r = 0;
	Frequences freqs;
	PredictedFrequences pred_freqs;
    FreqDBIO::load_frequences(freq_filename, freqs, pred_freqs, kmer_len, window_r);

	list<WindowBitsBothDirections> window_bits;
    int contig_number = 0;
    SeqLoader::for_every_clean_sequence_do(fasta_filename, [&](p_string line)
    {
        WindowBitsBothDirections bits;
//        const int SEARCH_WINDOW_R = 3;
        {
    	    for (int i=0; i <= line.len - kmer_len; i++)
	    	    compress(line, i, bits.direct, kmer_len, freqs, pred_freqs, window_r, false);
        }

//        cout << "rev complement" << endl;

        {
            string rev_complement(line.s, line.len);
            rev_complement = seq_transform_actg::apply_transformation(rev_complement, true, true);
    	    for (int i=0; i <= int(rev_complement.size()) - kmer_len; i++)
	    	    compress(p_string(&rev_complement[0], rev_complement.size()), i, bits.rev_complement, kmer_len, freqs, pred_freqs, window_r, true); //window_r);
        }

        window_bits.push_back(bits);

//        {
//	        auto direct_bits = bits.direct.sum_bits_with_frame_penalties();
//	        auto rev_complement_bits = bits.rev_complement.sum_bits_with_frame_penalties();
////            cout << contig_number << " contig bits: " << min(bits.direct.sum_bits()/bits.direct.windows.size(), bits.rev_complement.sum_bits()/bits.rev_complement.windows.size()) << endl;
//        }
        contig_number++;
    });

    double sum_bits = 0;
    size_t sum_windows = 0;
    int sum_intra_windows = 0;

    for (auto &bits : window_bits)
    {
	    auto direct_bits = bits.direct.sum_bits_with_frame_penalties();
	    auto rev_complement_bits = bits.rev_complement.sum_bits_with_frame_penalties();
        if (bits.direct.windows.size() != bits.rev_complement.windows.size())
            throw "bits.direct.windows.size() != bits.rev_complement.windows.size()";

//        cout << "direct         bits per window :" << direct_bits/bits.direct.windows.size() << endl;
//        cout << "rev complement bits per window :" << rev_complement_bits/bits.rev_complement.windows.size() << endl;

        sum_windows += bits.direct.windows.size();
        if (direct_bits < rev_complement_bits)
        {
            sum_bits += direct_bits;
            sum_intra_windows += bits.direct.sum_intra_windows();
        }
        else
        {
            sum_bits += rev_complement_bits;
            sum_intra_windows += bits.rev_complement.sum_intra_windows();
        }
    }

//	cout << "compressed" << endl;
//	cout << "sum bits: " << sum_bits << endl;
 //   if (sum_windows != 0)
//	    cout << "bits per window: " << sum_bits/sum_windows << endl;

    return CompressResult(sum_bits/sum_windows, double(sum_intra_windows)/sum_windows);
}

