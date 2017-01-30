#ifndef FREQ_DB_H_INCLUDED
#define FREQ_DB_H_INCLUDED

#include <map>
//#include <iostream>

typedef unsigned short hash_t; // kmers up to 16 nucleotides
typedef std::map<hash_t, double> OneKmerFreq;
typedef std::map<hash_t, OneKmerFreq > PredictedFrequences;

struct Frequences
{
	OneKmerFreq freqs;
	size_t total;
	Frequences() : total(0){}

	void add(hash_t s)
	{
		freqs[s]++;
		total++;
//        if ((total % (1024*1024)) == 0)
//            cerr << ".";
	}
};

struct Freq
{
    static void to_freq(OneKmerFreq &freq, size_t total)
    {
	    for (auto &it : freq)
		    it.second /= total;	
    }

    static size_t calculate_total(const OneKmerFreq &freq)
    {
	    size_t total = 0;
	    for (auto &it : freq)
		    total += it.second;

	    return total;
    }

    static void add(OneKmerFreq &a, const OneKmerFreq &b)
    {
//        total += b.total;
        for (auto &f : b)
            a[f.first] += f.second;
    }

    static void normalize(OneKmerFreq &a, int norm)
    {
        for (auto &f : a)
            f.second /= norm;
    }

    static void add(PredictedFrequences &a, const PredictedFrequences &b)
    {
        for (auto &kmer_it : b)
            add(a[kmer_it.first], kmer_it.second);
    }

    static void normalize(PredictedFrequences &a, int norm)
    {
        for (auto &kmer_it : a)
            normalize(kmer_it.second, norm);
    }
};

#endif