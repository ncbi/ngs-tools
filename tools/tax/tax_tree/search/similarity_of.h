#include <iostream>
#include <chrono>
#include <vector>
#include <thread>
#include <omp.h>
#include <algorithm>
#include "freq_db.h"
#include "freq_db_io.h"

using namespace std;
using namespace std::chrono;

typedef vector<double> Vector;

double norm_sq2(const Vector &v)
{
    double sum = 0;
    for (auto x : v)
        sum += x * x;
    
    return std::sqrt(sum);
}

double norm_sq2(const OneKmerFreq &v)
{
    double sum = 0;
    for (auto x : v)
        sum += x.second * x.second;
    
    return std::sqrt(sum);
}

double mul(const Vector &a, const Vector &b)
{
    if (a.size() != b.size())
        throw "diff a.size() != b.size()";

    double sum = 0;
    for (size_t i = 0; i < a.size(); i++)
        sum += a[i] * b[i];

    return sum;
}

double similarity_of(OneKmerFreq &a, OneKmerFreq &b) // todo: template
{
//    return 1.0;
    double x = 1.0;
    for (auto &a_it : a)
    {
//        va.push_back(a_it.second);
  //      vb.push_back(b[a_it.first]);
        x *= (1.0 - a_it.second) * (1.0 - b[a_it.first]);
//        cout << a_it.first << " " << a_it.second << " " << b[a_it.first] << " " << x << endl;
    }

    return 1.0 - x;
}

/*

double similarity_of(OneKmerFreq &a, OneKmerFreq &b) // todo: template
{
//    return 1.0;
    Vector va, vb;
    va.reserve(a.size());
    vb.reserve(a.size());

    for (auto &a_it : a)
    {
        cout << a_it.first << " " << a_it.second << " " << b[a_it.first] << endl;
        va.push_back(a_it.second);
        vb.push_back(b[a_it.first]);
    }

    auto na = norm_sq2(a);
    auto nb = norm_sq2(b);

    auto n = na * nb;

    auto sum = mul(va, vb);
    if (n > 0)
        return sum / n;

    if (sum == 0 && na == nb) // todo: remove sum == 0
        return 1.0;

    return 0;
}

double similarity_of(PredictedFrequences &a_freqs, PredictedFrequences &b_freqs, OneKmerFreq &fa, OneKmerFreq &fb)
{
    double sum = 0;
    for (auto &ita : fa)
    {
        auto kmer = ita.first;
        sum += ita.second * fb[kmer] * similarity_of(a_freqs[kmer], b_freqs[kmer]);
    }

    auto na = norm_sq2(fa);
    auto nb = norm_sq2(fb);

//    cout << "sum: " << sum << endl;
 //   cout << "norma a: " << na << endl;
  //  cout << "norma b: " << nb << endl;

    return sum/(na * nb);
}
*/

double bits_for(OneKmerFreq &a, OneKmerFreq &b) // todo: template
{
    //return 1.0;
    auto s = similarity_of(a, b);
    if (s == 0)
        return 20; // todo: tune

    return log2(1.0/s);
}


double similarity_of(PredictedFrequences &a_freqs, PredictedFrequences &b_freqs, OneKmerFreq &fa, OneKmerFreq &fb)
{
    double sum = 0;
    for (auto &ita : fa)
    {
        auto kmer = ita.first;
//        sum += ita.second * bits_for(a_freqs[kmer], b_freqs[kmer]); //  similarity_of(a_freqs[kmer], b_freqs[kmer]);
        double t = ita.second * similarity_of(a_freqs[kmer], b_freqs[kmer]);
        sum += t;
//        cout << sum << endl;
//        throw "";
    }

//    auto na = norm_sq2(fa);
//    auto nb = norm_sq2(fb);

//    cout << "sum: " << sum << endl;
 //   cout << "norma a: " << na << endl;
  //  cout << "norma b: " << nb << endl;

    return sum; ///fa.size();
}


bool file_exists(const string &filename)
{
    ifstream f(filename);
    return f.good();
}

double similarity_of(PredictedFrequences &a_freqs, OneKmerFreq &fa, const string &b_freq_filename)
{
    if (!file_exists(b_freq_filename)) // todo: remove
        return -1;

    int b_kmer_len = 0, b_window_r = 0;
	Frequences b_freqs;
	PredictedFrequences b_pred_freqs;
    FreqDBIO::load_frequences(b_freq_filename, b_freqs, b_pred_freqs, b_kmer_len, b_window_r);

    //if (a_kmer_len != b_kmer_len)
    //    throw "a_kmer_len != b_kmer_len";

    //if (a_window_r != b_window_r)
    //    throw "a_window_r != b_window_r";

    return similarity_of(a_freqs, b_pred_freqs, fa, b_freqs.freqs);
}


double similarity_of(const string &a_freq_filename, const string &b_freq_filename)
{
    if (!file_exists(a_freq_filename) || !file_exists(b_freq_filename)) // todo: remove
        return -1;

    int a_kmer_len = 0, a_window_r = 0;
	Frequences a_freqs;
	PredictedFrequences a_pred_freqs;
    FreqDBIO::load_frequences(a_freq_filename, a_freqs, a_pred_freqs, a_kmer_len, a_window_r);

    return similarity_of(a_pred_freqs, a_freqs.freqs, b_freq_filename);
}

