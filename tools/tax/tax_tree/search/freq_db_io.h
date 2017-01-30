#ifndef FREQ_DB_IO_H_INCLUDED
#define FREQ_DB_IO_H_INCLUDED

#include "freq_db.h"
#include <iostream>
#include <fstream>

struct FreqDBIO
{
    template <class ToStr>
    static void save_frequences_text(const std::string &filename, Frequences &freqs, PredictedFrequences &pred_freqs, ToStr &&to_str)
    {
        std::ofstream f(filename);
        f << freqs.total << std::endl;
        f << freqs.freqs.size() << std::endl;
        for (auto &fr : freqs.freqs)
            f << to_str(fr.first) << "\t" << fr.second << std::endl;

        f << pred_freqs.size() << std::endl;
        for (auto &kmer_it : pred_freqs)
        {
            f << to_str(kmer_it.first) << "\t" << kmer_it.second.size() << std::endl;
            for (auto &fr : kmer_it.second)
                f << to_str(fr.first) << "\t" << fr.second << std::endl;
        }
    }

    template <class X>
    static void write(std::ofstream &f, X x)
    {
		f.write((char*)&x, sizeof(x));
    }

    template <class X>
    static void read(std::ifstream &f, X &x)
    {
		f.read((char*)&x, sizeof(x));
    }

    static void save_frequences(const std::string &filename, Frequences &freqs, PredictedFrequences &pred_freqs, int kmer_len, int window_r)
    {
        std::ofstream f(filename, std::ios::out | std::ios::binary);
        {
            const int VERSION = 1;
            write(f, VERSION);
            write(f, kmer_len);
            write(f, window_r);
        }

        write(f, freqs.total);
        write(f, freqs.freqs.size());
        for (auto &fr : freqs.freqs)
        {
            write(f, fr.first);
            write(f, fr.second);
        }

        write(f, pred_freqs.size());

        for (auto &kmer_it : pred_freqs)
        {
            write(f, kmer_it.first);
            write(f, kmer_it.second.size());
            for (auto &fr : kmer_it.second)
            {
                write(f, fr.first);
                write(f, fr.second);
            }
        }

        if (f.fail())
            throw std::runtime_error("save_frequences failed");
    }

    static void load_frequences(const std::string &filename, Frequences &freqs, PredictedFrequences &pred_freqs, int &kmer_len, int &window_r)
    {
        std::ifstream f(filename, std::ios::in | std::ios::binary);
        if (!f.good())
            throw "cannot open freq file";

//        std::cerr << "load freq: " << filename << std::endl;

        // todo: remove?
        freqs = Frequences();
        pred_freqs = PredictedFrequences();

        {
            int version = 0;
            read(f, version);
            if (version != 1)
                throw "unsupported version";

            read(f, kmer_len);
            read(f, window_r);
        //    kmer_len = 6;
          //  window_r = 12;
        }

        {
            read(f, freqs.total);
            size_t freqs_size = 0;
            read(f, freqs_size);
//            std::cerr << freqs_size << " loading" << std::endl;
            for (size_t i = 0; i < freqs_size; i++)
            {
                hash_t kmer = 0;
                double freq = 0;
                read(f, kmer);
                read(f, freq);
                freqs.freqs[kmer] = freq;
            }

 //           std::cerr << freqs_size << " loaded" << std::endl;
        }


        {
            size_t pred_freqs_size;
            read(f, pred_freqs_size);

            for (size_t i = 0; i < pred_freqs_size; i++)
            {
                hash_t kmer = 0;
                read(f, kmer);
                size_t kmer_size = 0;
                read(f, kmer_size);
                auto &kmer_freqs = pred_freqs[kmer];

                for (size_t j = 0; j < kmer_size; j++)
                {
                    hash_t with_kmer = 0;
                    double freq = 0;
                    read(f, with_kmer);
                    read(f, freq);

                    kmer_freqs[with_kmer] = freq;
                }
            }
        }
    }

};

#endif