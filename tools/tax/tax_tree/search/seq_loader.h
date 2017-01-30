#ifndef SEQ_LOADER_H_INCLUDED
#define SEQ_LOADER_H_INCLUDED

#include "fasta.h"
#include "seq_cleaner.h"
#include <thread>
#include <iostream>
#include <string>

struct SeqLoader
{
    struct ReadySeq
    {
	    std::string seq;
	    SeqCleaner::p_strings clean_strings;
    };

    template <class Lambda>
    static size_t for_every_clean_sequence_do(const std::string &filename, Lambda &&lambda)
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
			    lambda(clean_string);

		    seq_index++;
		    if (seq_index % DOT_INTERVAL == 0)
			    std::cerr << ".";

		    loading_thread.join();
		    swap(processing_seq, loading_seq); // todo: move?
	    }

	    if (seq_index >= DOT_INTERVAL)
		    std::cerr << std::endl;

	    return total_size;
    }

private:
    static void swap(ReadySeq &a, ReadySeq &b)
    {
	    std::swap(a.seq, b.seq);
	    std::swap(a.clean_strings, b.clean_strings);
    }

    static void load_sequence(Fasta *_fasta, ReadySeq *_seq)
    {
	    Fasta &fasta = *_fasta;
	    ReadySeq &seq = *_seq;

	    seq.seq.clear(); // for better performance clear instead of constructor
	    seq.clean_strings.clear();

	    if (!fasta.get_next_sequence(seq.seq))
		    return;

	    SeqCleaner cleaner(seq.seq);
	    seq.clean_strings = move(cleaner.clean_strings);
    }
};

#endif