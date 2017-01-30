#include <iostream>
#include <chrono>
#include <thread>
#include <list>
#include <iomanip>
#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include "config_next_words.h"
#include "text_loader_mt.h"
#include "translation.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.1";

void print_current_time()
{
	auto t = std::time(nullptr);
	auto timeinfo = std::localtime(&t);
	const int BUFFER_SIZE = 256;
	char buffer[BUFFER_SIZE];
	strftime(buffer, BUFFER_SIZE, "%m/%d/%Y %H:%M:%S", timeinfo);
	cerr << "time is " << buffer << endl;
}

#define TRANSLATE 1

struct NextWords
{
//	typedef set<string> WordSet;
	typedef map<string, unsigned int> WordSetCount;
	typedef map<string, WordSetCount> KeyWords;

	Translation tr;

	KeyWords key_words;

	static bool has_stop(const string &s)
	{
		return s.find('*') != string::npos;
	}

	void add(const string &word, const string &next_word)
	{
#if !TRANSLATE
		key_words[word].insert(next_word); // todo: thread-safe
#else
		string word_tr;
		tr.translate(word, 0, word_tr);
		if (has_stop(word_tr))
			return;

		string next_word_tr;
		tr.translate(next_word, 0, next_word_tr);

		if (has_stop(next_word_tr))
			return;

//		key_words[word_tr].insert(next_word_tr); // todo: thread-safe
		key_words[word_tr][next_word_tr] ++;
#endif
	}
};

template <class Lambda>
void for_every_word(const string &s, Lambda &&lambda)
{
	istringstream iss(s);
	string word;
	while (iss)
	{
		iss >> word; // todo: optimize using istream_iterator
		lambda(word);
	}
}

void load_next_words(TextLoaderMTNoStore &seq_loader, NextWords &next_words)
{
	string seq;
	size_t seq_id = 0;
	int skipped = 0;
	size_t total_len = 0;
	while (seq_loader.load_next_sequence(seq, seq_id))
	{
		string prev_word;
		size_t bad_skips = 0, good_skips = 0;
		for_every_word(seq, [&](const string &word)
		{
			const int MIN_KEY_WORD_LEN = 12;
			if (word.length() < MIN_KEY_WORD_LEN)
			{
				skipped += word.length();
				return;
			}

			if (skipped > 50)
				bad_skips++;
			else
				good_skips++;

			skipped = 0;

			if (!prev_word.empty())
				next_words.add(prev_word, word);

			prev_word = word;
		});

		total_len+=seq.length();

		cerr << "skips: good " << good_skips << " bad " << bad_skips << " total " << (total_len + 500000)/1000000 << "M" << endl;
	}

//	cerr << ".";
}

void print_next_words(const NextWords &next_words)
{
	const int MAX_NEXT_WORDS = 100;

	for (auto &key_word : next_words.key_words)
	{
//		if (key_word.second.size() >= MAX_NEXT_WORDS)
//			continue;

		cout << key_word.first;
//		for (auto &word : key_word.second)
//			cout << " " << word;

		auto &word_set = key_word.second;
		for (auto &word_count : word_set)
		{
			auto count = word_count.second;
			if (count > 1)
				cout << " " << word_count.first;
		}


		cout << endl;
	}
}

//const int THREADS = 1;

int main(int argc, char const *argv[])
{
    try
    {
		ConfigNextWords config(argc, argv);
		cerr << "next_words version " << VERSION << endl;

		print_current_time();
		auto before = high_resolution_clock::now();

		TextLoaderMTNoStore seq_loader(config.input);

		NextWords next_words;

		//#pragma omp parallel num_threads(THREADS)
		{
			load_next_words(seq_loader, next_words);
		}

		cerr << endl << "printing results" << endl;
		print_next_words(next_words);

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
