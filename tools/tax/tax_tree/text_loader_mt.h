#ifndef TEXT_LOADER_MT_H_INCLUDED
#define TEXT_LOADER_MT_H_INCLUDED

#include <list>
#include <string>
#include <fstream>
#include <vector>
#include <mutex>
#include <thread>

struct TextLoaderMT
{
	std::ifstream f;
	std::mutex file_mutex;
	std::list< std::string > lines;

	TextLoaderMT(const std::string &filename) : f(filename)
	{
		if (f.fail())
			throw "cannot open input file";
	}

	std::string *load_next_sequence()
	{
		std::lock_guard<std::mutex> lock(file_mutex);
		if (!f.eof())
		{
			std::string line;
			std::getline(f, line);
			//if (line.empty())
			//{
			//	if (!f.eof())
			//		//std::cerr << "unexpected end of file" << std::endl;
			//	return nullptr;
			//}
			lines.push_back(line);
			return &*lines.rbegin();
		}

		return nullptr;
	}
};


struct TextLoaderMTNoStore
{
	std::ifstream f;
	std::mutex file_mutex;
	size_t seq_loaded;

	TextLoaderMTNoStore(const std::string &filename) : f(filename), seq_loaded(0)
	{
		if (f.fail())
			throw "cannot open input file";
	}

	bool load_next_sequence(std::string &s, size_t &seq_id)
	{
		std::lock_guard<std::mutex> lock(file_mutex);
//		if (f.eof())
//			return false;

//		std::string _s;
//		f.getline(s);
		s.clear();
		while (s.empty() && !f.eof())
			std::getline(f, s);
		
//		if (s.empty())
//			return load_next_sequence(s, seq_id);
//		s = _s; // todo: remove
//		f >> s;
//		seq_id = seq_loaded;
//		seq_loaded++;

		return !s.empty();
	}
};

// todo: traits
struct TextLoaderSTNoStore
{
	std::ifstream f;

	TextLoaderSTNoStore(const std::string &filename) : f(filename)
	{
		if (f.fail())
			throw "cannot open input file";
	}

	bool load_next_sequence(std::string &s)
	{
		s.clear();
		while (s.empty() && !f.eof())
			std::getline(f, s);

		return !s.empty();
	}
};

struct FastaWithTaxonomyLoader
{
	std::ifstream f;

	FastaWithTaxonomyLoader(const std::string &filename) : f(filename)
	{
		if (f.fail())
			throw "cannot open input file";
	}

	bool load_next_sequence(std::string &s, int &tax_id)
	{
		s.clear();
		tax_id = 0;
		if (f.fail() || f.eof())
			return false;

		f >> s;
		f >> tax_id;

		return !s.empty();
	}

};

#endif
