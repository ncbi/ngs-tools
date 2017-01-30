#ifndef KMERS_LOADING_H_INCLUDED
#define KMERS_LOADING_H_INCLUDED

#include <fstream>
#include <algorithm>
#include <list>
#include "kmers.h"
#include "filesystem.h"

struct KmerLoading
{
	struct Line
	{
		std::string s;
		size_t count;

		Line(const std::string &s = std::string(), size_t count = 0) : s(s), count(count) { }
		bool operator < (const Line &line) const { return count < line.count; }
	};

	std::list<Line> lines; // we need to store them here forever

	KmerLoading(Kmers &kmers, const std::string &file_mask)
	{
		auto filenames = FileSystem::get_filenames_by_mask(file_mask);
		for (auto &f : filenames)
			load_lines(f);

		load_kmers(kmers);
	}

	void load_lines(const std::string &filename)
	{
		std::ifstream f(filename);
		if (f.fail())
			throw std::string("cannot open file ") + filename;

//		double total_count = 0;

		while (!f.eof())
		{
			Line l;
			std::getline(f, l.s, '|');

			if (l.s.empty())
				break;

			std::string s_count;
			std::getline(f, s_count);
			l.count = std::stoll(s_count);
			lines.push_back(l);
		}
	}

	void load_kmers(Kmers &kmers)
	{
//		std::cout << "soring now" << std::endl;
//		lines.sort();
//		std::sort(lines.begin(), lines.end());
#if 1
		std::map<int, double> total_weight; // todo: replace map with vector
		for (auto &l: lines)
			total_weight[l.s.length()] += l.count;

		for (auto &l: lines)
			kmers.add(&l.s[0], l.s.length(), bits_for(l.count, total_weight[l.s.length()]));
#else
		double total_weight;
		for (auto &l: lines)
			total_weight += l.count;

		for (auto &l: lines)
			kmers.add(&l.s[0], l.s.length(), bits_for(l.count, total_weight));
#endif
	}

	static float bits_for(size_t count, double sum_weight)
	{
		const int MANY_BITS = 20;
		if (count == 0)
			return MANY_BITS;

		auto part = 1.0*count/sum_weight;
		return log2(float(1.0/part));
	}

};

#endif