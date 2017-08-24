#ifndef CONTIG_H_INCLUDED
#define CONTIG_H_INCLUDED

#include <string>
#include <list>

struct Contig
{
	const std::string seq;
//	const std::string desc;
	const double data_percent;
	double average_coverage;
	Contig(const std::string &seq, double data_percent, double average_coverage) : seq(seq), data_percent(data_percent), average_coverage(average_coverage){}

	bool operator < (const Contig &c2) const
	{
		return c2.data_percent < data_percent;
	}

//bool contig_less(const Contig &c1, const Contig &c2)
//{
//	return c2.data_percent < c1.data_percent;
//}

};

typedef std::list<Contig> Contigs;

static double percent_sum(const Contigs &contigs)
{
	double sum = 0;
	for (auto &c : contigs)
		sum += c.data_percent;

	return sum;
}

#endif
