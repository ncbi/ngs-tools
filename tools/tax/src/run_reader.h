#ifndef RUN_READER_H_INCLUDED
#define RUN_READER_H_INCLUDED

#include <string>
#include <map>
// todo: remove unnecessary includes
#include <ngs/ncbi/NGS.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>

struct RunReader
{
/*
	template <class Lambda>
	static void read_run(const ngs::String &acc, Lambda &&lambda)
	{
		ngs::ReadCollection run = ncbi::NGS::openReadCollection( acc );
		ngs::ReadIterator it = run.getReadRange( 1, run.getReadCount(), ngs::Read::all ); // todo: get biological only ?
		while ( it.nextRead() )
			while ( it.nextFragment() )
				lambda(it.getFragmentBases());
	}

	static size_t read_count(const ngs::String &acc, ngs::Read::ReadCategory cat)
	{
		ngs::ReadCollection run = ncbi::NGS::openReadCollection( acc );
		return run.getReadCount(cat);
	}
*/

	static long long int load_run(const ngs::String &acc, std::vector<std::string> &reads)
	{
		ngs::ReadCollection run = ncbi::NGS::openReadCollection( acc );
		std::cout << "open" << std::endl;
		auto cat = ngs::Read::unaligned;
//		ngs::ReadIterator it = run.getReads(cat); // run.getReadRange( 1, run.getReadCount(), cat); //ngs::Read::all ); // todo: get biological only ?
		ngs::ReadIterator it = run.getReadRange( 1, run.getReadCount(cat), cat);
		std::cout << "got reads" << std::endl;

		reads.clear();
//		reads.reserve(run.getReadCount(cat));
		long long int total_len = 0;
		size_t count = 0;

	    std::cout << "mem usage before nextread M " << mem_usage()/1000000 << std::endl;
		while ( it.nextRead() )
		{
			while ( it.nextFragment() )
			{
				auto s = it.getFragmentBases();
				total_len += s.size();
				reads.push_back(s.toString()); // todo: is it efficient enough ?
				if ((count % 1000000) == 0)
                {
					std::cout << count << " mem usage M " << mem_usage()/1000000 << std::endl;
                }
				count++;
			}
		}

		return total_len;
	}

};

#endif
