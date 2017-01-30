#include "kmers_loading.h"
#include <iostream>

using namespace std;

template <class A, class B>
void equal(A a, B b)
{
	if (a == b)
		cout << "ok" << endl;
	else
	{
		cout << "FAIL: " << a << " != " << b << endl;
//		exit(1);
	}
}

int main(int argc, char const *argv[])
{
	try
	{
		Kmers kmers;
		KmerLoading l(kmers, "./tests/kmer_loading_test_data/stat*.txt");

		auto kmer = kmers.get_kmer("NKHHH", 5);
		equal(int(kmer.bits), 13);
		equal(int(kmers.get_kmer("IL", 2).bits), 1);
		equal(int(kmers.get_kmer("GG", 2).bits), 0);
	}
	catch (const string &s)
	{
		cout << s << endl;
	}
	catch (const char *s)
	{
		cout << s << endl;
	}
}

