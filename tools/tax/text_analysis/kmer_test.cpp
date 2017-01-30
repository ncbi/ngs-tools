#include "kmers.h"
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
	{
		Kmers kmers;
		kmers.add("ABC", 3, 1);
		kmers.add("KOT", 3, 20);
		kmers.add("T", 1, 5);

		equal(kmers.max_kmer_len(), 3);
		equal(kmers.get_kmer("A", 1).empty(), true);
		equal(kmers.get_kmer("ABCCCC", 2).empty(), true);

		{
			auto kmer = kmers.get_kmer("ABCCCC", 3);
			equal(kmer.empty(), false);
			equal(kmer.bits, 1);
			equal(kmer.seq == nullptr, false);
		}

		{
			auto kmer = kmers.get_kmer("BCCCC", 4);
			equal(kmer.empty(), true);
		}

		{
			auto kmer = kmers.get_kmer("BCCCC", 40);
			equal(kmer.empty(), true);
		}

		{
			auto kmer = kmers.get_kmer("KOT", 3);
			equal(kmer.bits, 20);
			equal(kmer.empty(), false);
		}

		{
			auto kmer = kmers.get_kmer("T", 1);
			equal(kmer.bits, 5);
		}
	}
}

