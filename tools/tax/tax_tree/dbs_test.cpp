#include <cstdint>

typedef uint64_t hash_t;

#include "dbs.h"
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
		vector<int> kmers;
		kmers.push_back(10);
		kmers.push_back(20);
		kmers.push_back(30);

		DBS::save_dbs("./test/tmp.dbs", kmers, 22);
	}

	{
		vector<int> kmers;
		auto kmer_len = DBS::load_dbs("./test/tmp.dbs", kmers);
		equal(kmer_len, 22);
		equal(kmers.size(), 3);
		equal(kmers[0], 10);
		equal(kmers[1], 20);
		equal(kmers[2], 30);
	}
}
