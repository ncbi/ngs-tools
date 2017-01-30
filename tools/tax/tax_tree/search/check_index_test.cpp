#include "check_index.h"

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
		TaxIdTree tax_id_tree;
		TaxIdTreeLoader::load_tax_id_tree(tax_id_tree, "./test/tax.parents.small");

#if 1
		{
			Kmers kmers(tax_id_tree);
			load_kmers(kmers, "./test/kmers2");
			cerr << kmers.storage.size() << " kmers loaded" << endl;

//			load_kmers(kmers, "./test/kmers");
			check_kmers(kmers, "./test/tax1.fasta", 4, 32);
			check_kmers(kmers, "./test/tax2.fasta", 5, 32);
			check_kmers(kmers, "./test/tax3.fasta", 3, 32);

			print_kmers(kmers, 32);
			/* will be merged with: expect:
7			3					3
4			5					2
5			-					5
6			4 5					2
2			4 5 3				1 <- will not be printed

expect:

TTTCGGGTTTAAAAAGGTTCTCGAGGGCGCTC        2
ATTGGCTTTCGGGCGTTCTGGGGAAGGGTCAA        5
ACATTAATTTTCCCCGCCCCGAAAAGAATTCC        2
AAAAGAAAAGCCCCTTTTTACGGGTAATCGAC        3
			*/

			/*
AAAAGAAAAGCCCCTTTTTACGGGTAATCGAC        3       624     4096    4577
ACATTAATTTTCCCCGCCCCGAAAAGAATTCC        2       3218    6087
ATTGGCTTTCGGGCGTTCTGGGGAAGGGTCAA        100     101     3708    4096
TTTCGGGTTTAAAAAGGTTCTCGAGGGCGCTC        1       2       130583
TTTGGAATTTTTTTTTTTTTTTTTTCCGGCCC        1       2       3       5857    6945    7230
*/
		}

		{ // if 1 char variation is allowed
			Kmers kmers(tax_id_tree);
			load_kmers(kmers, "./test/kmers3");
			cerr << kmers.storage.size() << " kmers loaded" << endl;

			check_kmers(kmers, "./test/tax1.fasta", 4, 32);
			check_kmers(kmers, "./test/tax2.fasta", 5, 32);
			check_kmers(kmers, "./test/tax5.fasta", 3, 32);

			print_kmers(kmers, 32);
			/*
expect:

TTTGGAATTTTTTTTTTTTTTTTTTCCGGCCC        2
TTTCGGGTTTAAAAAGGTTCTCGAGGGCGCTC        2
ATTGGCTTTCTGGCGTTCTGGGGAAGGGTCAA        5
ACATTAATTTTCCCCGCCCCGAAAAGAATTCT        2
AAAAGAAAAGCCCCTTTTTACGGGTAATCGAC        3
*/
		}

#endif
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
