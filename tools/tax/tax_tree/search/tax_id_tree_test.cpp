#include "tax_id_tree.h"
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
		{
			TaxIdTree tree;
			TaxIdTreeLoader::load_tax_id_tree(tree, "./test/tax.parents.small");
			equal(tree.consensus_of(1, 1), 1);
			equal(tree.consensus_of(2, 3), 1);
			equal(tree.consensus_of(2, 1), 1);
			equal(tree.consensus_of(1, 2), 1);
			equal(tree.consensus_of(3, 7), 3);
			equal(tree.consensus_of(7, 3), 3);
			equal(tree.consensus_of(1, 7), 1);
			equal(tree.consensus_of(4, 6), 2);
			equal(tree.consensus_of(6, 4), 2);
			equal(tree.consensus_of(7, 4), 1);
		}

		{
			TaxIdTree tree;
			TaxIdTreeLoader::load_tax_id_tree(tree, "./test/tax.parents");
			equal(tree.consensus_of(470, 909768), 909768);
			equal(tree.consensus_of(1117, 909768), 2);
		}
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
