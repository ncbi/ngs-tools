#include "feature_tree.h"
#include "feature_tree_io.h"

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

Features fs(const string &s)
{
    return FeatureOperations::calculate_features(&s[0], s.length());
}

int main(int argc, char const *argv[])
{
    try
    {
        {
            Tree t;
            t.set_head(fs("good"));
            TreeIO::print_tree(t);
            cout << "-------------" << endl;

            t.insert_neighbour(t.head_id, fs("bad"));
            TreeIO::print_tree(t);
            cout << "-------------" << endl;

            t.insert_neighbour(1, fs("bod"));
            TreeIO::print_tree(t);
            cout << "-------------" << endl;

            t.insert_neighbour(1, fs("bad"));
            TreeIO::print_tree(t);
            cout << "-------------" << endl;

            TreeIO::save_tree(t, "./tree.save", 2);
            TreeIO::load_tree(t, "./tree.save");
            cout << "loaded:" << endl;
            TreeIO::print_tree(t);

            TreeIO::save_tree(t, "./tree.save", 2);
            TreeIO::load_tree(t, "./tree.save");
            cout << "loaded:" << endl;
            TreeIO::print_tree(t);
        }

        return 0;
    }
    catch ( exception & x )
    {
        cerr << x.what() << endl;
		cerr << "exit 3" << endl;
		return 3;
    }
    catch ( string & x )
    {
        cerr << x << endl;
		cerr << "exit 4" << endl;
		return 4;
    }
    catch ( const char * x )
    {
        cerr << x << endl;
		cerr << "exit 5" << endl;
		return 5;
    }
    catch ( ... )
    {
        cerr << "unknown exception" << endl;
//		cerr << "exit 6" << endl;
		return 6;
    }
}
