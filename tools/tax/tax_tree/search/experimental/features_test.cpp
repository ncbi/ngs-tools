#include "features.h"
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
            auto f1 = FeatureOperations::calculate_features("bad", 3);
            print_features(f1);
            cout << f1.norm << endl;
        }
        {
            auto f1 = FeatureOperations::calculate_features("baads", 5);
            print_features(f1);
            cout << f1.norm << endl;
        }
        {
            auto f1 = FeatureOperations::calculate_features("aaaoooo", 7);
            print_features(f1);
            cout << f1.norm << endl;
        }

        {
            auto f1 = FeatureOperations::calculate_features("aaa", 3);
            auto f2 = FeatureOperations::calculate_features("oooo", 4);
            FeatureOperations::add(f1, f2);
            print_features(f1);
            cout << f1.norm << endl;
        }

        {
            auto f1 = FeatureOperations::calculate_features("baads", 5);
            auto f2 = FeatureOperations::calculate_features("bad", 3);
            FeatureOperations::add(f1, f2);
            print_features(f1);
            cout << f1.norm << endl;
        }

        {
            auto f2 = FeatureOperations::calculate_features("baads", 5);
            auto f1 = FeatureOperations::calculate_features("bad", 3);
            FeatureOperations::add(f1, f2);
            print_features(f1);
            cout << f1.norm << endl;
        }

        {
            cout << FeatureOperations::get_similarity(FeatureOperations::calculate_features("baads", 5), FeatureOperations::calculate_features("bad", 3)) << endl;
            cout << FeatureOperations::get_similarity(FeatureOperations::calculate_features("bad", 3), FeatureOperations::calculate_features("bad", 3)) << endl;
            cout << FeatureOperations::get_similarity(FeatureOperations::calculate_features("boo", 3), FeatureOperations::calculate_features("bad", 3)) << endl;
            cout << FeatureOperations::get_similarity(FeatureOperations::calculate_features("boooooo", 7), FeatureOperations::calculate_features("bad", 3)) << endl;
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
