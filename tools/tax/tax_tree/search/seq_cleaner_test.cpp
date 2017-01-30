#include "seq_cleaner.h"

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
		    SeqCleaner sc("ACTG");
            equal(sc.clean_strings.size(), 1);
            auto s = *sc.clean_strings.begin();
            equal(string(s.s, s.len), "ACTG");
        }
        {
		    SeqCleaner sc("ACTGNAAG");
            equal(sc.clean_strings.size(), 2);
            auto it = sc.clean_strings.begin();
            equal(string(it->s, it->len), "ACTG");
            it++;
            equal(string(it->s, it->len), "AAG");
        }
        {
		    SeqCleaner sc("NNNWACTGNNNAAGYYN");
            equal(sc.clean_strings.size(), 2);
            auto it = sc.clean_strings.begin();
            equal(string(it->s, it->len), "ACTG");
            it++;
            equal(string(it->s, it->len), "AAG");
        }
        {
		    SeqCleaner sc("NNNWACTGNNNAAGYYNTTA");
            equal(sc.clean_strings.size(), 3);
            auto it = sc.clean_strings.begin();
            equal(string(it->s, it->len), "ACTG");
            it++;
            equal(string(it->s, it->len), "AAG");
            it++;
            equal(string(it->s, it->len), "TTA");
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
