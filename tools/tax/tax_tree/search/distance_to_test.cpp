#include "distance_to.h"

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
            OneKmerFreq a, b;
            {
                a[1] = 1.0;
                a[2] = 2.0;
            }
            {
                b[1] = 2.0;
                b[2] = 4.0;
            }

            equal(similarity_of(a, b), 1.0);
        }

        {
            OneKmerFreq a, b;
            {
                a[1] = 1.0;
                a[2] = 2.0;
            }
            {
                b[1] = 1.0;
                b[2] = 2.0;
            }

            equal(similarity_of(a, b), 1.0);
        }


        //{
        //    PredictedFrequences a, b;
        //    {
        //        OneKmerFreq x;
        //        x[1] = 0.5;
        //        x[2] = 0.5;
        //        a[9] = x;
        //    }

        //    {
        //        OneKmerFreq x;
        //        x[1] = 0.5;
        //        x[3] = 0.5;
        //        b[9] = x;
        //    }

        //    equal(distance_to(a, b), 1);
        //}

        //{
        //    PredictedFrequences a, b;
        //    {
        //        OneKmerFreq x;
        //        x[1] = 0.5;
        //        x[2] = 0.5;
        //        a[9] = x;
        //    }

        //    {
        //        OneKmerFreq x;
        //        x[1] = 0.5;
        //        x[2] = 0.5;
        //        b[10] = x;
        //    }

        //    equal(distance_to(a, b), 2);
        //}

        //{
        //    PredictedFrequences a, b;
        //    {
        //        OneKmerFreq x;
        //        x[1] = 0.5;
        //        x[2] = 0.5;
        //        a[9] = x;
        //    }

        //    //{
        //    //    OneKmerFreq x;
        //    //    x[1] = 0.5;
        //    //    x[2] = 0.5;
        //    //    b[10] = x;
        //    //}

        //    equal(distance_to(a, b), 1);
        //}

        //{
        //    PredictedFrequences a, b;
        //    {
        //        OneKmerFreq x;
        //        x[1] = 0.5;
        //        x[2] = 0.5;
        //        a[10] = x;
        //    }

        //    {
        //        OneKmerFreq x;
        //        x[1] = 0.5;
        //        x[2] = 0.5;
        //        b[10] = x;
        //    }

        //    equal(distance_to(a, b), 0);
        //}

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
