#include <ktst/unit_test.hpp> 

using namespace std;
using namespace ncbi::NK;

TEST_SUITE(SraSearchTestSuite);

int
main( int argc, char *argv [] )
{
    return SraSearchTestSuite(argc, argv);
}
