#include "filesystem.h"
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
	equal(FileSystem::dir_of("./"), "./");
	equal(FileSystem::dir_of("./tests/"), "./tests/");
	equal(FileSystem::dir_of("./tests/megamask*.txt"), "./tests/");

	equal(FileSystem::no_dir("./tests/"), "");
	equal(FileSystem::no_dir("./tests/megamask*.txt"), "megamask*.txt");
	equal(FileSystem::no_dir("./megamask*.txt"), "megamask*.txt");
	equal(FileSystem::no_dir("megamask*.txt"), "megamask*.txt");


	equal(FileSystem::matches_mask("./test", "*"), true);
	equal(FileSystem::matches_mask("translation", "trans*"), true);
	equal(FileSystem::matches_mask("translation", "trans*ion"), true);
	equal(FileSystem::matches_mask("translation", "*ion"), true);
	equal(FileSystem::matches_mask("translation", "*iov"), false);
	equal(FileSystem::matches_mask("translation", "tro*ion"), false);
	equal(FileSystem::matches_mask("translation", "tro*"), false);

	try
	{
		{
			auto files = FileSystem::get_filenames_by_mask("./*");
			equal(files.size(), 43);
		}
		{
			auto files = FileSystem::get_filenames_by_mask("./tests/kmer_loading_test_data/stat*.txt");
			equal(files.size(), 3);
			for (auto &f: files)
				cout << f << endl;
		}
	} 
	catch ( string &s)
	{
		cerr << s << endl;
	}

}

