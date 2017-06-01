/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

#include <ctime>
#include "seq_transform.h"
#include <iostream>
#include <chrono>
#include <assert.h>
#include <cstring>
#include <vector>

#include "hash.h"

using namespace std;
using namespace std::chrono;

#define ASSERT_EQUALS(a, b) { auto&& aa = a; auto&& bb = b; if (!(aa == bb)) { std::cerr << aa << " != " << bb << std::endl << std::flush; assert(a == b); } }
#define ASSERT(a) assert(a)
#define equal ASSERT_EQUALS // backwrad compatibility

// test manager
// do not use directly, use defines below
typedef void (*TestFunctionPtr)();
class TestManager {
public:
    static TestManager& get() {
        static TestManager instance;
        return instance;
    }

    void add(const char* name, TestFunctionPtr fn, bool disabled = false) {
        tests.push_back({name, fn, disabled});
    }

    int run(const char* name = nullptr) {
        int count = 0;
        for (auto& test: tests) {
            if (!name && test.disabled) {
                std::cerr << "WARNING: skipping disabled test " << test.name << std::endl;
                continue;
            }
            if ((!name && !test.disabled) || (name && !strcmp(name, test.name))) {
                std::cerr << "Running test " << test.name << std::endl;
                (*test.ptr)();
                ++count;
            }
        }
        return count;
    }

    int main(int argc, const char* const* argv) {
        if (argc > 1) {
            std::cerr << "Running " << (argc - 1) << " tests" << std::endl;
            for (int i = 1; i < argc; ++i) {
                if (!run(argv[i])) {
                    std::cerr << "Test not found: " << argv[i] << std::endl;
                    return 1;
                }
            }
        } else {
            if (tests.empty()) {
                std::cerr << "No tests found" << std::endl;
                return 1;
            } else {
                std::cerr << "Running " << tests.size() << " tests" << std::endl;
                run(nullptr);
            }
        }
        return 0;
    }

private:
    struct Test {
        const char* name;
        TestFunctionPtr ptr;
        bool disabled;
    };
    std::vector<Test> tests;
};
struct TestCase {
    TestCase(const char* name, TestFunctionPtr fn, bool disabled) { TestManager::get().add(name, fn, disabled); }
};

// define and register test
#define TEST(function) void function(); static TestCase _test_case_##function(#function, function, false); void function()
// define for temporarily disabled tests
#define DISABLED_TEST(function) void function(); static TestCase _test_case_##function(#function, function, true); void function()
// define standard main function
#define TEST_MAIN() int main(int argc, char** argv) { return TestManager::get().main(argc, argv); }


template <class TestKmerMap>
void build_test_map(const std::vector<string> &reads, TestKmerMap &kmer_map)
{
    for (auto& read: reads) {
        Hash<typename TestKmerMap::hash_t>::for_all_hashes_do(read, kmer_map.kmer_len, [&](typename TestKmerMap::hash_t hash) {
                kmer_map.add(hash);
                return true;
            });
    }
	kmer_map.optimize();
}
