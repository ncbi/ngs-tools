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
#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <chrono>
#include "omp_adapter.h"

typedef uint64_t hash_t;

#include "dbs.h"
#include "hash.h"
#include "config_seq_filter.h"
#include "seq_transform.h"
#include "p_string.h"


const std::string VERSION = "0.1";

using namespace std;
using namespace std::chrono;

struct LineReader
{
    std::istream *f = nullptr;
    LineReader()
    {
        f = &std::cin;
        std::ios_base::sync_with_stdio(false);
        std::cin.tie(nullptr);
    }

    bool readline(string &line)
    {
        std::getline(*f, line);
        // handling windows line endings
        if (!line.empty() && *line.rbegin() == '\r')
            line.erase(line.size() - 1);

        return (!f->eof());
    }
};

struct Filter
{
    typedef std::vector<hash_t> HashSortedArray;

	HashSortedArray hash_array;
    int kmer_len = 0;

    Filter(const string &db)
    {
	    kmer_len = DBSIO::load_dbs(db, hash_array);
    }

    bool in_db(hash_t hash)
    {
        hash = seq_transform<hash_t>::min_hash_variant(hash, kmer_len);
    	return std::binary_search(hash_array.begin(), hash_array.end(), hash);
    }

    bool fine_seq(const p_string seq)
    {
        bool found = false;
    	Hash<hash_t>::for_all_hashes_do(seq, kmer_len, [&](hash_t hash)
    	{
	        if (in_db(hash) > 0)
    		    found = true;
                return !found;
    	});

        return !found;
    }
};

bool is_nucl_string_char(char ch)
{
    return ch == 'A' || ch == 'C' || ch == 'T' || ch == 'G' || ch == 'N';
}

bool seems_line_nucl_string(p_string seq)
{
    return seq.len > 32 && is_nucl_string_char(seq.s[0]);
}

// limitations: accepts ACTG only, upper case, no N or other charachters
int main(int argc, char const *argv[])
{
    #ifdef __GLIBCXX__
    std::set_terminate(__gnu_cxx::__verbose_terminate_handler);
    #endif
    
    std::cerr << "seq_filter version " << VERSION << endl;
    auto before = high_resolution_clock::now();

	Config config(argc, argv);

    Filter filter(config.db);

    LineReader reader;

    struct Stat
    {
        size_t total = 0;
        size_t rejected = 0;
        size_t nucl_empty = 0;
        size_t nucl_one_char = 0;
        size_t nucl_valid = 0;
    } stat;

    string line;
    while (reader.readline(line))
    {
        auto const nucl_seq = seq_transform_actg::to_upper_inplace(line));
        if (nucl_seq.len == 0)
            stat.nucl_empty++;
        else if (nucl_seq.len == 1) // usually (*)
            stat.nucl_one_char++;
        else if (seems_line_nucl_string(nucl_seq))
            stat.nucl_valid++;

        if (filter.fine_seq(nucl_seq))
            cout << "pass" << endl;
        else
        {
            cout << "fail" << endl;
            stat.rejected++;
        }

        stat.total++;
    }

    cerr << "total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;
    std::cerr << "total sam lines processed: " << stat.total << endl;
    std::cerr << "sam lines rejected by filter: " << stat.rejected << endl;
    std::cerr << "sam lines having no sequence strings: " << stat.nucl_empty << endl;
    std::cerr << "sam lines having 1 char sequence strings: " << stat.nucl_one_char << endl;
    std::cerr << "sam lines having valid sequence strings: " << stat.nucl_valid << endl;
}
