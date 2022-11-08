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

#include "reader.h"

#include "tests.h"

#include "fasta_reader.h"
#include "mt_reader.h"
#include "aux_reader.h"

class DummyReader: public Reader {
    size_t idx;
    std::vector<std::string> reads;
public:
    DummyReader(const std::vector<std::string>& reads) : reads(reads), idx(0) {}

    virtual SourceStats stats() const { return SourceStats(reads.size(), reads.size()); }
    virtual float progress() const { return reads.empty() ? 1 : float(idx) / reads.size(); }

    bool read(Fragment* output) override final {
        if (idx >= reads.size()) {
            return false;
        }
        if (output) {
            output->spotid = std::to_string(idx+1);
            output->bases = reads[idx];
        }
        ++idx;
        return true;
    }
};

struct FragmentSort {
    bool operator() (const Reader::Fragment& f1, const Reader::Fragment& f2) const {
        if (f1.spotid != f2.spotid) {
            return f1.spotid < f2.spotid;
        } else {
            return f1.bases < f2.bases;
        }
    }
};

template <typename ReaderPtr>
static std::vector<Reader::Fragment> read_all(ReaderPtr reader) {
    std::vector<Reader::Fragment> result;
    Reader::Fragment fragment;
    while (reader->read(&fragment)) {
        result.push_back(fragment);
        ASSERT(reader->progress() >= 0 && reader->progress() <= 1);
    }
    ASSERT_EQUALS(reader->progress(), 1);
    return result;
}

template <typename ReaderPtr>
static std::vector<std::string> read_all_bases(ReaderPtr reader) {
    std::vector<std::string> result;
    Reader::Fragment fragment;
    while (reader->read(&fragment)) {
        result.push_back(fragment.bases);
        ASSERT(reader->progress() >= 0 && reader->progress() <= 1);
    }
    ASSERT_EQUALS(reader->progress(), 1);
    return result;
}

template <typename ReaderType>
struct Helper {
    template <typename ...Args>
    static std::vector<Reader::Fragment> read_all(Args... args) {
        ReaderType reader(args...);
        return ::read_all(&reader);
    }

    template <typename ...Args>
    static std::vector<std::string> read_all_bases(Args... args) {
        ReaderType reader(args...);
        return ::read_all_bases(&reader);
    }
};

TEST(fasta_reader) {
    FastaReader reader("./tests/data/SRR1068106.fasta");
    Reader::Fragment f;

    ASSERT(reader.read(&f));
    ASSERT_EQUALS(f.spotid, "SRR1068106.1");
    ASSERT_EQUALS(f.bases, "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCACTTCTGAA");

    ASSERT(reader.read(&f));
    ASSERT(reader.read(&f));

    ASSERT(reader.read(&f));
    ASSERT_EQUALS(f.spotid, "SRR1068106.4");
    ASSERT_EQUALS(f.bases, "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCAACTTCCTGAATGTTGACCTCG");

    size_t got_spots = 4;
    while (reader.read(&f))
        ++got_spots;

    ASSERT_EQUALS(got_spots, 226);
    ASSERT(reader.stats() == Reader::SourceStats(got_spots, got_spots));
    ASSERT_EQUALS(reader.progress(), 1);

    auto _unix = Helper<FastaReader>::read_all("./tests/data/SRR1068106.fasta");
    auto dos = Helper<FastaReader>::read_all("./tests/data/SRR1068106.fasta.dos");
    ASSERT(_unix == dos);
}

#if 0
template <typename ReaderType>
void test_mt_reader(const char* type, const char* path) {
    auto reference = Helper<ReaderType>::read_all(path);
    std::sort(reference.begin(), reference.end(), FragmentSort());
    for (int thread_count = 1; thread_count <= 16; thread_count <<= 2) {
        for (size_t chunk_size = 1; chunk_size <= 256; chunk_size <<= 2) {
            std::cout << "checking mt " << type << " reader with thread_count=" << thread_count << " chunk_size=" << chunk_size << " on path=" << path << std::endl;
            auto result = Helper<MTReader<ReaderType> >::read_all(thread_count, chunk_size, path);
            std::sort(result.begin(), result.end(), FragmentSort());
            ASSERT(reference == result);
        }
    }
}
TEST(mt_reader) {
    test_mt_reader<FastaReader>("fasta", "./tests/data/SRR1068106.fasta");
}
#endif

TEST(splitting_cutting_reader) {
    std::vector<std::string> source = {"ATGC",
                                       "ATNGC",
                                       "NATGCN",
                                       "",
                                       "NANTNGNCN",
                                       "NNNN"};
    std::vector<std::string> expected_split = {"ATGC",
                                               "AT", "GC",
                                               "ATGC",
                                               "", // okay to keep empty spots if they are in source, garbage in garbage out, but no crash
                                               "A", "T", "G", "C"};
    std::vector<std::string> expected_cut = {"ATGC",
                                             "AT"};
    std::vector<std::string> split_result = Helper<SplittingReader<DummyReader> >::read_all_bases(source);
    std::vector<std::string> cut_result = Helper<CuttingReader<DummyReader> >::read_all_bases(source);
    ASSERT(split_result == expected_split);
    ASSERT(cut_result == expected_cut);
}

void test_read(const char* path) {
    std::cerr << "Test reading: " << path << std::endl;
    auto reader = Reader::create(path, Reader::Params());
    std::vector<Reader::Fragment> chunk;
    int percent = -1;
    while (reader->read_many(chunk, 0)) {
        ASSERT(reader->progress() >= 0 && reader->progress() <= 1);
        int new_percent = int(100 * reader->progress());
        if (new_percent != percent) {
            percent = new_percent;
            std::cerr << percent << '%' << std::endl;
        }
    }
    ASSERT_EQUALS(reader->progress(), 1);
    std::cerr << "Done" << std::endl;
}

int main(int argc, char** argv) {
    if (argc == 3) {
        if (!strcmp(argv[1], "-read")) {
            test_read(argv[2]);
            return 0;
        }
    }
    return TestManager::get().main(argc, argv);
}
