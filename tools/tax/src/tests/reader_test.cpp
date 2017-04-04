#include "tests.h"

#include "reader.h"
#include "vdb_reader.h"
#include "fasta_reader.h"
#include "mt_reader.h"
#include "aux_reader.h"

class DummyReader: public Reader {
    size_t idx;
    std::vector<std::string> reads;
public:
    DummyReader(const std::vector<std::string>& reads) : reads(reads), idx(0) {}

    virtual SourceStats stats() const { return SourceStats(reads.size()); }
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

TEST(vdb_reader) {
    { // all
        VdbReader reader("tests/data/SRR1068106");
        Reader::Fragment f;
    
        ASSERT(reader.read(&f));
        ASSERT_EQUALS(f.spotid, "1");
        ASSERT_EQUALS(f.bases, "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCACTTCTGAA");
    
        ASSERT(reader.read(nullptr));
        ASSERT(reader.read(&f));
    
        ASSERT(reader.read(&f));
        ASSERT_EQUALS(f.spotid, "4");
        ASSERT_EQUALS(f.bases, "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCAACTTCCTGAATGTTGACCTCG");

        size_t got_spots = 4;
        while (reader.read(nullptr)) {
            ++got_spots;
        }
        ASSERT_EQUALS(got_spots, 226);
        ASSERT(reader.stats() == Reader::SourceStats(got_spots));
        ASSERT_EQUALS(reader.progress(), 1);
    }

    // unaligned only
    {
        VdbReader reader("tests/data/SRR1068106", false, true);
        Reader::Fragment f;
    
        ASSERT(reader.read(&f));
        ASSERT_EQUALS(f.spotid, "4");
        ASSERT_EQUALS(f.bases, "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCAACTTCCTGAATGTTGACCTCG");
    
        ASSERT(reader.read(nullptr));
        ASSERT(reader.read(&f));
    
        ASSERT(reader.read(&f));
        ASSERT_EQUALS(f.spotid, "33");
        ASSERT_EQUALS(f.bases, "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCAACTTCCTGAA");

        size_t got_spots = 4;
        while (reader.read(nullptr)) {
            ++got_spots;
        }
        ASSERT_EQUALS(got_spots, 12);
        auto stats = reader.stats();
        ASSERT_EQUALS(stats.spot_count, 12);
        ASSERT_EQUALS(stats.expected_spot_count, 12);
        ASSERT_EQUALS(stats.frags_per_spot, 1);
        ASSERT_EQUALS(reader.progress(), 1);
    }
}

TEST(vdb_reader_paired) {
    VdbReader reader("tests/data/ERR333883");
    auto fragments = read_all(&reader);
    ASSERT_EQUALS(fragments.size(), 6);
    ASSERT_EQUALS(fragments[0].spotid, "1");
    ASSERT_EQUALS(fragments[0].bases, "TTGGGACCTTAGCTGGCGGTCTGGGTTGTT");
    ASSERT_EQUALS(fragments[1].spotid, "1");
    ASSERT_EQUALS(fragments[1].bases, "GTGTCATGGGGAGCTGGTATCTCCAGGCGC");
    ASSERT_EQUALS(fragments[2].spotid, "2");
    ASSERT_EQUALS(fragments[2].bases, "CTCCTGACTGACCGATAGTGAACCAGTACC");
    ASSERT_EQUALS(fragments[3].spotid, "2");
    ASSERT_EQUALS(fragments[3].bases, "TAGTGTTACCCAACCTTCAACCTGCCCATG");
    ASSERT_EQUALS(fragments[4].spotid, "3");
    ASSERT_EQUALS(fragments[4].bases, "ATTAACGCAAACCGGCGACAGTTCTCGCAA");
    ASSERT_EQUALS(fragments[5].spotid, "3");
    ASSERT_EQUALS(fragments[5].bases, "ACGTCCACAGAAACCCGTTTCCCGCCGGTC");
    ASSERT(reader.stats() == Reader::SourceStats(3, 2));
}

TEST(fasta_reader) {
    FastaReader reader("tests/data/SRR1068106.fasta");
    Reader::Fragment f;
    
    ASSERT(reader.read(&f));
    ASSERT_EQUALS(f.spotid, "SRR1068106.1 1 length=499");
    ASSERT_EQUALS(f.bases, "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCACTTCTGAA");
    
    ASSERT(reader.read(nullptr));
    ASSERT(reader.read(&f));
    
    ASSERT(reader.read(&f));
    ASSERT_EQUALS(f.spotid, "SRR1068106.4 4 length=512");
    ASSERT_EQUALS(f.bases, "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCAACTTCCTGAATGTTGACCTCG");

    size_t got_spots = 4;
    while (reader.read(nullptr)) {
        ++got_spots;
    }
    ASSERT_EQUALS(got_spots, 226);
    ASSERT(reader.stats() == Reader::SourceStats(got_spots));
    ASSERT_EQUALS(reader.progress(), 1);

    auto unix = Helper<FastaReader>::read_all("tests/data/SRR1068106.fasta");
    auto dos = Helper<FastaReader>::read_all("tests/data/SRR1068106.fasta.dos");
    ASSERT(unix == dos);
}

TEST(vdb_fasta_equal) {
    {
        auto vdb = Helper<VdbReader>::read_all_bases("tests/data/SRR1068106");
        auto fasta = Helper<FastaReader>::read_all_bases("tests/data/SRR1068106.fasta");
        ASSERT(fasta == vdb);
    }
    {
        auto vdb = Helper<VdbReader>::read_all_bases("tests/data/SRR1068106", false, true);
        auto fasta = Helper<FastaReader>::read_all_bases("tests/data/SRR1068106.unaligned.fasta");
        ASSERT(fasta == vdb);
    }
}

TEST(vdb_qualities_reader) {
    auto no_qual = Helper<VdbReader>::read_all_bases("tests/data/SRR1068106", false);
    auto with_qual = Helper<VdbReader>::read_all_bases("tests/data/SRR1068106", true);
    ASSERT_EQUALS(no_qual.size(), with_qual.size());
    ASSERT(no_qual != with_qual);
    ASSERT_EQUALS(with_qual[4], "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACCTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGT!ACGTCGAGCTTCCATAGCGTAGTAGT");
    for (size_t i = 0; i < no_qual.size(); ++i) {
        auto nq = no_qual[i];
        auto wq = with_qual[i];
        std::replace( wq.begin(), wq.end(), '!', 'N');
        ASSERT_EQUALS(nq, wq);
    }
}

void test_aligned_vdb_reader(const char* path, bool read_qualities) {
    auto vdb = Helper<VdbReader>::read_all(path, read_qualities);
    auto aligned = Helper<AlignedVdbReader>::read_all(path, read_qualities);
    std::sort(vdb.begin(), vdb.end(), FragmentSort());
    std::sort(aligned.begin(), aligned.end(), FragmentSort());
    ASSERT(vdb == aligned);
    ASSERT(VdbReader(path, read_qualities).stats() == AlignedVdbReader(path, read_qualities).stats());
}
TEST(aligned_vdb_reader) {
    test_aligned_vdb_reader("tests/data/SRR1068106", false); // has aligned and unaligned single reads
    test_aligned_vdb_reader("tests/data/SRR1068106", true);
    test_aligned_vdb_reader("tests/data/ERR333883", false); // has aligned, unaligned and partially aligned _paired_ reads
    test_aligned_vdb_reader("tests/data/ERR333883", true); 
}

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
    test_mt_reader<FastaReader>("fasta", "tests/data/SRR1068106.fasta");
    test_mt_reader<VdbReader>("vdb", "tests/data/SRR1068106");
    test_mt_reader<AlignedVdbReader>("aligned vdb", "tests/data/SRR1068106");
}

TEST(filtering_reader) {
    std::vector<std::string> source = {"A", "C", "T", "G"};
    auto filter1 = [](const std::string& spotid) { return spotid != "2"; };
    auto filter2 = [](const std::string& spotid) { return spotid != "2" && spotid != "3"; };
    auto filter3 = [](const std::string& spotid) { return false; };
    auto filter4 = [](const std::string& spotid) { return true; };
    std::vector<std::string> expected1 = {"A", "T", "G"};
    std::vector<std::string> expected2 = {"A", "G"};
    std::vector<std::string> expected3 = {};
    std::vector<std::string> expected4 = {"A", "C", "T", "G"};
    std::vector<std::string> result1 = Helper<FilteringReader<DummyReader, CallableSpotFilter<decltype(filter1)> > >::read_all_bases(filter1, source);
    std::vector<std::string> result2 = Helper<FilteringReader<DummyReader, CallableSpotFilter<decltype(filter2)> > >::read_all_bases(filter2, source);
    std::vector<std::string> result3 = Helper<FilteringReader<DummyReader, CallableSpotFilter<decltype(filter3)> > >::read_all_bases(filter3, source);
    std::vector<std::string> result4 = Helper<FilteringReader<DummyReader, CallableSpotFilter<decltype(filter4)> > >::read_all_bases(filter4, source);
    ASSERT(result1 == expected1);
    ASSERT(result2 == expected2);
    ASSERT(result3 == expected3);
    ASSERT(result4 == expected4);

    {
        std::vector<Reader::Fragment> file_expected = {
            {"1", "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCACTTCTGAA"},
            {"20", "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCAACTTCT"},
            {"35", "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCGTTAAACCCCA"},
            {"79", "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATTGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACGCCG"},
            {"123", "AAGTCGTAACAAGGTCTCCGTTGGTGAACCAGCGGAGGGATCATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCACTTGTTGCCTCGGCGGATCAGCCCGCTCCCGGTAAACGGGACGGCCCGCCAGAGGACCCCTAAACTCTGTTTCTATATGTAACTTCTGAGTAAAACCATAAATAAATCAAAACTTTCAACAACGGATCTCTTGGTTCTGGCATCGATGAAGAACGCAGCAAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCGCCAGTATTCTGGCGGGCATGCCTGTTCGAGCGTCATTTCAACCCTCAAGCACAGCTTGGTGTTGGGACTCGCGTTAATTCGCGTTCCTCAAATTGATCGGCGGTCACGTCGAGCTTCCATAGCGTAGTAGTAAAACCCTCGTTACTGGTAATCGTCGCGGCCACG"},
        };
        IncludeFileSpotFilter file_filter("tests/data/filter.spots");
        ASSERT_EQUALS(file_filter.expected_spot_count(), 5);
        FilteringReader<VdbReader, IncludeFileSpotFilter> reader(file_filter, "tests/data/SRR1068106");
        auto stats = reader.stats();
        auto file_result = read_all(&reader);
        ASSERT(file_result == file_expected);
        ASSERT_EQUALS(stats.spot_count, 226);
        ASSERT_EQUALS(stats.expected_spot_count, 5);
    }

    {
        ExcludeFileSpotFilter file_filter("tests/data/filter.spots");
        FilteringReader<VdbReader, ExcludeFileSpotFilter> reader(file_filter, "tests/data/SRR1068106");
        auto file_result = read_all(&reader);
        ASSERT_EQUALS(file_result.size(), 226 - 5);
    }
    
}

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

TEST(reader_factory) {
    { // simple
        auto vdb = read_all_bases(Reader::create("tests/data/SRR1068106"));
        auto fasta = read_all_bases(Reader::create("tests/data/SRR1068106.fasta"));
        std::sort(vdb.begin(), vdb.end());
        std::sort(fasta.begin(), fasta.end());
        ASSERT(fasta == vdb);
    }
    { // filtering
        Reader::Params params;
        auto all = read_all_bases(Reader::create("tests/data/SRR1068106", params));
        params.filter_file = "tests/data/filter.spots";
        auto filtered = read_all_bases(Reader::create("tests/data/SRR1068106", params));
        ASSERT(all != filtered);
        ASSERT_EQUALS(filtered.size(), 5);
    }

    // no way to test qualities with this file, all of them are already Ns =(
    /*{ // qualities
        auto no_qual = read_all_bases(Reader::create("tests/data/SRR1068106", nullptr, false));
        auto with_qual = read_all_bases(Reader::create("tests/data/SRR1068106", nullptr, true));
        ASSERT_EQUALS(no_qual.size(), with_qual.size());
        ASSERT(no_qual != with_qual);
    }*/
    { // split
        Reader::Params params;
        auto cut = read_all_bases(Reader::create("tests/data/SRR1068106", params));
        params.split_non_atgc = true;
        auto split = read_all_bases(Reader::create("tests/data/SRR1068106", params));
        ASSERT(cut.size() < split.size());
    }
    { // unaligned
        Reader::Params params;
        auto full = read_all_bases(Reader::create("tests/data/SRR1068106", params));
        params.unaligned_only = true;
        auto unaligned = read_all_bases(Reader::create("tests/data/SRR1068106", params));
        auto unaligned2 = read_all_bases(Reader::create("tests/data/SRR1068106.unaligned.fasta"));
        ASSERT(full != unaligned);
        ASSERT(unaligned == unaligned2);
    }
}

void print_stats(const char* path) {
    ngs::ReadCollection run = ncbi::NGS::openReadCollection(path);
    auto treads = run.getReadCount(ngs::Read::all);
    auto areads = run.getReadCount(ngs::Read::aligned);
    auto freads = run.getReadCount(ngs::Read::fullyAligned);
    auto preads = run.getReadCount(ngs::Read::partiallyAligned);
    auto ureads = run.getReadCount(ngs::Read::unaligned);
    auto talign = run.getAlignmentCount();
    auto palign = run.getAlignmentCount(ngs::Alignment::primaryAlignment);
    auto salign = run.getAlignmentCount(ngs::Alignment::secondaryAlignment);
    auto it = run.getReads(ngs::Read::all);
    it.nextRead();
    auto fcount = it.getNumFragments();
    std::cerr << "total reads " << treads << std::endl;
    std::cerr << "aligned reads " << areads << std::endl;
    std::cerr << "fully aligned reads " << freads << std::endl;
    std::cerr << "partially aligned reads " << preads << std::endl;
    std::cerr << "unaligned reads " << ureads << std::endl;
    std::cerr << "fragments per read " << fcount << std::endl;
    std::cerr << "total alignments " << talign << std::endl;
    std::cerr << "primary alignments " << palign << std::endl;
    std::cerr << "secondary alignments " << salign << std::endl;
}

void test_read(const char* path) {
    std::cerr << "Test reading: " << path << std::endl;
    auto reader = Reader::create(path);
    std::vector<Reader::Fragment> chunk;
    int percent = -1;
    while (reader->read_many(chunk)) {
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
        if (!strcmp(argv[1], "-stat")) {
            print_stats(argv[2]);
            return 0;
        }
        if (!strcmp(argv[1], "-read")) {
            test_read(argv[2]);
            return 0;
        }
    }
    return TestManager::get().main(argc, argv);
}
