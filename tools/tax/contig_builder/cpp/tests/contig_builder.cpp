#include "tests.h"

TEST(contig_builder_1) {
    std::vector<string> reads;
    reads.push_back("CCCCCAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTGGGGGGAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCTTTTTTAGGACTAA");
    reads.push_back("CCCCCAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTGGGGGGAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCTTTTTTATTACGGGGGGGGG");
                reads.push_back("AAAAAAAAATTTTTTTTTTTTTTTTTTTGGGGGGAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCTTTTTTAGGACTGGGTTT");

    typedef KmerMap64 TestKmerMap;
    TestKmerMap kmer_map;
    build_test_map<TestKmerMap>(reads, kmer_map);

                            string seq = "TTTTTTTTTTTTTTTTTTTGGGGGGAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCTTTTTT";

    auto hash = Hash<TestKmerMap::hash_t>::hash_of(seq);
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'A');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'G');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'G');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'A');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 0);

    equal(ContigBuilder::get_next_contig<TestKmerMap>(kmer_map, Hash<TestKmerMap::hash_t>::hash_of(seq)), "CCCCCAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTGGGGGGAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCTTTTTTAGGACT");
}

TEST(contig_builder_2) {
    std::vector<string> reads;
    reads.push_back("CCCCCAAAAAAAAAAA");
    reads.push_back("CCCCCAAAAAAAAAAA");

    typedef KmerMap16 TestKmerMap;
    TestKmerMap kmer_map;
    build_test_map<TestKmerMap>(reads, kmer_map);

    string seq = "CCCCCAAAAAAAAAAA";
    auto hash = Hash<TestKmerMap::hash_t>::hash_of(seq);
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 0);

    equal(ContigBuilder::get_next_contig<TestKmerMap>(kmer_map, Hash<TestKmerMap::hash_t>::hash_of(seq)), "CCCCCAAAAAAAAAAA");
}

TEST(contig_builder_3) {
    std::vector<string> reads;
    reads.push_back("CCCCCAAAAAAAAAAAT");
    reads.push_back("CCCCCAAAAAAAAAAAT");

    typedef KmerMap16 TestKmerMap;
    TestKmerMap kmer_map;
    build_test_map<TestKmerMap>(reads, kmer_map);

    string seq = "CCCCCAAAAAAAAAAA";
    auto hash = Hash<TestKmerMap::hash_t>::hash_of(seq);
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 0);
    equal(ContigBuilder::get_next_contig<TestKmerMap>(kmer_map, Hash<TestKmerMap::hash_t>::hash_of(seq)), "CCCCCAAAAAAAAAAAT");
}

TEST(contig_builder_4)
{
    std::vector<string> reads;
    reads.push_back("CCCCCAAGAAAAAAAATTTTTTTGTTTTTTATGGGGCTGGGGGGGGGGCCCCCCCCCCCCCCC");
    reads.push_back("CCCCCAAGAAAAAAAATTTATTTGTTTTTTATGGGGCTGGGAAAAAAACCCC");
    reads.push_back("GCCCCAAGAAAAAAAATTTTTTTGTTTTTTATGGGGCTGGGAAAAAAACCCCCCCCCCCCCCAC");
    reads.push_back("GGGGCTGGGAAAAAAACCCCCCCCCCCCCCCC");

    typedef KmerMap16 TestKmerMap;
    TestKmerMap kmer_map;
    build_test_map<TestKmerMap>(reads, kmer_map);

    string seq = "CCCCCAAGAAAAAAAA";
    auto hash = Hash<TestKmerMap::hash_t>::hash_of(seq);
    cout << "t block" << endl;
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'G');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'A');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');

    cout << "g block" << endl;
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'G');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'G');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'G');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'G');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'T');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'G');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'G');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'G');

    cout << "a block" << endl;
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'A');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'A');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'A');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'A');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'A');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'A');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'A');

    cout << "c block" << endl;
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 'C');

    equal(ContigBuilder::get_next_contig<TestKmerMap>(kmer_map, Hash<TestKmerMap::hash_t>::hash_of("AAAAAATTTTTTTGTTT")), "CCCCCAAGAAAAAAAATTTTTTTGTTTTTTATGGGGCTGGGAAAAAAACCCCCCCCCCCCCC");
    equal(ContigBuilder::choose_next_letter<TestKmerMap>(kmer_map, &hash), 0);
}

TEST_MAIN();
