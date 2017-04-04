#include "tests.h"

TEST(coverage) {
    std::vector<string> reads;
    reads.push_back("CCCCCAAGAAAAAAAATTTTTTAGTTTTTTATGGGGCTGGGGGGGGGGCCCCCCCCCCTAACC");
    reads.push_back("CCCCCAAGAAAAAAAATTTATTTGTTTTTTATGGGGCTGGGAAAAAAACCCC");
    reads.push_back("GCCCCAAGAAAAAAAATTTTTTTGTTTTTTATGGGGCTGGGAAAAAAACCCCCCCCCCTAACAC");
    reads.push_back("GGGGCTGGGAAAAAAACCCCCCCCCCTAACCC");

    typedef KmerMap16 TestKmerMap;
    TestKmerMap kmer_map;
    build_test_map<TestKmerMap>(reads, kmer_map);
    auto cov = Coverage<TestKmerMap>("CCCCCAAGAAAAAAAATTTTTTTGTTTTTTATGGGGCTGGGAAAAAAACCCCCCCCCCTAACG", kmer_map);
    equal(cov.size(), 48);
    equal(cov[0], 2);
    equal(cov[1], 3);
    equal(cov[2], 3);
    equal(cov[3], 3);
    equal(cov[4], 2);
    equal(cov[5], 2);
    equal(cov[6], 2);
    equal(cov[7], 0);
    equal(cov[8], 0);
    equal(cov[9], 0);
    equal(cov[10], 0);
    equal(cov[11], 0);
    equal(cov[12], 0);
    equal(cov[13], 0);
    equal(cov[14], 0);
    equal(cov[15], 0);
    equal(cov[16], 0);
    equal(cov[17], 0);
    equal(cov[18], 0);
    equal(cov[19], 0);
    equal(cov[20], 2);
    equal(cov[21], 2);
    equal(cov[22], 2);
    equal(cov[23], 3);
    equal(cov[24], 3);
    equal(cov[25], 3);
    equal(cov[26], 2);
    equal(cov[27], 2);
    equal(cov[28], 2);
    equal(cov[29], 2);
    equal(cov[30], 2);
    equal(cov[31], 2);
    equal(cov[31], 2);
    equal(cov[33], 3);
    equal(cov[34], 3);
    equal(cov[35], 3);
    equal(cov[36], 3);
    equal(cov[37], 2);
    equal(cov[38], 2);
    equal(cov[39], 2);
    equal(cov[40], 2);
    equal(cov[41], 2);
    equal(cov[42], 2);
    equal(cov[43], 2);
    equal(cov[44], 2);
    equal(cov[45], 2);
    equal(cov[46], 2);
    equal(cov[47], 0);
}

TEST_MAIN();
