#include "tests.h"

TEST(kmer_map_short) {
    const int KMER_LEN = 6;
    const int KMER_BUCKETS = 1;

    typedef KmerMap<unsigned int, KMER_LEN, KMER_BUCKETS> TestKmerMap;
    TestKmerMap kmer_map;

    kmer_map.add(Hash<unsigned int>::hash_of("ACTGAC"));
    kmer_map.add(Hash<unsigned int>::hash_of("GTCAGT"));

    kmer_map.add(Hash<unsigned int>::hash_of("AAAATT"));
    kmer_map.add(Hash<unsigned int>::hash_of("AAAATT"));
    kmer_map.add(Hash<unsigned int>::hash_of("AATTTT"));

    kmer_map.add(Hash<unsigned int>::hash_of("AGGGGG"));

    Begins<TestKmerMap> begins(kmer_map);
    unsigned int hash = 0;
    equal(begins.next(&hash), true);
    equal(Hash<unsigned int>::str_from_hash(hash, KMER_LEN), "AAAATT");
    equal(begins.next(&hash), true);
    equal(Hash<unsigned int>::str_from_hash(hash, KMER_LEN), "ACTGAC");
    equal(begins.next(&hash), true);
    equal(Hash<unsigned int>::str_from_hash(hash, KMER_LEN), "AGGGGG");
    equal(begins.next(&hash), false);
}

TEST(kmer_map_long) {
    typedef KmerMap64 TestKmerMap;
    TestKmerMap kmer_map;

    kmer_map.add(Hash<TestKmerMap::hash_t>::hash_of("CTATACGATCGAGGTCATCGACCTGATGAAGGACCCGGCCTTGGCGCAGCGCGACCAGATCGTC"));
    kmer_map.add(Hash<TestKmerMap::hash_t>::hash_of("ATATGCTAGCTCCAGTAGCTGGACTACTTCCTGGGCCGGAACCGCGTCGCGCTGGTCTAGCAGG"));
    kmer_map.add(Hash<TestKmerMap::hash_t>::hash_of("ATATGCTAGCTCCAGTAGCTGGACTACTTCCTGGGCCGGAACCGCGTCGCGCTGGTCTAGCAGG"));

    Begins<TestKmerMap> begins(kmer_map);
    TestKmerMap::hash_t hash = 0;
    equal(begins.next(&hash), true);
    equal(Hash<TestKmerMap::hash_t>::str_from_hash(hash, 64), "ATATGCTAGCTCCAGTAGCTGGACTACTTCCTGGGCCGGAACCGCGTCGCGCTGGTCTAGCAGG");
    equal(begins.next(&hash), true);
    equal(Hash<TestKmerMap::hash_t>::str_from_hash(hash, 64), "CTATACGATCGAGGTCATCGACCTGATGAAGGACCCGGCCTTGGCGCAGCGCGACCAGATCGTC");
    equal(begins.next(&hash), false);
}

TEST_MAIN();
