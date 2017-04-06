#include "tests.h"

TEST(kmer_map) {
    const int KMER_LEN = 12;
    const int KMER_BUCKETS = 256;

    typedef KmerMap<unsigned int, KMER_LEN, KMER_BUCKETS> TestKmerMap;
    TestKmerMap kmer_map;

    int hash = Hash<unsigned int>::hash_of("ACTGACTGACTG");
    kmer_map.add(hash);

    equal(kmer_map.coverage_of(hash), 1);
    equal(kmer_map.size(), 1);
    equal(kmer_map.originally_complement(hash), false);
    equal(kmer_map.originally_reverse(hash), false);

    //for (int i=0; i<TestKmerMap::Count::MAX_COUNT + 100; i++)
    //	kmer_map.add(hash);

    //equal(kmer_map.coverage_of(hash), TestKmerMap::Count::MAX_COUNT);
    //equal(kmer_map.originally_complement(hash), false);
    //equal(kmer_map.originally_reverse(hash), false);

    kmer_map.add(Hash<unsigned int>::hash_of("CAGTCAGTCAGT"));
    equal(kmer_map.coverage_of(hash), 2);
    equal(kmer_map.size(), 1);
    equal(kmer_map.originally_complement(Hash<unsigned int>::hash_of("CAGTCAGTCAGT")), true);
    equal(kmer_map.originally_reverse(Hash<unsigned int>::hash_of("CAGTCAGTCAGT")), true);

    kmer_map.add(Hash<unsigned int>::hash_of("AAAAAAAAAAAC"));
    equal(kmer_map.size(), 2);

    kmer_map.add(Hash<unsigned int>::hash_of("AAAAAAAAAACC"));
    equal(kmer_map.size(), 3);

    kmer_map.add(Hash<unsigned int>::hash_of("AAATAAAAAACC"));
    equal(kmer_map.size(), 4);

    kmer_map.add(Hash<unsigned int>::hash_of("GTTTTTTTTTTT"));
    equal(kmer_map.size(), 4);
    kmer_map.optimize();
    equal(kmer_map.size(), 2);
    equal(kmer_map.coverage_of(Hash<unsigned int>::hash_of("ACTGACTGACTG")), 2);
    kmer_map.remove(Hash<unsigned int>::hash_of("ACTGACTGACTG"));
    equal(kmer_map.coverage_of(Hash<unsigned int>::hash_of("ACTGACTGACTG")), 0);
    kmer_map.restore(Hash<unsigned int>::hash_of("ACTGACTGACTG"));
    equal(kmer_map.coverage_of(Hash<unsigned int>::hash_of("ACTGACTGACTG")), 2);

    equal(kmer_map.coverage_of(Hash<unsigned int>::hash_of("GTTTTTTTTTTT")), 2);
    kmer_map.remove(Hash<unsigned int>::hash_of("GTTTTTTTTTTT"));
    equal(kmer_map.coverage_of(Hash<unsigned int>::hash_of("GTTTTTTTTTTT")), 0);
    kmer_map.restore(Hash<unsigned int>::hash_of("GTTTTTTTTTTT"));
    equal(kmer_map.coverage_of(Hash<unsigned int>::hash_of("GTTTTTTTTTTT")), 2);

    equal(kmer_map.originally_complement(Hash<unsigned int>::hash_of("AAAAAAAAAAAC")), false);
    equal(kmer_map.originally_reverse(Hash<unsigned int>::hash_of("AAAAAAAAAAAC")), false);

    equal(kmer_map.originally_complement(Hash<unsigned int>::hash_of("GTTTTTTTTTTT")), true);
    equal(kmer_map.originally_reverse(Hash<unsigned int>::hash_of("GTTTTTTTTTTT")), true);

    equal(kmer_map.originally_complement(Hash<unsigned int>::hash_of("TTTTTTTTTTTG")), true);
    equal(kmer_map.originally_reverse(Hash<unsigned int>::hash_of("TTTTTTTTTTTG")), true);
}

TEST(kmer_map64) {
    KmerMap64 kmers64;
    auto hash = Hash<KmerMap64::hash_t>::hash_of("CTATACGATCGAGGTCATCGACCTGATGAAGGACCCGGCCTTGGCGCAGCGCGACCAGATCGTC");
    kmers64.add(hash);
    equal(kmers64.coverage_of(hash), 1);
    equal(kmers64.size(), 1);
}

TEST_MAIN();
