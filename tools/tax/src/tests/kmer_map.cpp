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
