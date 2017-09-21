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

#include <map>
#include <memory>

#include "tests.h"


TEST(hash_1) {
    string seq = "TCTCCGAGCCCACGAGAC";
    auto hash = Hash<uint64_t>::hash_of(&seq[0], seq.length());
    ASSERT_EQUALS(hash, 41191556913);

    string s = Hash<uint64_t>::str_from_hash(hash, 18);
    ASSERT_EQUALS(s, seq);
}

TEST(hash_left_right) {
    string seq_left  = "TCTCCGAGCCCACGAG";
    string seq_right = "GTCAGTCAGTCAAAAA";

    string seq_full = seq_left + seq_right;

    auto hash_left = Hash<unsigned int>::hash_of(seq_left);
    auto hash_right = Hash<unsigned int>::hash_of(seq_right);
    auto hash_full = Hash<uint64_t>::hash_of(seq_full);

    ASSERT_EQUALS(Hash<uint64_t>::left<unsigned int>(hash_full), hash_left);
    ASSERT_EQUALS(Hash<uint64_t>::right<unsigned int>(hash_full), hash_right);
}

TEST(hash_left_right_long) {
#if defined ( __GNUC__ )
    string seq_left  = "TATACGATCGAGGTCATCGACCTGATGAAGGA";
    string seq_right = "CCCGGCCTTGGCGCAGCGCGACCAGATCGTCG";

    string seq_full = seq_left + seq_right;

    auto hash_left = Hash<uint64_t>::hash_of(seq_left);
    auto hash_right = Hash<uint64_t>::hash_of(seq_right);
    auto hash_full = Hash<__uint128_t>::hash_of(seq_full);

    ASSERT_EQUALS(Hash<__uint128_t>::left<uint64_t>(hash_full), hash_left);
    ASSERT_EQUALS(Hash<__uint128_t>::right<uint64_t>(hash_full), hash_right);
#endif
}

TEST(hash_next) {
#if defined ( __GNUC__ )
    string read = "CTATACGATCGAGGTCATCGACCTGATGAAGGACCCGGCCTTGGCGCAGCGCGACCAGATCGTCGCGATCCCGACGCTG";
    auto hash0 = Hash<__uint128_t>::hash_of(&read[0], 64);
    auto hash1 = Hash<__uint128_t>::hash_of(&read[1], 64);
    auto hash2 = Hash<__uint128_t>::hash_of(&read[2], 64);
    auto hash3 = Hash<__uint128_t>::hash_of(&read[3], 64);
    auto hash4 = Hash<__uint128_t>::hash_of(&read[4], 64);
    auto hash5 = Hash<__uint128_t>::hash_of(&read[5], 64);
    auto hash = Hash<__uint128_t>::hash_next(&read[1], hash0, 64);
    ASSERT_EQUALS(hash == hash1, true);

    hash = Hash<__uint128_t>::hash_next(&read[2], hash, 64);
    ASSERT_EQUALS(hash == hash2, true);

    hash = Hash<__uint128_t>::hash_next(&read[3], hash, 64);
    ASSERT_EQUALS(hash == hash3, true);

    hash = Hash<__uint128_t>::hash_next(&read[4], hash, 64);
    ASSERT_EQUALS(hash == hash4, true);

    hash = Hash<__uint128_t>::hash_next(&read[5], hash, 64);
    ASSERT_EQUALS(hash == hash5, true);
#endif
}

TEST(hash_to_and_from) 	{
    int KMER_LEN = 32;
    typedef uint64_t hash_t;
    auto hash = Hash<hash_t>::hash_of("AAACTCTCGAGCACCTGCCGCTCGGGGAGGCC");
    ASSERT_EQUALS(hash, 115348763461549301ULL);
    ASSERT_EQUALS(Hash<hash_t>::str_from_hash(hash, KMER_LEN), string("AAACTCTCGAGCACCTGCCGCTCGGGGAGGCC"));
    auto hash_rev_compl = Hash<hash_t>::hash_of("GGCCTCCCCGAGCGGCAGGTGCTCGAGAGTTT");
    ASSERT_EQUALS(hash_rev_compl, 17696177292584799466ULL);
    ASSERT_EQUALS(Hash<hash_t>::str_from_hash(hash_rev_compl, KMER_LEN), string("GGCCTCCCCGAGCGGCAGGTGCTCGAGAGTTT"));
}

TEST(hash_for_all) {
    cout << "for_all_hashes_do test" << endl;
    string seq = "TCTCCGAGCCCACGAGAC";
    std::map<unsigned int, int> counter;

    Hash<unsigned int>::for_all_hashes_do(seq, 8, [&](unsigned int hash)
        {
            counter[hash]++;
            return true;
        });

    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("TCTCCGAG")], 1);
    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("CTCCGAGC")], 1);
    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("TCCGAGCC")], 1);
    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("CCGAGCCC")], 1);
    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("CGAGCCCA")], 1);
    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("GAGCCCAC")], 1);
    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("GAGCCCAC")], 1);
    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("AGCCCACG")], 1);
    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("GCCCACGA")], 1);
    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("CCCACGAG")], 1);
    ASSERT_EQUALS(counter[Hash<unsigned int>::hash_of("CCACGAGA")], 1);
}

TEST_MAIN();
