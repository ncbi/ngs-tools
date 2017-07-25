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

#if 0
TEST(seq_transform_1) {
    int KMER_LEN = 32;
    typedef uint64_t hash_t;
    auto hash = Hash<hash_t>::hash_of("AAACTCTCGAGCACCTGCCGCTCGGGGAGGCC");
    seq_transform<hash_t>::for_all_1_char_variations_do(hash, KMER_LEN, [&](hash_t hash)
                                                        {
                                                            cout << Hash<hash_t>::str_from_hash(hash, KMER_LEN) << endl;
                                                            return true;
                                                        });
}
#endif

TEST(seq_transform_2) 	{
    int KMER_LEN = 32;
    typedef uint64_t hash_t;
    auto hash = Hash<hash_t>::hash_of("AAACTCTCGAGCACCTGCCGCTCGGGGAGGCC");
    ASSERT_EQUALS(hash, 115348763461549301ULL);
    ASSERT_EQUALS(Hash<hash_t>::str_from_hash(hash, KMER_LEN), string("AAACTCTCGAGCACCTGCCGCTCGGGGAGGCC"));
    auto hash_rev_compl = Hash<hash_t>::hash_of("GGCCTCCCCGAGCGGCAGGTGCTCGAGAGTTT");
    ASSERT_EQUALS(hash_rev_compl, 17696177292584799466ULL);
    ASSERT_EQUALS(Hash<hash_t>::str_from_hash(hash_rev_compl, KMER_LEN), string("GGCCTCCCCGAGCGGCAGGTGCTCGAGAGTTT"));

    ASSERT_EQUALS(hash_rev_compl, seq_transform<hash_t>::to_rev_complement(hash, KMER_LEN));
    ASSERT_EQUALS(seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN), hash);
    ASSERT_EQUALS(seq_transform<hash_t>::min_hash_variant(hash_rev_compl, KMER_LEN), hash);
}


TEST(seq_transform_3) {
    //	static hash_t bit_reverse(hash_t hash, int kmer_len)
	ASSERT_EQUALS(seq_transform<uint64_t>::bit_reverse(0, 1), 0);
	ASSERT_EQUALS(seq_transform<uint64_t>::bit_reverse(1, 1), 1);
	ASSERT_EQUALS(seq_transform<uint64_t>::bit_reverse(1, 2), 1 << 2);
	ASSERT_EQUALS(seq_transform<uint64_t>::bit_reverse(44, 3), 14); // 101100

	ASSERT_EQUALS(seq_transform<uint64_t>::bit_reverse(41191556913, 18), 20620627302);
	// 100110010111001101010100011100110001
	// 010011001101000101011100110101100110

	ASSERT_EQUALS(seq_transform<uint64_t>::to_rev_complement(0, 1), 2);
	ASSERT_EQUALS(seq_transform<uint64_t>::to_rev_complement(1, 1), 3);
	ASSERT_EQUALS(seq_transform<uint64_t>::to_rev_complement(1, 2), 14);
	ASSERT_EQUALS(seq_transform<uint64_t>::to_rev_complement(41191556913, 18), 61870139340);
	// 100110010111001101010100011100110001
	// 010011001101000101011100110101100110 - rev
	//  3 2 1 2 1 3 2 3 3 3 1 2 1 3 3 0 3 0
	// 111001100111101111110110011111001100


	auto hash1 = Hash<unsigned int>::hash_of("ACTG");
	ASSERT_EQUALS(hash1, 27);
	auto hash2 = Hash<unsigned int>::hash_of("TGAC"); // compl
	ASSERT_EQUALS(hash2, 177);
	auto hash3 = Hash<unsigned int>::hash_of("CAGT"); // rev compl
	ASSERT_EQUALS(hash3, 78);
	auto hash4 = Hash<unsigned int>::hash_of("GTCA"); // rev
	ASSERT_EQUALS(hash4, 228);
	ASSERT_EQUALS(seq_transform<unsigned int>::min_hash_variant(hash1, 4), 27);
	ASSERT_EQUALS(seq_transform<unsigned int>::min_hash_variant(hash2, 4), 177);
	ASSERT_EQUALS(seq_transform<unsigned int>::min_hash_variant(hash3, 4), 27);
	ASSERT_EQUALS(seq_transform<unsigned int>::min_hash_variant(hash4, 4), 177);

	ASSERT_EQUALS(seq_transform<uint64_t>::apply_transformation(Hash<uint64_t>::hash_of("TAAAAAAAAACTGGGG"), 16, false, false), Hash<uint64_t>::hash_of("TAAAAAAAAACTGGGG"));
	ASSERT_EQUALS(seq_transform<uint64_t>::apply_transformation(Hash<uint64_t>::hash_of("TAAAAAAAAACTGGGG"), 16, true, false), Hash<uint64_t>::hash_of("GGGGTCAAAAAAAAAT"));
	ASSERT_EQUALS(seq_transform<uint64_t>::apply_transformation(Hash<uint64_t>::hash_of("TAAAAAAAAACTGGGG"), 16, true, true), Hash<uint64_t>::hash_of("CCCCAGTTTTTTTTTA"));
	ASSERT_EQUALS(seq_transform<uint64_t>::apply_transformation(Hash<uint64_t>::hash_of("TAAAAAAAAACTGGGG"), 16, false, true), Hash<uint64_t>::hash_of("ATTTTTTTTTGACCCC"));

    {
        string s1 = "TAAAAAAAAACTGGGG";
        seq_transform_actg::apply_transformation(s1, false, false);
	    ASSERT_EQUALS(s1, "TAAAAAAAAACTGGGG");
    }

    {
        string s2 = "TAAAAAAAAACTGGGG";
        seq_transform_actg::apply_transformation(s2, true, false);
	    ASSERT_EQUALS(s2, "GGGGTCAAAAAAAAAT");
    }

    {
        string s3 = "TAAAAAAAAACTGGGG";
        seq_transform_actg::apply_transformation(s3, true, true);
	    ASSERT_EQUALS(s3, "CCCCAGTTTTTTTTTA");
    }

    {
        string s4 = "TAAAAAAAAACTGGGG";
        seq_transform_actg::apply_transformation(s4, false, true);
	    ASSERT_EQUALS(s4, "ATTTTTTTTTGACCCC");
    }

//	ASSERT_EQUALS(seq_transform_actg::apply_transformation("TAAAAAAAAACTGGGG", true, false), "GGGGTCAAAAAAAAAT");
//	ASSERT_EQUALS(seq_transform_actg::apply_transformation("TAAAAAAAAACTGGGG", true, true), "CCCCAGTTTTTTTTTA");
//	ASSERT_EQUALS(seq_transform_actg::apply_transformation("TAAAAAAAAACTGGGG", false, true), "ATTTTTTTTTGACCCC");
}

TEST_MAIN();
