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
    auto hash_rev_compl = Hash<hash_t>::hash_of("GGCCTCCCCGAGCGGCAGGTGCTCGAGAGTTT");

    equal(seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN), hash);
    equal(seq_transform<hash_t>::min_hash_variant(hash_rev_compl, KMER_LEN), hash);
}


TEST(seq_transform_3) {
    //	static hash_t bit_reverse(hash_t hash, int kmer_len)
	equal(seq_transform<uint64_t>::bit_reverse(0, 1), 0);
	equal(seq_transform<uint64_t>::bit_reverse(1, 1), 1);
	equal(seq_transform<uint64_t>::bit_reverse(1, 2), 1 << 2);
	equal(seq_transform<uint64_t>::bit_reverse(44, 3), 14); // 101100
	
	equal(seq_transform<uint64_t>::bit_reverse(41191556913, 18), 20620627302);
	// 100110010111001101010100011100110001
	// 010011001101000101011100110101100110

	equal(seq_transform<uint64_t>::to_rev_complement(0, 1), 2);
	equal(seq_transform<uint64_t>::to_rev_complement(1, 1), 3);
	equal(seq_transform<uint64_t>::to_rev_complement(1, 2), 14);
	equal(seq_transform<uint64_t>::to_rev_complement(41191556913, 18), 61870139340);
	// 100110010111001101010100011100110001
	// 010011001101000101011100110101100110 - rev
	//  3 2 1 2 1 3 2 3 3 3 1 2 1 3 3 0 3 0
	// 111001100111101111110110011111001100


	auto hash1 = Hash<unsigned int>::hash_of("ACTG");
	equal(hash1, 27);
	auto hash2 = Hash<unsigned int>::hash_of("TGAC"); // compl
	equal(hash2, 177);
	auto hash3 = Hash<unsigned int>::hash_of("CAGT"); // rev compl
	equal(hash3, 78);
	auto hash4 = Hash<unsigned int>::hash_of("GTCA"); // rev
	equal(hash4, 228);
	equal(seq_transform<unsigned int>::min_hash_variant(hash1, 4), 27);
	equal(seq_transform<unsigned int>::min_hash_variant(hash2, 4), 177);
	equal(seq_transform<unsigned int>::min_hash_variant(hash3, 4), 27);
	equal(seq_transform<unsigned int>::min_hash_variant(hash4, 4), 177);

	equal(seq_transform<uint64_t>::apply_transformation(Hash<uint64_t>::hash_of("TAAAAAAAAACTGGGG"), 16, false, false), Hash<uint64_t>::hash_of("TAAAAAAAAACTGGGG"));
	equal(seq_transform<uint64_t>::apply_transformation(Hash<uint64_t>::hash_of("TAAAAAAAAACTGGGG"), 16, true, false), Hash<uint64_t>::hash_of("GGGGTCAAAAAAAAAT"));
	equal(seq_transform<uint64_t>::apply_transformation(Hash<uint64_t>::hash_of("TAAAAAAAAACTGGGG"), 16, true, true), Hash<uint64_t>::hash_of("CCCCAGTTTTTTTTTA"));
	equal(seq_transform<uint64_t>::apply_transformation(Hash<uint64_t>::hash_of("TAAAAAAAAACTGGGG"), 16, false, true), Hash<uint64_t>::hash_of("ATTTTTTTTTGACCCC"));

    {
        string s1 = "TAAAAAAAAACTGGGG";
        seq_transform_actg::apply_transformation(s1, false, false);
	    equal(s1, "TAAAAAAAAACTGGGG");
    }

    {
        string s2 = "TAAAAAAAAACTGGGG";
        seq_transform_actg::apply_transformation(s2, true, false);
	    equal(s2, "GGGGTCAAAAAAAAAT");
    }

    {
        string s3 = "TAAAAAAAAACTGGGG";
        seq_transform_actg::apply_transformation(s3, true, true);
	    equal(s3, "CCCCAGTTTTTTTTTA");
    }

    {
        string s4 = "TAAAAAAAAACTGGGG";
        seq_transform_actg::apply_transformation(s4, false, true);
	    equal(s4, "ATTTTTTTTTGACCCC");
    }

//	equal(seq_transform_actg::apply_transformation("TAAAAAAAAACTGGGG", true, false), "GGGGTCAAAAAAAAAT");
//	equal(seq_transform_actg::apply_transformation("TAAAAAAAAACTGGGG", true, true), "CCCCAGTTTTTTTTTA");
//	equal(seq_transform_actg::apply_transformation("TAAAAAAAAACTGGGG", false, true), "ATTTTTTTTTGACCCC");
}

TEST_MAIN();
