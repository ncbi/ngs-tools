#ifndef SEQ_TRANSFORM_H_INCLUDED
#define SEQ_TRANSFORM_H_INCLUDED

#include <type_traits>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <assert.h>
//#include "hash.h"

class byte_map 
{
    uint8_t map[256];
public:
    template <typename Lambda>
    byte_map(Lambda fn)
    {
        for (int i = 0; i < 256; ++i)
            map[i] = fn(i);
    }

    uint8_t operator[] (uint8_t key) const { return map[key]; }
};

uint8_t reverse_byte_kmers(uint8_t src) {
    uint8_t dst;
    dst |= (src & 3);
	src >>= 2;
    dst <<= 2;
    dst |= (src & 3);
	src >>= 2;
    dst <<= 2;
    dst |= (src & 3);
	src >>= 2;
    dst <<= 2;
    dst |= (src & 3);
    return dst;
};

static const byte_map reverse_map(reverse_byte_kmers);

template <class hash_t>
struct seq_transform
{
	static hash_t min_hash_variant(hash_t hash, int kmer_len, bool *complement, bool *reverse)
	{
		auto min_hash = hash;

		hash = to_rev_complement(hash, kmer_len);
		if (hash < min_hash)
        {
            *complement = *reverse = true;
			return hash;
        }

        *complement = *reverse = false;
		return min_hash;
	}
    
	static hash_t min_hash_variant(hash_t hash, int kmer_len)
	{
		auto min_hash = hash;

		hash = to_rev_complement(hash, kmer_len);
		if (hash < min_hash)
			return hash;

		return min_hash;
	}
    
	static hash_t bit_reverse(hash_t hash, int kmer_len);

    static hash_t to_rev_complement(hash_t hash, int kmer_len);

	static hash_t apply_transformation(hash_t hash, int kmer_len, bool need_reverse, bool need_compl)
	{
		if (need_reverse && need_compl)
			return to_rev_complement(hash, kmer_len);

		if (!need_reverse && !need_compl)
			return hash;

		if (need_reverse)
			return bit_reverse(hash, kmer_len);

		hash = to_rev_complement(hash, kmer_len);
		return bit_reverse(hash, kmer_len);
	}

/*
	static hash_t change_hash_bits_for_char(hash_t hash, int i, int kmer_len, int bits)
	{
		// i = kmer-len - i - 1; ?
		hash &= ~(hash_t(3) << i*2);
		hash |= (hash_t(bits) << i*2);

		return hash;
	}


	template <class Lambda>
	static void for_all_1_char_variations_do(hash_t hash, int kmer_len, Lambda &&lambda)
	{
		for (int i=0; i < kmer_len; i++)
		for (int bits = 0; bits < 4; bits++)
		{
			auto mod_hash = change_hash_bits_for_char(hash, i, kmer_len, bits);
			if (mod_hash != hash)
				if (!lambda(mod_hash))
					return;
		}
	}

	template <class Lambda>
	static void for_all_2_char_variations_do(hash_t hash, int kmer_len, Lambda &&lambda)
	{
		for (int i = 0; i < kmer_len; i++)
		for (int j = i + 1; j < kmer_len; j++)
		for (int bitsi = 0; bitsi < 4; bitsi++)
		for (int bitsj = 0; bitsj < 4; bitsj++)
		{
			auto mod_hash = change_hash_bits_for_char(hash, i, kmer_len, bitsi);
			mod_hash = change_hash_bits_for_char(mod_hash, j, kmer_len, bitsj);
			if (mod_hash != hash)
				if (!lambda(mod_hash))
					return;
		}
	}

	template <class Lambda>
	static void for_all_3_char_variations_do(hash_t hash, int kmer_len, Lambda &&lambda)
	{
		for (int i = 0; i < kmer_len; i++)
		for (int j = i + 1; j < kmer_len; j++)
		for (int k = j + 1; k < kmer_len; k++)
		for (int bitsi = 0; bitsi < 4; bitsi++)
		for (int bitsj = 0; bitsj < 4; bitsj++)
		for (int bitsk = 0; bitsk < 4; bitsk++)
		{
			auto mod_hash = change_hash_bits_for_char(hash, i, kmer_len, bitsi);
			mod_hash = change_hash_bits_for_char(mod_hash, j, kmer_len, bitsj);
			mod_hash = change_hash_bits_for_char(mod_hash, k, kmer_len, bitsk);
			if (mod_hash != hash)
				if (!lambda(mod_hash))
					return;
		}
	}
*/

	//static void to_complement(std::string &seq)
	//{
	//	for (size_t i=0; i<seq.length(); i++)
	//		seq[i] = complement_letter(seq[i]);
	//}

//private:
//

	//template <class Lambda>
	//static void for_every_hash_variant_do(hash_t hash, int kmer_len, Lambda &&lambda) // todo: optimize ?
	//{
	//	lambda(hash);
	//	hash = bit_reverse(hash, kmer_len);
	//	lambda(hash);
	//	hash = to_rev_complement(hash, kmer_len);
	//	lambda(hash);
	//	hash = bit_reverse(hash, kmer_len);
	//	lambda(hash);
	//}
};

template <>
uint64_t seq_transform<uint64_t>::bit_reverse(uint64_t hash, int kmer_len)
{
    uint64_t dst =
        (uint64_t)reverse_map[hash      ] << 56 |
        (uint64_t)reverse_map[hash >>  8] << 48 |
        (uint64_t)reverse_map[hash >> 16] << 40 |
        (uint64_t)reverse_map[hash >> 24] << 32 |
        (uint64_t)reverse_map[hash >> 32] << 24 |
        (uint64_t)reverse_map[hash >> 40] << 16 |
        (uint64_t)reverse_map[hash >> 48] <<  8 |
        (uint64_t)reverse_map[hash >> 56];
    dst >>= sizeof(uint64_t) * 8 - kmer_len * 2;
    return dst;
}
template <typename hash_t> 
hash_t seq_transform<hash_t>::bit_reverse(hash_t hash, int kmer_len)
{
    hash_t dst = 0;
    for (int i=0; i<kmer_len; i++)
    {
        dst <<= 2;
        dst |= (hash & 3);
        hash >>= 2;
    }
    return dst;
}

template <>
uint64_t seq_transform<uint64_t>::to_rev_complement(uint64_t hash, int kmer_len)
{
    uint64_t comp_xor = 0xAAAAAAAAAAAAAAAAull;
    uint64_t comp = hash ^ comp_xor;
    return bit_reverse(comp, kmer_len);
}
template <typename hash_t>
hash_t seq_transform<hash_t>::to_rev_complement(hash_t hash, int kmer_len)
{
    hash_t dst = 0;
    for (int i=0; i<kmer_len; i++)
    {
        dst <<= 2;
        dst |= (hash + 2) & 3;
        hash >>= 2;
    }
    return dst;
}

struct seq_transform_actg
{
	static void apply_transformation(std::string &s, bool need_reverse, bool need_compl)
	{
		if (need_reverse)
			std::reverse(s.begin(), s.end());

		if (need_compl)
			for (auto &c : s)
				c = complement_letter(c);
	}

    static void to_rev_complement(std::string &s)
    {
        apply_transformation(s, true, true);
    }

private:
	static char complement_letter(char ch) // todo: optimize ?
	{
		switch (ch)
		{
			case 'A': return 'T';
			case 'C': return 'G';
			case 'T': return 'A';
			case 'G': return 'C';
		};

		return ch; // ? todo: think
	}
};

#endif
