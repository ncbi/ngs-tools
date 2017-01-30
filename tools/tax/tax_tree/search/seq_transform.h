#ifndef SEQ_TRANSFORM_H_INCLUDED
#define SEQ_TRANSFORM_H_INCLUDED

#include <string>
#include <algorithm>
//#include "hash.h"

template <class hash_t>
struct seq_transform
{
	static hash_t min_hash_variant(hash_t hash, int kmer_len, bool *complement = nullptr, bool *reverse = nullptr)
	{
		bool _compl = false, _rev = false;
		if (!complement)
			complement = &_compl;
		if (!reverse)
			reverse = &_rev;

		auto min_hash = hash;
		*complement = false;
		*reverse = false;

		hash = bit_reverse(hash, kmer_len);
		if (hash < min_hash)
		{
			min_hash = hash;
			*complement = false;
			*reverse = true;
		}

		hash = to_rev_complement(hash, kmer_len);
		if (hash < min_hash)
		{
			min_hash = hash;
			*complement = true;
			*reverse = false;
		}

		hash = bit_reverse(hash, kmer_len);
		if (hash < min_hash)
		{
			min_hash = hash;
			*complement = true;
			*reverse = true;
		}

		return min_hash;
	}

	static hash_t min_hash_variant2(hash_t hash, int kmer_len)
	{
		auto min_hash = hash;

		hash = to_rev_complement(hash, kmer_len);
		if (hash < min_hash)
			return hash;

		return min_hash;
	}

	static hash_t bit_reverse(hash_t hash, int kmer_len)
	{
		hash_t dst;
		for (int i=0; i<kmer_len; i++)
		{
			dst <<= 2;
			dst |= (hash & 3);
			hash >>= 2;
		}

		return dst;
	}

	static hash_t to_rev_complement(hash_t hash, int kmer_len)
	{
		const char COMPLEMENT[4] = {2, 3, 0, 1};
		hash_t dst;
		for (int i=0; i<kmer_len; i++)
		{
			dst <<= 2;
			dst |= COMPLEMENT[hash & 3];
			hash >>= 2;
		}

		return dst;
	}

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

struct seq_transform_actg
{
	static std::string apply_transformation(std::string s, bool need_reverse, bool need_compl)
	{
		if (need_reverse)
			std::reverse(s.begin(), s.end());

		if (need_compl)
			for (auto &c : s)
				c = complement_letter(c);
			
		return s;
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
