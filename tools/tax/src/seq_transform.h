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

};

/*
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
*/

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
