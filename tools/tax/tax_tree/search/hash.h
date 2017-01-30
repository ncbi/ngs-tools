#ifndef HASH_H_INCLUDED
#define HASH_H_INCLUDED

#include <algorithm>
//typedef uint64_t hash_t;

template <class hash_t>
struct Hash
{
/*
	static const int MAX_HASH_BITS = 16;

	static size_t hash_bits(int kmer_len)
	{
		return std::min(MAX_HASH_BITS, kmer_len);
	}

	static hash_t hash_mask(int kmer_len)
	{
		return hash_t(1) << (kmer_len * 2); 
	}

	*/
	// todo: works only for 16, 32, 64
	template <class To>
	static To left(hash_t x)
	{
		union 
		{
			hash_t x_big;
			struct
			{
				To lo;
				To hi;
			} x_small;
		} u;

		u.x_big = x;

		return u.x_small.hi;

		//return To(x >> (sizeof(x)*8/2));
	}

	template <class To>
	static To right(hash_t x)
	{
		union 
		{
			hash_t x_big;
			struct
			{
				To lo;
				To hi;
			} x_small;
		} u;

		u.x_big = x;

		return u.x_small.lo;

//		return To(x);
	}

	//static hash_t hash_values(int kmer_len)
	//{
	//	return hash_t(1) << hash_bits(kmer_len);
	//}

	static unsigned int hash_bits(int kmer_len)
	{
		return 2 * kmer_len;
	}

	static hash_t hash_of(const std::string &s)
	{
		return hash_of(&s[0], s.length());
	}

	static hash_t hash_of(const char *s, int kmer_len)
	{
		hash_t hash = 0;
		for (int i=0; i < kmer_len; i++)
			hash = update_hash(s[i], hash);

		return hash;
	}

    static hash_t hash_of_amino_one(const char *s) // standard code. todo: optimize
    {
        switch (s[0])
        {
            case 'T':
                switch (s[1])
                {
                    case 'T':
                        switch (s[2])
                        {
                            case 'T': return 0; // F
                            case 'C': return 0; // F
                            case 'A': return 1; // L
                            case 'G': return 1; // L
                        }
                    case 'C':
                        return 2; // S

                    case 'A':
                        switch (s[2])
                        {
                            case 'T': return 3; // Y
                            case 'C': return 3; // Y
                            case 'A': return 4; // *
                            case 'G': return 4; // *
                        }
                    case 'G':
                        switch (s[2])
                        {
                            case 'T': return 5; // C
                            case 'C': return 5; // C
                            case 'A': return 4; // *
                            case 'G': return 6; // W
                        }
                }
            case 'C':
                switch (s[1])
                {
                    case 'T': return 1; // L
                    case 'C': return 7; // P
                    case 'A':
                        switch (s[2])
                        {
                            case 'T': return 8; // H
                            case 'C': return 8; // H
                            case 'A': return 9; // Q
                            case 'G': return 9; // Q
                        }
                    case 'G': return 10; // R
                }

            case 'A':
                switch (s[1])
                {
                    case 'T':
                        switch (s[2])
                        {
                            case 'T': return 11; // I
                            case 'C': return 11; // I
                            case 'A': return 11; // I
                            case 'G': return 12; // M
                        }

                    case 'C':
                        return 13; // T

                    case 'A':
                        switch (s[2])
                        {
                            case 'T': return 14; // N
                            case 'C': return 14; // N
                            case 'A': return 15; // K
                            case 'G': return 15; // K
                        }

                    case 'G':
                        switch (s[2])
                        {
                            case 'T': return 2; // S
                            case 'C': return 2; // S
                            case 'A': return 10; // R
                            case 'G': return 10; // R
                        }

                }

            case 'G':
                switch (s[1])
                {
                    case 'T': return 16; // V
                    case 'C': return 17; // A
                    case 'A':
                        switch (s[2])
                        {
                            case 'T': return 18; // D
                            case 'C': return 18; // D
                            case 'A': return 19; // E
                            case 'G': return 19; // E
                        }
                    case 'G': return 20; // G
                }
        };

        throw "hash_of_amino_one - bad bases";
    }

    static const int AMINO_COUNT = 21;

	static hash_t hash_of_amino(const char *s, int kmer_len)
	{
        if ((kmer_len % 3) != 0)
            throw "(kmer_len % 3) != 0";

		hash_t hash = 0;
		for (int i=0; i < kmer_len; i+=3)
        {
            hash *= AMINO_COUNT;
			hash += hash_of_amino_one(&s[i]);
        }

		return hash;
	}

	static hash_t update_hash(char ch, hash_t hash)
	{
		hash <<= 2;
		hash_t new_bit1 = (ch & 2) >> 1;
		hash_t new_bit2 = (ch & 4) >> 2;
		hash |= ( new_bit1 | (new_bit2 << 1)); // todo: can remove one extra shift
		return hash;
	}

	static hash_t hash_next(const char *s, hash_t hash, int kmer_len)
	{
		return hash_next(s[kmer_len -1], hash, kmer_len);
	}

	static hash_t hash_next(char ch, hash_t hash, int kmer_len)
	{
		hash &= ( (hash_t(1) << (kmer_len*2 -2)) -1  );
		hash = update_hash(ch, hash);
		return hash;
	}

	static std::string str_from_hash(hash_t hash, int kmer_len)
	{
		std::string s;
		s.reserve(kmer_len);
		for (int i=0; i<kmer_len; i++, hash >>=2)
			s += hash_char(hash & 3);

		std::reverse(s.begin(), s.end()); // todo: optimize. get rid of reverse. 
		return s;
	}

	static char hash_char(hash_t hash)
	{
		switch (hash)
		{
			case 0: return 'A';
			case 1: return 'C';
			case 2: return 'T';
			case 3: return 'G';
		};

		return 'N';
	}

	template <class Lambda>
	static void for_all_hashes_do(const std::string &s, int kmer_len, Lambda &&lambda)
	{
		if (s.length() < kmer_len)
			return;

		auto hash = hash_of(&s[0], kmer_len);
		for (int i=0; i <= s.size() - kmer_len; i++, hash = hash_next(&s[i], hash, kmer_len))
			if (!lambda(hash))
				break;
	}
};

namespace std
{
	template<>
	struct hash<__uint128_t>
	{
		size_t operator() (__uint128_t x) const
		{
			return Hash<__uint128_t>::left<uint64_t>(x) ^ Hash<__uint128_t>::right<uint64_t>(x);
				//  size_t(x); // todo: left xor right
		}
	};
}



#endif
