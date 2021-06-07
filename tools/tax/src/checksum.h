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

#pragma once

#define AUX_READER_NO_BASES 0
#define AUX_READER_LOG_CHECKSUM 0

#if AUX_READER_LOG_CHECKSUM
#define AUX_READER_LOG_PERIODIC_CHECKSUM 0

struct CheckSum {
    /**
     * These checksums are insensitive to changes in the order of bases.
     * Additionally, the first one is insensitive to changes in strand.
     *
     * This comes from the commutativity of modular multiplication, i.e.:
     * ab mod N = a mod N * b mod N mod N = b mod N * a mod N mod N = ba mod N
     */
    uint64_t cksumGC = 1;
    uint64_t cksumACGT = 1;

    // This one is not order-insensitive.
    uint64_t fnv1a_hash = 0xcbf29ce484222325; ///< the FNV-1a initial value
    uint64_t readCount = 0;

    void update(std::string const &bases) {
        /**
         * A set of 17-bit primes.
         * These happen to be in a sequence but that's not important. It is not
         * important that they are primes, only that they're coprime, so that
         * there's no "mixing" between them in the arithmetic.
         *
         * The final result is:
         *   pA^A * pC^C * pG^G * pT^T * pN^N mod M
         *   pAT^AT * pCG^CG * pN^N mod M
         * where:
         *   A is the count of 'A's, C is the count of 'C's, etc. and N is the count of not ACGT.
         * It's easy to see how the order of the bases doesn't matter in the final result.
         */
        auto const pN = 130981ull;
        auto const pA = 130987ull;
        auto const pC = 131009ull;
        auto const pG = 131011ull;
        auto const pT = 131023ull;
        auto const pAT = 131041ull;
        auto const pCG = 131059ull;
        auto const modulus = 0x100000001B3ull; ///< the FNV-1a prime

        for (auto base : bases) {
            auto baseValue = pN;
            auto baseUnstrandedValue = pN;

            switch (base) {
            case 'A':
                baseValue = pA;
                baseUnstrandedValue = pAT;
                break;
            case 'T':
                baseValue = pT;
                baseUnstrandedValue = pAT;
                break;
            case 'C':
                baseValue = pC;
                baseUnstrandedValue = pCG;
                break;
            case 'G':
                baseValue = pG;
                baseUnstrandedValue = pCG;
                break;
            }
            cksumGC = (cksumGC * baseUnstrandedValue) % modulus;
            cksumACGT = (cksumACGT * baseValue) % modulus;
            fnv1a_hash = (fnv1a_hash ^ baseValue) * modulus;
        }
        readCount += 1;
#if AUX_READER_LOG_PERIODIC_CHECKSUM
        if ((readCount & (0x100000ull - 1)) == 0)
            log(true);
#endif
    }
    void log(bool shortlog = false) const {
        char buf[32];
        auto &&to_hex = [&](uint64_t const &value) -> char const * {
            auto x = value;
            char *cp = buf + sizeof(buf);

            *--cp = '\0';
            do {
                *--cp = "0123456789ABCDEF"[x & 0xF];
                x >>= 4;
            } while (x != 0);
            return cp;
        };
#if AUX_READER_LOG_PERIODIC_CHECKSUM
        if (shortlog) {
            LOG("checksums: " << to_hex(fnv1a_hash) << ' ' << to_hex(cksumACGT) << ' ' << to_hex(cksumGC));
            return;
        }
#else
        ((void)shortlog); // suppress unused parameter warning
#endif
        LOG("checksum: " << to_hex(fnv1a_hash));
        LOG("unordered checksum: " << to_hex(cksumACGT));
        LOG("unordered+unstranded checksum: " << to_hex(cksumGC));
    }
};
#endif

