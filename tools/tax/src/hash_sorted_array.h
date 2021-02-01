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

#include <vector>
#include "dbs.h"

typedef uint64_t hash_t;

struct KmerTax : public DBS::KmerTax
{
    KmerTax(hash_t kmer = 0, int tax_id = 0) : DBS::KmerTax(kmer, tax_id){} // todo: remove constructor from KmerTax for faster loading ?

    bool operator < (const KmerTax &x) const // for binary search by hash
    {
        return kmer < x.kmer;
    }
};

typedef int tax_t;
typedef std::vector<KmerTax> HashSortedArray;

// todo: remove duplicates
static tax_t find_hash(hash_t hash, tax_t default_value, HashSortedArray &hash_array)
{
    auto first = hash_array.begin();
    auto last = hash_array.end();
    first = std::lower_bound(first, last, KmerTax(hash, 0));
    return ((first == last) || (hash < first->kmer) ) ? default_value : first->tax_id;
}

