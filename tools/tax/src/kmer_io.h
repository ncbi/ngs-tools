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

#include <iostream>
#include "hash.h"
#include "kmers.h"
#include "kmers_multi.h"
#include "seq_transform.h"
#include "tax_id_tree.h"

//test tc

struct KmerIO
{
    static hash_t kmer_from(const char *s, int from, int len)
    {
	    hash_t kmer = Hash<hash_t>::hash_of(s + from, len);
	    return seq_transform<hash_t>::min_hash_variant(kmer, len);
    }

    static std::string str_kmer(hash_t kmer, int kmer_len)
    {
	    return Hash<hash_t>::str_from_hash(kmer, kmer_len);
    }

    static void print_kmers(const Kmers &kmers, int kmer_len)
    {
	    for (auto &kmer : kmers.storage)
	    {
		    auto tax_id = kmer.second;
		    if (tax_id != TaxIdTree::ROOT)
			    std::cout << str_kmer(kmer.first, kmer_len) << '\t' << tax_id << std::endl;
	    }
    }

    static void print_kmers(const KmersMulti &kmers, int kmer_len)
    {
	    for (auto &kmer : kmers.storage)
        {
            std::cout << str_kmer(kmer.first, kmer_len);

		    for (auto tax_id : kmer.second)
			     std::cout << '\t' << tax_id;
                 
            std::cout << std::endl;
        }
    }

};

