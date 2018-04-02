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

#ifndef KMERS_MULTI_H_INCLUDED
#define KMERS_MULTI_H_INCLUDED

#include <set>
#include <string>
#include <map>
#include <unordered_map>
#include "tax_id_tree.h"
#include "kmer_hash.h"

//typedef uint64_t hash_t;

struct KmersMulti
{
    typedef std::set<tax_id_t> TaxIds;

	std::unordered_map<hash_t, TaxIds> storage;

	bool has_kmer(hash_t kmer) const
	{
		auto it = storage.find(kmer);
		return it != storage.end();
	}

	void add_kmer(hash_t kmer, tax_id_t tax_id)
	{
		auto &at = storage[kmer];
        at.insert(tax_id);
	}

	bool has_kmer_but_not_tax(hash_t kmer, tax_id_t tax_id) const
	{
		auto it = storage.find(kmer);
		if (it == storage.end())
			return false;

		if (it->second.find(tax_id) == it->second.end())
			return true;

		return false;
	}

};

#endif
