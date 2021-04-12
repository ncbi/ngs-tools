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

#ifndef KMERS_H_INCLUDED
#define KMERS_H_INCLUDED

#include <set>
#include <string>
#include <map>
#include <unordered_map>
#include "tax_id_tree.h"

typedef uint64_t hash_t;

struct Kmers
{
	const TaxIdTree &tax_id_tree;

	std::unordered_map<hash_t, tax_id_t> storage;

	Kmers(const TaxIdTree &tax_id_tree) : tax_id_tree(tax_id_tree)
	{
		storage.reserve(128*1024*1024); // todo: tune
	}

	bool has_kmer(hash_t kmer) const
	{
		auto it = storage.find(kmer);
		return it != storage.end();
	}

	void add_kmer(hash_t kmer, tax_id_t tax_id)
	{
		auto &at = storage[kmer];
		if (at == tax_id) // todo: remove?
			return;

		if (at == 0)
			at = tax_id;
		else
			at = tax_id_tree.consensus_of(tax_id, at); // todo: avoid writing if match?
	}

	bool has_kmer_but_not_tax(hash_t kmer, tax_id_t tax_id) const
	{
		auto it = storage.find(kmer);
		if (it == storage.end())
			return false;

		tax_id_t stored_tax_id = it->second;
		if (!stored_tax_id)
			return true;

		return !tax_id_tree.a_sub_b(tax_id, stored_tax_id);
	}

	// obsolete and inefficient. left only for old check_index implementation todo: remove
	void update_kmer(hash_t kmer, tax_id_t tax_id)
	{
		if (storage.find(kmer) == storage.end())
			return;

		auto &at = storage[kmer];
		if (at == tax_id) // todo: remove?
			return;

		if (at == 0)
			at = tax_id;
		else
			at = tax_id_tree.consensus_of(tax_id, at); // todo: avoid writing if match?
	}

};

#endif
