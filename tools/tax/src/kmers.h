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

//	std::map<hash_t, tax_id_t> storage;
	std::unordered_map<hash_t, tax_id_t> storage;

	Kmers(const TaxIdTree &tax_id_tree) : tax_id_tree(tax_id_tree)
	{
		storage.reserve(100*1000*1000);
	}

	bool has_kmer(hash_t kmer) const
	{
		auto it = storage.find(kmer);
		return it != storage.end();
	}

	void add_kmer(hash_t kmer, tax_id_t tax_id)
	{
//		std::cout << "add kmer " << kmer << " " << tax_id << std::endl;
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

		if (!tax_id)
		{
			std::cerr << "!tax_id in has_kmer_but .. test" << std::endl;
			throw 1;
		}

		return !tax_id_tree.a_sub_b(tax_id, stored_tax_id);
	}

	//void add_kmers(hash_t kmer, const TaxIds &tax_ids)
	//{
	//	auto &at = storage[kmer];
	//	at.insert(tax_ids.begin(), tax_ids.end()); // todo: swap? move?
	//}

	//void add_tax_for_existing(hash_t kmer, tax_id_t tax_id)
	//{
	//	auto it = storage.find(kmer);
	//	if (it == storage.end())
	//		return;

	//	auto &at = it->second;

	//	if (at.size() < MAX_TAXES)
	//		at.insert(tax_id);
	//}
};

#endif
