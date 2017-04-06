#ifndef MERGE_H_INCLUDED
#define MERGE_H_INCLUDED

#include <vector>

template <class MatchId>
void merge_vectors(const std::vector < vector<MatchId> > &matched_ids, std::vector<MatchId> &ids)
{
	ids.clear();
	int ids_count = 0;
	for (int i = 0; i < matched_ids.size(); i++)
		ids_count += matched_ids[i].size();

	ids.reserve(ids_count);
	for (int i = 0; i < matched_ids.size(); i++)
	for (auto seq_id : matched_ids[i])
		ids.push_back(seq_id);
}

#endif