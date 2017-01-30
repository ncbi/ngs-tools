#ifndef REDUCE_H_INCLUDED
#define REDUCE_H_INCLUDED

#include <vector>
#include <algorithm>
#include <numeric>

struct reduce
{
	static int median_of(std::vector<int> &v)
	{
		std::sort(v.begin(), v.end());
		return v[v.size()/2];
	}

	static int average_of(const std::vector<int> &v)
	{
		if (v.empty())
			return 0;

		int sum = std::accumulate(v.begin(), v.end(), 0);
		return sum/v.size();
	}

	static int median_of_const(const std::vector<int> &_v)
	{
		std::vector<int> v(_v);
		return median_of(v);
	}

	static std::vector<int> median_filter(const std::vector<int> &v, int radius) // todo: make it right
	{
		std::vector<int> result;
		result.reserve(v.size());

		for (int pos =0; pos < radius; pos++)
			result.push_back(v[pos]);

		std::vector<int> rs;
		auto diameter = radius * 2 + 1;
		rs.resize(diameter);
		for (int pos = radius; pos < v.size() - radius; pos++)
		{
			for (int dx = -radius; dx <= radius; dx++) // todo: deque or queue
				rs[radius + dx] = v[pos + dx];

			result.push_back(median_of(rs));
		}

		while (result.size() < v.size())
			result.push_back(v[result.size()]);

		return result;
	}
};

#endif
