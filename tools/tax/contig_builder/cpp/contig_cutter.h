#include <vector>
#include <string>
#include <list>

struct ContigCutter
{
public:
	static const int DEFAULT_SENSITIVITY = 6;
	static std::list<std::string> cut(const std::string &seq, const std::vector<int> &cov, int kmer_len, int sensitivity = DEFAULT_SENSITIVITY)
	{
		auto cut_points = identify_cut_points(cov, kmer_len, sensitivity);
		return cut_by_points(seq, cut_points);
	}

	static std::list<std::string> cut_by_points(const std::string &seq, const std::list<int> &cut_points)
	{
		std::list<std::string> contigs;
		//if (cut_points.empty())
		//{
		//	contigs.push_back(seq);
		//	return contigs;
		//}

		int from = 0;
		for (auto to : cut_points)
		{
			if (to > seq.size()) 
			{
				throw std::runtime_error("ContigCutter::cut to > seq.size()"); // todo: remove ?
				to = seq.size();
			}

			contigs.push_back(seq.substr(from, to - from));
			from = to;
		}

//		std::cerr << "from " << from << " len " << seq.size() - from << std::endl;
		contigs.push_back(seq.substr(from, seq.size() - from));
		return contigs;
	}

public:

	//static int left_pos(int right_pos, int distance, int left_bound)
	//{
	//	int pos = right_pos - distance;
	//	return pos < left_bound ? left_bound : pos;
	//}

	//static int current_piece_coverage(const std::vector<int> &cov, int right_pos, int left_bound)
	//{
	//	const int LEFT_POS_DISTANCE = 10;

	//}


	static std::list<int> identify_cut_points(const std::vector<int> &cov, int kmer_len, int sensitivity = DEFAULT_SENSITIVITY)
	{
		std::list<int> points;
		if (cov.empty())
			return points;

		struct CurrentContig
		{
			int contig_start;

			CurrentContig(int contig_start): contig_start(contig_start){}
			int coverage(const std::vector<int> &cov, int right_pos) const 
			{ 
				int from = std::max(contig_start, right_pos - 40);
				int to = std::max(contig_start, right_pos - 20);
				int sum_cov = 0;
				for (int i=from; i<=to; i++)
					sum_cov += cov[i];

				return sum_cov/(to - from + 1);
			}

		} current_contig(0);
	
		for (int right_pos = 1; right_pos < cov.size();)
		{
			auto left_cov = current_contig.coverage(cov, right_pos); //current_piece_coverage(cov, right_pos, left_bound); //cov[left_pos(right_pos, LEFT_POS_DISTANCE, left_bound)];
			auto right_cov = cov[right_pos];
//			std::cout << left_cov << " " << right_cov << " at right pos " << right_pos << std::endl;
			if (jumps_down(left_cov, right_cov, sensitivity))
			{
//				std::cout << "jumps down " << left_cov << " " << right_cov << " at right pos " << right_pos << std::endl;
				right_pos += kmer_len;
				if (right_pos < cov.size())
					points.push_back(right_pos);
				current_contig = CurrentContig(right_pos);
			}
			else if (jumps_up(left_cov, right_cov, sensitivity))
			{
//				std::cout << "jumps up " << left_cov << " " << right_cov << " at right pos " << right_pos << std::endl;
				points.push_back(right_pos);
				current_contig = CurrentContig(right_pos);
			}

			right_pos ++;
		}

		return points;
	}

private:
	static const int MIN_SENSITIVITY = 200;

	static bool jumps_down(int prev, int now, int sensitivity)
	{
		if (prev < MIN_SENSITIVITY && now < MIN_SENSITIVITY)
			return false;

		return now < prev / sensitivity;
	}

	static bool jumps_up(int prev, int now, int sensitivity)
	{
		if (prev < MIN_SENSITIVITY && now < MIN_SENSITIVITY)
			return false;

		return now > prev * sensitivity;
	}


};


