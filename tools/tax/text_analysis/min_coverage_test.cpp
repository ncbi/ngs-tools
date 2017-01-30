#include "min_coverage.h"
#include <iostream>

using namespace std;

template <class A, class B>
void equal(A a, B b)
{
	if (a == b)
		cout << "ok" << endl;
	else
	{
		cout << "FAIL: " << a << " != " << b << endl;
//		exit(1);
	}
}

void print_it(Sequence *seq, size_t pos, const Kmer &kmer, const KmerPositions &kmers)
{
	cout << pos << " " << kmers.coverage_of(kmer) << " ";
	print_string(kmer.s, kmer.len);
	cout << endl;

	print_used(seq->used);
}


int main(int argc, char const *argv[])
{
	try
	{
		{
			Sequence s;
			s.seq = "12345";
			s.used.resize(s.seq.size());

			equal(s.used_mark(0), false);
			equal(s.used_mark(1), false);
			equal(s.used_mark(2), false);
			equal(s.used_mark(3), false);
			equal(s.used_mark(4), false);

			equal(find_first_used_kmer_mark(1, 10, &s), 0);
			equal(find_last_used_kmer_mark(1, 10, &s), 0);

			s.mark_used(1);

			equal(s.used_mark(0), false);
			equal(s.used_mark(1), true);
			equal(s.used_mark(2), false);

			equal(find_first_used_kmer_mark(1, 10, &s), 1);
			equal(find_last_used_kmer_mark(1, 10, &s), 1);

			equal(find_first_used_kmer_mark(2, 10, &s), 0);
			equal(find_last_used_kmer_mark(2, 10, &s), 0);

			s.mark_used(2);
			equal(find_first_used_kmer_mark(1, 10, &s), 1);
			equal(find_last_used_kmer_mark(1, 10, &s), 2);
		}
		{
			Sequence s;
			s.seq = "12345";
			s.used.resize(s.seq.size());
			mark_as_used(&s, 1, 3);

			equal(s.used_mark(0), false);
			equal(s.used_mark(1), true);
			equal(s.used_mark(2), true);
			equal(s.used_mark(3), true);
			equal(s.used_mark(4), false);
		}

		{
			Sequence s;
			s.seq = "12345";
			s.used.resize(s.seq.size());
			mark_as_used(&s, 2, 30);

			equal(s.used_mark(0), false);
			equal(s.used_mark(1), false);
			equal(s.used_mark(2), true);
			equal(s.used_mark(3), true);
			equal(s.used_mark(4), true);
		}

		{
			Sequence s1;
			s1.seq = "ABCDEF";
			s1.used.resize(s1.seq.size());

			Sequence s2;
			s2.seq = "BCDEFG";
			s2.used.resize(s2.seq.size());

			KmerPositions kmers;
			kmers.add(&s1, 0, 4);
			equal(kmers.coverage_of(p_string("ABCD", 4)), 1);
			kmers.add(&s1, 1, 4);
			equal(kmers.coverage_of(p_string("BCDE", 4)), 1);
			kmers.add(&s2, 0, 4);
			equal(kmers.coverage_of(p_string("BCDE", 4)), 2);
			equal(kmers.coverage_of(p_string("CDEF", 4)), 0);
			equal(kmers.coverage_of(p_string("ABCD", 4)), 1);

			equal(kmers.get_positions(p_string("ABCD", 4)).size(), 1);
			for (auto p : kmers.get_positions(p_string("ABCD", 4)))
			{
				equal(p.seq, &s1);
				equal(p.pos, 0);
			}

			{
				equal(kmers.get_positions(p_string("BCDE", 4)).size(), 2);
				int i = 0;
				Sequence* test_seq[2] = {&s2, &s1};
				int test_pos[2] = {0, 1};
				for (auto p : kmers.get_positions(p_string("BCDE", 4)))
				{
					equal(p.seq, test_seq[i]);
					equal(p.pos, test_pos[i]);
					i++;
				}
			}

			{
				cout << "storage test" << endl;
				TextLoaderSTNoStore seq_loader("./1.fasta");
				KmerPositions kmers;
				KmerBuilder kmer_builder(seq_loader, kmers, 4);
				equal(kmer_builder.sequences.size(), 1);
				equal(kmers.storage.size(), 7);
				{
					int i = 0;
					int test_pos[2] = {0, 1};
					equal(kmers.get_positions(p_string("AAAA", 4)).size(), 2);
					for (auto p : kmers.get_positions(p_string("AAAA", 4)))
					{
						equal(p.pos, test_pos[i]);
						i++;
					}
				}
				{
					int i = 0;
					int test_pos[2] = {4, 8};
					cout << "check ACTG" << endl;
					equal(kmers.get_positions(p_string("ACTG", 4)).size(), 2);
					equal(kmers.coverage_of(p_string("ACTG", 4)), 2);
					for (auto p : kmers.get_positions(p_string("ACTG", 4)))
					{
						equal(p.pos, test_pos[i]);
						i++;
					}
				}

				Sequence *seq = &*kmer_builder.sequences.begin();
				equal(find_last_used_kmer_mark(1, 4+4, seq), 0);
				equal(find_next_free_pos(1, 4, seq, kmers), 4);
				equal(kmers.coverage_of(p_string("ACTG", 4)), 2);

				compress(kmer_builder, kmers, 4, print_it); 
			}

			{
				cout << "storage test 2" << endl;
				TextLoaderSTNoStore seq_loader("./2.fasta");
				KmerPositions kmers;
				KmerBuilder kmer_builder(seq_loader, kmers, 4);
				compress(kmer_builder, kmers, 4, print_it); 
			}

		}

	}
	catch (const string &s)
	{
		cout << s << endl;
	}
	catch (const char *s)
	{
		cout << s << endl;
	}
}

