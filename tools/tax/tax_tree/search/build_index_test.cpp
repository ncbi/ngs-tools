#include "build_index.h"

//using namespace std;

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

int main(int argc, char const *argv[])
{
    try
    {
		{
			Kmers kmers;
//			add_kmers(kmers, "./test/seq.fasta", 10, 100);
			add_kmers(kmers, "/panfs/pan1.be-md.ncbi.nlm.nih.gov/tax_analysis/human_asm.38.5.fasta", 10, 8000);

			print_kmers(kmers);
			return 0;
		}

		{
			const char *seq = "JJJACTGXXAGGX";
			equal(string(&seq[3], 4), string("ACTG"));

			SeqCleaner cleaner(seq);
			equal(cleaner.clean_strings.size(), 2);
			auto it = cleaner.clean_strings.begin();
//			cout << "prin str" << endl;
//			print_string(it->s, it->len);
//			cout << endl;
			auto s1 = string(it->s, it->len);
//			cout << "s1 is" << s1;
			equal(s1, string("ACTG"));
			equal(it->len, 4);
			it++;
			equal(string(it->s, it->len), string("AGG"));
			equal(it->len, 3);
		}

		{
			SeqCleaner cleaner("ACTG");
			equal(cleaner.clean_strings.size(), 1);
			auto it = cleaner.clean_strings.begin();
			equal(string(it->s), string("ACTG"));
			equal(it->len, 4);
		}

		{
			SeqCleaner cleaner("ACTGNA");
			equal(cleaner.clean_strings.size(), 2);
			auto it = cleaner.clean_strings.begin();
			equal(string(it->s), string("ACTGNA"));
			equal(it->len, 4);
			it++;
			equal(string(it->s), string("A"));
			equal(it->len, 1);
		}

		Kmers kmers;
		equal(kmers.has_kmer("AA"), false);
		kmers.add_kmer("AA", 10);
		equal(kmers.has_kmer("AA"), true);
		equal(kmers.has_kmer("AC"), false);
		kmers.add_kmer("AC", 11);
		kmers.add_kmer("AA", 14);
		equal(kmers.has_kmer("AA"), true);
		equal(kmers.has_kmer("AC"), true);

//		print_kmers(kmers);

		FileListLoader l("./test/files.list.small");
		equal(l.files.size(), 3);
		auto it = l.files.begin();

		equal(it->filesize, 6517412686);
		equal(it->filename, "./tree/Bacteria/Firmicutes/Bacilli/Lactobacillales/Streptococcaceae/Streptococcus/Streptococcus pneumoniae/1313.fasta");
		it++;
		equal(it->filesize, 5408935484);
		equal(it->filename, "./tree/Bacteria/Proteobacteria/Gammaproteobacteria/Pseudomonadales/Pseudomonadaceae/Pseudomonas/Pseudomonas aeruginosa/287.fasta");
		it++;
		equal(it->filesize, 5370001716);
		equal(it->filename, "./tree/Bacteria/Proteobacteria/Gammaproteobacteria/Enterobacteriales/Enterobacteriaceae/Salmonella/Salmonella enterica subsp. enterica serovar Typhi/90370.fasta");

		equal(tax_id_from("./tree/Bacteria/Firmicutes/Bacilli/Lactobacillales/Streptococcaceae/Streptococcus/Streptococcus pneumoniae/1313.fasta"), 1313);
    }
    catch ( exception & x )
    {
        cerr << x.what() << endl;
//		cerr << "exit 3" << endl;
		return 3;
    }
    catch ( string & x )
    {
        cerr << x << endl;
//		cerr << "exit 4" << endl;
		return 4;
    }
    catch ( const char * x )
    {
        cerr << x << endl;
//		cerr << "exit 5" << endl;
		return 5;
    }
    catch ( ... )
    {
        cerr << "unknown exception" << endl;
//		cerr << "exit 6" << endl;
		return 6;
    }
}
