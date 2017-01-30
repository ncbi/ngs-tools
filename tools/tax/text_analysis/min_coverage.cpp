#include "min_coverage.h"

int main(int argc, char const *argv[])
{
    try
    {
		ConfigMinCoverage config(argc, argv);
		cerr << "min_coverage version " << VERSION << endl;

		print_current_time();
		auto before = high_resolution_clock::now();

		TextLoaderSTNoStore seq_loader(config.reference);

		KmerPositions kmers;
		KmerBuilder kmer_builder(seq_loader, kmers, config.kmer_len);
		cerr << "kmers loaded at (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;

		compress(kmer_builder, kmers, config.kmer_len, config.lookup_len, 
			[](Sequence *seq, size_t pos, const Kmer &kmer, const KmerPositions &kmers)
			{
//				cout << pos << " " << kmers.coverage_of(kmer) << " ";
				print_string(kmer.s, kmer.len);
				cout << endl;
//				print_used(seq->used);
			},
			[](Sequence *seq)
			{
//				print_used(seq->used);
				cerr << ".";
			}
		);

		cerr << "total time (sec) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;

		exit(0); // dont want to wait for KmerMaps destructors
        return 0;
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
