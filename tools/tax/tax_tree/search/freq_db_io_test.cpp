#include "freq_db_io.h"

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

int main(int argc, char const *argv[])
{
    try
    {
        {
            Frequences freqs;
            PredictedFrequences pred_freqs;
            FreqDBIO::load_frequences("./test/tax5.fasta.out", freqs, pred_freqs);
            cout << "loaded" << endl;
            FreqDBIO::save_frequences("./test/tax5.fasta.out_reconstructed", freqs, pred_freqs);
        }

        {
            Frequences freqs;
            PredictedFrequences pred_freqs;
            FreqDBIO::load_frequences("/panfs/pan1.be-md.ncbi.nlm.nih.gov/tax_analysis/tax_tree_2_june/tree/Viruses/Retro-transcribing viruses/Retroviridae/Orthoretrovirinae/Lentivirus/Primate lentivirus group/Human immunodeficiency virus 1/11676.fasta.out", freqs, pred_freqs);
            FreqDBIO::save_frequences("/panfs/pan1.be-md.ncbi.nlm.nih.gov/tax_analysis/tax_tree_2_june/tree/Viruses/Retro-transcribing viruses/Retroviridae/Orthoretrovirinae/Lentivirus/Primate lentivirus group/Human immunodeficiency virus 1/11676.fasta.out.reconstructed", freqs, pred_freqs);
//            FreqDBIO::load_frequences("/panfs/pan1.be-md.ncbi.nlm.nih.gov/tax_analysis/tax_tree_2_june/tree/Viruses/ssRNA viruses/ssRNA negative-strand viruses/Bunyaviridae/Orthobunyavirus/Simbu virus/35306.fasta.out", freqs, pred_freqs);
 //           FreqDBIO::save_frequences("/panfs/pan1.be-md.ncbi.nlm.nih.gov/tax_analysis/tax_tree_2_june/tree/Viruses/ssRNA viruses/ssRNA negative-strand viruses/Bunyaviridae/Orthobunyavirus/Simbu virus/35306.fasta.out.reconstructed", freqs, pred_freqs);
//            FreqDBIO::load_frequences("./test/human.out", freqs, pred_freqs);
  //          FreqDBIO::save_frequences("./test/human.out_reconstructed", freqs, pred_freqs);
        }
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
