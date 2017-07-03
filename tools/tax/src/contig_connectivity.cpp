#include <iostream>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <unordered_map>
#include <vector>
#include "seq_transform.h"
#include "fasta.h"
#include "hash.h"
#include "reader.h"
#include "config_contig_connectivity.h"

using namespace std;
using namespace std::chrono;

const string VERSION = "0.10";

typedef uint64_t hash_t;
const int KMER_LEN = 32;

struct Contig
{
    string desc, seq;
    Contig(const string &desc, const string &seq) : desc(desc), seq(seq){}
};

typedef std::vector<Contig> Contigs;

Contigs load_contigs(const string &filename)
{
    Contigs contigs;
    Fasta fasta(filename);

    string seq;
    while (fasta.get_next_sequence(seq))
        contigs.push_back(Contig(fasta.sequence_description(), seq));

    return contigs;
}

template <class Lambda>
void for_all_reads_do(const string &accession, Lambda &&lambda)
{
    Reader::Params reader_params;
    reader_params.unaligned_only = true;
    reader_params.read_qualities = false;

    auto reader = Reader::create(accession, reader_params);

#if MULTITHREADED
    const int THREADS = 4;
    #pragma omp parallel num_threads(THREADS)
#endif
    {
        std::vector<Reader::Fragment> chunk;
        bool done = false;
        while (!done) 
        {
#if MULTITHREADED
            #pragma omp critical (read)
#endif
            {
                done = !reader->read_many(chunk);   
            }

            for (auto& frag: chunk) 
                lambda(frag.bases);
        }
    }
}

struct ContigPos
{
    int contig, pos;
    ContigPos(int contig, int pos) : contig(contig), pos(pos){}
};

typedef vector<ContigPos> ContigPoss;
typedef std::unordered_map<hash_t, ContigPoss> ContigMap;


void load_contig_map(const Contigs &contigs, ContigMap &contig_map)
{
    for (int contig = 0; contig < contigs.size(); contig++)
    {
        int pos = 0;
        Hash<hash_t>::for_all_hashes_do(contigs[contig].seq, KMER_LEN, [&](hash_t hash)
        {
            hash = seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN);
            contig_map[hash].push_back(ContigPos(contig, pos));
            pos++;
            return true;
        });
    }
}

int pos_distance(const ContigPos &a, const ContigPos &b)
{
    return std::abs(a.pos - b.pos);
}

int about_the_same(int a, int b)
{
    return abs(a - b) <= 2; // <= len / 20; // todo: tune
}

int read_distance(int direct_pos, int rev_compl_pos, int read_len)
{
    int d = read_len - direct_pos - rev_compl_pos - KMER_LEN;
    if (d < 0)
        throw std::runtime_error("read_distance < 0");

    return d;
}

struct Conn
{
    int len = 0;
    Conn() = default;
};

typedef vector<Conn> Connectivity;
typedef vector<Connectivity> Connectivities;

void connect(Connectivity &conn, const ContigPos &start_pos, const ContigPos &end_pos)
{
    if (conn.size() <= start_pos.pos || conn.size() <= end_pos.pos)
        throw std::runtime_error("conn.size() <= start_pos.pos || conn.size() <= end_pos.pos");

    if (end_pos.pos >= start_pos.pos)
        conn[start_pos.pos].len = std::max(conn[start_pos.pos].len, end_pos.pos - start_pos.pos + KMER_LEN);
    else
        connect(conn, end_pos, start_pos);        
}

bool update_connectivity(const ContigPos &start_pos, int read_start_pos, const string &rev_complement, Connectivities &conns, const ContigMap &contig_map)
{
    int pos = 0;
    bool connected = false;

    Hash<hash_t>::for_all_hashes_do(rev_complement, KMER_LEN, [&](hash_t hash)
    {
        hash = seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN);
        auto contig_poss = contig_map.find(hash);
        if (contig_poss != contig_map.end())
            for (auto &end_pos : contig_poss->second)
                if (end_pos.contig == start_pos.contig && about_the_same(pos_distance(end_pos, start_pos), read_distance(read_start_pos, pos, rev_complement.size())))
                {
                    connect(conns[start_pos.contig], start_pos, end_pos);
                    connected = true;
                    return false;
                }

        pos++;
        return pos < rev_complement.size()/2 - KMER_LEN;
    });

    return connected;
}

void update_connectivity(const string &bases, const string &rev_complement, Connectivities &conns, const ContigMap &contig_map)
{
    int pos = 0;

    Hash<hash_t>::for_all_hashes_do(bases, KMER_LEN, [&](hash_t hash)
    {
        hash = seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN);
        auto contig_poss = contig_map.find(hash);
        if (contig_poss != contig_map.end())
            for (auto &start_pos : contig_poss->second)
                if (update_connectivity(start_pos, pos, rev_complement, conns, contig_map)) 
                    return false; // todo: think

        pos++;
        return pos < bases.size()/2;
    });
}

Connectivities get_contig_connectivity(const Contigs &contigs, const string &accession)
{
    Connectivities conns(contigs.size());
    for (int i = 0; i < contigs.size(); i++)
        conns[i].resize(contigs[i].seq.size());

    ContigMap contig_map;
    load_contig_map(contigs, contig_map);

    size_t counter = 0;
    for_all_reads_do(accession, [&](const string &bases)
    {
        string rev_complement = bases;
        seq_transform_actg::to_rev_complement(rev_complement);
        update_connectivity(bases, rev_complement, conns, contig_map);
        update_connectivity(rev_complement, bases, conns, contig_map);
        counter ++;
        if (counter % 1024 == 0)
            cerr << ".";
    });

    return conns;
}


void print_connectivity(const Connectivities &conns)
{
    for (auto &conn : conns)
    {
        int len = 0;
        for (auto &c : conn)
        {
            len = std::max(c.len, len);
            cout << len << '\t';
            len--;
        }

        cout << endl;
    }
}
/*

void print_connectivity(const Connectivities &conns)
{
    auto &conn = conns[0];
    {
        int len = 0;
        for (auto &c : conn)
        {
            len = std::max(c.len, len);
            cout << len << endl;
            len--;
        }

        cout << endl;
    }
}
*/

int main(int argc, char const *argv[])
{
	Config config(argc, argv);

	auto before = high_resolution_clock::now();

	auto contigs = load_contigs(config.fasta_filename);
    if (contigs.empty())
        return 0;

    auto conns = get_contig_connectivity(contigs, config.accession);
    print_connectivity(conns);

	cerr << "total time (s) " << std::chrono::duration_cast<std::chrono::seconds>( high_resolution_clock::now() - before ).count() << endl;

    return 0;
}
