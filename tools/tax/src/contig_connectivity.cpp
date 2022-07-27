#include <iostream>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <omp.h>
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

const int THREADS = 16;
#define MULTITHREADED 1


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
                done = !reader->read_many(chunk, 0);   
            }

#if 0
            for (auto& frag: chunk) 
                lambda(frag.bases);
#else
            {
                string spotid;
                vector<string> spot;

                for (auto& frag: chunk) 
                    if (frag.spotid == spotid)
                        spot.push_back(frag.bases);
                    else
                        {
                            if (!spot.empty())
                                lambda(spot);

                            spot.clear();
                            spot.push_back(frag.bases);
                            spotid = frag.spotid;
                        }

                if (!spot.empty())
                    lambda(spot);
            }

#endif
        }
    }
}

struct ContigPos
{
    int contig = 0, pos = 0;

    ContigPos() = default;
    ContigPos(int contig, int pos) : contig(contig), pos(pos){}
};

typedef vector<ContigPos> ContigPoss;
typedef std::unordered_map<hash_t, ContigPoss> ContigMap;


struct Connect
{
    ContigPos contig_pos;
    int len = 0;

    Connect() = default;
    Connect(const ContigPos &contig_pos, int len) : contig_pos(contig_pos), len(len){}
};

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

struct ConnectivitiesMT
{
    vector<Connectivities> conn_per_thread;

    ConnectivitiesMT(const Contigs &contigs) : conn_per_thread(THREADS)
    {
        // todo: if (max threads is 0, resize to 1)
        for (auto &conns : conn_per_thread)
        {
            conns.resize(contigs.size());
            for (int i = 0; i < contigs.size(); i++)
                conns[i].resize(contigs[i].seq.size());
        }
    }

    Connectivities reduce()
    {
        for (int thread_i = 1; thread_i < conn_per_thread.size(); thread_i++)
            for (int seq = 0; seq < conn_per_thread[thread_i].size(); seq++)
                for (int pos = 0; pos < conn_per_thread[thread_i][seq].size(); pos++)
                    conn_per_thread[0][seq][pos].len = std::max(conn_per_thread[0][seq][pos].len, conn_per_thread[thread_i][seq][pos].len);
        
        return conn_per_thread[0];
    }

    Connectivities &get()
    {
        int thread_id = omp_get_thread_num();
        if (thread_id < 0 || thread_id >= conn_per_thread.size())
            throw std::runtime_error("thread_id < 0 || thread_id >= conn_per_thread.size()");

        return conn_per_thread[thread_id];
    }

};

Connect connect(Connectivity &conn, const ContigPos &start_pos, const ContigPos &end_pos)
{
    if (conn.size() <= start_pos.pos || conn.size() <= end_pos.pos)
        throw std::runtime_error("conn.size() <= start_pos.pos || conn.size() <= end_pos.pos");

    if (end_pos.pos >= start_pos.pos)
    {
        conn[start_pos.pos].len = std::max(conn[start_pos.pos].len, end_pos.pos - start_pos.pos + KMER_LEN);
        return Connect(start_pos, conn[start_pos.pos].len);
    }
    else
        return connect(conn, end_pos, start_pos);        
}

bool looks_like_paired_read_distance(int a, int b) // todo: use statistics
{
    return std::abs(a - b) < 1000; // todo: think
}

void connect_spot(ConnectivitiesMT &conns, const Connect &read1_connect, const Connect &read2_connect)
{
    if (read1_connect.len <= 0 || read2_connect.len <= 0)
        return;

    if (read1_connect.contig_pos.contig != read2_connect.contig_pos.contig)
        return;        

    if (!looks_like_paired_read_distance(read2_connect.contig_pos.pos, read1_connect.contig_pos.pos))
        return;

    if (read1_connect.contig_pos.contig < 0 || read1_connect.contig_pos.contig >= conns.get().size())
        throw std::runtime_error("read1_connect.contig_pos.contig < 0 || read1_connect.contig_pos.contig >= conns.get().size()");

    auto &conn = conns.get()[read1_connect.contig_pos.contig];

    if (read2_connect.contig_pos.pos >= read1_connect.contig_pos.pos)
        conn[read1_connect.contig_pos.pos].len = std::max(conn[read1_connect.contig_pos.pos].len, read2_connect.contig_pos.pos - read1_connect.contig_pos.pos + read2_connect.len);
    else
        connect_spot(conns, read2_connect, read1_connect);
}

Connect update_connectivity(const ContigPos &start_pos, int read_start_pos, const string &rev_complement, ConnectivitiesMT &conns, const ContigMap &contig_map)
{
    int pos = 0;
    Connect connect_result;

    Hash<hash_t>::for_all_hashes_do(rev_complement, KMER_LEN, [&](hash_t hash)
    {
        hash = seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN);
        auto contig_poss = contig_map.find(hash);
        if (contig_poss != contig_map.end())
            for (auto &end_pos : contig_poss->second)
                if (end_pos.contig == start_pos.contig && about_the_same(pos_distance(end_pos, start_pos), read_distance(read_start_pos, pos, rev_complement.size())))
                {
                    connect_result = connect(conns.get()[start_pos.contig], start_pos, end_pos);
                    return false;
                }

        pos++;
        return pos < rev_complement.size()/2 - KMER_LEN;
    });

    return connect_result;
}

Connect update_connectivity(const string &bases, const string &rev_complement, ConnectivitiesMT &conns, const ContigMap &contig_map)
{
    Connect connect;
    int pos = 0;

    Hash<hash_t>::for_all_hashes_do(bases, KMER_LEN, [&](hash_t hash)
    {
        hash = seq_transform<hash_t>::min_hash_variant(hash, KMER_LEN);
        auto contig_poss = contig_map.find(hash);
        if (contig_poss != contig_map.end())
            for (auto &start_pos : contig_poss->second)
            {
                connect = update_connectivity(start_pos, pos, rev_complement, conns, contig_map);
                if (connect.len > 0)
                    return false; // todo: think
            }

        pos++;
        return pos < bases.size()/2;
    });

    return connect;
}

Connect process_contig_connectivity_per_read(const string &bases, ConnectivitiesMT &conns, const ContigMap &contig_map)
{
    string rev_complement = bases;
    seq_transform_actg::to_rev_complement(rev_complement);
    auto connect_direct = update_connectivity(bases, rev_complement, conns, contig_map);
    auto connect_rev_complement = update_connectivity(rev_complement, bases, conns, contig_map);

    return connect_direct.len > connect_rev_complement.len ? connect_direct : connect_rev_complement;
}    

Connectivities get_contig_connectivity(const Contigs &contigs, const string &accession)
{
    ConnectivitiesMT conns(contigs);

    ContigMap contig_map;
    load_contig_map(contigs, contig_map);

    size_t counter = 0;
    for_all_reads_do(accession, [&](const vector<string> &spot)
    {
        if (spot.size() == 2)
            connect_spot(conns, process_contig_connectivity_per_read(spot[0], conns, contig_map), process_contig_connectivity_per_read(spot[1], conns, contig_map));
        else
            for (auto &bases : spot)
                process_contig_connectivity_per_read(bases, conns, contig_map);

        counter ++;
        if (counter % 1024 == 0)
            cerr << ".";
    });

    return conns.reduce();
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
