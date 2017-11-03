/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================
*
*/

#include <boost/program_options.hpp>

#include "readsgetter.hpp"
#include "concurrenthash.hpp"
#include "graphdigger.hpp"

using namespace boost::program_options;
using namespace DeBruijn;

int main(int argc, const char* argv[]) {
    for(int n = 0; n < argc; ++n)
        cerr << argv[n] << " ";
    cerr << endl << endl;

    int kmer;
    int ncores;
    double vector_percent;
    int min_count;
    vector<string> sra_list;
    vector<string> fasta_list;
    size_t estimated_kmer_num;
    bool usepairedends;

    options_description all("Options");
    all.add_options()
        ("help,h", "Produce help message")
        ("version,v", "Print version")
        ("sra_run", value<vector<string>>(), "Input sra run accession (could be used multiple times for different runs) [string]")
        ("fasta", value<vector<string>>(), "Input fasta file(s) (could be used multiple times for different runs) [string]")
        ("use_paired_ends", "Use pairing information from paired reads in fasta input [flag]")
        ("kmer", value<int>()->default_value(21), "Kmer length [integer]")
        ("min_count", value<int>()->default_value(2), "Minimal count for kmers retained for comparing alternate choices [integer]")
        ("vector_percent ", value<double>()->default_value(0.05, "0.05"), "Count for  vectors as a fraction of the read number [float [0,1)]")

        ("hash_count", "Use hash counter [flag]")
        ("estimated_kmers", value<int>()->default_value(100), "Estimated number of unique kmers for bloom filter (M) only for hash count [integer]")
        ("skip_bloom_filter", "Don't do bloom filter; use --estimated_kmers as the hash table size (only for hash count) [flag]")
        ("memory", value<int>()->default_value(400), "Memory available (GB) only for sorted count [integer]")

        ("dbg_out", value<string>(), "De Bruijn graph output")
        ("text_out", value<string>(), "Text kmer output")
        ("hist", value<string>(), "File for histogram [string]")

        ("cores", value<int>()->default_value(0), "Number of cores to use (default all) [integer]");

    try {
        variables_map argm;                                // boost arguments
        store(parse_command_line(argc, argv, all), argm);
        notify(argm);    

        if(argm.count("help")) {
#ifdef SVN_REV
            cerr << "SVN revision:" << SVN_REV << endl << endl;
#endif
            cerr << all << "\n";
            return 1;
        }

        if(argm.count("version")) {
            cerr << "kmercounter v.2.0" << endl;
#ifdef SVN_REV
            cerr << "SVN revision:" << SVN_REV << endl << endl;
#endif
            return 0;
        }

        if(!argm.count("sra_run") && !argm.count("fasta")) {
            cerr << "Provide some input reads" << endl;
            cerr << all << "\n";
            return 1;
        }

        if(argm.count("sra_run")) {
            sra_list = argm["sra_run"].as<vector<string>>();
            unsigned num = sra_list.size();
            sort(sra_list.begin(), sra_list.end());
            sra_list.erase(unique(sra_list.begin(),sra_list.end()), sra_list.end());
            if(sra_list.size() != num)
                cerr << "WARNING: duplicate input entries were removed from SRA run list" << endl; 
            usepairedends = true;
        }
        if(argm.count("fasta")) {
            fasta_list = argm["fasta"].as<vector<string>>();
            unsigned num = fasta_list.size();
            sort(fasta_list.begin(), fasta_list.end());
            fasta_list.erase(unique(fasta_list.begin(),fasta_list.end()), fasta_list.end());
            if(fasta_list.size() != num)
                cerr << "WARNING: duplicate input entries were removed from fasta file list" << endl; 
            usepairedends = argm.count("use_paired_ends");
        }

        ncores = thread::hardware_concurrency();
        if(argm["cores"].as<int>()) {
            int nc = argm["cores"].as<int>();
            if(nc < 0) {
                cerr << "Value of --cores must be >= 0" << endl;
                exit(1);
            } else if(nc > ncores) {
                cerr << "WARNING: number of cores was reduced to the hardware limit of " << ncores << " cores" << endl;
            } else if(nc > 0) {
                ncores = nc;
            }
        }

        min_count = argm["min_count"].as<int>();
        if(min_count <= 0) {
            cerr << "Value of --min_count must be > 0" << endl;
            exit(1);
        }

        estimated_kmer_num =  argm["estimated_kmers"].as<int>();
        kmer = argm["kmer"].as<int>();

        CReadsGetter readsgetter(sra_list, fasta_list, vector<string>(), ncores, usepairedends, false);
        vector_percent = argm["vector_percent "].as<double>();
 
        if(argm.count("hash_count")) {
            readsgetter.ClipAdaptersFromReads_HashCounter(vector_percent, estimated_kmer_num, argm.count("skip_bloom_filter"));
            readsgetter.PrintAdapters();

            size_t MB = 1000000;
            CKmerHashCounter counter(readsgetter.Reads(), kmer, min_count, estimated_kmer_num*MB, true, ncores, argm.count("skip_bloom_filter"));

            if(argm.count("text_out")) {
                ofstream out(argm["text_out"].as<string>());
                if(!out.is_open()) {
                    cerr << "Can't open file " << argm["text_out"].as<string>() << endl;
                    return 1;             
                }

                CKmerHashCount& hash = counter.Kmers();
                for(auto index = hash.Begin(); index != hash.End(); ++index) {
                    auto rslt = index.GetElement();
                    out << rslt.first.toString(kmer) << "\t" << rslt.second->Count() << "\t" << (rslt.second->m_data.Load() >> 32) << endl;
                }
            }

            if(argm.count("hist")) {
                ofstream out(argm["hist"].as<string>());
                if(!out.is_open()) {
                    cerr << "Can't open file " << argm["hist"].as<string>() << endl;
                    return 1;             
                }
                TBins bins = counter.Kmers().GetBins();
                for(auto& bin : bins)
                    out << bin.first << '\t' << bin.second << endl;
            }

            if(argm.count("dbg_out")) {
                counter.GetBranches();
                CDBHashGraph graph(move(counter.Kmers()), true);
                ofstream dbg_out(argm["dbg_out"].as<string>());
                if(!dbg_out.is_open()) {
                    cerr << "Can't open file " << argm["dbg_out"].as<string>() << endl;
                    exit(1);
                }
                graph.Save(dbg_out);
            }                
        } else {
            int64_t memory = argm["memory"].as<int>();

            readsgetter.ClipAdaptersFromReads_SortedCounter(vector_percent, memory);
            readsgetter.PrintAdapters();

            int64_t GB = 1000000000;            
            CKmerCounter counter(readsgetter.Reads(), kmer, min_count, true, GB*memory, ncores);
            CKmerCount& kmers = counter.Kmers();

            if(argm.count("text_out")) {
                ofstream out(argm["text_out"].as<string>());
                if(!out.is_open()) {
                    cerr << "Can't open file " << argm["text_out"].as<string>() << endl;
                    return 1;             
                }

                for(size_t index = 0; index < kmers.Size(); ++index) {
                    pair<TKmer,size_t> rslt = kmers.GetKmerCount(index);
                    out << rslt.first.toString(kmer) << "\t" << (uint32_t)rslt.second << "\t" << (rslt.second >> 32) << endl;
                }
            }

            if(argm.count("hist")) {
                ofstream out(argm["hist"].as<string>());
                if(!out.is_open()) {
                    cerr << "Can't open file " << argm["hist"].as<string>() << endl;
                    return 1;             
                }
                map<int,size_t> hist;
                for(size_t index = 0; index < kmers.Size(); ++index) {
                    ++hist[kmers.GetCount(index)];                  // count clipped to integer automatically
                }
                TBins bins(hist.begin(), hist.end());
                for(auto& bin : bins)
                    out << bin.first << '\t' << bin.second << endl;
            }

            if(argm.count("dbg_out")) {
                counter.GetBranches();
                map<int,size_t> hist;
                for(size_t index = 0; index < kmers.Size(); ++index) {
                    ++hist[kmers.GetCount(index)];                  // count clipped to integer automatically
                }
                TBins bins(hist.begin(), hist.end());
                CDBGraph graph(move(kmers), move(bins), true);
                ofstream dbg_out(argm["dbg_out"].as<string>());
                if(!dbg_out.is_open()) {
                    cerr << "Can't open file " << argm["dbg_out"].as<string>() << endl;
                    exit(1);
                }
                graph.Save(dbg_out);
            }
        }
        
        cerr << "DONE" << endl;
        exit(0);
    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        exit(1);
    }

}

