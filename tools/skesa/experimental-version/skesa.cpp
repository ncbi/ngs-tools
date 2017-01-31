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
#include "assembler.hpp"

using namespace boost::program_options;
using namespace DeBruijn;

int main(int argc, const char* argv[]) {
    for(int n = 0; n < argc; ++n)
        cerr << argv[n] << " ";
    cerr << endl << endl;

    int ncores;
    int steps;
    double fraction;
    int jump;
    int low_count;
    int min_count;
    int min_kmer;
    bool usepairedends;
    int maxkmercount;
    ofstream contigs_out;
    ofstream all_out;
    ofstream hist_out;
    ofstream connected_reads_out;
    ofstream dbg_out;
    int memory;
    int max_kmer_paired = 0;
    vector<string> sra_list;
    vector<string> fasta_list;
    vector<string> fastq_list;
    bool gzipped;
    int mincontig;

    options_description general("General options");
    general.add_options()
        ("help,h", "Produce help message")
        ("memory", value<int>()->default_value(32), "Memory available (GB) [integer]")
        ("cores", value<int>()->default_value(0), "Number of cores to use (default all) [integer]");

    options_description input("Input/output options : at least one input providing reads for assembly must be specified");
    input.add_options()
        ("fasta", value<vector<string>>(), "Input fasta file(s) (could be used multiple times for different runs) [string]")
        ("fastq", value<vector<string>>(), "Input fastq file(s) (could be used multiple times for different runs) [string]")
        ("sra_run", value<vector<string>>(), "Input sra run accession (could be used multiple times for different runs) [string]")
        ("gz", "Input files are gzipped [flag]")
        ("contigs_out", value<string>(), "Output file for contigs (stdout if not specified) [string]");

    options_description assembly("Assembly options");
    assembly.add_options()
        ("kmer", value<int>()->default_value(21), "Minimal kmer length for assembly [integer]")
        ("min_count", value<int>()->default_value(2), "Minimal count for kmers retained for comparing alternate choices [integer]")
        ("use_paired_ends", "Use pairing information from paired reads in input [flag]")
        ("insert_size", value<int>(), "Expected insert size for paired reads (if not provided, it will be estimated) [integer]")
        ("steps", value<int>()->default_value(11), "Number of assembly iterations from minimal to maximal kmer length in reads [integer]")
        ("max_kmer_count", value<int>()->default_value(10), "Minimum acceptable average count for estimating the maximal kmer length in reads [integer]")
        ("fraction", value<double>()->default_value(0.1, "0.1"), "Maximum noise to signal ratio acceptable for extension [float [0,1)]")
        ("min_dead_end", value<int>()->default_value(50), "Ignore dead end paths shorter than this when comparing alternate extensions [integer]")
        ("low_count", value<int>()->default_value(6), "Minimal count for kmers used in assembly [integer]")
        ("min_contig", value<int>()->default_value(200), "Minimal contig length reported in output [integer]");

    options_description debug("Debugging options");
    debug.add_options()
        ("all", value<string>(), "Output fasta for each iteration [string]")
        ("dbg_out", value<string>(), "Output kmer file [string]")
        ("hist", value<string>(), "File for histogram [string]")
        ("connected_reads", value<string>(), "File for connected paired reads [string]");

    options_description all("");
    all.add(general).add(input).add(assembly).add(debug); 

    try {
        variables_map argm;                                // boost arguments
        store(parse_command_line(argc, argv, all), argm);
        notify(argm);    

        if(argm.count("help")) {
            cerr << all << "\n";
            return 1;
        }

        if(!argm.count("fasta") && !argm.count("fastq") && !argm.count("sra_run")) {
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
        }
        if(argm.count("fasta")) {
            fasta_list = argm["fasta"].as<vector<string>>();
            unsigned num = fasta_list.size();
            sort(fasta_list.begin(), fasta_list.end());
            fasta_list.erase(unique(fasta_list.begin(),fasta_list.end()), fasta_list.end());
            if(fasta_list.size() != num)
                cerr << "WARNING: duplicate input entries were removed from fasta file list" << endl; 
        }
        if(argm.count("fastq")) {
            fastq_list = argm["fastq"].as<vector<string>>();
            unsigned num = fastq_list.size();
            sort(fastq_list.begin(), fastq_list.end());
            fastq_list.erase(unique(fastq_list.begin(),fastq_list.end()), fastq_list.end());
            if(fastq_list.size() != num)
                cerr << "WARNING: duplicate input entries were removed from fastq file list" << endl; 
        }
        gzipped = argm.count("gz");
   
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

        steps = argm["steps"].as<int>();
        if(steps <= 0) {
            cerr << "Value of --steps must be > 0" << endl;
            exit(1);
        }
        fraction = argm["fraction"].as<double>();
        if(fraction >= 1.) {
            cerr << "Value of --fraction must be < 1 (more than 0.25 is not recommended)" << endl;
            exit(1);
        }
        if(fraction < 0.) {
            cerr << "Value of --fraction must be >= 0" << endl;
            exit(1);
        }
        jump = argm["min_dead_end"].as<int>();
        if(jump < 0) {
            cerr << "Value of --min_dead_end must be >= 0" << endl;
            exit(1);
        }
        low_count = argm["low_count"].as<int>();
        if(low_count <= 0) {
            cerr << "Value of --low_count must be > 0" << endl;
            exit(1);
        }
        if(argm.count("insert_size"))
            max_kmer_paired = argm["insert_size"].as<int>();
        if(max_kmer_paired < 0) {
            cerr << "Value of --insert_size must be >= 0" << endl;
            exit(1);
        }
        mincontig = argm["min_contig"].as<int>();
        if(mincontig <= 0) {
            cerr << "Value of --min_contig must be > 0" << endl;
            exit(1);
        }
        min_count = argm["min_count"].as<int>();
        if(min_count <= 0) {
            cerr << "Value of --min_count must be > 0" << endl;
            exit(1);
        }
        if(low_count < min_count) {
            cerr << "WARNING: --low_count changed from " << low_count << " to " << min_count << " as that is the minimum count retained" << endl;
            low_count = min_count;
        }
        min_kmer = argm["kmer"].as<int>();
        if(min_kmer < 21 || min_kmer%2 ==0) {
            cerr << "Kmer must be an odd number >= 21" << endl;
            return 1;
        }
        usepairedends = argm.count("use_paired_ends");
        maxkmercount = argm["max_kmer_count"].as<int>();
        if(maxkmercount <= 0) {
            cerr << "Value of --max_kmer_count must be > 0" << endl;
            exit(1);
        }

        memory = argm["memory"].as<int>();
        if(memory <= 0) {
            cerr << "Value of --memory must be > 0" << endl;
            exit(1);
        }

        if(argm.count("contigs_out")) {
            contigs_out.open(argm["contigs_out"].as<string>());
            if(!contigs_out.is_open()) {
                cerr << "Can't open file " << argm["contigs_out"].as<string>() << endl;
                exit(1);
            }
        }
    
        if(argm.count("all")) {
            all_out.open(argm["all"].as<string>());
            if(!all_out.is_open()) {
                cerr << "Can't open file " << argm["all"].as<string>() << endl;
                exit(1);
            }
        }


        if(argm.count("hist")) {
            hist_out.open(argm["hist"].as<string>());
            if(!hist_out.is_open()) {
                cerr << "Can't open file " << argm["hist"].as<string>() << endl;
                exit(1);
            }
        }

        if(argm.count("connected_reads")) {
            connected_reads_out.open(argm["connected_reads"].as<string>());
            if(!connected_reads_out.is_open()) {
                cerr << "Can't open file " << argm["connected_reads"].as<string>() << endl;
                exit(1);
            }
        }

        if(argm.count("dbg_out")) {
            dbg_out.open(argm["dbg_out"].as<string>());
            if(!dbg_out.is_open()) {
                cerr << "   Can't open file " << argm["dbg_out"].as<string>() << endl;
                exit(1);
            }
        }

        CReadsGetter readsgetter(sra_list, fasta_list, fastq_list, ncores, usepairedends, gzipped);
        CDBGAssembler assembler(fraction, jump, low_count, steps, min_count, min_kmer, usepairedends, max_kmer_paired, maxkmercount, memory, ncores, readsgetter.Reads()); 

        CDBGraph& first_graph = *assembler.Graphs().begin()->second;
        int num = 0; 
        ostream& out = contigs_out.is_open() ? contigs_out : cout;
        for(string& contig : assembler.Contigs()) {
            if((int)contig.size() >= mincontig) {
                CReadHolder rh(false);
                rh.PushBack(contig);
                double abundance = 0; // average count of kmers in contig
                for(CReadHolder::kmer_iterator itk = rh.kbegin(first_graph.KmerLen()); itk != rh.kend(); ++itk) {
                    CDBGraph::Node node = first_graph.GetNode(*itk);
                    abundance += first_graph.Abundance(node);
                }
                abundance /= contig.size()-first_graph.KmerLen()+1;
                out << ">Contig_" << ++num << "_" << abundance << "\n" << contig << endl;;
            } 
        }  
        
        if(all_out.is_open()) {
            auto graphp = assembler.Graphs().begin();
            for(auto& contigs : assembler.AllIterations()) {
                int nn = 0;
                for(auto& contig : contigs) 
                    all_out << ">kmer" << graphp->first << "_" << ++nn << endl << contig << endl;
                ++graphp;
            }
        } 

        if(hist_out.is_open()) {
            for(auto& gr : assembler.Graphs()) {
                const TBins& bins = gr.second->GetBins();
                for(auto& bin : bins)
                    hist_out << gr.first << '\t' << bin.first << '\t' << bin.second << endl;
            }
        }

        if(connected_reads_out.is_open()) {
            CReadHolder connected_reads = assembler.ConnectedReads();
            int num = 0;
            for(CReadHolder::string_iterator is = connected_reads.sbegin(); is != connected_reads.send(); ++is) {
                string s = *is;
                connected_reads_out << ">ConnectedRead_" << ++num << "\n" << s << "\n";
            }
        }

        if(dbg_out.is_open()) {
            for(auto& gr : assembler.Graphs())
                gr.second->Save(dbg_out);
        }

    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        exit(1);
    }

    cerr << "DONE" << endl;

    return 0;
}
