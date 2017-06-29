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
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/seek.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "DBGraph.hpp"
#include "graphdigger.hpp"

using namespace boost::program_options;
using namespace DeBruijn;

vector<pair<double,char>> FilteredNodeSuccessors(CDBGraph::Node node, CDBGraph& graph, CDBGraphDigger& graph_digger, bool reverse) {
    vector<pair<double,char>> rslt;
    if(!node)
        return rslt;

    vector<CDBGraph::Successor> successors = graph.GetNodeSuccessors(node);
    sort(successors.begin(), successors.end(), [&](const CDBGraph::Successor& a, const CDBGraph::Successor& b) {return graph.Abundance(a.m_node) > graph.Abundance(b.m_node);});
    for(auto& suc : successors)
        rslt.push_back(make_pair(graph.Abundance(suc.m_node),suc.m_nt));
    graph_digger.FilterNeighbors(successors, true);
    for(auto& r :  rslt) {
        char nt = r.second;
        if(reverse)
            r.second = Complement(r.second);
        if(find_if(successors.begin(),successors.end(),[nt](const CDBGraph::Successor& s){return s.m_nt == nt;}) == successors.end())
            r.second = tolower(r.second);
    }

    return rslt;
}

int main(int argc, const char* argv[])
{
    options_description all("Test options");
    all.add_options()
        ("help,h", "Produce help message")
        ("version,v", "Print version")
        ("dbg", value<string>(), "de Bruijn graph")
        ("fraction", value<double>()->default_value(0.1, "0.1"), "Threshold for extension")
        ("mindeadend", value<int>()->default_value(50), "Ignore dead ends shorter than this")
        ("lowcount", value<int>()->default_value(6), "Minimal count for filtering")
        ("genome", value<string>(), "Assembled genome")
        ("kmer", value<int>(), "Kmer length for testing")
        ("test_kmer", value<string>(), "Build contig for test kmer");

    string dbg;
    string genome_file;
    double fraction;
    int jump;
    int low_count;
    int kmer = 0;
    string test_kmer;
    variables_map argm;                                // boost arguments

    try {
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
            cerr << "dbgtester v.1.0" << endl;
#ifdef SVN_REV
            cerr << "SVN revision:" << SVN_REV << endl << endl;
#endif
            return 0;
        }

        if(argm.count("dbg"))
            dbg = argm["dbg"].as<string>();
        else {
            cerr << "Provide de Bruijn graph" << endl;
            cerr << all << "\n";
            return 1;
        }

        fraction = argm["fraction"].as<double>();
        jump = argm["mindeadend"].as<int>();
        low_count = argm["lowcount"].as<int>();
        if(argm.count("genome"))
            genome_file = argm["genome"].as<string>();
        if(argm.count("kmer"))
            kmer = argm["kmer"].as<int>();
        if(argm.count("test_kmer"))
            test_kmer = argm["test_kmer"].as<string>();

    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        cerr << all << "\n";
        return 1;
    }

    map<int,CDBGraph*> graphs;
    ifstream file(dbg);
    file.seekg (0, file.end);
    streampos file_length = file.tellg();
    file.seekg (0, file.beg);
    while(file.tellg() != file_length) {
        CDBGraph* graphp = new CDBGraph(file);
        graphs[graphp->KmerLen()] = graphp;
        cerr << "Loaded kmer: " << graphp->KmerLen() << endl;
    }    
    
    if(!genome_file.empty()) {
        ifstream fasta(genome_file);
        if(!fasta.is_open()) {
            cerr << "Can't open file " << genome_file << endl;
            return 1;             
        }
        map<string,string> genome;
        string line;
        string* contigp;
        int genome_size = 0;        
        while(getline(fasta,line)) {
            if(line[0] == '>') {
                size_t l = line.find_first_of(" \t");
                contigp = &(genome[line.substr(1,l-1)]);
            } else {
                for (auto & c: line) c = toupper(c);
                *contigp += line;
                genome_size += line.size();
            }
        }
        cerr << "Genome size: " << genome_size << " Number of contigs: " << genome.size() << endl;
        CDBGraph* graphp = graphs.begin()->second;
        if(kmer > 0) {
            graphp = graphs[kmer];
            if(graphp == 0) {
                cerr << "Wrong kmer length" << endl;
                return 1;
            }
        }

        CDBGraph& graph(*graphp);
        int kmer_len = graph.KmerLen();
        CDBGraphDigger graph_digger(graph, fraction, jump, low_count); 

        for(auto& nc : genome) {
            const string& acc = nc.first;
            string& contig = nc.second;
            map<int,int> bins;
            int not_in_graph = 0;
            for(int pos = 0; pos < (int)contig.size(); ++pos) {
                CDBGraph::Node lnode = 0;
                int labundance = 0;
                double lplusf = 0;
                if(pos >= kmer_len-1) {
                    lnode = graph.GetNode(contig.substr(pos-kmer_len+1,kmer_len));
                    if(lnode) {
                        labundance = graph.Abundance(lnode);
                        lplusf = graph.PlusFraction(lnode);
                        ++bins[labundance];      
                    } else  {
                    ++not_in_graph;
                    }
                }
                CDBGraph::Node rnode = 0;
                int rabundance = 0;
                double rplusf = 0;
                if(pos <= (int)contig.size()-kmer_len) {
                    rnode = graph.GetNode(contig.substr(pos,kmer_len));
                    if(rnode) {
                        rabundance = graph.Abundance(rnode);
                        rplusf = graph.PlusFraction(rnode);
                    }
                }                    

                cout << acc << "\t" << pos+1 << "\t" << contig[pos] << "\t";
                if(labundance > 0)
                    cout << labundance << "\t" << lplusf << "/" << 1-lplusf << "\t";
                else
                    cout << "-\t-\t";
                if(rabundance > 0)
                    cout << rabundance << "\t" << rplusf << "/" << 1-rplusf << "\t";
                else
                    cout << "-\t-\t";
        
                vector<pair<double,char>> sinfo;
                if(pos >= kmer_len)
                    sinfo = FilteredNodeSuccessors(graph.GetNode(contig.substr(pos-kmer_len,kmer_len)), graph, graph_digger, false);
                if(sinfo.empty()) {
                    cout << "-";
                } else {
                    for(int i = 0; i < (int)sinfo.size(); ++i) {
                        if(i > 0)
                            cout << ";";
                        cout << sinfo[i].second << "/" << sinfo[i].first;
                    }
                }
                cout << "\t";
                sinfo.clear();
                if(pos < (int)contig.size()-kmer_len)
                    sinfo = FilteredNodeSuccessors(CDBGraph::ReverseComplement(graph.GetNode(contig.substr(pos+1,kmer_len))), graph, graph_digger, true);
                if(sinfo.empty()) {
                    cout << "-";
                } else {
                    for(int i = 0; i < (int)sinfo.size(); ++i) {
                        if(i > 0)
                            cout << ";";
                        cout << sinfo[i].second << "/" << sinfo[i].first;
                    }
                }
                cout << endl;
            }
        }

        return 0;

    } else if(!test_kmer.empty()) {
        CDBGraph* graphp = graphs[test_kmer.size()];
        if(!graphp) {
            cerr << "Wrong test_kmer size" << endl;
            return 1;
        }
        CDBGraphDigger graph_digger(*graphp, fraction, jump+test_kmer.size(), low_count, true);         
        //        CDBGraphDigger graph_digger(*graphp, fraction, jump, low_count, false);         
        
        CDBGraph::Node node = graphp->GetNode(test_kmer);
        if(!node) {
            cerr << "Kmer: " << test_kmer << " not in graph" << endl;
            return 1;
        }
        SContig scontig = graph_digger.GetContigForKmer(node, 0);

        if(scontig.LenMax() == 0) {
            cerr << "No assebly" << endl;
            return 0;
        }

        /*
        TContigSequenceList contigs(1, scontig.m_seq);
        CombineSimilarContigs(contigs, 0, test_kmer.size()-1);
        auto& contig = contigs.front();
        */

        auto& contig = scontig.m_seq;

        /*
        int num = 0;
        for(auto& chunk : contig) {
            cerr << "Chunk" << ++num << endl;
            for(auto& seq : chunk) {
                for(char c : seq)
                    cerr << c;
                cerr << endl;
            }
        }
        */

        
        string first_variant;
        for(auto& lst : contig) {
            for(char c : lst.front()) {
                if(c != '-')
                    first_variant.push_back(c);
            }
        }
        CReadHolder rh(false);
        rh.PushBack(first_variant);
        double abundance = 0; // average count of kmers in contig
        for(CReadHolder::kmer_iterator itk = rh.kbegin(graphp->KmerLen()); itk != rh.kend(); ++itk) {
            CDBGraph::Node node = graphp->GetNode(*itk);
            abundance += graphp->Abundance(node);
        }
        abundance /= first_variant.size()-graphp->KmerLen()+1;
        cout << ">Contig_" << abundance << "\n" << first_variant << endl;

        int pos = 0;
        for(unsigned chunk = 0; chunk < contig.size()-1; ++chunk) { //output variants
            int chunk_len = 0;
            for(char c : contig[chunk].front()) {
                if(c != '-')
                    ++chunk_len;
            }
            if(contig.VariableChunk(chunk)) {
                int left = min(100,(int)contig[chunk-1].front().size());
                int right = min(100,(int)contig[chunk+1].front().size());
                int var = 0;
                auto it = contig[chunk].begin();
                for(++it; it != contig[chunk].end(); ++it) {
                    auto& variant = *it;
                    cout << ">Variant_" << ++var << "_for_Contig:" << pos-left+1 << "_" << pos+chunk_len+right << "\n";
                    for(int l = left ; l > 0; --l)
                        cout << *(contig[chunk-1].front().end()-l);
                    for(char c : variant) {
                        if(c != '-')
                            cout << c;
                    }
                    for(int r = 0; r < right; ++r)
                        cout << contig[chunk+1].front()[r];
                    cout << endl;
                }
            }
            pos += chunk_len;
        }
            
    } else {
        cerr << "Provide genome or test_kmer" << endl;
        cerr << all << "\n";
        return 1;
    }
    

    return 0;
}

