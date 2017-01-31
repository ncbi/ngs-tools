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

#ifndef _DBGAssembler_
#define _DBGAssembler_

#include <random>
#include "DBGraph.hpp"
#include "counter.hpp"
#include "graphdigger.hpp"

#include "concurrenthash.hpp"

namespace DeBruijn {

    /******************************
    General description

    CDBGAssembler implements the SKESA assembling algorithm.

    1. It uses the counts for kmers with the minimal kmer length specified (default 21 bp) to estimate the maximal kmer length
       (starting from average mate length) that has sufficient coverage requested in maxkmercount. If reads are paired and
       insert size isn't specified, it estimates the insert size by assembling between mates for a sample of the reads.

    2. It assembles iteratively starting from minimal to maximal kmer length in a specified number of steps. Each step builds a
       de Bruijn graph for the kmer size for that iteration and uses it to improve previously assembled contigs. After each
       assembly iteration, the reads already used in the contigs are removed from further consideration.

    3. If reads are paired, it uses the reads that are not marked as used and the set of de Bruijn graphs built in 2) to connect
       the mate pairs. 

    4. Using the paired reads connected in 3), it performs three additional assembly iterations with the kmer size up
       to the insert size.
    *******************************/

    class CDBGAssembler {
    public:
        // fraction - Maximal noise to signal ratio of counts acceptable for extension
        // jump - minimal length of accepted dead ends; i.e. dead ends shorter than this length are ignored
        // low_count - minimal count for kmers in a contig
        // steps - number of assembly iterations from minimal to maximal kmer size in reads
        // min_count - minimal kmer count to be included in a de Bruijn graph
        // min_kmer - the minimal kmer size for the main steps
        // usepairedends - whether or not to use paired ends
        // max_kmer_paired - insert size (0 if not known)
        // maxkmercount - the minimal average count for estimating the maximal kmer
        // memory - the upper bound for memory use (GB)
        // ncores - number of threads
        // raw_reads - reads (for effective multithreading, number of elements in the list should be >= ncores)
        
        CDBGAssembler(double fraction, int jump, int low_count, int steps, int min_count, int min_kmer, bool usepairedends, 
                      int max_kmer_paired, int maxkmercount, int memory, int ncores, list<array<CReadHolder,2>>& raw_reads) : 
            m_fraction(fraction), m_jump(jump), m_low_count(low_count), m_steps(steps), m_min_count(min_count), m_min_kmer(min_kmer), m_usepairedends(usepairedends),
            m_max_kmer_paired(max_kmer_paired), m_maxkmercount(maxkmercount), m_memory(memory), m_ncores(ncores), m_raw_reads(raw_reads) {

            m_scan_window = 50; // the size-1 of the contig's flank area used for extensions and connections
            m_max_kmer = m_min_kmer;
            m_insert_size = 0;

            for(auto& reads : m_raw_reads) {
                m_raw_pairs.push_back({reads[0], CReadHolder(false)});
            }    
            m_connected_reads.resize(m_raw_reads.size(), {CReadHolder(false), CReadHolder(true)});

            //graph for minimal kmer
            double average_count = GetGraph(m_min_kmer, m_raw_reads, true);
            if(average_count == 0)
                throw runtime_error("Reads are too short for selected minimal kmer length");

            // estimate genome
            double total_seq = 0;
            size_t total_reads = 0;
            for(auto& reads : m_raw_reads) {
                total_seq += reads[0].TotalSeq()+reads[1].TotalSeq();
                total_reads += reads[0].ReadNum()+reads[1].ReadNum();
            }
            int read_len = total_seq/total_reads+0.5;
            cerr << endl << "Average read length: " << read_len << endl;
            size_t genome_size = m_graphs[m_min_kmer]->GenomeSize();
            cerr << "Genome size estimate: " << genome_size << endl << endl;

            {// first iteration
                ImproveContigs(min_kmer);
                if(m_contigs.back().empty())
                    throw runtime_error("Was not able to assemble anything");
            }

            //estimate max_kmer
            if(m_steps > 1 && average_count > maxkmercount) {
                m_max_kmer = read_len+1-double(maxkmercount)/average_count*(read_len-min_kmer+1);
                m_max_kmer = min(TKmer::MaxKmer(), m_max_kmer);
                while(m_max_kmer > min_kmer) {
                    m_max_kmer -= 1-m_max_kmer%2;           // odd kmers desired
                    CKmerCounter kmer_counter(m_raw_reads, m_max_kmer, m_min_count, true, AvailableMemory(), m_ncores);
                    if(kmer_counter.Kmers().Size() < 100) { // find a kmer length with at least 100 distinct kmers at that length
                        m_max_kmer -= read_len/25;          // reduce maximal kmer length by a small amount based on read length
                        continue;
                    }
                    double average_count_for_max_kmer = kmer_counter.AverageCount();
                    if(average_count_for_max_kmer >= maxkmercount)
                        break;
                    else 
                        m_max_kmer -= read_len/25;                                    
                }
                m_max_kmer = max(m_max_kmer, min_kmer);
                cerr << endl << "Average count: " << average_count << " Max kmer: " << m_max_kmer << endl;
            }
            
            //estimate insert size
            if(steps > 1 || m_usepairedends) {
                if(m_max_kmer_paired == 0) {
                    size_t mates = 0;
                    for(auto& rh : m_raw_reads)
                        mates += rh[0].ReadNum();
                    unsigned sample_size = 10000; // use 10000 reads for connecting to estimate insert size
                    unordered_set<size_t> selection;
                    if(mates/2 > 2*sample_size) {  // make random choice for reads
                        default_random_engine generator;
                        uniform_int_distribution<size_t> distribution(0,mates/2-1);
                        for(unsigned s = 0; s < sample_size; ) {
                            if(selection.insert(distribution(generator)).second) 
                                ++s;
                        }
                    } else if(mates/2 > 0) { // too few paired reads so using all : may be > sample_size but <= twice that size
                        for(size_t i = 0; i <= mates/2-1; ++i)
                            selection.insert(i);
                    }

                    if(!selection.empty()) {
                        CStopWatch timer;
                        timer.Restart();

                        list<array<CReadHolder,2>> mate_pairs;
                        size_t mp = 0;
                        int sub_sample = sample_size/m_ncores;
                        size_t num = 0;
                        for(auto& reads : m_raw_reads) {
                            for(CReadHolder::string_iterator is = reads[0].sbegin(); is != reads[0].send(); ++is, ++mp) {
                                if(selection.count(mp)) {
                                    if((num++)%sub_sample == 0)
                                        mate_pairs.push_back({CReadHolder(true), CReadHolder(false)});
                                    mate_pairs.back()[0].PushBack(is);
                                    mate_pairs.back()[0].PushBack(++is);
                                } else {
                                    ++is;
                                }
                            }
                        }
                
                        int long_insert_size = 2000; // we don't expect inserts to be longer than 2000 bp for this program
                        CDBGraphDigger graph_digger(*m_graphs[min_kmer], m_fraction, m_jump, m_low_count);
                        list<array<CReadHolder,2>> connected_mate_pairs = graph_digger.ConnectPairs(mate_pairs, long_insert_size, m_ncores);
                        CReadHolder connected_mates(false);
                        for(auto& mp : connected_mate_pairs) {
                            for(CReadHolder::string_iterator is = mp[0].sbegin(); is != mp[0].send(); ++is)
                                connected_mates.PushBack(is);
                        }

 
                        m_max_kmer_paired = connected_mates.N50();
                        cerr << endl << "N50 for inserts: " << m_max_kmer_paired << endl << endl;

                    }
                }  
                m_max_kmer_paired = min(m_max_kmer_paired,TKmer::MaxKmer());
                m_insert_size = 3*m_max_kmer_paired; // we don't expect spread of histogram to go beyond three times expected insert

                CleanReads();               
            }
                
            //main iterations
            if(m_steps > 1) {
                if(m_max_kmer > 1.5*m_min_kmer) {
                    double alpha = double(m_max_kmer-m_min_kmer)/(steps-1); // find desired distance between consecutive kmers
                    for(int step = 1; step < m_steps; ++step) {
                        int kmer_len = min_kmer+step*alpha+0.5;             // round to integer
                        kmer_len -= 1-kmer_len%2;                           // get odd kmer
                        if(GetGraph(kmer_len, m_raw_reads, true) == 0) {
                            cerr << "Empty graph for kmer length: " << kmer_len << " skipping this and longer kmers" << endl;
                            break;
                        }
                        ImproveContigs(kmer_len);
                        CleanReads();
                    }
                } else {
                    cerr << "WARNING: iterations are disabled" << endl;
                }
            }
            
            // three additional iterations with kmers (usually) longer than read length and upto insert size
            if(m_usepairedends && m_insert_size > 0 && m_max_kmer_paired > 1.5*m_max_kmer) {
                ConnectPairsIteratively();

                array<int,3> long_kmers;
                long_kmers[0] = 1.25*m_max_kmer;
                long_kmers[2] = m_max_kmer_paired;
                long_kmers[1] = (long_kmers[0]+long_kmers[2])/2;
                    
                for(int kmer_len : long_kmers) {
                    kmer_len -= 1-kmer_len%2;
                    if(GetGraph(kmer_len, m_connected_reads, false) == 0) {
                        cerr << "Empty graph for kmer length: " << kmer_len << " skipping this and longer kmers" << endl;
                        break;
                    }
                    ImproveContigs(kmer_len);
                }
            }                                              
        }        

        map<int,CDBGraph*>& Graphs() { return m_graphs; }
        TStrList& Contigs() { return m_contigs.back(); }
        vector<TStrList>& AllIterations() { return m_contigs; }
        CReadHolder ConnectedReads() const {
            CReadHolder connected_reads(false);
            for(const auto& cr : m_connected_reads) {
                for(CReadHolder::string_iterator is = cr[0].sbegin(); is != cr[0].send(); ++is)
                    connected_reads.PushBack(is);
            }
            return connected_reads;
        }

        virtual ~CDBGAssembler() {
            for(auto& graph : m_graphs)
                delete graph.second;    
        }

    private:
        // connects paired reads using all constructed de Bruijn graphs 
        void ConnectPairsIteratively() {
            for(auto& gr : m_graphs) {
                int kmer_len = gr.first;
                cerr << endl << "Connecting mate pairs using kmer length: " << kmer_len << endl;
                CDBGraphDigger graph_digger(*gr.second, m_fraction, m_jump, m_low_count);
                list<array<CReadHolder,2>> connected_reads_temp = graph_digger.ConnectPairs(m_raw_pairs, m_insert_size, m_ncores);
                list<array<CReadHolder,2>>::iterator pairedi = m_connected_reads.begin();
                list<array<CReadHolder,2>>::iterator rawi = m_raw_pairs.begin();
                for(auto& pr : connected_reads_temp) {
                    swap((*rawi)[0], pr[1]);                                                         // keep still not connected            
                    for(CReadHolder::string_iterator is = pr[0].sbegin() ;is != pr[0].send(); ++is)  // add new connected reads             
                        (*pairedi)[0].PushBack(*is);
                    ++rawi;
                    ++pairedi;
                }
            }

            size_t connected = 0;
            for(auto& rh : m_connected_reads)
                connected += rh[0].ReadNum();        
            cerr << "Totally connected: " << connected << endl;
        }

        // scans kmers for all assembled contigs and creates a map 
        // the key is the smaller of two possible kmer directions
        // the value is a tupe:
        //     int -     position on contig
        //     bool -    the same as the key or reverse complemented
        //     string* - pointer to the contig
        typedef TKmerMap<tuple<int, bool, const string*>> TKmerToContig;
        TKmerToContig GetAssembledKmers() const {
            int min_len = max(m_max_kmer_paired, m_max_kmer)+2*m_scan_window;
            int kmer_len = m_graphs.rbegin()->first;
            TKmerToContig assembled_kmers(kmer_len);
            size_t knum = 0;
            for(const string& contig : m_contigs.back()) {
                if((int)contig.size() >= min_len)
                    knum += contig.size()-kmer_len+1;
            }
            assembled_kmers.Reserve(knum);
            for(const string& contig : m_contigs.back()) {
                if((int)contig.size() >= min_len) {
                    CReadHolder rh(false);
                    rh.PushBack(contig);
                    int pos = contig.size()-kmer_len;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) { // iteration from last kmer to first  
                        TKmer kmer = *ik;
                        TKmer rkmer = revcomp(kmer, kmer_len);
                        if(kmer < rkmer)
                            assembled_kmers[kmer] = make_tuple(pos, true, &contig);
                        else
                            assembled_kmers[rkmer] = make_tuple(pos, false, &contig);
                    }
                }
            }

            return assembled_kmers;
        }

        // finds if a read belongs to any of the contigs
        // return tuple:
        //   int -     position on the contig (-1 if not found)
        //   int -     +1 if in positive strand; -1 if in negative strand
        //   string* - pointer to the contig 
        static tuple<int, int, const string*> FindMatchForRead(const CReadHolder::string_iterator& is, TKmerToContig& assembled_kmers) {
            int rlen = is.ReadLen();
            int kmer_len = assembled_kmers.KmerLen();

            int plus = 1;
            tuple<int, bool, const string*>* rsltp = 0;
            int knum = rlen-kmer_len+1;
            for(CReadHolder::kmer_iterator ik = is.KmersForRead(kmer_len); rsltp == 0 && knum > 0; --knum, ++ik) {
                TKmer kmer = *ik;
                TKmer rkmer = revcomp(kmer, kmer_len);
                TKmer* kmerp = &kmer;
                plus = 1;
                if(rkmer < kmer) {
                    kmerp = &rkmer;
                    plus = -plus;
                }
                rsltp = assembled_kmers.Find(*kmerp);
            }

            int pos = -1; // position on contig of the 'outer' read end (aka insert end)    
            const string* sp = 0;
            if(rsltp != 0) {
                sp = get<2>(*rsltp); // pointer to the contig
                if(!get<1>(*rsltp))
                    plus = -plus;
                if(plus > 0)
                    pos = get<0>(*rsltp)-knum;
                else
                    pos = get<0>(*rsltp)+kmer_len-1+knum;
            }
            
            return make_tuple(pos, plus, sp);
        }
    
        // removes reads if they belong to already assembled contigs
        // using contig sequence creates artificial connected pairs when both mates are placed
        //
        // assembled_kmers - a map of all kmers in already assembled contigs
        // margin - the minimal distance from an edge of a contig for a read to be removed
        // insert_size - the upper limit for insert size
        // raw_reads - reads
        // connected_reads - pointer to connected reads (nullp if not used)
        static void RemoveUsedReadsJob(TKmerToContig& assembled_kmers, int margin, int insert_size, array<CReadHolder,2>& raw_reads, CReadHolder* connected_reads) {
            int kmer_len = assembled_kmers.KmerLen();

            {
                CReadHolder cleaned_reads(true);
                CReadHolder::string_iterator is1 = raw_reads[0].sbegin();
                CReadHolder::string_iterator is2 = raw_reads[0].sbegin();
                ++is2;
                for( ; is2 != raw_reads[0].send(); ++is1, ++is1, ++is2, ++is2) {
                    if((int)min(is1.ReadLen(), is2.ReadLen()) < kmer_len) {
                        if(connected_reads) {    // keep short pairs for connection         
                            cleaned_reads.PushBack(is1);
                            cleaned_reads.PushBack(is2);                                 
                        } else {                 // give chance to be used as unpaired  
                            raw_reads[1].PushBack(is1);
                            raw_reads[1].PushBack(is2);
                        }
                        continue;
                    }

                    tuple<int, int, const string*> rslt1 = FindMatchForRead(is1, assembled_kmers);
                    int pos1 = get<0>(rslt1);
                    int plus1 = get<1>(rslt1);
                    const string* sp1 = get<2>(rslt1);
                    if(pos1 >= 0) {
                        int clen = sp1->size();
                        if((plus1 > 0 && pos1 >= margin && pos1+insert_size-1 < clen-margin) || (plus1 < 0 && pos1-insert_size+1 >= margin && pos1 < clen-margin))
                            continue;
                    }

                    // check for second mate in case first mate was of bad quality and not found in contigs 
                    tuple<int, int, const string*> rslt2 = FindMatchForRead(is2, assembled_kmers);
                    int pos2 = get<0>(rslt2);
                    int plus2 = get<1>(rslt2);
                    const string* sp2 = get<2>(rslt2);
                    if(pos2 >= 0) {
                        int clen = sp2->size();
                        if((plus2 > 0 && pos2 >= margin && pos2+insert_size-1 < clen-margin) || (plus2 < 0 && pos2-insert_size+1 >= margin && pos2 < clen-margin))
                            continue;
                    }

                    if(pos1 >= 0 && pos2 >= 0 && sp1 == sp2 && plus1 != plus2) { // same contig, different strands      
                        int clen = sp1->size();
                        if((plus1 > 0 && pos1 >= margin && pos2 < clen-margin) || (plus1 < 0 && pos2 >= margin && pos1 < clen-margin)) { // deep inside     
                            continue; 
                        } else if(connected_reads) {
                            if((plus1 > 0 && pos1 >= 0 && pos2 < clen) || (plus1 < 0 && pos2 >= 0 && pos1 < clen)) {                     // inside but not deep     
                                int a = min(pos1,pos2);
                                int b = max(pos1,pos2);
                                connected_reads->PushBack(sp1->substr(a,b-a+1));
                                continue;
                            }
                        }
                    }

                    cleaned_reads.PushBack(is1);
                    cleaned_reads.PushBack(is2);                                 
                }
                cleaned_reads.Swap(raw_reads[0]);
            }

            if(!connected_reads) {          
                CReadHolder cleaned_reads(false);
                for(CReadHolder::string_iterator is = raw_reads[1].sbegin() ;is != raw_reads[1].send(); ++is) {
                    int rlen = is.ReadLen();
                    if(rlen < kmer_len)
                        continue;        
            
                    tuple<int, int, const string*> rslt = FindMatchForRead(is, assembled_kmers);
                    int pos = get<0>(rslt);
                    int plus = get<1>(rslt);
                    const string* sp = get<2>(rslt);
                    if(pos >= 0) {
                        int clen = sp->size();
                        if((plus > 0 && pos >= margin && pos+rlen-1 < clen-margin) || (plus < 0 && pos-rlen+1 >= margin && pos < clen-margin))
                            continue;
                    }

                    cleaned_reads.PushBack(is);
                }            
                cleaned_reads.Swap(raw_reads[1]);
            }
        }

        // removes used reads from the read set used for de Bruijn graphs
        // assembled_kmers - a map of all kmers in already assembled contigs
        // margin - the minimal distance from an edge of a contig for a read to be removed
        // insert_size - the upper limit for insert size
        // ncores - number of threads
        // raw_reads - reads
        static void RemoveUsedReads(TKmerToContig& assembled_kmers, int margin, int insert_size, int ncores, list<array<CReadHolder,2>>& raw_reads) {
            list<function<void()>> jobs;
            for(auto& job_input : raw_reads) {
                jobs.push_back(bind(RemoveUsedReadsJob, ref(assembled_kmers), margin, insert_size, ref(job_input), (CReadHolder*)0));                
            }
            RunThreads(ncores, jobs);
        }

        // removes used reads from the read set used for pair connection and from already connected (by contig sequence) reads
        // assembled_kmers - a map of all kmers in already assembled contigs
        // margin - the minimal distance from an edge of a contig for a read to be removed
        // insert_size - the upper limit for insert size
        // ncores - number of threads
        // raw_reads - reads
        // connected_reads - already connected by contig sequence reads
        static void RemoveUsedPairs(TKmerToContig& assembled_kmers, int margin, int insert_size, int ncores, list<array<CReadHolder,2>>& raw_reads, list<array<CReadHolder,2>>& connected_reads) {
            list<function<void()>> jobs;
            auto icr = connected_reads.begin();
            for(auto& job_input : raw_reads) {
                jobs.push_back(bind(RemoveUsedReadsJob, ref(assembled_kmers), margin, insert_size, ref(job_input), &(*icr++)[1]));                
            }
            RunThreads(ncores, jobs);
        }

        // removes used reads from the read set used for de Bruijn graphs and from the read set used for pair connection 
        // removes paired reads not needed as they are already connected by contig sequence reads
        // creates new set of reads to use
        void CleanReads() {
            CStopWatch timer;
            timer.Restart();
            TKmerToContig assembled_kmers = GetAssembledKmers();
            cerr << "Contigs: " << m_contigs.back().size() << " Assembled kmers: " << assembled_kmers.Size() << endl;

            RemoveUsedReads(assembled_kmers, m_max_kmer+m_jump+m_scan_window, m_insert_size, m_ncores, m_raw_reads);
            RemoveUsedReads(assembled_kmers, m_jump+m_scan_window, m_insert_size, m_ncores, m_connected_reads);
            RemoveUsedPairs(assembled_kmers, m_jump+m_scan_window, m_insert_size, m_ncores, m_raw_pairs, m_connected_reads);

            size_t reads = 0;
            for(auto& rh : m_raw_reads)
                reads += rh[0].ReadNum()+rh[1].ReadNum();
            cerr << "Cleaned reads: " << reads << endl;
            reads = 0;
            for(auto& rh : m_raw_pairs)
                reads += rh[0].ReadNum()+rh[1].ReadNum();
            cerr << "Reads for connection: " << reads << endl;
            reads = 0;
            for(auto& rh : m_connected_reads)
                reads += rh[0].ReadNum()+rh[1].ReadNum();
            cerr << "Internal reads: " << reads << endl;
            cerr << "Reads cleaned in " << timer.Elapsed();                                               
        }

        // improves previously assembled contigs using a longer kmer
        void ImproveContigs (int kmer_len) {
            CDBGraph& graph = *m_graphs[kmer_len];
            cerr << "Kmer: " << kmer_len << " Graph size: " << graph.GraphSize() << " Contigs in: " << (m_contigs.empty() ? 0 : m_contigs.back().size()) << endl;
            CStopWatch total;
            total.Restart();
            
            CStopWatch timer;
            timer.Restart();
            //convert strings to SContig and mark visited kmers 
            TContigList scontigs = ConverToSContigAndMarkVisited(); 
            cerr << endl << "Mark used kmers in " << timer.Elapsed();

            timer.Restart();
            //create new contigs using not yet included kmers
            CDBGraphDigger graph_digger_no_jump(graph, m_fraction, 0, m_low_count);
            unsigned min_len_for_new_seeds = 3*kmer_len;      // short ones are likely to be noise
            TContigList new_seeds = graph_digger_no_jump.GenerateNewSeeds(min_len_for_new_seeds, m_ncores);
            cerr << "New seeds in " << timer.Elapsed();
            cerr << "New seeds: " << new_seeds.size() << endl;

            //add new seeds
            scontigs.splice(scontigs.end(), new_seeds);
            for(auto& contig : scontigs)
                contig.SelectMinDirection();
            scontigs.sort();

            timer.Restart();
            CDBGraphDigger graph_digger(graph, m_fraction, m_jump, m_low_count);
            graph_digger.ConnectAndExtenContigs(scontigs, m_scan_window, m_ncores); 

            // convert back to strings  
            m_contigs.push_back(TStrList());
            for(auto& contig : scontigs) {
                string s(contig.m_seq.begin(), contig.m_seq.end());
                m_contigs.back().push_back(s);
            }
            cerr << "Connections and extensions in " << timer.Elapsed();

            vector<size_t> contigs_len;
            size_t genome_len = 0;
            for(auto& contig : m_contigs.back()) {
                contigs_len.push_back(contig.size());
                genome_len += contigs_len.back();                            
            }
            sort(contigs_len.begin(), contigs_len.end());
            size_t n50 = 0;
            int l50 = 0;
            size_t len = 0;
            for(int j = (int)contigs_len.size()-1; j >= 0 && len < 0.5*genome_len; --j) {
                ++l50;
                n50 = contigs_len[j];
                len += contigs_len[j];
            }
            cerr << "Contigs out: " << contigs_len.size() << " Genome: " << genome_len << " N50: " << n50 << " L50: " << l50 << endl; 
            cerr << "Assembled in " << total.Elapsed() << endl; 
        }


        // converts contigs from the previous iteration into SContig and marks visited the nodes in the graph
        TContigList ConverToSContigAndMarkVisited() {
            TContigList scontigs;
            if(m_contigs.empty())
                return scontigs;

            vector<pair<const string*, SAtomic<uint8_t>>> contig_is_taken;
            for(const auto& contig : m_contigs.back())
                contig_is_taken.push_back(make_pair(&contig,SAtomic<uint8_t>(0)));
            vector<TContigList> scontigs_for_threads(m_ncores);
            list<function<void()>> jobs;
            for(auto& sc : scontigs_for_threads) 
                jobs.push_back(bind(&CDBGAssembler::ConverToSContigAndMarkVisitedJob, this, ref(contig_is_taken), ref(sc)));
            RunThreads(m_ncores, jobs);

            for(auto& sc : scontigs_for_threads)
                scontigs.splice(scontigs.end(), sc);
    
            return scontigs;
        }

        // one-thread worker for ConverToSContigAndMarkVisited()
        void ConverToSContigAndMarkVisitedJob(vector<pair<const string*, SAtomic<uint8_t>>>& contig_is_taken, TContigList& scontigs) {
            CDBGraph& graph = *m_graphs.rbegin()->second; // last graph
            int kmer_len = graph.KmerLen();
            for(auto& cnt : contig_is_taken) {
                if(!cnt.second.Set(1))
                    continue;
                const string& contig = *cnt.first;
                if((int)contig.size() >= kmer_len+2*m_scan_window)
                    scontigs.push_back(SContig(contig, graph)); // constructor sets visited in graph           
            }
        }

        // estimates available memory
        int64_t AvailableMemory() const {
            int64_t GB = 1000000000;
            int64_t mem_available = GB*m_memory;
            int64_t mem_used = 0;
            for(const auto& reads : m_raw_reads)
                mem_used += reads[0].MemoryFootprint()+reads[1].MemoryFootprint();
            for(const auto& reads : m_raw_pairs)
                mem_used += reads[0].MemoryFootprint()+reads[1].MemoryFootprint();
            for(const auto& reads : m_connected_reads)
                mem_used += reads[0].MemoryFootprint()+reads[1].MemoryFootprint();
            for(auto& graph : m_graphs)
                mem_used += graph.second->MemoryFootprint();

            return mem_available-mem_used;

        }

        // counts kmers and build a de Bruijn graph; returns average count of kmers in the graph
        // kmer_len - the size of the kmer
        // reads - reads from input or connected internally
        // is_stranded - whether or not stranded information is meaningful
        double GetGraph(int kmer_len, const list<array<CReadHolder,2>>& reads, bool is_stranded) {
            CKmerCounter kmer_counter(reads, kmer_len, m_min_count, is_stranded, AvailableMemory(), m_ncores);
            if(kmer_counter.Kmers().Size() == 0)
                throw runtime_error("Insufficient coverage");
                
            double average_count = kmer_counter.AverageCount();
            kmer_counter.GetBranches();
            TKmerCount& sorted_kmers =  kmer_counter.Kmers();
            if(sorted_kmers.Size() == 0)
                return 0;

            map<int,size_t> hist;
            for(size_t index = 0; index < sorted_kmers.Size(); ++index) {
                ++hist[sorted_kmers.GetCount(index)];                  // count clipped to integer automatically
            }
            TBins bins(hist.begin(), hist.end());
            m_graphs[kmer_len] = new CDBGraph(move(sorted_kmers), move(bins), is_stranded);

            return average_count;
        }


        double m_fraction;                                   // Maximal noise to signal ratio of counts acceptable for extension
        int m_jump;                                          // minimal length of accepted dead ends
        int m_low_count;                                     // minimal kmer count to be included in a contig
        int m_steps;                                         // number of main steps
        int m_min_count;                                     // minimal kmer count to be included in a de Bruijn graph
        int m_min_kmer;                                      // the minimal kmer size for the main steps
        bool m_usepairedends;                                // whether or not to use paired ends
        int m_max_kmer_paired;                               // insert size
        int m_insert_size;                                   // upper bound for the insert size
        int m_maxkmercount;                                  // the minimal average count for estimating the maximal kmer
        int m_memory;                                        // the upper bound for memory use (GB)
        int m_ncores;                                        // number of threads

        int m_scan_window;                                   // the size-1 of the contig's flank area used for extensions and connections
        int m_max_kmer;                                      // maximal kmer size for the main steps

        list<array<CReadHolder,2>>& m_raw_reads;             // original reads - will be reduced gradually
        list<array<CReadHolder,2>> m_raw_pairs;              // paired original reads for connection - will be reduced gradually
        list<array<CReadHolder,2>> m_connected_reads;        // connected pairs (long reads)
        map<int,CDBGraph*> m_graphs;                         // De Bruijn graphs for mutiple kmers
        vector<TStrList> m_contigs;                          // assembled contigs for each iteration
    };


}; // namespace
#endif /* _DBGAssembler_ */
