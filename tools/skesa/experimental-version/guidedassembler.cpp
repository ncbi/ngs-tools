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

#include <unordered_map>
#include <math.h>
#include "readsgetter.hpp"
#include "counter.hpp"
#include "guidedpath.hpp"

using namespace boost::program_options;
using namespace DeBruijn;
using namespace std;


class CGuidedAssembler {
public:
    enum {eRedundant = 1, eLeftDeadEnd = 2, eRightDeadEnd = 4, eCircular = 8}; 
    typedef list<tuple<int, TRange, string>> TargetInfo;  // score, contig range, acc
    struct SContigData {
        int m_status = 0;
        TargetInfo m_tinfo;

        TRange TargetRange() const {
            TRange target_range = get<1>(m_tinfo.front());
            for(auto& inf : m_tinfo) {
                target_range.first = min(target_range.first, get<1>(inf).first);
                target_range.second = max(target_range.second, get<1>(inf).second);
            }
            return target_range;
        }
        string FastaDefLine() const {
            string def;
            for(auto& info : m_tinfo) {
                int score = get<0>(info);
                TRange range = get<1>(info);
                string acc = get<2>(info);
                if(info != m_tinfo.front())
                    def += " ";
                def += acc+":"+to_string(range.first)+":"+to_string(range.second)+":"+to_string(score);
            }
            if(m_status&eCircular)
                def += " Cir";
            if(m_status&eLeftDeadEnd)
                def += " LD";
            if(m_status&eRightDeadEnd)
                def += " RD";
            return def;
        }
        bool operator<(const SContigData& a) const { return m_tinfo < a.m_tinfo; }
    };
    typedef map<string, SContigData> TContigs;

    CGuidedAssembler(int kmer_len, int min_count, double fraction, int memory, int match, int mismatch, int gap_open, int gap_extend, int drop_off, int ncores, list<array<CReadHolder,2>>& raw_reads) : 
        m_kmer_len(kmer_len), m_match(match), m_mismatch(mismatch), m_gap_open(gap_open), m_gap_extend(gap_extend), m_drop_off(drop_off), m_delta(match, mismatch), 
        m_ncores(ncores), m_raw_reads(raw_reads) {

        
        //read fasta
        {
            char c;
            if(!(cin >> c) || c != '>')
                throw runtime_error("Invalid fasta file format for targets");
            string record;
            while(getline(cin, record, '>')) {
                size_t first_ret = min(record.size(),record.find('\n'));
                if(first_ret == string::npos)
                    throw runtime_error("Invalid fasta file format for targets");
                string acc = record.substr(0, first_ret);
                acc = acc.substr(0, acc.find_first_of(" \t"));
                string target = record.substr(first_ret+1);
                target.erase(remove(target.begin(),target.end(),'\n'),target.end());
                for(char& c : target) c = toupper(c);
                m_targets.emplace_back(target,acc,0);
            }
        }

        int64_t GB = 1000000000;
        int jump = 50;
        int low_count = min_count;
        {
            CKmerCounter kmer_counter(m_raw_reads, m_kmer_len, min_count, true, GB*memory, m_ncores);
            if(kmer_counter.Kmers().Size() == 0)
                throw runtime_error("Insufficient coverage");
            m_average_count = kmer_counter.AverageCount();
            cerr << "Average count: " << m_average_count << endl;
            kmer_counter.GetBranches();
            TKmerCount& sorted_kmers = kmer_counter.Kmers();        
            map<int,size_t> hist;
            for(size_t index = 0; index < sorted_kmers.Size(); ++index) {
                ++hist[sorted_kmers.GetCount(index)];                  // count clipped to integer automatically        
            }
            TBins bins(hist.begin(), hist.end());
            m_graphp.reset(new CDBGraph(move(sorted_kmers), move(bins), true));
            m_graphdiggerp.reset(new CDBGraphDigger(*m_graphp, fraction, jump, low_count));
        }

        CStopWatch timer;
        timer.Restart();
        { // hash for kmers
            m_kmer_hash.resize(size_t(numeric_limits<uint16_t>::max())+1);

            for(size_t index = 0; index < m_graphp->GraphSize(); ++index) {
                CDBGraph::Node node = 2*(index+1);
                const uint16_t* wordp = (uint16_t*)m_graphp->getPointer(node);  // points to the last 8 sybols of kmer stored in graph
                m_kmer_hash[*wordp].emplace_front(node);  // LAST 8 symbols for forward kmer
                uint16_t rev_word;
                uint16_t* rev_wordp = &rev_word;
                *(uint8_t*)rev_wordp = revcomp_4NT [*((uint8_t*)wordp+1)];
                *((uint8_t*)rev_wordp+1) = revcomp_4NT [*(uint8_t*)wordp];
                m_kmer_hash[*rev_wordp].emplace_front(m_graphp->ReverseComplement(node));  // FIRST 8 symbols for reversed kmer
            }    
        }
        cerr << "Graph hash in " << timer.Elapsed();

        timer.Restart();
        { // assemble
            vector<TContigs> thread_rslts(ncores); // [seq] score, blen, acc
            
            list<function<void()>> jobs;
            for(int thr = 0; thr < m_ncores; ++thr) {
                jobs.push_back(bind(&CGuidedAssembler::AssemblerJob, this, ref(thread_rslts[thr])));
            }
            RunThreads(m_ncores, jobs);

            for(auto& trslt : thread_rslts) {
                for(auto& rslt : trslt) {
                    SContigData& cdata = m_rslts[rslt.first];
                    cdata.m_status |= rslt.second.m_status;
                    auto& lst = cdata.m_tinfo;
                    AddInfo(lst, rslt.first.size(), rslt.second.m_tinfo);
                }
            }
        }
        cerr << "Assembling in " << timer.Elapsed();

        //replace remaining SNPs by ambiguous symbols
        timer.Restart();
        {
            map<int, list<pair<deque<char>,SContigData>>> len_contigs;
            for(auto& rslt : m_rslts) {
                const string& contig = rslt.first;
                SContigData& info = rslt.second;
                len_contigs[contig.size()].emplace_back(deque<char>(contig.begin(),contig.end()),info);
            }
            m_rslts.clear();
            
            list<function<void()>> jobs;
            for(auto& len_group : len_contigs)
                jobs.push_back(bind(&CGuidedAssembler::FindSNPsJob, this, ref(len_group.second)));
            RunThreads(ncores, jobs);

            for(auto& len_group : len_contigs) {
                for(auto& mem : len_group.second) {
                    string contig(mem.first.begin(), mem.first.end());
                    m_rslts[contig] = mem.second;
                }
            }                
        }
        cerr << "Additional SNP detection in in " << timer.Elapsed();

        timer.Restart();
        {
            vector<SAtomic<uint8_t>> centinel;
            size_t total_reads = 0;
            list<function<void()>> jobs;
            for(auto& job_input : m_raw_reads) {
                total_reads += job_input[0].ReadNum()+job_input[1].ReadNum();
                jobs.push_back(bind(&CGuidedAssembler::HashForReadsJob, this, ref(centinel), ref(job_input)));
            }
            m_reads_hash_table.resize(total_reads);
            centinel.resize(total_reads);
            RunThreads(ncores, jobs);
        }
        cerr << "Reads hash in " << timer.Elapsed();
               
        timer.Restart();
        {
            list<TContigs> clean_rslts_from_threads;
            list<function<void()>> jobs;
            for(TContigs::iterator contigp = m_rslts.begin(); contigp != m_rslts.end(); ++contigp) {
                clean_rslts_from_threads.push_back(TContigs());
                jobs.push_back(bind(&CGuidedAssembler::CleanContigJob, this, contigp, ref(clean_rslts_from_threads.back())));
            }
            RunThreads(ncores, jobs);

            m_rslts.clear();
            for(auto& from_thread : clean_rslts_from_threads) {
                for(auto& rslt : from_thread) {
                    m_rslts[rslt.first] = rslt.second;
                }
            }
        }
        cerr << "Clean contigs in " << timer.Elapsed();                
        
        timer.Restart();
        {
            TUtilMap util;
            m_kmer_for_duplicates = MaxPrec*32;
            for(auto& rslt : m_rslts) {
                const string& contig = rslt.first;
                util.insert({contig, TUtilInfo()}); // centinel = 0, dupl = false, deque empty, inheritance empty
                m_kmer_for_duplicates = min(m_kmer_for_duplicates, (int)floor((double)contig.size()/3)); // at least 3 kmers for contig
            }

            cerr << "Kmer for duplicates: " << m_kmer_for_duplicates << endl;

            int hash_table_size = 100000;
            vector<TContigShifts> hash_contigs(hash_table_size);
            {
                // kmer hash values 
                list<function<void()>> jobs;
                for(int thr = 0; thr < ncores; ++thr) {
                    jobs.push_back(bind(&CGuidedAssembler::KmerHashValuesJob, this, ref(util), hash_table_size));
                }
                RunThreads(ncores, jobs);
            }

            // hash map
            for(auto it = m_rslts.begin(); it != m_rslts.end(); ++it) {
                deque<int>& hash_valuesp = get<1>(util[it->first])[0];
                for(unsigned pos = 0; pos < hash_valuesp.size(); ++pos)
                    hash_contigs[hash_valuesp[pos]].emplace_front(it,pos);                
            }            

            { // mark and remove identical substrings
                for(auto& util_info : util)
                    get<0>(util_info.second) = 0;   // clear centinel
                
                list<function<void()>> jobs;
                for(int thr = 0; thr < ncores; ++thr) {
                    jobs.push_back(bind(&CGuidedAssembler::MarkDuplicatesJob, this, ref(util), ref(hash_contigs)));
                }
                RunThreads(ncores, jobs);
                
                //transfer info from duplicates
                bool keep_doing = true;
                while(keep_doing) {
                    keep_doing = false;
                    for(auto& util_info : util) {
                        TContigShifts& inherit = get<2>(util_info.second);
                        for(auto infp : inherit) {
                            TargetInfo lst = m_rslts[util_info.first].m_tinfo;
                            int shift = infp.second;
                            for(auto& info : lst) {
                                TRange& range = get<1>(info);
                                if(shift >= 0) {
                                    range.first += shift;
                                    range.second += shift;
                                } else {
                                    swap(range.first, range.second);
                                    range.first = -shift-range.first;
                                    range.second = -shift-range.second;
                                }
                            }
                            keep_doing = (AddInfo(infp.first->second.m_tinfo, infp.first->first.size(), lst) || keep_doing);
                        }
                    }                
                }

                //remove duplicates
                for(auto it_loop = m_rslts.begin(); it_loop != m_rslts.end(); ) {
                    auto it = it_loop++;
                    if(it->second.m_status&eRedundant)
                        m_rslts.erase(it);                   
                }
            }
        }            
        cerr << "Remove duplicates in " << timer.Elapsed();
    }

    TContigs Contigs() const { return m_rslts; }
    CDBGraph& Graph() { return *m_graphp; }

private:

    typedef forward_list<pair<TContigs::iterator, int>> TContigShifts;
    typedef tuple<SAtomic<uint8_t>, deque<int>[2], TContigShifts> TUtilInfo; // centinel, kmers hash values(+/-), inheritance
    typedef map<string, TUtilInfo> TUtilMap;

    typedef forward_list<CReadHolder::string_iterator> TSimpleHashElement;
    typedef vector<TSimpleHashElement> TSimpleHash;

    enum {eLeftComplete = 1, eRightComplete = 2, eBothComplete = 3}; // completeness direction relative contig
    typedef pair<CReadHolder::string_iterator, int> TKeyForHits;  // key - [string_iterator to read, strand]
    typedef tuple<int, int , int, forward_list<char>, int> THitPosition;      // contig start/stop, read start or read end, list of ambigous choiced, remaining length
    typedef forward_list<THitPosition> TValueForHits;                        // value - list of {contig start/stop, read start or read end, list of ambigous choiced}
    struct SHashForHits { size_t operator() (const TKeyForHits& read_strand) const { return read_strand.first.Hash()^hash<int>()(read_strand.second); }};
    typedef unordered_map<TKeyForHits, TValueForHits, SHashForHits> TReadHits;

    static bool AddInfo(TargetInfo& to, int to_len, const TargetInfo& from) {
        bool changed = false;
        for(auto& inf_from : from) {
            bool included = false;
            TRange range_from = get<1>(inf_from);
            if(range_from.first < 0) {
                range_from.first += to_len;
                range_from.second += to_len;
            }
            for(auto& inf_to : to) {
                TRange range_to = get<1>(inf_to);
                if(range_to.first < 0) {
                    range_to.first += to_len;
                    range_to.second += to_len;
                }
                if(range_to == range_from) { // same interval   
                    if(get<2>(inf_from) > get<2>(inf_to)) {
                        inf_to = inf_from;
                        changed = true;
                    }
                    included = true;                    
                } else if(range_to.first <= range_from.first && range_to.second >= range_from.second) { // to larger   
                    included = true;
                } else if(range_to.first >= range_from.first && range_to.second <= get<1>(inf_from).second) { // from larger   
                    inf_to = inf_from;
                    changed = true;
                    included = true;
                }                
            }
            if(!included) {
                to.push_back(inf_from);
                changed = true;
            }
        }
        if(changed) {
            to.sort();
            to.erase(unique(to.begin(),to.end()), to.end());
            reverse(to.begin(), to.end());        
        }
        return changed;
    }

    void FindSNPsJob(list<pair<deque<char>, SContigData>>& group) {
        // find potentially circular
        list<pair<deque<char>, SContigData>> circulars;
        for(auto it_loop = group.begin(); it_loop != group.end(); ) {
            auto it = it_loop++;
            auto& contig = it->first;
            if(equal(contig.begin(), contig.begin()+m_kmer_len-1, contig.end()-m_kmer_len+1))
                circulars.splice(circulars.end(), group, it);
        }

        //rotate circular candidates
        circulars.sort([](const pair<deque<char>,SContigData>& a, pair<deque<char>,SContigData>& b){ return a.second.m_tinfo > b.second.m_tinfo; });
        for(auto it_loop = circulars.begin(); it_loop != circulars.end(); ) {
            auto it = it_loop++;
            auto& contigi = it->first;
            contigi.erase(contigi.end()-m_kmer_len+1, contigi.end());
            it->second.m_status |= eCircular;
            int circular_length = contigi.size();

            auto RotateRange = [=](TRange& range, int shift) {
                range.first += shift;
                range.second += shift;
                if(range.first < 0)
                   range.first += circular_length; 
                else if(range.first >= circular_length)
                    range.first -= circular_length;
                if(range.second < 0)
                   range.second += circular_length; 
                if(range.second >= circular_length)
                    range.second -= circular_length;
                if(range.first > range.second)
                    range.first -= circular_length; // negative
            };

            // recalculate range because of clip
            TRange& rangei = get<1>(it->second.m_tinfo.front());
            RotateRange(rangei, 0);

            group.splice(group.end(), circulars, it); // contigi is still valid
            vector<char> query(contigi.begin(), contigi.end());
            for(auto jt_loop = it_loop; jt_loop != circulars.end(); ) {
                auto jt = jt_loop++;
                auto& contigj = jt->first;
                vector<char> subject(contigj.begin(), contigj.end());
                int gap_open = 10*m_match*circular_length; // we don't want indels
                int gap_extend = gap_open;
                CCigar cigar = LclAlign(&query[0], circular_length, &subject[0], circular_length, gap_open, gap_extend, m_delta.matrix);
                if(cigar.QueryRange().second-cigar.QueryRange().first+1 > circular_length/3) { // reasonably good to try
                    int qleft = cigar.QueryRange().first;
                    int sleft = cigar.SubjectRange().first;
                    int snew_origin = sleft-qleft;
                    if(snew_origin < 0)
                        snew_origin = circular_length+snew_origin;
                    deque<char> new_contigj(circular_length);
                    rotate_copy(contigj.begin(), contigj.begin()+snew_origin, contigj.end()-m_kmer_len+1, new_contigj.begin());
                    int diff = 0;
                    int ambig_diff = 0;
                    for(int i = 0; i < circular_length; ++i) {
                        if(contigi[i] != new_contigj[i]) {
                            ++diff;
                            if(!MatchWithAmbiguousDNA(contigi[i], new_contigj[i]))
                                ++ambig_diff;
                        }
                    }
                    if(ambig_diff < 0.1*circular_length) {
                        // recalculate range because of rotation
                        TRange& rangej = get<1>(jt->second.m_tinfo.front());
                        RotateRange(rangej, -snew_origin);

                        if(diff > 0) {
                            group.emplace_back(new_contigj, jt->second);
                        } else {
                            AddInfo(group.back().second.m_tinfo, circular_length, jt->second.m_tinfo);
                            group.back().second.m_status |= jt->second.m_status;
                        }
                        
                        if(jt == it_loop)
                            ++it_loop;
                        circulars.erase(jt);
                    }
                }
            }
        }

        bool keep_doing = (group.size() > 1);
        int len = group.front().first.size();
        while(keep_doing) {
            keep_doing = false;
            for(int l = 0; l < len; ++l) {
                //rotate    
                for(auto& mem : group) {
                    auto& contig = mem.first;
                    contig.push_back(contig.front());
                    contig.pop_front();
                }
                //sort  
                group.sort();
                //combine seqs with equal prefix
                for(auto it = group.begin(); it != group.end(); ++it) {
                    auto& contigi = it->first;
                    set<char> last_symbols;
                    for(auto jt_loop = it; jt_loop != group.end(); ) {
                        auto jt = jt_loop++;
                        if(jt != it) {
                            auto& contigj = jt->first;
                            if(equal(contigi.begin(), contigi.end()-1, contigj.begin())) { // only last symbol different    
                                string ls = FromAmbiguousIUPAC[contigj.back()];
                                last_symbols.insert(ls.begin(),ls.end());
                                AddInfo(it->second.m_tinfo, contigi.size(), jt->second.m_tinfo);
                                it->second.m_status |= jt->second.m_status;
                                group.erase(jt);
                                keep_doing = true;
                            } else {
                                break;
                            }
                        }
                    }
                    if(!last_symbols.empty()) {
                        string ls = FromAmbiguousIUPAC[contigi.back()];
                        last_symbols.insert(ls.begin(),ls.end());
                        string final_ls(last_symbols.begin(),last_symbols.end());
                        contigi.back() = ToAmbiguousIUPAC[final_ls];
                    }
                }
            }
        }
    }

    int AlignReads(const string& contig, const vector<int>& ambig_positions, TReadHits& hits) const {
        string variant = contig;
        // scan contig and find all m_kmer_len hits to reads
        int contig_len = contig.size();
        for(int p = 0; p <= contig_len-m_kmer_len; ++p) {
            int kmer_left = p;
            int kmer_right = p+m_kmer_len-1;

            stack<pair<int,char>> remaining_choices;
            for(int ap : ambig_positions) {
                if(ap < kmer_left)
                    continue;
                if(ap > kmer_right)
                    break;
                string ambig = FromAmbiguousIUPAC[contig[ap]];
                variant[ap] = ambig[0];
                for(unsigned i = 1; i < ambig.size(); ++i)
                    remaining_choices.emplace(ap, ambig[i]);
            }

            auto FindHits = [&](const TKmer& kmer, int strand, int kmer_len) {
                for(CReadHolder::string_iterator is : m_reads_hash_table[kmer.oahash()%m_reads_hash_table.size()]) {
                    int rlen = is.ReadLen();
                    CReadHolder::kmer_iterator ik = is.KmersForRead(kmer_len);
                    TKeyForHits key(is, strand);
                    if(*ik == kmer)                                // read end
                        hits[key].emplace_front(kmer_left, kmer_right, strand > 0 ? eRightComplete : eLeftComplete, forward_list<char>(), rlen-m_kmer_len);
                    else if(*(ik += rlen-kmer_len) == kmer)       // read start
                        hits[key].emplace_front(kmer_left, kmer_right, strand > 0 ? eLeftComplete : eRightComplete, forward_list<char>(), rlen-m_kmer_len);                
                }
            };

            TKmer first_kmer(variant.begin()+kmer_left, variant.begin()+kmer_right+1);
            FindHits(first_kmer, 1, m_kmer_len);
            FindHits(revcomp(first_kmer, m_kmer_len), -1, m_kmer_len);

            while(!remaining_choices.empty()) {
                int current_ap = remaining_choices.top().first;
                variant[current_ap] = remaining_choices.top().second;
                remaining_choices.pop();
                for(int ap : ambig_positions) {
                    if(ap <= current_ap)
                        continue;
                    if(ap > kmer_right)
                        break;
                    string ambig = FromAmbiguousIUPAC[contig[ap]];
                    variant[ap] = ambig[0];
                    for(unsigned i = 1; i < ambig.size(); ++i)
                        remaining_choices.emplace(ap, ambig[i]);
                }
                TKmer second_kmer(variant.begin()+kmer_left, variant.begin()+kmer_right+1);
                FindHits(second_kmer, 1, m_kmer_len);
                FindHits(revcomp(second_kmer, m_kmer_len), -1, m_kmer_len);
            }
        }

        int max_read_len = 0;
        // expand hits to first mismatch
        for(auto& hit : hits) {
            const CReadHolder::string_iterator& is = hit.first.first;
            int rlen = is.ReadLen();
            max_read_len = max(max_read_len, rlen);
            int strand = hit.first.second;
            auto& positions = hit.second;
            positions.sort();
            string read = *is;
            if(strand < 0)
                ReverseComplementSeq(read.begin(), read.end());

            for(THitPosition& pos : positions) {
                int& hit_left = get<0>(pos);
                int& hit_right = get<1>(pos);
                int& etype = get<2>(pos);
                forward_list<char>& ambig_choices = get<3>(pos);
                int& remaining_len = get<4>(pos);
                if(etype == eLeftComplete) {  // scan to right
                    int r = 0;
                    int p = hit_left;
                    for( ; r < rlen && p < contig_len; ++r, ++p) {
                        string ambig = FromAmbiguousIUPAC[contig[p]];
                        if(ambig.find(read[r]) == string::npos)
                            break;
                        hit_right = p;
                        if(ambig.size() > 1)
                            ambig_choices.emplace_front(read[r]);  // reverse order!                                                          
                    }
                    ambig_choices.reverse();
                    if(r == rlen) { // full align - delete matching hit
                        etype = eBothComplete;
                        remaining_len = 0;
                        positions.remove_if([=](const THitPosition& hp) { return  get<2>(hp) == eRightComplete && get<1>(hp) == hit_right; });
                    } else {
                        TKmer rkmer(read.begin()+r-m_kmer_len, read.begin()+r);
                        CDBGraph::Node rnode = m_graphp->GetNode(rkmer);
                        
                        remaining_len = 0;
                        while(remaining_len < rlen-r) {
                            vector<CDBGraph::Successor> successors = m_graphdiggerp->GetReversibleNodeSuccessors(rnode);
                            char next_symbol = read[r+remaining_len];
                            auto rslt = find_if(successors.begin(), successors.end(), [=](const CDBGraph::Successor& suc){return suc.m_nt == next_symbol;});
                            if(rslt != successors.end()) {
                                ++remaining_len;
                                rnode = rslt->m_node;
                            } else {
                                break;
                            }
                        }
                    }
                } else {  // scan to left - no matching hit (partial)
                    int p = hit_right;
                    for(int r = rlen-1; r >= 0 && p >= 0; --r, --p) {
                        string ambig = FromAmbiguousIUPAC[contig[p]];
                        if(ambig.find(read[r]) == string::npos) {
                            TKmer lkmer(read.begin()+r+1, read.begin()+r+m_kmer_len+1);
                            CDBGraph::Node lnode = m_graphp->ReverseComplement(m_graphp->GetNode(lkmer));

                            remaining_len = 0;
                            while(remaining_len < r+1) {
                                vector<CDBGraph::Successor> successors = m_graphdiggerp->GetReversibleNodeSuccessors(lnode);
                                char next_symbol = Complement(read[r-remaining_len]);
                                auto rslt = find_if(successors.begin(), successors.end(), [=](const CDBGraph::Successor& suc){return suc.m_nt == next_symbol;});
                                if(rslt != successors.end()) {
                                    ++remaining_len;
                                    lnode = rslt->m_node;
                                } else {
                                    break;
                                }
                            }

                            break;
                        }
                        hit_left = p;                        
                        if(ambig.size() > 1)
                            ambig_choices.emplace_front(read[r]); // direct order!                        
                    }
                }
            }
        }

        return max_read_len;
    }

    typedef tuple<int,int,int,int> THitExtend;  // position, extend (hit length-1), expected extend, read/pair length
    void SelectHitsLtoR(const string& variant, const vector<int>& ambig_positions, const TReadHits& hits, vector<THitExtend>& hits_for_variant) const {
        auto ClipRight = [&] (const THitPosition& pos, int left, int& right) {
            auto ambig = get<3>(pos).begin();
            if(ambig != get<3>(pos).end()) {
                for(int p : ambig_positions) {
                    if(p < left)
                        continue;
                    if(p > right)
                        break;
                    if(*(ambig++) != variant[p]) {
                        right = p-1;
                        break;
                    }
                }
            }
        };

        hits_for_variant.clear();
        hits_for_variant.reserve(hits.size());
        for(auto& hit : hits) {
            const auto& positions = hit.second;
            for(const THitPosition& pos : positions) {
                int left = get<0>(pos);
                if(left > 0 && !(get<2>(pos)&eLeftComplete))
                    continue;
                int right = get<1>(pos);
                int remaining_len = get<4>(pos);
                int expected_right = (get<2>(pos)&eRightComplete) ? right : right+remaining_len;   // one of the ends is always complete
                int rlen = right-left+1+remaining_len;
                ClipRight(pos, left, right);
                if(left > right)
                    continue;
                hits_for_variant.emplace_back(left, right-left, expected_right-left, rlen);
            }
        }    
        sort(hits_for_variant.begin(), hits_for_variant.end(), [](const THitExtend& a, const THitExtend& b) { return get<0>(a) < get<0>(b); });
    }

    void SelectHitsRtoL(const string& variant, const vector<int>& ambig_positions, const TReadHits& hits, vector<THitExtend>& hits_for_variant) const {
        auto ClipLeft = [&] (const THitPosition& pos, int& left, int right) {
            auto ambig = get<3>(pos).begin();
            if(ambig != get<3>(pos).end()) {
                for(int p : ambig_positions) {
                    if(p < left)
                        continue;
                    if(p > right)
                        break;
                    if(*(ambig++) != variant[p])
                        left = p+1;
                }
            }
        };

        hits_for_variant.clear();
        hits_for_variant.reserve(hits.size());
        for(auto& hit : hits) {
            const auto& positions = hit.second;
            for(const THitPosition& pos : positions) {
                int right = get<1>(pos);
                if(right < (int)variant.size()-1 && !(get<2>(pos)&eRightComplete))
                    continue;
                int left = get<0>(pos);
                int remaining_len = get<4>(pos);
                int expected_left = (get<2>(pos)&eLeftComplete) ? left : left-remaining_len; // one of the ends is always complete
                int rlen = right-left+1+remaining_len;
                ClipLeft(pos, left, right);
                if(left > right)
                    continue;
                hits_for_variant.emplace_back(right, right-left, right-expected_left, rlen);
            }
        }    
        sort(hits_for_variant.begin(), hits_for_variant.end(), [](const THitExtend& a, const THitExtend& b) { return get<0>(a) > get<0>(b); });
    }

    int SelectPairedHitsLtoR(const string& variant, const vector<int>& ambig_positions, const TReadHits& hits, vector<THitExtend>& hits_for_variant) {
        auto ClipRight = [&] (const THitPosition& pos, int left, int& right) {
            auto ambig = get<3>(pos).begin();
            if(ambig != get<3>(pos).end()) {
                for(int p : ambig_positions) {
                    if(p < left)
                        continue;
                    if(p > right)
                        break;
                    if(*(ambig++) != variant[p]) {
                        right = p-1;
                        return true;
                    }
                }
            }
            return false;
        };

        hits_for_variant.clear();
        hits_for_variant.reserve(hits.size());
        vector<int> paired_lengths;
        unordered_set<TKeyForHits, SHashForHits> included_mates;
        for(auto hitp_loop = hits.begin(); hitp_loop != hits.end(); ++hitp_loop) {
            auto hitap = hitp_loop;
            CReadHolder::string_iterator isa = hitap->first.first;
            int strand = hitap->first.second;
            if(!included_mates.emplace(isa, strand).second || !isa.HasMate())
                continue;

            CReadHolder::string_iterator isb = isa.GetMate();
            auto hitbp = hits.find(make_pair(isb,-strand));
            if(hitbp == hits.end())
                continue;

            included_mates.emplace(isb, -strand);
            if(strand < 0) {
                swap(hitap, hitbp);
                strand = -strand;
            }

            const auto& positionsa = hitap->second;
            for(const THitPosition& posa : positionsa) {
                int lefta = get<0>(posa);
                if(lefta > 0 && !(get<2>(posa)&eLeftComplete))
                    continue;

                int righta = get<1>(posa);
                int remaining_lena = get<4>(posa);
                int expected_lefta = (get<2>(posa)&eLeftComplete) ? lefta : lefta-remaining_lena; // one of the ends is always complete
                if(!ClipRight(posa, lefta, righta) && lefta <= righta) {
                    auto& hitb = *hitbp;
                    const auto& positionsb = hitb.second;
                    for(const THitPosition& posb : positionsb) {
                        int leftb = get<0>(posb);
                        int rightb = get<1>(posb);
                        int remaining_lenb = get<4>(posb);
                        int expected_rightb = (get<2>(posb)&eRightComplete) ? rightb : rightb+remaining_lenb;
                        ClipRight(posb, leftb, rightb);
                        if(leftb > rightb)
                            continue;
                        if(rightb > righta && (leftb <= righta || (get<2>(posb)&eLeftComplete))) {
                            int insert = expected_rightb-expected_lefta+1;
                            hits_for_variant.emplace_back(lefta, rightb-lefta, expected_rightb-lefta, insert);
                            paired_lengths.push_back(insert);
                        }
                    }                
                }
            }
        }

        sort(hits_for_variant.begin(), hits_for_variant.end(), [](const THitExtend& a, const THitExtend& b) { return get<0>(a) < get<0>(b); });
        sort(paired_lengths.begin(), paired_lengths.end()); 

        if(paired_lengths.empty())
            return 0;
        else
            return paired_lengths[0.75*paired_lengths.size()];
    }

    int SelectPairedHitsRtoL(const string& variant, const vector<int>& ambig_positions, const TReadHits& hits, vector<THitExtend>& hits_for_variant) {
        auto ClipLeft = [&] (const THitPosition& pos, int& left, int right) {
            auto ambig = get<3>(pos).begin();
            bool clip = false;
            if(ambig != get<3>(pos).end()) {
                for(int p : ambig_positions) {
                    if(p < left)
                        continue;
                    if(p > right)
                        break;
                    if(*(ambig++) != variant[p]) {
                        left = p+1;
                        clip = true;
                    }
                }
            }
            return clip;
        };

        hits_for_variant.clear();
        hits_for_variant.reserve(hits.size());
        vector<int> paired_lengths;
        unordered_set<TKeyForHits, SHashForHits> included_mates;
        for(auto hitp_loop = hits.begin(); hitp_loop != hits.end(); ++hitp_loop) {
            auto hitap = hitp_loop;
            CReadHolder::string_iterator isa = hitap->first.first;
            int strand = hitap->first.second;
            if(!included_mates.emplace(isa, strand).second || !isa.HasMate())
                continue;

            CReadHolder::string_iterator isb = isa.GetMate();
            auto hitbp = hits.find(make_pair(isb,-strand));
            if(hitbp == hits.end())
                continue;

            included_mates.emplace(isb, -strand);
            if(strand > 0) {
                swap(hitap, hitbp);
                strand = -strand;
            }

            const auto& positionsa = hitap->second;
            for(const THitPosition& posa : positionsa) {
                int righta = get<1>(posa);
                if(righta < (int)variant.size()-1 && !(get<2>(posa)&eRightComplete))
                    continue;

                int lefta = get<0>(posa);
                int remaining_lena = get<4>(posa);
                int expected_righta = (get<2>(posa)&eRightComplete) ? righta : righta+remaining_lena; // one of the ends is always complete
                if(!ClipLeft(posa, lefta, righta) && lefta <= righta) {
                    auto& hitb = *hitbp;
                    const auto& positionsb = hitb.second;
                    for(const THitPosition& posb : positionsb) {
                        int leftb = get<0>(posb);
                        int rightb = get<1>(posb);
                        int remaining_lenb = get<4>(posb);
                        int expected_leftb = (get<2>(posb)&eLeftComplete) ? leftb : leftb-remaining_lenb; // one of the ends is always complete
                        ClipLeft(posb, leftb, rightb);
                        if(leftb > rightb)
                            continue;
                        if(leftb < lefta && (rightb >= lefta || (get<2>(posb)&eRightComplete))) {
                            int insert = expected_righta-expected_leftb+1;
                            hits_for_variant.emplace_back(righta, righta-leftb, righta-expected_leftb, insert);
                            paired_lengths.push_back(insert);
                        }
                    }
                }
            }
        }    

        sort(hits_for_variant.begin(), hits_for_variant.end(), [](const THitExtend& a, const THitExtend& b) { return get<0>(a) > get<0>(b); });
        sort(paired_lengths.begin(), paired_lengths.end());

        if(paired_lengths.empty())
            return 0;
        else
            return paired_lengths[0.75*paired_lengths.size()];
    }

       
    void CleanContigJob(TContigs::iterator contigp, TContigs& rslt) {
        const string& contig = contigp->first;
        if((int)contig.size() < m_kmer_len) {
            rslt[contig] = contigp->second;
            return;
        }

        int contig_len = contig.size();
        vector<int> ambig_positions;
        for(int i = 0; i < contig_len; ++i) {
            if(FromAmbiguousIUPAC[contig[i]].size() > 1)
                ambig_positions.push_back(i);
        }

        TReadHits hits;
        int max_read_len = AlignReads(contig, ambig_positions, hits);

        list<deque<char>> selected_snps(1);
        string variant = contig;
        stack<pair<int,char>> remaining_choices;
        for(int ap : ambig_positions) {
            string ambig = FromAmbiguousIUPAC[contig[ap]];
            variant[ap] = ambig[0];
            selected_snps.front().push_back(variant[ap]);
            for(unsigned i = 1; i < ambig.size(); ++i)
                remaining_choices.emplace(ap, ambig[i]);
        }  


        TRange target_range = contigp->second.TargetRange();
        TRange clip_range{0, contig_len-1};

        while(true) { // loop through all ambiguous symbols combinations
            CReadHolder rh(false);
            rh.PushBack(variant);
            deque<CDBGraph::Node> nodes;
            for(CReadHolder::kmer_iterator itk = rh.kbegin(m_kmer_len); itk != rh.kend(); ++itk)  // gives kmers in reverse order! 
                nodes.push_front(m_graphp->GetNode(*itk));

            vector<double> labundance(contig_len, 0);
            vector<double> rabundance(contig_len, 0);
            for(int i = m_kmer_len-1; i < contig_len; ++i)
                labundance[i] = m_graphp->Abundance(nodes[i-m_kmer_len+1]);
            for(int i = 0; i <= contig_len-m_kmer_len; ++i)
                rabundance[i] = m_graphp->Abundance(nodes[i]);

            int critical_count = 5;  //10;
            int critical_extend = 0.25*max_read_len+0.5;
            double threshold = 3.e-4;

            vector<THitExtend> hits_for_variant;

            auto ScanLeftRight = [&](int insert, char skip_type) {
                if(insert == 0)
                   return make_pair(false, 0);

                size_t hit_num = 0;
                int max_extend = 0;
                int clip_counts = 0;
                int clip_counts_from_uniq = 0;
                double accum_read_len = 0;
                int accum_read_num = 0;
                int buff = (skip_type == 'A' ? max_read_len/2 : max_read_len);
                deque<tuple<int,double>> prev_reads(buff);
                for(int from = 0; from < contig_len; ++from) {
                    --max_extend;

                    accum_read_num -= get<0>(prev_reads.front());
                    accum_read_len -= get<1>(prev_reads.front());
                    prev_reads.pop_front();
                    prev_reads.emplace_back(0, 0.);

                    double abundance_from = labundance[from];
                    if(abundance_from == 0)
                        abundance_from = rabundance[from];

                    for( ; hit_num < hits_for_variant.size() && get<0>(hits_for_variant[hit_num]) == from; ++hit_num) {
                        int hit_insert = get<3>(hits_for_variant[hit_num]);
                        if(hit_insert > insert)
                            continue;
                        int extend = get<1>(hits_for_variant[hit_num]);
                        if(extend > max_extend) {
                            max_extend = extend;
                            clip_counts = 0;
                            clip_counts_from_uniq = 0;
                        }
                        if(extend == max_extend) {
                            if(get<2>(hits_for_variant[hit_num]) > extend)
                                ++clip_counts;                       
                            if(abundance_from <= 1.5*m_average_count && get<2>(hits_for_variant[hit_num]) >= extend+critical_extend)
                                ++clip_counts_from_uniq;
                        }

                        ++accum_read_num;
                        accum_read_len += hit_insert;
                        ++get<0>(prev_reads.back());
                        get<1>(prev_reads.back()) += hit_insert;
                    }
                    if(max_extend < 1)
                        continue;
                    int to = from+max_extend;
                    if(to < m_kmer_len-1)
                        continue;
                    if(to >= contig_len-1)
                        break;

                    auto Successors = [&]() {
                        TKmer kmer(variant.begin()+to-m_kmer_len+1, variant.begin()+to+1);
                        CDBGraph::Node node = m_graphp->GetNode(kmer); 
                        vector<CDBGraph::Successor> successors = m_graphdiggerp->GetReversibleNodeSuccessors(node);
                        return successors;
                    };

                    if(clip_counts_from_uniq >= critical_count && Successors().size() > 1) 
                        return make_pair(true, to);

                    if(accum_read_num > 10 &&  max_read_len >= m_kmer_len && to >= insert-1) {  
                        double expected_count = m_average_count;
                        if(skip_type == 'C')
                            expected_count /= 2;                                                                         // there are 2x less pairs than reads
                        double density = min(double(accum_read_num)/buff, expected_count/(max_read_len-m_kmer_len+1));   // number of reads/pairs per genome base (coverage/rlen)
                        double rlen = accum_read_len/accum_read_num;                                                     // local average read/pairs lengts
                        double margin = -log(threshold)/density;                                                         // read length distribution is delta func around rlen                  
                        if(max_extend < rlen-margin && Successors().size() > 1)                                                                 
                            return make_pair(true, to);
                    } 

                    if(max_extend > 0 && clip_counts > 5 && to >= insert-1) {
                        auto successors = Successors();
                        if(successors.empty()) {
                            return make_pair(true, to);
                        } else if(successors.size() > 1) {
                            double fraction = 0;
                            char next_symbol = variant[to+1];
                            auto rslt = find_if(successors.begin(), successors.end(), [=](const CDBGraph::Successor& suc){return suc.m_nt == next_symbol;});
                            if(rslt != successors.end()) {
                                fraction = m_graphp->Abundance(rslt->m_node);
                                double total = 0;
                                for(auto& suc: successors)
                                    total += m_graphp->Abundance(suc.m_node);
                                fraction /= total;                            
                            }
 
                            if(fraction == 0 || pow(1.-fraction, clip_counts) < threshold)
                                return make_pair(true, to);                        
                        }
                    }
                }

                return make_pair(false, 0);
            };            

            auto ScanRightLeft = [&](int insert, char skip_type) {
                if(insert == 0)
                    return make_pair(false, 0);

                size_t hit_num = 0;
                int max_extend = 0;
                int clip_counts = 0;
                int clip_counts_from_uniq = 0;
                double accum_read_len = 0;
                int accum_read_num = 0;
                int buff = (skip_type == 'B' ? max_read_len/2 : max_read_len);
                deque<tuple<int,double>> prev_reads(buff);
                for(int from = contig_len-1; from >=0; --from) {
                    --max_extend;

                    accum_read_num -= get<0>(prev_reads.front());
                    accum_read_len -= get<1>(prev_reads.front());
                    prev_reads.pop_front();
                    prev_reads.emplace_back(0, 0.);

                    double abundance_from = rabundance[from];
                    if(abundance_from == 0)
                        abundance_from = labundance[from];

                    for( ; hit_num < hits_for_variant.size() && get<0>(hits_for_variant[hit_num]) == from; ++hit_num) {
                        int hit_insert = get<3>(hits_for_variant[hit_num]);
                        if(hit_insert > insert)
                            continue;
                        int extend = get<1>(hits_for_variant[hit_num]);
                        if(extend > max_extend) {
                            max_extend = extend;
                            clip_counts = 0;
                            clip_counts_from_uniq = 0;
                        }
                        if(extend == max_extend) {
                            if(get<2>(hits_for_variant[hit_num]) > extend)
                                ++clip_counts;                       
                            if(abundance_from <= 1.5*m_average_count && get<2>(hits_for_variant[hit_num]) >= extend+critical_extend)
                                ++clip_counts_from_uniq;
                        }

                        ++accum_read_num;
                        accum_read_len += hit_insert;
                        ++get<0>(prev_reads.back());
                        get<1>(prev_reads.back()) += hit_insert;
                    }
                    if(max_extend < 1)
                        continue;
                    int to = from-max_extend;
                    if(to > contig_len-m_kmer_len)
                        continue;
                    if(to <= 0)
                        break;

                    auto Successors = [&]() {
                        TKmer kmer(variant.begin()+to, variant.begin()+to+m_kmer_len);
                        CDBGraph::Node node = m_graphp->ReverseComplement(m_graphp->GetNode(kmer)); 
                        vector<CDBGraph::Successor> successors = m_graphdiggerp->GetReversibleNodeSuccessors(node);
                        return successors;
                    };

                    if(clip_counts_from_uniq >= critical_count && Successors().size() > 1) 
                        return make_pair(true, to);

                    if(accum_read_num > 10 &&  max_read_len >= m_kmer_len && to <= contig_len-insert) {
                        double expected_count = m_average_count;
                        if(skip_type == 'D')
                            expected_count /= 2;                                                                         // there are 2x less pairs than reads
                        double density = min(double(accum_read_num)/buff, expected_count/(max_read_len-m_kmer_len+1));   // number of reads/pairs per genome base (coverage/rlen)
                        double rlen = accum_read_len/accum_read_num;                                                     // local average read/pairs lengts
                        double margin = -log(threshold)/density;                                                         // read length distribution is delta func around rlen
                        if(max_extend < rlen-margin && Successors().size() > 1)                                                               
                            return make_pair(true, to);                        
                    }

                    if(max_extend > 0 && clip_counts > 5 && to <= contig_len-insert) {
                        auto successors = Successors();
                        if(successors.empty()) {
                            return make_pair(true, to); 
                        } else if(successors.size() > 1) {
                            double fraction = 0;
                            char next_symbol = Complement(variant[to-1]);
                            auto rslt = find_if(successors.begin(), successors.end(), [=](const CDBGraph::Successor& suc){return suc.m_nt == next_symbol;});
                            if(rslt != successors.end()) {
                                fraction = m_graphp->Abundance(rslt->m_node);
                                double total = 0;
                                for(auto& suc: successors)
                                    total += m_graphp->Abundance(suc.m_node);
                                fraction /= total;                            
                            }

                            if(fraction == 0 || pow(1.-fraction, clip_counts) < threshold)
                                return make_pair(true, to);                                                                              
                        }                                                
                    }  
                }

                return make_pair(false, 0);
            };


            char skip_type;
            pair<bool, int> skip;
            {
               skip_type = 'A';
               SelectHitsLtoR(variant, ambig_positions, hits, hits_for_variant);            
               skip = ScanLeftRight(max_read_len, skip_type);            
            }
            if(skip.first) {
                if(skip.second < target_range.first) {
                    skip.first = false;
                    clip_range.first = max(clip_range.first, skip.second+1);
                } else if(skip.second > target_range.second) {
                    skip.first = false;
                    clip_range.second = min(clip_range.second, skip.second);
                }
            }
            
            if(!skip.first) {
                skip_type = 'B';
                SelectHitsRtoL(variant, ambig_positions, hits, hits_for_variant);
                skip = ScanRightLeft(max_read_len, skip_type);
            }
            
            if(skip.first) {
                if(skip.second < target_range.first) {
                    skip.first = false;
                    clip_range.first = max(clip_range.first, skip.second);
                } else if(skip.second > target_range.second) {
                    skip.first = false;
                    clip_range.second = min(clip_range.second, skip.second-1);
                }
            }            
            
            if(!skip.first) {
                skip_type = 'C';
                int insert = SelectPairedHitsLtoR(variant, ambig_positions, hits, hits_for_variant);
                skip = ScanLeftRight(insert, skip_type);            
            }  
                     
            if(skip.first) {
                if(skip.second < target_range.first) {
                    skip.first = false;
                    clip_range.first = max(clip_range.first, skip.second+1);
                } else if(skip.second > target_range.second) {
                    skip.first = false;
                    clip_range.second = min(clip_range.second, skip.second);
                }
            }
               
            if(!skip.first) {
                skip_type = 'D';
                int insert = SelectPairedHitsRtoL(variant, ambig_positions, hits, hits_for_variant);
                skip = ScanRightLeft(insert, skip_type);
            } 
                       
            if(skip.first) {
                if(skip.second < target_range.first) {
                    skip.first = false;
                    clip_range.first = max(clip_range.first, skip.second);
                } else if(skip.second > target_range.second) {
                    skip.first = false;
                    clip_range.second = min(clip_range.second, skip.second-1);
                }
            }
                
            if(skip.first) {
                lock_guard<mutex> guard(m_out_mutex);
                cerr << ">Skip variant " << skip_type << ": " << skip.second << " " << contigp->second.FastaDefLine() << "\n" << variant << endl;
                selected_snps.pop_back();
            }
            
            if(remaining_choices.empty()) {
                break;
            } else {
                int current_ap = remaining_choices.top().first;
                variant[current_ap] = remaining_choices.top().second;
                remaining_choices.pop();
                selected_snps.push_back(deque<char>());
                for(int ap : ambig_positions) {
                    if(ap > current_ap) {
                        string ambig = FromAmbiguousIUPAC[contig[ap]];
                        variant[ap] = ambig[0];
                        for(unsigned i = 1; i < ambig.size(); ++i)
                            remaining_choices.emplace(ap, ambig[i]);
                    }
                    selected_snps.back().push_back(variant[ap]);
                }
            }
        }

        if(selected_snps.empty())
            return;

        for(int n = 0; n < (int)ambig_positions.size(); ++n) {
            //rotate
            for(auto& snps : selected_snps) {
                snps.push_back(snps.front());
                snps.pop_front();
            }
            //sort
            selected_snps.sort();
            //combine snps with equal prefix
            for(auto it = selected_snps.begin(); it != selected_snps.end(); ++it) {
                string last_symbols;
                for(auto jt_loop = it; jt_loop != selected_snps.end(); ) {
                    auto jt = jt_loop++;
                    if(jt != it) {
                        if(equal(it->begin(), it->end()-1, jt->begin())) {
                            last_symbols.push_back(jt->back());
                            selected_snps.erase(jt);
                        } else {
                            break;
                        }
                    }
                }
                if(!last_symbols.empty()) {
                    last_symbols.push_back(it->back());
                    it->back() = ToAmbiguousIUPAC[AmbiguousString(last_symbols)];
                }
            }
        }

        SContigData cdata = contigp->second;
        if(clip_range.first > 0) {
            cdata.m_status &= ~eCircular;
            cdata.m_status &= ~eLeftDeadEnd;
            for(auto& inf : cdata.m_tinfo) {
                get<1>(inf).first -= clip_range.first;
                get<1>(inf).second -= clip_range.first;
            }
        }
        if(clip_range.second < contig_len-1) {
            cdata.m_status &= ~eCircular;
            cdata.m_status &= ~eRightDeadEnd;
        }
        variant = variant.substr(clip_range.first, clip_range.second-clip_range.first+1);
    
        for(auto& snps : selected_snps) {
            for(unsigned i = 0; i < ambig_positions.size(); ++i) {
                variant[ambig_positions[i]-clip_range.first] = snps[i];
            }
            rslt[variant] = cdata;            
        }
    }   

    void HashForReadsJob(vector<SAtomic<uint8_t>>& centinel, array<CReadHolder,2>& reads) {
        for(int p = 0; p < 2; ++p) {
            for(CReadHolder::string_iterator is = reads[p].sbegin(); is != reads[p].send(); ++is) {
                int rlen = is.ReadLen();
                if(rlen >= m_kmer_len) {
                    CReadHolder::kmer_iterator ik = is.KmersForRead(m_kmer_len);
                    // last read kmer
                    size_t index = (*ik).oahash()%m_reads_hash_table.size();
                    //grub bucket   
                    while(!centinel[index].Set(1, 0));
                    m_reads_hash_table[index].emplace_front(is);
                    //release bucket    
                    centinel[index] = 0;
                    // skip to first kmer (kmers are in inverse order) 
                    index = (*(ik += rlen-m_kmer_len)).oahash()%m_reads_hash_table.size();
                    //grub bucket   
                    while(!centinel[index].Set(1, 0));
                    m_reads_hash_table[index].emplace_front(is);
                    //release bucket    
                    centinel[index] = 0;
                }
            }
        }
    }

    void MarkDuplicatesJob(TUtilMap& util, vector<TContigShifts>& hash_contigs) {

        int dead_end_zone = 2*m_kmer_len;
        
        for(auto& util_info : util) {
            if(!get<0>(util_info.second).Set(1,0))
                continue;

            string contig = util_info.first;
            int contig_len = contig.size();
            TContigs::iterator contig_it = m_rslts.find(contig);
            SContigData& cdata = contig_it->second;
            int& status = cdata.m_status;
            TContigShifts& inherit = get<2>(util_info.second);
            for(int dir = 0; dir < 2 && !(status&eRedundant); ++dir) {
                deque<int>& hash_values = get<1>(util_info.second)[dir];
                bool left_dead_end = status&eLeftDeadEnd;
                bool right_dead_end = status&eRightDeadEnd;
                TRange target_range = cdata.TargetRange();
                if(dir > 0) {
                    ReverseComplementSeq(contig.begin(), contig.end());
                    swap(left_dead_end, right_dead_end);
                    swap(target_range.first, target_range.second);
                    target_range.first = contig_len-1-target_range.first;
                    target_range.second = contig_len-1-target_range.second;
                }

                array<int,2> match_points = {0, contig_len-m_kmer_for_duplicates};
                if(left_dead_end)
                    match_points[0] += m_kmer_for_duplicates;
                if(right_dead_end)
                    match_points[1] -= m_kmer_for_duplicates;
                for(int match_point : match_points) {
                    if(status&eRedundant)
                        break;
                    int index = hash_values[match_point];
                    for(auto& match : hash_contigs[index]) {
                        const string& other_contig = match.first->first;
                        int other_contig_len = other_contig.size();
                        int start = match.second-match_point;

                        if(contig_it == match.first)                                   // same contig
                            continue;
                        if(start < 0 || start+contig_len > other_contig_len)           // won't fit
                            continue;
                        if(other_contig_len == contig_len && other_contig < contig_it->first) // only one of two equal size (compare to not reversed)
                            continue;

                        SContigData& other_cdata = match.first->second;
                        TRange other_target_range = other_cdata.TargetRange();
                        deque<int>& other_hash_values = get<1>(util[other_contig])[0];
                        int diff = 0;
                        for(int pos = match_points[0]; pos <= match_points[1]; pos += m_kmer_for_duplicates) {
                            if(hash_values[pos] != other_hash_values[start+pos])
                                ++diff;
                        }
                        if(diff <= 1) {
                            bool target_included = other_target_range.first <= target_range.first+start && other_target_range.second >= target_range.second+start;
                            int left_zone = 0;
                            int right_zone = 0;
                            if(left_dead_end && target_included && start > 0 && contig_len >= dead_end_zone) {
                                left_zone = dead_end_zone;
                                int zone_diff = 0;
                                for(int i = 0; i < dead_end_zone; ++i) {
                                    if(!MatchWithAmbiguousDNA(contig[i], other_contig[i+start]))
                                        ++zone_diff;
                                }
                                if(zone_diff > 5)
                                    continue;
                            } 
                            if(right_dead_end && target_included && start+contig_len < other_contig_len && contig_len >= dead_end_zone) {
                                right_zone = dead_end_zone;
                                int zone_diff = 0;
                                for(int i = contig_len-dead_end_zone; i < contig_len; ++i) {
                                    if(!MatchWithAmbiguousDNA(contig[i], other_contig[i+start]))
                                        ++zone_diff;
                                }
                                if(zone_diff > 5)
                                    continue;
                            } 

                            if(left_zone+right_zone < contig_len && !equal(contig.begin()+left_zone, contig.end()-right_zone, other_contig.begin()+start+left_zone, MatchWithAmbiguousDNA))
                                continue;                            

                            int shift = (dir == 0 ? start : -(start+contig_len-1));  // diff between contig starts; negative for minus strand
                            inherit.emplace_front(match.first, shift);
                            status |= eRedundant;
                        }
                    }
                }
            }
        }
    }

    void KmerHashValuesJob(TUtilMap& util, int hash_table_size) {

        for(auto& util_info : util) {
            if(!get<0>(util_info.second).Set(1,0))
                continue;

            string contig = util_info.first;
            for(char& c : contig)
                c = FromAmbiguousIUPAC[c].front();
            deque<int>& hash_valuesp = get<1>(util_info.second)[0];
            CReadHolder rhp(false);
            rhp.PushBack(contig);
            for(CReadHolder::kmer_iterator itk = rhp.kbegin(m_kmer_for_duplicates); itk != rhp.kend(); ++itk)  { // gives kmers in reverse order! 
                TKmer k = *itk;
                hash_valuesp.push_front(k.oahash()%hash_table_size);
            }
        
            ReverseComplementSeq(contig.begin(), contig.end());
            deque<int>& hash_valuesm = get<1>(util_info.second)[1];
            CReadHolder rhm(false);
            rhm.PushBack(contig);
            for(CReadHolder::kmer_iterator itk = rhm.kbegin(m_kmer_for_duplicates); itk != rhm.kend(); ++itk)  { // gives kmers in reverse order!     
                TKmer k = *itk;
                hash_valuesm.push_front(k.oahash()%hash_table_size);
            }
        }
    }

    pair<string, bool> ExtendToFirstFork(CDBGraph::Node node, unordered_set<CDBGraph::Node>& used_nodes, int direction) const {
        string s;
        while(true) {
            vector<CDBGraph::Successor> successors = m_graphdiggerp->GetReversibleNodeSuccessors(node);
            //            vector<CDBGraph::Successor> successors = m_graphp->GetNodeSuccessors(node);
            //            m_graphdiggerp->FilterNeighbors(successors);
            if(successors.empty())
                return make_pair(s, true);                            
            if(successors.size() != 1)
                return make_pair(s, false);                        
            node = successors[0].m_node;
            if(used_nodes.insert(direction > 0 ? node : m_graphp->ReverseComplement(node)).second)
                s.push_back(successors[0].m_nt);
            else
                return make_pair(s, false);            
        }
    }
    void AssemblerJob(TContigs& rslts) {
        
        for(auto& item : m_targets) {
            if(!get<2>(item).Set(1,0))
                continue;

            string& target = get<0>(item);
            string& acc = get<1>(item);
            unordered_map<CDBGraph::Node, int> initial_target_kmers; // [node] position on target

            for(int p = 0; p <= (int)target.size()-8; ++p) {
                string seed = target.substr(p, 8); // last 8 symbols
                if(seed.find_first_not_of("ACGT") != string::npos)
                    continue;

                uint16_t word = 0;
                for(char c : seed) {
                    word = word << 2;
                    word += (find(bin2NT.begin(), bin2NT.end(), c) - bin2NT.begin());
                }

                array<int,2> pos = {p+8-m_kmer_len, p};
                array<unique_ptr<TKmer>,2> tkmerp;
                if(pos[0] >= 0) {
                    string tkmer = target.substr(pos[0], m_kmer_len);
                    if(tkmer.find_first_not_of("ACGT") == string::npos) {
                        tkmerp[0].reset(new TKmer(tkmer));                        
                    }
                }
                if(pos[1] <= (int)target.size()-m_kmer_len) {
                    string tkmer = target.substr(p, m_kmer_len);
                    if(tkmer.find_first_not_of("ACGT") == string::npos) {
                        ReverseComplementSeq(tkmer.begin(), tkmer.end());
                        tkmerp[1].reset(new TKmer(tkmer));
                    }
                }

                auto CountMatches = [](uint64_t a, uint64_t b) {
                    uint64_t w = ~(a^b);       // each match will produce 11                                           
                    w &= (w << 1);             // upper bit is 1 only if both bits are 1 (only for matches)
                    w &= 0xAAAAAAAAAAAAAAAAUL; // each match has 10; all mismatches 00
                    return _mm_popcnt_u64(w);  // count number set of bits == matches 
                };

                for(CDBGraph::Node node : m_kmer_hash[word]) {
                    int dir = node%2;
                    if(tkmerp[dir]) {                        
                        const uint64_t* targetp = tkmerp[dir]->getPointer();
                        const uint64_t* graphp = m_graphp->getPointer(node);
                        int prec = tkmerp[dir]->getSize()/64;// number of 8-byte blocks
                        int matches = -(prec*32-m_kmer_len); // empty positions will match  
                        for(int i = 0; i < prec; ++i)
                            matches += CountMatches(*(targetp+i), *(graphp+i));

                        if(matches > 0.8*m_kmer_len)
                            initial_target_kmers[node] = pos[dir];                                        
                    }
                }
            }            

            bool skip_target = false;
            TContigs rslts_for_target;
            int all_max_score = 0;

            for(int dir = 0; dir < 2 && !skip_target; ++dir) {
                unordered_map<CDBGraph::Node, int> target_kmers(initial_target_kmers);
                while(!target_kmers.empty()) {
                    unordered_map<CDBGraph::Node, int>::iterator first_kmerp;
                    typedef unordered_map<CDBGraph::Node, int>::value_type elem_t;
                    if(dir == 0)
                        first_kmerp = min_element(target_kmers.begin(), target_kmers.end(), [](const elem_t& a, const elem_t& b) { return a.second < b.second; });
                    else
                        first_kmerp = min_element(target_kmers.begin(), target_kmers.end(), [](const elem_t& a, const elem_t& b) { return a.second > b.second; }); 
                    
                    CDBGraph::Node initial_node = first_kmerp->first;
                    string kmer_seq = m_graphp->GetNodeSeq(initial_node);
                    int first_matching_kmer = first_kmerp->second;
                    int kmer_score = m_match*m_kmer_len;
                    for(int i = 0; i < m_kmer_len; ++i) {
                        if(kmer_seq[i] != target[first_matching_kmer+i])
                            kmer_score -= m_match+m_mismatch;
                    }
                    target_kmers.erase(first_kmerp);
                    list<tuple<int,int,string>> left_extensions;
                    if(first_matching_kmer > 0) {
                        string left_target = target.substr(0, first_matching_kmer);
                        ReverseComplementSeq(left_target.begin(), left_target.end());
                        CGudedPath extender(m_graphp->ReverseComplement(initial_node), left_target, *m_graphdiggerp, m_delta, m_gap_open, m_gap_extend, m_drop_off);
                        while(extender.ProcessNextEdge()) {                    
                            if(dir == 1) {
                                for(CDBGraph::Node n : extender.LastStepNodes())
                                    target_kmers.erase(m_graphp->ReverseComplement(n));
                            }
                            if(extender.PathEnd()) {
                                CGudedPath::SPathChunk rslt = extender.GetBestPart();
                                ReverseComplementSeq(rslt.m_seq.begin(), rslt.m_seq.end());
                                left_extensions.emplace_back(rslt.m_score, rslt.m_tlen, rslt.m_seq);
                                if(left_extensions.size() > m_max_variants) {
                                    skip_target = true;
                                    break;
                                }                                
                            }
                        }
                        if(skip_target)
                            break;
                    }
                    if(left_extensions.empty())
                        left_extensions.emplace_back(0, 0, "");
                    list<tuple<int,int,string>> right_extensions;
                    if(first_matching_kmer < (int)target.size()-m_kmer_len) {
                        string right_target =  target.substr(first_matching_kmer+m_kmer_len);  
                        CGudedPath extender(initial_node, target.substr(first_matching_kmer+m_kmer_len), *m_graphdiggerp, m_delta, m_gap_open, m_gap_extend, m_drop_off);
                        while(extender.ProcessNextEdge()) {
                            if(dir == 0) {
                                for(CDBGraph::Node n : extender.LastStepNodes())
                                    target_kmers.erase(n);
                            }
                            if(extender.PathEnd()) {
                                CGudedPath::SPathChunk rslt = extender.GetBestPart();
                                right_extensions.emplace_back(rslt.m_score, rslt.m_tlen, rslt.m_seq);
                                if(right_extensions.size() > m_max_variants) {
                                    skip_target = true;
                                    break;
                                }                                
                            }                      
                        }
                        if(skip_target)
                            break;
                    }
                    if(right_extensions.empty())
                        right_extensions.emplace_back(0, 0, "");
                    left_extensions.sort();
                    right_extensions.sort();
                    int max_score = get<0>(left_extensions.back())+kmer_score+get<0>(right_extensions.back());

                    for(auto& left : left_extensions) {
                        for(auto& right : right_extensions) {
                            int score = get<0>(left)+kmer_score+get<0>(right);
                            int tlen = get<1>(left)+m_kmer_len+get<1>(right);
                            if(score <= max_score-m_drop_off || tlen <= 0.5*target.size())
                                continue;
                            
                            string extension = get<2>(left)+kmer_seq+get<2>(right);                            
                            CCigar cigar = LclAlign(extension.c_str(), extension.size(), target.c_str(), target.size(), m_gap_open, m_gap_extend, m_delta.matrix);
                            if(cigar.SubjectRange().second-cigar.SubjectRange().first+1 <= 0.5*target.size())
                                continue;
                            score = cigar.Score(extension.c_str(), target.c_str(), m_gap_open, m_gap_extend, m_delta.matrix);
                            all_max_score = max(all_max_score, score);
                            extension = extension.substr(cigar.QueryRange().first, cigar.QueryRange().second-cigar.QueryRange().first+1);
                            
                            unordered_set<CDBGraph::Node> used_nodes;
                            {
                                string ext_first(extension);
                                for(char& c : ext_first)
                                    c = FromAmbiguousIUPAC[c].front();
                                CReadHolder rh(false);
                                rh.PushBack(ext_first);
                                for(CReadHolder::kmer_iterator itk = rh.kbegin(m_kmer_len); itk != rh.kend(); ++itk)
                                    used_nodes.insert(m_graphp->GetNode(*itk));
                            }                                

                            int left_pos = 0;
                            int right_pos = extension.size()-1;

                            string right_kmer = extension.substr(extension.size()-m_kmer_len);
                            for(char& c : right_kmer)
                                c = FromAmbiguousIUPAC[c].front();
                            CDBGraph::Node right_node = m_graphp->GetNode(right_kmer);
                            auto rslt = ExtendToFirstFork(right_node, used_nodes, 1);
                            extension += rslt.first;
                            bool right_dead_end = rslt.second;

                            string left_kmer = extension.substr(0, m_kmer_len);
                            for(char& c : left_kmer)
                                c = FromAmbiguousIUPAC[c].front();
                            CDBGraph::Node left_node = m_graphp->GetNode(left_kmer);
                            rslt = ExtendToFirstFork(m_graphp->ReverseComplement(left_node), used_nodes, -1);
                            bool left_dead_end = rslt.second;
                            string left_extra = rslt.first;
                            ReverseComplementSeq(left_extra.begin(), left_extra.end());
                            extension = left_extra+extension;

                            left_pos += left_extra.size();
                            right_pos += left_extra.size();

                            SContigData& cdata = rslts_for_target[extension];
                            if(left_dead_end)
                                cdata.m_status |= eLeftDeadEnd;
                            if(right_dead_end)
                                cdata.m_status |= eRightDeadEnd;
                            auto& lst = cdata.m_tinfo;
                            if(lst.empty())
                                lst.emplace_back(score, TRange(left_pos,right_pos), acc);
                            else if(score > get<0>(lst.front()))
                                lst.front() = make_tuple(score, TRange(left_pos,right_pos), acc); // only one info for each sequence
                        }                    
                    }
                }            
            }

            if(skip_target) {
                lock_guard<mutex> guard(m_out_mutex);
                cerr << "Too many variants for target: " << acc << endl;
            } else {
                for(auto& rslt : rslts_for_target) {
                    int score = get<0>(rslt.second.m_tinfo.front());
                    if(score <= all_max_score-m_drop_off)
                        continue;

                    SContigData& cdata = rslts[rslt.first];
                    cdata.m_status |= rslt.second.m_status;
                    auto& lst = cdata.m_tinfo;
                    AddInfo(lst, rslt.first.size(), rslt.second.m_tinfo); 
                }
            }
        }               
    }

    vector<forward_list<CDBGraph::Node>> m_kmer_hash;

    TSimpleHash m_reads_hash_table;

    unsigned m_max_variants = 10000;
    int m_kmer_for_duplicates = MaxPrec*32;

    list<tuple<string,string,SAtomic<uint8_t>>> m_targets; // seq, accsession, centinel
    TContigs m_rslts;
    unique_ptr<CDBGraphDigger> m_graphdiggerp;
    unique_ptr<CDBGraph> m_graphp;
    int m_kmer_len;
    int m_kmer_len_for_seeds;
    int m_match;
    int m_mismatch;
    int m_gap_open;
    int m_gap_extend;
    int m_drop_off;
    SMatrix m_delta;
    int m_ncores;
    list<array<CReadHolder,2>>& m_raw_reads;
    double  m_average_count;
    mutex m_out_mutex;

};

int main(int argc, const char* argv[]) {
    for(int n = 0; n < argc; ++n)
        cerr << argv[n] << " ";
    cerr << endl << endl;

    int ncores;
    double fraction;
    int min_count;
    int kmer;
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
    int drop_off;
    vector<string> sra_list;
    vector<string> fasta_list;
    vector<string> fastq_list;
    bool gzipped;
    int memory;

    ofstream contigs_out;
    ofstream profile_out;

    options_description general("General options");
    general.add_options()
        ("help,h", "Produce help message")
        ("memory", value<int>()->default_value(32), "Memory available (GB) [integer]")
        ("cores", value<int>()->default_value(0), "Number of cores to use (default all) [integer]")
        ("contigs_out", value<string>(), "Output file for contigs (stdout if not specified) [string]")
        ("profile_out", value<string>(), "Output file for coverage profile (optional) [string]");

    options_description input("Input/output options : at least one input providing reads for assembly must be specified");
    input.add_options()
        ("fasta", value<vector<string>>(), "Input fasta file(s) (could be used multiple times for different runs) [string]")
        ("fastq", value<vector<string>>(), "Input fastq file(s) (could be used multiple times for different runs) [string]")
        ("sra_run", value<vector<string>>(), "Input sra run accession (could be used multiple times for different runs) [string]")
        ("gz", "Input files are gzipped [flag]");

    options_description assembly("Assembly options");
    assembly.add_options()
        ("kmer", value<int>()->default_value(21), "Minimal kmer length for assembly [integer]")
        ("min_count", value<int>()->default_value(2), "Minimal count for kmers retained for comparing alternate choices [integer]")
        ("fraction", value<double>()->default_value(0.05, "0.05"), "Maximum noise to signal ratio acceptable for extension [float [0,1)]")
        ("match", value<int>()->default_value(1), "Bonus for match")
        ("mismatch", value<int>()->default_value(2), "Penalty for mismatch")
        ("gap_open", value<int>()->default_value(5), "Penalty for gap opening")
        ("gap_extend", value<int>()->default_value(2), "Penalty for gap extension")
        ("drop_off", value<int>()->default_value(25), "Maximal decrease of score");

    options_description all("");
    all.add(general).add(input).add(assembly); 

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

        fraction = argm["fraction"].as<double>();
        if(fraction >= 1.) {
            cerr << "Value of --fraction must be < 1 (more than 0.25 is not recommended)" << endl;
            exit(1);
        }
        if(fraction < 0.) {
            cerr << "Value of --fraction must be >= 0" << endl;
            exit(1);
        }
        min_count = argm["min_count"].as<int>();
        if(min_count <= 0) {
            cerr << "Value of --min_count must be > 0" << endl;
            exit(1);
        }
        kmer = argm["kmer"].as<int>();
        if(kmer%2 ==0) {
            cerr << "Kmer must be an odd number" << endl;
            return 1;
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

        if(argm.count("profile_out")) {
            profile_out.open(argm["profile_out"].as<string>());
            if(!profile_out.is_open()) {
                cerr << "Can't open file " << argm["profile_out"].as<string>() << endl;
                exit(1);
            }
        }

        match = argm["match"].as<int>();
        mismatch = argm["mismatch"].as<int>();
        gap_open = argm["gap_open"].as<int>();
        gap_extend = argm["gap_extend"].as<int>();
        drop_off = argm["drop_off"].as<int>();

        CReadsGetter readsgetter(sra_list, fasta_list, fastq_list, ncores, true, gzipped);
        CGuidedAssembler gassembler(kmer, min_count, fraction, memory, match, mismatch, gap_open, gap_extend, drop_off, ncores, readsgetter.Reads());

        CGuidedAssembler::TContigs rslts = gassembler.Contigs();

        CStopWatch timer;
        timer.Restart();
        int count = 0;
        ostream& out = contigs_out.is_open() ? contigs_out : cout;
        for(auto& rslt : rslts) {
            const string& contig = rslt.first;
            string acc = "Contig"+to_string(++count)+"_kmer"+to_string(kmer)+"_improved";
            out << ">" << acc << " " << rslt.second.FastaDefLine() << "\n" << contig << endl;

            if(profile_out.is_open()) {
                auto AbundanceForAllChoices = [&](const string& kmer) {
                    auto& graph = gassembler.Graph(); 
                    map<int,string> ambig_positions;
                    for(int i = 0; i < (int)kmer.size(); ++i) {
                        string ambig = FromAmbiguousIUPAC[kmer[i]];
                        if(ambig.size() > 1)
                            ambig_positions[i] = ambig;
                    }

                    string variant = kmer;
                    stack<pair<int,char>> remaining_choices;
                    for(auto& pos_sympols : ambig_positions) {
                        int ap = pos_sympols.first;
                        string& ambig = pos_sympols.second;
                        variant[ap] = ambig[0];
                        for(unsigned i = 1; i < ambig.size(); ++i)
                            remaining_choices.emplace(ap, ambig[i]);
                    }

                    int abundance = graph.Abundance(graph.GetNode(variant));
                    while(!remaining_choices.empty()) {
                        int current_ap = remaining_choices.top().first;
                        variant[current_ap] = remaining_choices.top().second;
                        remaining_choices.pop();
                        for(auto& pos_sympols : ambig_positions) {
                            int ap = pos_sympols.first;
                            if(ap <= current_ap)
                                continue;
                            string& ambig = pos_sympols.second;
                            variant[ap] = ambig[0];
                            for(unsigned i = 1; i < ambig.size(); ++i)
                                remaining_choices.emplace(ap, ambig[i]);
                        }
                        abundance += graph.Abundance(graph.GetNode(variant));
                    }
                    
                    return abundance;
                };
               
                int contig_len = contig.size();
                for(int p = 0; p < contig_len; ++p) {
                    profile_out << acc << '\t' << p+1 << '\t';
                    string symb = FromAmbiguousIUPAC[contig[p]];
                    for(unsigned i = 0; i < symb.size(); ++i) {
                        if(i > 0)
                            profile_out << '/';
                        profile_out << symb[i];
                    }
                    profile_out << '\t';
                    for(unsigned i = 0; i < symb.size(); ++i) {
                        int labundance = 0;
                        int rabundance = 0;
                        if(p >= kmer-1) {
                            string kmer_seq = contig.substr(p-kmer+1, kmer);
                            kmer_seq.back() = symb[i];
                            labundance = AbundanceForAllChoices(kmer_seq);
                        }
                        if(p <= contig_len-kmer) {
                            string kmer_seq = contig.substr(p, kmer);
                            kmer_seq.front() = symb[i];
                            rabundance = AbundanceForAllChoices(kmer_seq);
                        }                        
                        if(i > 0)
                            profile_out << '/';
                        profile_out << max(labundance,rabundance);
                    }
                    profile_out << endl;
                }
            }
        }
        cerr << "Output in " << timer.Elapsed();
    
    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        exit(1);
    }

    cerr << "DONE" << endl;

    return 0;
}
