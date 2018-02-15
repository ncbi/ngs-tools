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

#include <math.h>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <atomic>
#include <future>
#include <thread>
#include <string>
#include <deque>
#include <forward_list>

#include "glb_align.hpp"
#include "concurrenthash.hpp"
#include "ngs_includes.hpp"

#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/seek.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/copy.hpp>

#include <assert.h> 

using namespace boost::program_options;
using namespace boost::iostreams;
using namespace DeBruijn;
using namespace std;
using namespace ngs;

typedef vector<char> TSeqV;
typedef pair<string, TSeqV> TSeq;
typedef list<TSeq> TSeqHolder;

typedef pair<string, SAtomic<unsigned>> TCnt;
typedef CForwardList<TCnt> TCounter;
typedef vector<TCounter> TCounterV;
typedef tuple<string, TSeqV, TCounterV> TGenomeSeq;
typedef list<TGenomeSeq> TGenomeHolder;
typedef pair<int, TGenomeSeq*> TKmerPos;
typedef CForwardList<TKmerPos> TKmerLocations;

class CKmerLocationsHash : public CKmerHashMap<TKmerLocations, 8> {
public:
    CKmerLocationsHash(int kmer_len, size_t size) : CKmerHashMap(kmer_len, size) {}
    void Insert(const TKmer& kmer, int pos, TGenomeSeq* contigp) {
        size_t index = kmer.oahash()%m_table_size;
        FindOrInsertInBucket(kmer, index)->PushFront(TKmerPos(pos,contigp));
    }
};

struct SAlignParams{
    SAlignParams(variables_map& argm) {
        m_match = argm["match"].as<int>();
        m_mismatch = argm["mismatch"].as<int>();
        m_gapopen = argm["gapopen"].as<int>();
        m_gapextend = argm["gapextend"].as<int>();
        m_word_size = argm["word"].as<int>();
        m_min_compart = (int)(m_word_size*argm["min_compart"].as<double>()+0.5);
        m_entropy_level = argm["entropy"].as<double>();
        m_coverage = argm["coverage"].as<double>();
        if(argm.count("sra_run"))
            m_sra_run = argm["sra_run"].as<string>();
    }

    int m_word_size;
    int m_min_compart;
    double m_entropy_level;
    int m_match;
    int m_mismatch;
    int m_gapopen;
    int m_gapextend;
    double m_coverage;    
    string m_sra_run;
};

struct SAlignRslt {
    vector<char> m_out;
    size_t m_hits = 0;
    size_t m_alignments = 0;
    size_t m_attempts = 0;
};

struct SHit {
    SHit(int lng, int shrt, int l) : m_long(lng), m_short(shrt), m_l(l) {};
    int m_long;
    int m_short;
    int m_l;
    int m_score = 0;
};

struct SHitP : public SHit {
    SHitP(int lng, int shrt, int l) : SHit(lng, shrt, l) {};
    bool operator < (const SHitP& h) const { return m_long <  h.m_long; }
    void MarkUsed() {
        for(auto left = m_left; left != nullptr; left = left->m_left)
            left->m_score = -1;
    }

    SHitP* m_left = nullptr;
};

class CCompartment : public list<SHit> {
public:
    CCompartment(const SHitP& hit) {
        push_front(hit);
        for(auto left = hit.m_left; left != nullptr; left = left->m_left)
            push_front(*left);        
    }
    void ConnectHits(const TSeqV& read_seq, const TSeqV& contig_seq, int match, int mismatch, int gapcost) {

        auto MisMatches = [](TSeqV::const_iterator i1, TSeqV::const_iterator i2, int len) {
            int mismatches = 0;
            for( ; len > 0; --len, ++i1, ++i2) {
                if(*i1 != *i2)
                    ++mismatches;
            }
            return mismatches;
        };

        auto next = begin();
        next->m_score = next->m_l*match;
        for(auto prev = next++; next != end(); ) {
            if(next->m_long-prev->m_long == next->m_short-prev->m_short) { // same diag
                int gap_len = next->m_long-(prev->m_long+prev->m_l);
                int dist = MisMatches(read_seq.begin()+prev->m_short+prev->m_l, contig_seq.begin()+prev->m_long+prev->m_l, gap_len);
                int gap_score = gap_len*match-dist*mismatch;
                if(dist*mismatch < gapcost && gap_score+min(prev->m_score, next->m_score) > 0) {
                    prev->m_score += gap_score+next->m_score;
                    prev->m_l += gap_len+next->m_l;
                    next = erase(next);
                } else {
                    prev = next++;
                }
            } else {
                prev = next++;
            }
        }
    }
    bool SimpleCompartment() const {
        auto it = begin();
        return it != end() && ++it == end(); // only one hit
    }            
};

class CHitsHolder : public forward_list<SHitP> {
public:
    bool ExtendHit(int lng, int shrt, int l) {
        for(auto& hit :  *this) {
            if(lng-hit.m_long == shrt-hit.m_short) {  // same diag
                int overlap = hit.m_long+hit.m_l-lng;
                if(overlap >= 0) {
                    hit.m_l += l-overlap; // because we scan read forward this is garanteed to be > 0
                    return true;                
                }
            }
        }

        return false;
    }
    void AddHit(int lng, int shrt, int l) {
        if(!ExtendHit(lng, shrt, l)) 
           emplace_front(lng, shrt, l); 
    }
    void ScoreCompartments(int match, int mismatch, int gapopen, int gapextend) {
        sort();
        for(auto i = begin(); i != end(); ++i) {
            i->m_score = match*i->m_l;
            int ilong = i->m_long+i->m_l;
            int ishort = i->m_short+i->m_l;
            for(auto j = begin(); j != i; ++j) {
                int jlong = j->m_long+j->m_l;
                int jshort = j->m_short+j->m_l;
                if(j->m_long < i->m_long && j->m_short < i->m_short && jlong < ilong && jshort < ishort) {
                    int score = i->m_score+j->m_score;                    
                    int max_overlap = max(jlong - i->m_long, jshort - i->m_short);
                    if(max_overlap > 0)                                           // seeds do overlap
                        score -= match*max_overlap;                    
                    int diag_diff = abs(ilong-jlong-ishort+jshort);
                    if(diag_diff > 0)                                            // not same diag
                        score -= gapopen+gapextend*diag_diff;                    
                    if(score > i->m_score) {
                        i->m_score = score;
                        i->m_left = &(*j);
                    }
                }
            }
        }
    }
};
       

class CHits {
public:
    CHits() : m_best(m_hfr.end()) {};
    CHits(CHits&& other) {
        m_hfr.splice_after(m_hfr.before_begin(), other.m_hfr);
        m_best = other.m_best;
    }

    struct SHit;
    typedef forward_list<SHit> THitsForRead;

    struct SHit {
        SHit(int64_t lng, int shrt, int l) : m_longp(lng), m_shortp(shrt), m_l(l), m_score(0) {};
        bool operator < (const SHit& h) const { return m_longp <  h.m_longp; }

        THitsForRead::iterator m_left;
        int64_t m_longp;
        int m_shortp;
        int m_l;
        int m_score;
    };

    THitsForRead& HitList() { return m_hfr; }

    bool ExtendHit(int64_t lng, int shrt, int l) {
        for(auto& hit :  m_hfr) {
            if(lng-hit.m_longp == shrt-hit.m_shortp) {  // same diag
                int overlap = hit.m_longp+hit.m_l-lng;
                if(overlap >= 0) {
                    hit.m_l += l-overlap; // because we scan read forward this is garanteed to be > 0
                    return true;                
                }
            }
        }

        return false;
    }

    void AddHit(int lng, int shrt, int l) {
        if(!ExtendHit(lng, shrt, l)) 
           m_hfr.emplace_front(lng, shrt, l); 
    }

    void Sort() { m_hfr.sort(); }

    int CalculateCoverage() {
        for(auto ih = m_hfr.begin(); ih != m_hfr.end(); ++ih) {
            ih->m_score = ih->m_l;
            ih->m_left = m_hfr.end();
            int iqr = ih->m_longp+ih->m_l;
            int isr = ih->m_shortp+ih->m_l;
            for(THitsForRead::iterator jh = m_hfr.begin(); jh != ih; ++jh) {
                int jqr = jh->m_longp+jh->m_l;
                int jsr = jh->m_shortp+jh->m_l;            
                if(iqr > jqr && isr > jsr) {
                    int score = ih->m_l+jh->m_score;
                    int max_overlap = max(int(jqr-ih->m_longp),jsr-ih->m_shortp);
                    if(max_overlap > 0)                                           // seeds do overlap
                        score -= max_overlap;

                    if(score > ih->m_score) {
                        ih->m_score = score;
                        ih->m_left = jh;
                    }
                }
            }
        }

        m_best = m_hfr.begin();
        for(auto ih = m_hfr.begin(); ih != m_hfr.end(); ++ih) {
            if(ih->m_score > m_best->m_score)
                m_best = ih;
        }

        return m_best->m_score;
    }

    void FinalizeCompartment() {
        THitsForRead::iterator right = m_hfr.end();
        for(THitsForRead::iterator br = m_best ; br != m_hfr.end(); br = br->m_left) {
            m_hfr.erase_after(br,right);
            right = br;
        }
        m_hfr.erase_after(m_hfr.before_begin(),right);
    }

    void RemoveFirst() {
        int p = m_hfr.front().m_longp;
        while(!m_hfr.empty() && m_hfr.front().m_longp == p)
            m_hfr.pop_front(); 
    }

    THitsForRead::iterator FirstHit() { return m_hfr.begin(); }
    THitsForRead::iterator LastHit() { return m_best; }

    void ConnectHits(const TSeqV& read_seq, const TSeqV& contig_seq, int match, int mismatch, int gapcost) {

        auto MisMatches = [](TSeqV::const_iterator i1, TSeqV::const_iterator i2, int len) {
            int mismatches = 0;
            for( ; len > 0; --len, ++i1, ++i2) {
                if(*i1 != *i2)
                    ++mismatches;
            }
            return mismatches;
        };

        auto next = m_hfr.begin();
        next->m_score = next->m_l*match;
        for(auto prev = next++; next != m_hfr.end(); ) {
            next->m_score = next->m_l*match;
            if(next->m_longp-prev->m_longp == next->m_shortp-prev->m_shortp) { // same diag
                int gap_len = next->m_longp-(prev->m_longp+prev->m_l);
                int dist = MisMatches(read_seq.begin()+prev->m_shortp+prev->m_l, contig_seq.begin()+prev->m_longp+prev->m_l, gap_len);
                int gap_score = gap_len*match-dist*mismatch;
                if(dist*mismatch < gapcost && gap_score+min(prev->m_score, next->m_score) > 0) {
                    prev->m_score += gap_score+next->m_score;
                    prev->m_l += gap_len+next->m_l;
                    if(next == m_best)
                        m_best = prev;
                    next = m_hfr.erase_after(prev);
                } else {
                    prev = next++;
                }
            } else {
                prev = next++;
            }
        }
    }

    bool SimpleCompartment() const {
        auto it = m_hfr.begin();
        return it != m_hfr.end() && ++it == m_hfr.end(); // only one hit
    }

private:
    CHits(const CHits&);            // copy will invalidate m_best
    CHits& operator=(const CHits&); // copy will invalidate m_best
    THitsForRead m_hfr;
    THitsForRead::iterator m_best;
};

void GetReadsFromSRAJob(const string& sra_run, size_t from, size_t to, TSeqHolder& reads) {
    ReadCollection run = ncbi::NGS::openReadCollection (sra_run);
    ReadIterator it = run.getReadRange (from+1, to-from+1, Read::all);

    while(it.nextRead()) {
        int fragments = it.getNumFragments ();
        if(fragments == 2) { // paired read
            it.nextFragment();
            reads.push_back(TSeq());
            StringRef s1 = it.getFragmentBases();
            int read_length1 = s1.size();
            reads.back().second.insert(reads.back().second.end(), s1.data(), s1.data()+read_length1);
            StringRef id1 = it.getReadName();
            reads.back().first.insert(reads.back().first.end(), id1.data(), id1.data()+id1.size());
            reads.back().first += "_1";
                
            it.nextFragment();
            reads.push_back(TSeq());
            StringRef s2 = it.getFragmentBases();
            int read_length2 = s2.size();
            reads.back().second.insert(reads.back().second.end(), s2.data(), s2.data()+read_length2);
            StringRef id2 = it.getReadName();
            reads.back().first.insert(reads.back().first.end(), id2.data(), id2.data()+id2.size());
            reads.back().first += "_2";
        } else {             // unpaired read   
            while(it.nextFragment()) {
                reads.push_back(TSeq());
                StringRef s = it.getFragmentBases();
                int read_length = s.size();
                reads.back().second.insert(reads.back().second.end(), s.data(), s.data()+read_length);
                StringRef id = it.getReadName();
                reads.back().first.insert(reads.back().first.end(), id.data(), id.data()+id.size());
            } 
        }               
    }
}

template <typename Iterator>
Iterator FirstACGT(Iterator b, Iterator e) { return  find_if(b, e, [](char c){return c=='A' || c== 'C' || c=='G' || c=='T';}); }

template <typename Iterator>
Iterator FirstNotACGT(Iterator b, Iterator e) { return  find_if(b, e, [](char c){return c!='A' && c!= 'C' && c!='G' && c!='T';}); }

void IncrementCount(TCounterV& counts, const TCharAlign& align, int cstart) {
    for(int i = 0; i < (int)align.first.size(); ++cstart) {
        string s(1,align.first[i++]);
        for( ; i < (int)align.first.size() && align.second[i] == '-'; ++i)
            s.push_back(align.first[i]);                                        
        
        TCounter& lst = counts[cstart];                                
        auto existing_head = lst.Head();
        bool found = false;
        for(auto p = existing_head; p != nullptr; p = p->m_next) {
            if(p->m_data.first == s) {
                ++p->m_data.second.m_atomic;
                found = true;
            }
        }
        if(!found) {
            TCounter::SNode* nodep = new TCounter::SNode;
            nodep->m_data = make_pair(s, 1);
            nodep->m_next = existing_head;
            while(!found && !lst.TryPushFront(nodep)) {
                for(auto p = nodep->m_next; !found && p != existing_head; p = p->m_next) {
                    if(p->m_data.first == s) {
                        delete nodep;
                        ++p->m_data.second.m_atomic;
                        found = true;
                    }
                }
            }
        }
    }
}

void GetAlignsJob(const TSeqHolder& reads, CKmerLocationsHash& contig_hash, const SAlignParams& align_params, SAlignRslt& align_rslts) {
    int word_size = align_params.m_word_size;
    int min_compart = align_params.m_min_compart;
    int match = align_params.m_match;
    int mismatch = align_params.m_mismatch;
    int gapopen = align_params.m_gapopen;
    int gapextend = align_params.m_gapextend;
    int gapcost = gapopen+gapextend;
    double coverage = align_params.m_coverage;    
    string sra_run = align_params.m_sra_run;
    SMatrix delta(match, mismatch);

    back_insert_device<vector<char>> snk{align_rslts.m_out};
    filtering_ostream out;
    out.push(gzip_compressor());
    out.push(snk);

    for(auto& read : reads) {
        string read_acc = (sra_run.empty() ? "" : sra_run+".")+read.first;
        auto read_seq = read.second; // copy
        int read_length = read_seq.size();
        for(int strand = 1; strand >= -1; strand -= 2) {
            unordered_map<TGenomeSeq*, CHitsHolder> contig_hits;
            for(auto left = FirstACGT(read_seq.begin(), read_seq.end()); left != read_seq.end(); ) {
                auto right = FirstNotACGT(left+1, read_seq.end());
                unsigned len = right-left;
                if(len >= (unsigned)word_size) {
                    CReadHolder rh(false);
                    rh.PushBack(left, len);
                    unsigned read_pos = left-read_seq.begin();
                    CReadHolder::kmer_iterator ik = rh.kbegin(word_size);
                    unsigned knum = right-left-word_size+1;
                    for(ik += knum-1; read_pos <= right-read_seq.begin()-word_size; ++read_pos, ik += -1) { // forward iteration
                    
                        auto rslt = contig_hash.Find(*ik);
                        if(rslt != nullptr) {
                            for(auto& hit : *rslt) {
                                contig_hits[hit.second].AddHit(hit.first, read_pos, word_size);
                            }
                        }
                    }
                }
                left = FirstACGT(right, read_seq.end());
            }   
                        
            for(auto& ch : contig_hits) {
                auto& contig = ch.first;
                auto& contig_acc = get<0>(*contig);
                auto& contig_seq = get<1>(*contig);
                auto& counts = get<2>(*contig);
                int contig_length = contig_seq.size();

                list<CCompartment> compartments;
                {
                    auto& hits = ch.second;
                    hits.ScoreCompartments(match, mismatch, gapopen, gapextend);

                    forward_list<SHitP*> compartment_ptrs;
                    for(auto& hit : hits)
                        compartment_ptrs.push_front(&hit);
                    compartment_ptrs.sort([](SHitP* a, SHitP* b){ return a->m_score > b->m_score; });
                    for(auto& ptr : compartment_ptrs) {
                        if(ptr->m_score > 0)
                            ptr->MarkUsed();
                        if(ptr->m_score >= match*min_compart) {
                            compartments.emplace_back(*ptr);
                            compartments.back().ConnectHits(read_seq, contig_seq, match, mismatch, gapcost);
                            ++align_rslts.m_attempts;
                            align_rslts.m_hits += compartments.back().size();                
                        }
                    }
                }
                                
                set<tuple<int,int,int,int>> previous_aligns;
                for(auto& compartment : compartments) {
                    if(compartment.SimpleCompartment()) { // only one hit
                        auto hit = compartment.front();
                        int rstart = hit.m_short;
                        int rstop = hit.m_short+hit.m_l-1;
                        if(hit.m_score > match*max(rstart, read_length-rstop-1)) { // 'main' portion
                            int cstart = hit.m_long;
                            int cstop = hit.m_long+hit.m_l-1;
                        
                            CCigar lcigar;
                            if(rstart > 0) {
                                int max_gap = max(0,(match*rstart-gapopen)/max(match, gapextend));
                                int clen = min(rstart+max_gap, cstart);
                                lcigar = LclAlign(read_seq.data(), rstart, contig_seq.data()+cstart-clen, clen, gapopen, gapextend, false, true, delta.matrix);
                                int ext = lcigar.QueryRange().second-lcigar.QueryRange().first+1;
                                if(ext > 0) {
                                    rstart -= ext;
                                    cstart -= lcigar.SubjectRange().second-lcigar.SubjectRange().first+1;
                                }
                            }
                            
                            CCigar rcigar;
                            if(rstop < read_length-1) {
                                int max_gap = max(0,(match*(read_length-rstop-1)-gapopen)/max(match, gapextend));
                                int clen = min(read_length-rstop-1+max_gap, contig_length-cstop-1);
                                rcigar = LclAlign(read_seq.data()+rstop+1, read_length-rstop-1, contig_seq.data()+cstop+1, clen, gapopen, gapextend, true, false, delta.matrix);
                                int ext = rcigar.QueryRange().second-rcigar.QueryRange().first+1;
                                if(ext > 0) {
                                    rstop += ext;
                                    cstop += rcigar.SubjectRange().second-rcigar.SubjectRange().first+1;
                                }
                            }                            

                            if(rstop-rstart+1 >= coverage*read_length && previous_aligns.emplace(rstart,rstop,cstart,cstop).second) {
                                CCigar cigar(rstop, cstop);
                                cigar.PushFront(rcigar);
                                cigar.PushFront(CCigar::SElement(hit.m_l, 'M'));
                                cigar.PushFront(lcigar);
                                int dist = cigar.Distance(read_seq.data(), contig_seq.data());
                                out << read_acc << '\t' << (strand > 0 ? 0 : 16) << '\t' << contig_acc << '\t' << cstart+1 << '\t' << 255 << '\t' << cigar.CigarString(0, read_length) << "\t*\t0\t0\t";
                                for(char c : read_seq)
                                    out << c;                                        
                                out << "\t*\tNM:i:" << dist << endl; 
                                ++align_rslts.m_alignments;

                                if(!counts.empty())
                                    IncrementCount(counts, cigar.ToAlign(read_seq.data(), contig_seq.data()), cstart);
                            }

                            continue;
                        }
                    }
                                        
                    auto first_hit = compartment.front();
                    auto last_hit = compartment.back();
                    
                    int rstart = first_hit.m_short;
                    int rstop = last_hit.m_short+last_hit.m_l-1;
                    int cstart = first_hit.m_long;
                    int cstop = last_hit.m_long+last_hit.m_l-1;

                    CCigar cigar;
                    while(true) {
                        int flank_left = min(cstart, max(0,(match*rstart-gapopen)/max(match, gapextend))+rstart);
                        int flank_right = min(contig_length-cstop-1, max(0,(match*(read_length-rstop-1)-gapopen)/max(match, gapextend))+read_length-rstop-1);
                        cstart = cstart-flank_left;
                        cstop = cstop+flank_right;
                        int llen = cstop-cstart+1;
                        CCigar lcl_cigar = LclAlign(read_seq.data(), read_length, contig_seq.data()+cstart, llen, gapopen, gapextend, delta.matrix); 
                        TRange read_range = lcl_cigar.QueryRange();
                        TRange contig_range = lcl_cigar.SubjectRange();
                        rstart = read_range.first;
                        rstop = read_range.second;
                        cstop = cstart+contig_range.second;
                        cstart = cstart+contig_range.first;
                        
                        cigar = CCigar(rstop, cstop);
                        cigar.PushFront(lcl_cigar);
                        
                        flank_left = min(cstart, max(0,(match*rstart-gapopen)/max(match, gapextend))+rstart);
                        flank_right = min(contig_length-cstop-1, max(0,(match*(read_length-rstop-1)-gapopen)/max(match, gapextend))+read_length-rstop-1);
                        if(flank_left <= contig_range.first && flank_right <= llen-contig_range.second-1)
                            break;
                    }


                    if(rstop-rstart+1 >= coverage*read_length && previous_aligns.emplace(rstart,rstop,cstart,cstop).second) {
                        ++align_rslts.m_alignments;
                        int dist = cigar.Distance(read_seq.data(), contig_seq.data()); 
                        out << read_acc << '\t' << (strand > 0 ? 0 : 16) << '\t' << contig_acc << '\t' << cstart+1 << '\t' << 255 << '\t' << cigar.CigarString(0, read_length) << "\t*\t0\t0\t";
                        for(char c : read_seq)
                            out << c;                                        
                        out << "\t*\tNM:i:" << dist << endl;  
                                                           
                        if(!counts.empty())
                            IncrementCount(counts, cigar.ToAlign(read_seq.data(), contig_seq.data()), cstart);
                    }
                }                                    
            }                                                

            ReverseComplementSeq(read_seq.begin(), read_seq.end());
        }
    }

    out << flush;
    out.pop();
}

typedef tuple<TGenomeSeq*,size_t,size_t> TContigChunk;
typedef list<TContigChunk> THashJob;

void ContigHashInputs(int word, size_t job_length,  TGenomeHolder& contigs, list<THashJob>& job_inputs) {
    size_t assigned_length_from_contig = 0;
    for(auto it = contigs.begin(); it != contigs.end(); ) {
        job_inputs.emplace_back();
        size_t current_job_length = 0;
        while(current_job_length < job_length && it != contigs.end()) {
            size_t max_chunk = get<1>(*it).size() - assigned_length_from_contig;
            if(current_job_length + max_chunk <= job_length) {  // the rest of the file could be assigned           
                job_inputs.back().emplace_back(&(*it), assigned_length_from_contig, get<1>(*it).size()-1);
                current_job_length += max_chunk;
                ++it;
                assigned_length_from_contig = 0;
            } else {   // something left for another job            
                size_t chunk = job_length - current_job_length;
                size_t right_with_extra = min(get<1>(*it).size()-1, assigned_length_from_contig+chunk-1+word-1);
                job_inputs.back().emplace_back( &(*it), assigned_length_from_contig, right_with_extra);
                assigned_length_from_contig += chunk;
                current_job_length = job_length;
            }
        }
    }
}

void ContigHashJob(CKmerLocationsHash& kmer_locations, TGenomeHolder& contigs, THashJob& job, int word_size, double entropy_level) {
    for(auto& chunk : job) {
        auto& contig = *get<0>(chunk);
        auto& seq = get<1>(contig);
        auto end = seq.begin()+get<2>(chunk)+1;
        for(auto left = FirstACGT(seq.begin()+get<1>(chunk), end); left < end; ) {
            auto right = FirstNotACGT(left+1, end);
            int len = right-left;
            if(len >= word_size) {
                CReadHolder rh(false);
                rh.PushBack(left, len);
                int pos = right-seq.begin()-word_size;
                for(CReadHolder::kmer_iterator ik = rh.kbegin(word_size) ; ik != rh.kend(); ++ik, --pos) { // iteration from last kmer to first
                    double entropy = Entropy(seq.begin()+pos, word_size);
                    if(entropy > entropy_level)
                        kmer_locations.Insert(*ik, pos, &contig);
                }
            }
            left = FirstACGT(right, end);
        }
    }    
}




int main(int argc, const char* argv[]) {
    for(int n = 0; n < argc; ++n)
        cerr << argv[n] << " ";
    cerr << endl << endl;

    options_description general("General options");
    general.add_options()
        ("help,h", "Produce help message")
        ("version,v", "Print version")
        ("cores", value<int>()->default_value(0), "Number of cores to use (default all) [integer]");

    options_description input("Input/output");
    input.add_options()
        ("sra_run", value<string>(), "Input sra run accession [string]")
        ("fasta", value<string>(), "Input fasta for reads [string]")
        ("contigs", value<string>(), "Input fasta for contigs [string]")
        ("out", value<string>(), "Output file (will be gzipped, suffix gz added) [string]")
        ("consensus", value<string>(), "Output file consensus data [string]");

    options_description word_params("Parameters for finding words");
    word_params.add_options()
        ("word", value<int>()->default_value(16), "Word size [integer]")
        ("min_compart", value<double>()->default_value(2, "2.0"), "Minimal compartment in words [float]")
        ("entropy", value<double>()->default_value(0.51, "0.51"), "Minimal entropy for a word [float]");

    options_description sw_params("Parameters for Smith-Waterman alignment");
    sw_params.add_options()
        ("match", value<int>()->default_value(1), "Match penalty (>0) [integer]")
        ("mismatch", value<int>()->default_value(2), "Mismatch penalty (>0) [integer]")
        ("gapopen", value<int>()->default_value(5), "Penalty for opening gap (>=0) [integer]")
        ("gapextend", value<int>()->default_value(2), "Penalty for extending gap gap (>0) [integer]")
        ("coverage", value<double>()->default_value(0.5, "0.5"), "Threshold coverage for alignment output [float]");

    options_description all("");
    all.add(general).add(input).add(word_params).add(sw_params); 

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
            cerr << "k_aligner v.1.0" << endl;
#ifdef SVN_REV
            cerr << "SVN revision:" << SVN_REV << endl << endl;
#endif
            return 0;
        }

        if(!argm.count("sra_run") && !argm.count("fasta")) {
            cerr << "Provide reads" << endl;
            cerr << all << "\n";
            return 1;
        }

        if(!argm.count("contigs")) {
            cerr << "Provide some input contigs" << endl;
            cerr << all << "\n";
            return 1;
        }

        if(!argm.count("out")) {
            cerr << "Provide output file" << endl;
            cerr << all << "\n";
            return 1;
        }

        int ncores = thread::hardware_concurrency();
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

        SAlignParams align_params(argm);

        list<TSeqHolder> read_chunks;

        //read fasta reads
        if(argm.count("fasta")) {
            cerr << "Getting reads from fasta..." << endl;
            CStopWatch sw;
            sw.Restart();
            TSeqHolder reads;
            ifstream fasta(argm["fasta"].as<string>());
            if(!fasta.is_open()) {
                cerr << "Can't open file " << argm["fasta"].as<string>() << endl;
                return 1;             
            }
            char c;
            if(!(fasta >> c) || c != '>')
                throw runtime_error("Invalid fasta file format for reads");
            string record;
            while(getline(fasta, record, '>')) {
                size_t first_ret = min(record.size(),record.find('\n'));
                if(first_ret == string::npos)
                    throw runtime_error("Invalid fasta file format for reads");
                string acc = record.substr(0, first_ret);
                acc = acc.substr(0, acc.find_first_of(" \t"));
                string read = record.substr(first_ret+1);
                read.erase(remove(read.begin(),read.end(),'\n'),read.end());
                for(char& c : read) c = toupper(c);
                reads.push_back(TSeq(acc, TSeqV(read.begin(), read.end())));
            }
            size_t total = reads.size();
            size_t job_length = total/ncores+1;
            size_t remaining = total;
            while(!reads.empty()) {
                read_chunks.emplace_back();
                size_t chunk = min(job_length, remaining);
                read_chunks.back().splice(read_chunks.back().end(), reads, reads.begin(), next(reads.begin(), chunk));
                remaining -= chunk;;
            }

            cerr << "Reads from fasta acquired in " << sw.Elapsed() << endl;
        }

        
        if(argm.count("sra_run")) {
            cerr << "Getting reads from SRA..." << endl;
            CStopWatch sw;
            sw.Restart();
            ReadCollection run = ncbi::NGS::openReadCollection(align_params.m_sra_run);
            size_t count = run.getReadCount();
            size_t job_length = count/ncores+1;

            list<function<void()>> jobs;
            size_t from = 0;
            while(from < count) {
                size_t to = min(count-1, from+job_length-1);
                read_chunks.push_back(TSeqHolder());
                jobs.push_back(bind(GetReadsFromSRAJob, ref(align_params.m_sra_run), from, to, ref(read_chunks.back())));
                from = to+1;
            }
            RunThreads(ncores, jobs);

            cerr << "Reads from SRA acquired in " << sw.Elapsed() << endl;
        }
                                
        size_t total = 0;
        for(auto& lst : read_chunks)
            total += lst.size();            
        cerr << "Reads: " << total << endl;

        TGenomeHolder contigs;
        size_t total_len = 0;

        //read fasta
        {
            CStopWatch sw;
            sw.Restart();
            ifstream fasta(argm["contigs"].as<string>());
            if(!fasta.is_open()) {
                cerr << "Can't open file " << argm["contigs"].as<string>() << endl;
                return 1;             
            }
            char c;
            if(!(fasta >> c) || c != '>')
                throw runtime_error("Invalid fasta file format for contigs");
            string record;
            while(getline(fasta, record, '>')) {
                size_t first_ret = min(record.size(),record.find('\n'));
                if(first_ret == string::npos)
                    throw runtime_error("Invalid fasta file format for contigs");
                string acc = record.substr(0, first_ret);
                acc = acc.substr(0, acc.find_first_of(" \t"));
                string contig = record.substr(first_ret+1);
                contig.erase(remove(contig.begin(),contig.end(),'\n'),contig.end());
                for(char& c : contig) c = toupper(c);
                //                contigs.push_back(make_tuple(acc, TSeqV(contig.begin(), contig.end())), TCounterV());
                contigs.emplace_back();
                get<0>(contigs.back()) = acc;
                get<1>(contigs.back()).insert(get<1>(contigs.back()).end(), contig.begin(), contig.end());
                if(argm.count("consensus"))
                    get<2>(contigs.back()).resize(contig.size());
                total_len += contig.size();
            }
            cerr << "Contigs: " << contigs.size() << endl;
            cerr << "Contigs input in " << sw.Elapsed() << endl;
        }

        CKmerLocationsHash kmer_locations(align_params.m_word_size, total_len);
        //create kmer hash
        {
            CStopWatch sw;
            sw.Restart();

            size_t job_length = total_len/ncores+1;
            list<THashJob> job_inputs;
            ContigHashInputs(align_params.m_word_size, job_length,  contigs, job_inputs);

            list<function<void()>> jobs;
            for(auto& job_input : job_inputs) {
                jobs.push_back(bind(ContigHashJob, ref(kmer_locations), ref(contigs), ref(job_input), align_params.m_word_size, align_params.m_entropy_level)); 
            }
                
            RunThreads(ncores, jobs);
            cerr << "Hash table: " << kmer_locations.TableSize() << endl;
            cerr << "Contig hash in " << sw.Elapsed() << endl;            
        }

                
        list<SAlignRslt> align_rslts;
        {
            CStopWatch sw;
            sw.Restart();

            list<function<void()>> jobs;
            for(auto& reads : read_chunks) {
                align_rslts.push_back(SAlignRslt());
                jobs.push_back(bind(GetAlignsJob, ref(reads), ref(kmer_locations), ref(align_params), ref(align_rslts.back())));
            }
            RunThreads(ncores, jobs);

            size_t hits = 0;
            size_t alignments = 0;
            size_t attempts = 0;
            for(auto& rslt : align_rslts) {
                hits += rslt.m_hits;
                alignments += rslt.m_alignments;
                attempts += rslt.m_attempts;                            
            }
            cerr << "Hits: " << hits << " Alignment attempts: " << attempts << " Alignments: " << alignments << endl;
            cerr << "Alignments in " << sw.Elapsed() << endl;

            if(alignments > 0) {
                sw.Restart();
                string file = argm["out"].as<string>()+".gz";
                ofstream out(file);
                if(!out.is_open()) {
                    cerr << "Can't open file " << file << endl;
                    return 1;             
                }
                out.close();
                for(auto& rslt : align_rslts) {
                    if(rslt.m_alignments > 0) {
                        array_source source{rslt.m_out.data(), rslt.m_out.size()};
                        filtering_istream is;
                        is.push(source);
                        out.open(file, ofstream::out | ofstream::app);
                        copy(is, out); // closes streams               
                    }
                }
                if(!out) {
                    cerr << "Can't write to file " << file << endl;
                    exit(1);
                }

                if(argm.count("consensus")) {
                    ofstream cons_out(argm["consensus"].as<string>());
                    if(!cons_out.is_open()) {
                        cerr << "Can't open file " << argm["consensus"].as<string>() << endl;
                        return 1;             
                    }
                    for(auto& contig : contigs) {
                        string& acc = get<0>(contig);
                        auto& seq = get<1>(contig);
                        auto& counts = get<2>(contig);
                        for(unsigned pos = 0; pos < seq.size(); ++pos) {
                            cons_out << acc << '\t' << pos+1 << '\t' << seq[pos];
                            vector<TCnt> sorted(counts[pos].begin(), counts[pos].end());
                            sort(sorted.begin(), sorted.end(), [](const TCnt& a, const TCnt b) { if(a.second == b.second) return a.first < b.first; else return a.second > b.second; });
                            for(auto& cnt : sorted)
                                cons_out << '\t' << cnt.first << ':' << cnt.second.Load();
                            cons_out << "\n";
                        }
                    }
                    cons_out.close();
                    if(!cons_out) {
                        cerr << "Can't write to file " << argm["consensus"].as<string>() << endl;
                        exit(1);
                    }
                }

                cerr << "Output in " << sw.Elapsed() << endl;
            } else {
                cerr << "No alignments found" << endl;
            }
        }
                
    
        cerr << "DONE" << endl;
        exit(0);
    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        exit(1);
    }
}
