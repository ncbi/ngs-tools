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

#ifndef _DeBruijn_Graph_
#define _DeBruijn_Graph_

#include <iostream>
#include <bitset>
#include <unordered_map>
#include <unordered_set>
#include <deque>
#include <forward_list>
#include <atomic>
#include <future>
#include <thread>
#include <boost/timer/timer.hpp>
#include <cmath>

#include "Integer.hpp"

// This file contains classes which facilitate basic operation of storing reads, counting kmers,
// and creating and traversing a de Bruijn graph

using namespace std;
namespace DeBruijn {

    // Wraps around atomic<> to make it possible to use in containers
    // IMPORTANT: don't concurrently create or modify containers of SAtomic!!!!!
    template <typename T>
    struct SAtomic {
        typedef T Type;
        SAtomic(T t = 0) { m_atomic.store(t); }
        SAtomic(const atomic<T> &a) { m_atomic.store(a.load()); }
        SAtomic(const SAtomic &other) { m_atomic.store(other.m_atomic.load()); }
        SAtomic& operator=(const SAtomic &other) { 
            m_atomic.store(other.m_atomic.load());
            return *this;
        }        
        SAtomic& operator=(T t) { 
            m_atomic.store(t); 
            return *this;
        } 
        bool Set(T value, T expected = 0) { return m_atomic.compare_exchange_strong(expected, value); }
        operator T() const { return m_atomic.load(); }

        atomic<T> m_atomic;
    };

    class CStopWatch : public boost::timer::cpu_timer {
    public:
        void Restart() { start(); }
        string Elapsed() const { return format(); }
    };


    // runs ncores threads until all jobs are exhausted
    void RunThreads(int ncores, list<function<void()>>& jobs) {
        typedef list<future<void>> ThreadsStatus;
        ThreadsStatus active_threads_status;

        //        int total_jobs = jobs.size();
        //        cerr << "Remaining " << total_jobs << " jobs from " << total_jobs << endl;

        //create ncores threads
        for(int i = 0; i < ncores && !jobs.empty(); ++i) {
            active_threads_status.push_front(async(launch::async, jobs.front()));
            jobs.pop_front();
        }

        //for each finished thread create a new one until done  
        chrono::milliseconds span (1);
        while(!active_threads_status.empty()) {
            for(auto iloop = active_threads_status.begin(); iloop != active_threads_status.end(); ) {
                auto done = iloop++;
                if(done->wait_for(span) == future_status::timeout)   // not ready    
                    continue;                

                done->get();
                active_threads_status.erase(done);
                if(!jobs.empty()) {
                    active_threads_status.push_front(async(launch::async, jobs.front()));
                    jobs.pop_front();
                }
                //                cerr << "Remaining jobs " << jobs.size()+active_threads_status.size() << " from " << total_jobs << endl;    
            }
        }
        //        cerr << endl; 
    }


    class CKmerCount {
    // Class for kmer counting and searching implemented using a boost::variant of vector<pair<LargeInt<N>,size_t>>
    // Currently, maximum N defined in config.hpp is 16 that allows kmers of length at most 512 to be stored.
    // Only smaller (in the bit encoding) of a kmer and its reverse complement is stored
    //
    // When vector is sorted, binary search on the first element of the pair that represents a kmer can be used for retrieving
    // the information stored for the kmer in the second element of the pair.
    // First 32 bits of the second element stores total count for kmer (self and reverse complement)
    // Remaining 32 bits store count for kmer for self only during the counting operations but are modified to additionally store
    // branching information when used inside CDBGraph

    public:
        typedef TKmerCountN Type;

        CKmerCount(int kmer_len = 0) : m_kmer_len(kmer_len) {
            if(m_kmer_len > 0)
                m_container = CreateVariant<TKmerCountN, TLargeIntVec>((m_kmer_len+31)/32);
        }
        size_t Size() const { return apply_visitor(container_size(), m_container); }         // number of elements in the container
        void Reserve(size_t rsrv) { apply_visitor(reserve(rsrv), m_container); }             // reserves memory for rsrv elements
        void Clear() { apply_visitor(clear(), m_container); }                                // clears container (doesn't release memory)
        size_t Capacity() const { return apply_visitor(container_capacity(), m_container); } // tells how many elements could be stored in reserved memory
        size_t ElementSize() const { return apply_visitor(element_size(), m_container); }    // size of one vector element in bytes
        size_t MemoryFootprint() const { return Capacity()*ElementSize(); }                  // reserved memory in bytes
        void PushBack(const TKmer& kmer, size_t count) {                                     // push back one element
            if(m_kmer_len == 0)
                throw runtime_error("Can't insert in uninitialized container");
            apply_visitor(push_back(kmer, count), m_container); 
        }
        void PushBackElementsFrom(const CKmerCount& other) {         // push back elements from other container
            if(m_kmer_len == 0)
                throw runtime_error("Can't insert in uninitialized container");
            apply_visitor(push_back_elements(), m_container, other.m_container); 
        }
        size_t Find(const TKmer& kmer) const { return apply_visitor(find_kmer(kmer), m_container); } // finds index for a kmer (returns Size() if not found)
        void UpdateCount(size_t count, size_t index) { apply_visitor(update_count(count, index), m_container); }          // updates count at the index position
        size_t GetCount(size_t index) const { return apply_visitor(get_count(index), m_container); }                      // gets count at the index position
        pair<TKmer,size_t> GetKmerCount(size_t index) const { return apply_visitor(get_kmer_count(index), m_container); } // gets kmer and count at the index position
        const uint64_t* getPointer(size_t index) { return apply_visitor(get_pointer(index), m_container); }               // gets access to binary kmer sequence
        int KmerLen() const { return m_kmer_len; }
        void Sort() { apply_visitor(container_sort(), m_container); }
        void SortAndExtractUniq(int min_count, CKmerCount& uniq) {  // sorts container, aggregates counts, copies elements with count >= min_count into uniq
            uniq = CKmerCount(m_kmer_len); // init
            Sort();
            apply_visitor(extract_uniq(min_count), m_container, uniq.m_container);
        }
        void SortAndUniq(int min_count) { // sorts container, aggregate counts, keeps elements with count >= min_count
            Sort();
            apply_visitor(uniq(min_count), m_container);
        }
        void MergeTwoSorted(const CKmerCount& other) { // merges with other assuming both sorted
            if(m_kmer_len != other.KmerLen())
                throw runtime_error("Can't merge kmers of different lengths");
            apply_visitor(merge_sorted(), m_container, other.m_container);
        }
        void Swap(CKmerCount& other) { // swaps with other
            swap(m_kmer_len, other.m_kmer_len);
            apply_visitor(swap_with_other(), m_container, other.m_container);    
        }
        void Save(ostream& out) const { 
            out.write(reinterpret_cast<const char*>(&m_kmer_len), sizeof(m_kmer_len));
            apply_visitor(save(out), m_container);
        }
        void Load(istream& in) {
            in.read(reinterpret_cast<char*>(&m_kmer_len), sizeof(m_kmer_len));
            m_container = CreateVariant<TKmerCountN, TLargeIntVec>((m_kmer_len+31)/32);
            apply_visitor(load(in), m_container);
        }

    private:

        struct find_kmer : public boost::static_visitor<size_t> { 
            find_kmer(const TKmer& k) : kmer(k) {}
            template <typename T> size_t operator()(const T& v) const { 
                typedef typename T::value_type pair_t;
                typedef typename pair_t::first_type large_t;
                auto it = lower_bound(v.begin(), v.end(), kmer.get<large_t>(), [](const pair_t& element, const large_t& target){ return element.first < target; });
                if(it == v.end() || it->first != kmer.get<large_t>())
                    return v.size();
                else
                    return it-v.begin();
            } 
            const TKmer& kmer;
        };
        struct reserve : public boost::static_visitor<> { 
            reserve(size_t r) : rsrv(r) {}
            template <typename T> void operator() (T& v) const { v.reserve(rsrv); }
            size_t rsrv;
        };
        struct container_size : public boost::static_visitor<size_t> { template <typename T> size_t operator()(const T& v) const { return v.size();} };
        struct container_capacity : public boost::static_visitor<size_t> { template <typename T> size_t operator()(const T& v) const { return v.capacity();} };
        struct element_size : public boost::static_visitor<size_t> { template <typename T> size_t operator()(const T& v) const { return sizeof(typename  T::value_type);} };
        struct clear : public boost::static_visitor<> { template <typename T> void operator()(T& v) const { v.clear();} };
        struct push_back : public boost::static_visitor<> { 
            push_back(const TKmer& k, size_t c) : kmer(k), count(c) {}
            template <typename T> void operator() (T& v) const {
                typedef typename  T::value_type::first_type large_t;
                v.push_back(make_pair(kmer.get<large_t>(), count)); 
            }
            const TKmer& kmer;
            size_t count;
        };
        struct push_back_elements : public boost::static_visitor<> {
            template <typename T> void operator() (T& a, const T& b) const { a.insert(a.end(), b.begin(), b.end()); }
            template <typename T, typename U> void operator() (T& a, const U& b) const { throw runtime_error("Can't copy from different type container"); }
        };
        struct merge_sorted : public boost::static_visitor<> {
            template <typename T> void operator() (T& a, const T& b) const { 
                T merged;
                merged.reserve(a.size()+b.size());
                merge(a.begin(), a.end(), b.begin(), b.end(), back_inserter(merged));
                merged.swap(a); 
            }
            template <typename T, typename U> void operator() (T& a, const U& b) const { throw runtime_error("Can't merge different type containers"); }
        };
        struct update_count : public boost::static_visitor<> {
            update_count(size_t c, size_t i) : count(c), index(i) {}
            template <typename T> void operator() (T& v) const { v[index].second = count; }  
            size_t count;
            size_t index;
        };
        struct get_count : public boost::static_visitor<size_t> {
            get_count(size_t i) : index(i) {}
            template <typename T> size_t operator() (T& v) const { return v[index].second; }        
            size_t index;
        };
        struct get_kmer_count : public boost::static_visitor<pair<TKmer,size_t>> {
            get_kmer_count(size_t i) : index(i) {}
            template <typename T> pair<TKmer,size_t> operator() (T& v) const { return make_pair(TKmer(v[index].first), v[index].second); }        
            size_t index;
        };
        struct get_pointer : public boost::static_visitor<const uint64_t*> {
            get_pointer(size_t i) : index(i) {}
            template <typename T> const uint64_t* operator() (T& v) const { return v[index].first.getPointer(); }        
            size_t index;
        };
        struct container_sort : public boost::static_visitor<> { template <typename T> void operator() (T& v) const { sort(v.begin(), v.end()); }};
        struct swap_with_other : public boost::static_visitor<> { 
            template <typename T> void operator() (T& a, T& b) const { a.swap(b); }
            template <typename T, typename U> void operator() (T& a, U& b)  const { throw runtime_error("Can't swap different type containers"); }
        };
        struct uniq : public boost::static_visitor<> {
            uniq(int mc) : min_count(mc) {}
            template <typename T> void operator() (T& v) const {
                typedef typename T::iterator iter_t;
                iter_t nextp = v.begin();
                for(iter_t ip = v.begin(); ip != v.end(); ) {
                    iter_t workp = ip;
                    while(++ip != v.end() && workp->first == ip->first)
                        workp->second += ip->second;   // accumulate all 8 bytes; we assume that count will not spill into higher half
                    if((uint32_t)workp->second >= min_count)
                        *nextp++ = *workp;
                }
                v.erase(nextp, v.end());
            }

            unsigned min_count;
        };
        struct extract_uniq : public boost::static_visitor<> {
            extract_uniq(int mc) : min_count(mc) {}
            template <typename T> void operator() (T& a, T& b) const {
                if(a.empty()) return;
                size_t num = 1;
                uint32_t count = a[0].second;  // count only 4 bytes!!!!!!
                for(size_t i = 1; i < a.size(); ++i) {
                    if(a[i-1].first < a[i].first) {
                        if(count >= min_count)
                            ++num;
                        count = a[i].second;
                    } else {
                        count += a[i].second;
                    }
                }
                if(count < min_count)
                    --num;
                b.reserve(num+1);
                b.push_back(a[0]);
                for(size_t i = 1; i < a.size(); ++i) {
                    if(b.back().first < a[i].first) {
                        if((uint32_t)b.back().second < min_count)
                            b.pop_back();
                        b.push_back(a[i]);
                    } else {
                        b.back().second += a[i].second;  // accumulate all 8 bytes; we assume that count will not spill into higher half
                    }
                }
                if((uint32_t)b.back().second < min_count)
                    b.pop_back();
            }
            template <typename T, typename U> void operator() (T& a, U& b) const { throw runtime_error("Can't extract into different type container"); }

            unsigned min_count;
        };
        struct save : public boost::static_visitor<> {
            save(ostream& out) : os(out) {}
            template <typename T> void operator() (T& v) const {
                size_t num = v.size();
                os.write(reinterpret_cast<const char*>(&num), sizeof num); 
                if(num > 0)
                    os.write(reinterpret_cast<const char*>(&v[0]), num*sizeof(v[0]));
            }
            ostream& os;
        };
        struct load : public boost::static_visitor<> {
            load(istream& in) : is(in) {}
            template <typename T> void operator() (T& v) const {
                size_t num;
                is.read(reinterpret_cast<char*>(&num), sizeof num);
                if(num > 0) {
                    v.resize(num);
                    is.read(reinterpret_cast<char*>(&v[0]), num*sizeof(v[0]));
                }
            }
            istream& is;
        };

        Type m_container;
        int m_kmer_len;
    };
    typedef CKmerCount TKmerCount; // for compatibility with previous code
    
    typedef vector<pair<int,size_t>> TBins; // pair of position,count
    
    template <typename V> class CKmerMap {
    // A hash with kmer as a key
    // Implemented using a boost::variant of unordered_map<<LargeInt<N>,V> with maximal N = 16 which allows kmer size up to 512

    public:
        typedef V MappedType;
        typedef TKmerMapN<V> Type;

        CKmerMap(int kmer_len = 0) : m_kmer_len(kmer_len) {
            if(m_kmer_len > 0)
                m_container = CreateVariant<TKmerMapN<V>, TLargeIntMap, V>((m_kmer_len+31)/32);
        }
        size_t Size() const { return apply_visitor(container_size(), m_container); }  // number of elements in the container
        void Reserve(size_t rsrv) { apply_visitor(reserve(rsrv), m_container); }      // reserves hash table for rsrv elements
        void Clear() { apply_visitor(clear(), m_container); }                         // clear hash table 
        V& operator[] (const TKmer& kmer) { 
            if(m_kmer_len == 0) 
                throw runtime_error("Can't insert in uninitialized container");
            return apply_visitor(mapper(kmer), m_container);
        }
        V* Find(const TKmer& kmer) { return apply_visitor(find(kmer), m_container); } // returns nullptr if not found
        int KmerLen() const { return m_kmer_len; }

        template <typename Prob> 
        void GetInfo(Prob& prob) { apply_visitor(get_info<Prob>(prob), m_container); } // scans the containier and calls prob(v) for each mapped element

    private:

        template <typename Prob>
        struct get_info : public boost::static_visitor<> {
            get_info(Prob& p) : prob(p) {}
            template <typename T> void operator()(T& v) const {
                for(auto& val : v)
                    prob(val.second);
            }

            Prob& prob;
        };
        struct container_size : public boost::static_visitor<size_t> { template <typename T> size_t operator()(const T& v) const { return v.size();} };
        struct clear : public boost::static_visitor<> { template <typename T> void operator()(const T& v) const { v.clear();} };
        struct reserve : public boost::static_visitor<> { 
            reserve(size_t r) : rsrv(r) {}
            template <typename T> void operator() (T& v) const { v.reserve(rsrv); }
            size_t rsrv;
        };
        struct mapper : public boost::static_visitor<V&> { 
            mapper(const TKmer& k) : kmer(k) {}
            template <typename T> V& operator()(T& v) const { 
                typedef typename  T::key_type large_t;
                return v[kmer.get<large_t>()];
            } 
            const TKmer& kmer;
        };
        struct find : public boost::static_visitor<V*> {
            find(const TKmer& k) : kmer(k) {}
            template <typename T> V* operator()(T& v) const { 
                typedef typename  T::key_type large_t;
                typename T::iterator it = v.find(kmer.get<large_t>());
                if(it != v.end())
                    return &(it->second);
                else
                    return 0;
            } 
            const TKmer& kmer;
        };

        Type m_container;
        int m_kmer_len;
    };
    template <typename V>
    using  TKmerMap = CKmerMap<V>; // for compatibility with previous code

    // simple heuristic to find a valley/peak in a histogram
    int FindValleyAndPeak(const TBins& bins, int rlimit) {
        int SLOPE_LEN = 5;
        int peak = min(rlimit,(int)bins.size()-SLOPE_LEN-1);
        while(peak >= SLOPE_LEN) {
            bool maxim = true;
            for(int i = 1; i <= SLOPE_LEN && maxim; ++i) 
                maxim = bins[peak+i].second < bins[peak].second;
            for(int i = 1; i <= SLOPE_LEN && maxim; ++i) 
                maxim = bins[peak-i].second < bins[peak].second;
            if(maxim)
                break;
            --peak;
        }

        if(peak < SLOPE_LEN)
            return -1;

        int valley = 0;
        for(int i = 1; i <= peak; ++i) {
            if(bins[i].second < bins[valley].second)
                valley = i;
        }
        if(valley == peak)
            return -1;

        for(int i = valley; i < (int)bins.size(); ++i) {
            if(bins[i].second > bins[peak].second)
                peak = i;
        }
        
        if(bins[valley].second < 0.7*bins[peak].second)
            return valley;
        else            
            return -1;
    }

    // a simple heuristic to find main range in a histogram
    pair<int,int> HistogramRange(const TBins& bins) {  // returns <valley,rlimit>; valley == -1 if not found
        unsigned MIN_NUM = 100;
        size_t gsize = 0;
        for(auto& bin : bins) {
            if(bin.second >= MIN_NUM) 
                gsize += bin.first*bin.second;
        }

        int rl = 0;
        size_t gs = 0;
        for(auto& bin : bins) {
            gs += bin.first*bin.second;
            ++rl;
            if(gs > 0.8*gsize)
                break;
        }

        int valley = -1;
        int rlimit = rl;
        size_t genome = 0;

        while(true) {
            int v = FindValleyAndPeak(bins, rl);

            size_t g = 0;
            for(int i = max(0, v); i <= rl; ++i)
                g += bins[i].second;

            if((v >= 0 && g > genome) || g > 10*genome) {
                valley = v;
                rlimit = rl;
                genome = g;
                //                cerr << valley << " " << rlimit << " " << genome << endl;
            }

            if(v < 0)
                break;
            rl = v;
        }
        
        return make_pair(valley, rlimit);
    }

    // Implementation of de Bruijn graph based on TKmerCount which stores kmer (smaller in the bit encoding of self and its reverse
    // complement), its count, fraction of times the stored kmer was seen as self, and information for presence/absence in graph
    // for each of the eight possible extensions to which this kmer can be connected
    // Allows basic traversing operations such as find kmer and its abundance (count) or find successors for a kmer
    // We use a node-centric definition of de Bruijn graph in which nodes of the graph are kmers
    class CDBGraph {
    public:

        // Construct graph from counted kmers and histogram
        // is_stranded indicates if count include reliable direction information (PlusFraction() and MinusFraction() could be used)
        CDBGraph(const TKmerCount& kmers, const TBins& bins, bool is_stranded) : m_graph_kmers(kmers.KmerLen()), m_bins(bins), m_is_stranded(is_stranded) {
            m_graph_kmers.PushBackElementsFrom(kmers);
            string max_kmer(m_graph_kmers.KmerLen(), bin2NT[3]);
            m_max_kmer = TKmer(max_kmer);
            m_visited.resize(GraphSize(), 0);
        }

        // Construct graph from temporary containers
        CDBGraph(TKmerCount&& kmers, TBins&& bins, bool is_stranded) :  m_graph_kmers(kmers.KmerLen()), m_is_stranded(is_stranded) {
            m_graph_kmers.Swap(kmers);
            m_bins.swap(bins);
            string max_kmer(m_graph_kmers.KmerLen(), bin2NT[3]);
            m_max_kmer = TKmer(max_kmer);
            m_visited.resize(GraphSize(), 0);
        }

        // Load from a file
        CDBGraph(istream& in) {
            m_graph_kmers.Load(in);
            string max_kmer(m_graph_kmers.KmerLen(), bin2NT[3]);
            m_max_kmer = TKmer(max_kmer);

            int bin_num;
            in.read(reinterpret_cast<char*>(&bin_num), sizeof bin_num);
            for(int i = 0; i < bin_num; ++i) {
                pair<int, size_t> bin;
                in.read(reinterpret_cast<char*>(&bin), sizeof bin);
                m_bins.push_back(bin);
            }

            in.read(reinterpret_cast<char*>(&m_is_stranded), sizeof m_is_stranded);
            m_visited.resize(GraphSize(), 0);
        }

        // Save in a file
        void Save(ostream& out) const {
            m_graph_kmers.Save(out);
            int bin_num = m_bins.size();
            out.write(reinterpret_cast<const char*>(&bin_num), sizeof bin_num);         
            out.write(reinterpret_cast<const char*>(&m_bins[0]), bin_num*(sizeof m_bins[0])); 
            out.write(reinterpret_cast<const char*>(&m_is_stranded), sizeof m_is_stranded);
        }

        // These two functions map kmers to integer indexes which could be used to retrieve kmer properties
        // 0 is returned for kmers not present in the graph
        // positive even numbers are for stored kmers
        // positive odd numbers are for reverse complement of stored kmers 
        typedef size_t Node;
        Node GetNode(const TKmer& kmer) const {   // finds kmer in graph
            TKmer rkmer = revcomp(kmer, KmerLen());
            if(kmer < rkmer) {
                size_t index = m_graph_kmers.Find(kmer);
                return (index == GraphSize() ? 0 : 2*(index+1));
            } else {
                size_t index = m_graph_kmers.Find(rkmer);
                return (index == GraphSize() ? 0 : 2*(index+1)+1);
            }
        }
        Node GetNode(const string& kmer_seq) const {   // finds kmer in graph
            if(kmer_seq.find_first_not_of("ACGT") != string::npos || (int)kmer_seq.size() != KmerLen())   // invalid kmer
                return 0;
            TKmer kmer(kmer_seq);
            return GetNode(kmer);
        }

        // for all access with Node there is NO check that node is in range !!!!!!!!
        int Abundance(const Node& node) const { // total count for a kmer
            if(node == 0)
                return 0;
            else
                return m_graph_kmers.GetKmerCount(node/2-1).second;  // automatically clips out branching information!
        }
        // 32 bit count; 8 bit branching; 8 bit not used yet; 16 bit +/-
        double MinusFraction(const Node& node) const {  // fraction of the times kmer was seen in - direction
            double plusf = PlusFraction(node);
            return min(plusf,1-plusf);
        }
        double PlusFraction(const Node& node) const {  // fraction of the times kmer was seen in + direction
            double plusf = double(m_graph_kmers.GetKmerCount(node/2-1).second >> 48)/numeric_limits<uint16_t>::max();
            if(node%2)
                plusf = 1-plusf;
            return plusf;
        }
        TKmer GetNodeKmer(const Node& node) const {  // returns kmer as TKmer
            if(node%2 == 0) 
                return m_graph_kmers.GetKmerCount(node/2-1).first;
            else
                return revcomp(m_graph_kmers.GetKmerCount(node/2-1).first, KmerLen());
        }
        string GetNodeSeq(const Node& node) const { // returnd kmer as string
            return GetNodeKmer(node).toString(KmerLen());
        }
        const uint64_t* getPointer(const Node& node) { return m_graph_kmers.getPointer(node/2-1); }        

        // multithread safe way to set visited value; returns true if value was as expected before and has been successfully changed
        // 1 is used for permanent holding; 2 is used for temporary holding
        bool SetVisited(const Node& node, uint8_t value=1, uint8_t expected=0) {
            return m_visited[node/2-1].Set(value, expected);
        }

        bool ClearVisited(const Node& node) { // multithread safe way to clear visited value; returns true if value was set before
            if(m_visited[node/2-1].Set(0, 1))
                return true;
            else
                return m_visited[node/2-1].Set(0, 2);
        }

        uint8_t IsVisited(const Node& node) const { // returns visited value
            return m_visited[node/2-1];
        }

        void ClearHoldings() { // clears temporary holdings
            for(auto& v : m_visited) 
                if(v == 2) v = 0;            
        }
        
        struct Successor {
            Successor(const Node& node, char c) : m_node(node), m_nt(c) {}
            Node m_node;
            char m_nt;
        };
        
        // Returns successors of a node 
        // These are nodes representing kmers produced by extending the right end of the kmer for
        // this node by one base and removing the leftmost base of the kmer
        // Each successor stores the successor's node and the extra base
        // Finding predecessors is done by finding successors of reverse complement of the kmer for the node
        vector<Successor> GetNodeSuccessors(const Node& node) const {
            vector<Successor> successors;
            if(!node)
                return successors;

            uint8_t branch_info = (m_graph_kmers.GetCount(node/2-1) >> 32);
            bitset<4> branches(node%2 ? (branch_info >> 4) : branch_info);
            if(branches.count()) {
                TKmer shifted_kmer = (GetNodeKmer(node) << 2) & m_max_kmer;
                for(int nt = 0; nt < 4; ++nt) {
                    if(branches[nt]) {
                        Node successor = GetNode(shifted_kmer + TKmer(KmerLen(), nt));
                        successors.push_back(Successor(successor, bin2NT[nt]));
                    }
                }
            }

            return successors;            
        }

        // Revese complement node
        static Node ReverseComplement(Node node) {
            if(node != 0)
                node = (node%2 == 0 ? node+1 : node-1);

            return node;
        }

        int KmerLen() const { return m_graph_kmers.KmerLen(); }             // returns kmer length
        size_t GraphSize() const { return m_graph_kmers.Size(); }           // returns total number of elements
        size_t ElementSize() const { return m_graph_kmers.ElementSize(); }  // element size in bytes
        size_t MemoryFootprint() const {                                    // reserved memory in bytes
            return m_graph_kmers.MemoryFootprint()+m_visited.capacity()+sizeof(TBins::value_type)*m_bins.capacity(); 
        }
        bool GraphIsStranded() const { return m_is_stranded; }              // indicates if graph contains stranded information

        // returns minimum position for stored histogram
        int HistogramMinimum() const {
            pair<int,int> r = HistogramRange(m_bins);
            if(r.first < 0)
                return 0;
            else 
                return m_bins[r.first].first;
        }

        // useus simple heuristic to evaluate the genome size
        int GenomeSize() const {
            pair<int,int> grange =  HistogramRange(m_bins);
            if(grange.first < 0)
                grange.first = 0;
            int genome = 0;
            for(int i = grange.first; i <= grange.second; ++i)
                genome += m_bins[i].second;            
        
            return genome;
        }

        // returns histogram
        const TBins& GetBins() const { return m_bins; }

    private:

        TKmerCount m_graph_kmers;     // only the minimal kmers are stored  
        TKmer m_max_kmer;             // contains 1 in all kmer_len bit positions  
        TBins m_bins;
        vector<SAtomic<uint8_t>> m_visited;
        bool m_is_stranded;
    };


    // Stores DNA sequences using 4 letter alphabet
    // The sequences and kmers could be accessed sequentially using iterator-type classes
    // 
    class CReadHolder {
    public:
        CReadHolder(bool contains_paired) :  m_total_seq(0), m_front_shift(0), m_contains_paired(contains_paired) {};

        // inserts read at the end 
        void PushBack(const string& read) {
            int shift = (m_total_seq*2 + m_front_shift)%64;
            for(int i = (int)read.size()-1; i >= 0; --i) {  // put backward for kmer compatibility
                if(shift == 0)
                    m_storage.push_back(0);
                m_storage.back() += ((find(bin2NT.begin(), bin2NT.end(),  read[i]) - bin2NT.begin()) << shift);
                shift = (shift+2)%64;
            }
            m_read_length.push_back(read.size());
            m_total_seq += read.size();
        }

        // insert sequence from other container
        class string_iterator;
        void PushBack(const string_iterator& is) {
            size_t read_len = is.ReadLen();
            m_read_length.push_back(read_len);
            size_t destination_first_bit = m_front_shift+2*m_total_seq;
            m_total_seq += read_len;
            m_storage.resize((m_front_shift+2*m_total_seq+63)/64);

            const CReadHolder& other_holder = *is.m_readholderp;
            size_t bit_from = is.m_readholderp->m_front_shift+is.m_position;
            size_t bit_to = bit_from+2*read_len;
            other_holder.CopyBits(bit_from, bit_to, m_storage, destination_first_bit, m_storage.size());
        }

        // removes first sequence
        void PopFront() {
            m_total_seq -= m_read_length.front();
            if(m_total_seq == 0) {
                Clear();
            } else {
                int nextp = m_front_shift+2*m_read_length.front();
                m_read_length.pop_front();
                m_front_shift = nextp%64;
                for(int num = nextp/64; num > 0; --num)
                    m_storage.pop_front();
            }
        }

        // swaps contents with other
        void Swap(CReadHolder& other) {
            swap(m_storage, other.m_storage);
            swap(m_read_length, other.m_read_length);
            swap(m_total_seq, other.m_total_seq);
            swap(m_front_shift, other. m_front_shift);
        }

        // deletes all sequences and releases  memory
        void Clear() { CReadHolder(m_contains_paired).Swap(*this); }

        // Total nucleotide count of the sequnce
        size_t TotalSeq() const { return m_total_seq; }

        // Maximal length of included sequences
        size_t MaxLength() const { 
            if(m_read_length.empty())
                return 0;
            else
                return *max_element(m_read_length.begin(), m_read_length.end()); 
        }

        // the number of kmers of give length that could be generated
        size_t KmerNum(unsigned kmer_len) const {
            size_t num = 0;
            if(m_read_length.empty())
                return num;
            for(auto l : m_read_length) {
                if(l >= kmer_len)
                    num += l-kmer_len+1;
            }
            return num;
        }

        // total number of sequences
        size_t ReadNum() const { return m_read_length.size(); }

        size_t MemoryFootprint() const { return 8*m_storage.size()+4*m_read_length.size(); }  // memory in bytes

        // shortest sequence length at xx% of total length
        size_t NXX(double xx) const {
            vector<uint32_t> read_length(m_read_length.begin(), m_read_length.end());
            sort(read_length.begin(), read_length.end());
            size_t nxx = 0;
            size_t len = 0;
            for(int j = (int)read_length.size()-1; j >= 0 && len < xx*m_total_seq; --j) {
                nxx = read_length[j];
                len += read_length[j];
            }
            
            return nxx;
        }
        // shortest sequence length at 50% of total length
        size_t N50() const { return NXX(0.5); }

        // iterator-type clas to access kmers
        class kmer_iterator;
        kmer_iterator kend() const { return kmer_iterator(0, *this, 2*m_total_seq); }
        kmer_iterator kbegin(int kmer_len) const { return kmer_iterator(kmer_len, *this); }

        class kmer_iterator {
        public:

            // dereference operator; returns value!
            TKmer operator*() const {
                TKmer kmer(m_kmer_len, 0);
                uint64_t* guts = kmer.getPointer();
                size_t bit_from = m_readholderp->m_front_shift+m_position;
                size_t bit_to = bit_from+2*m_kmer_len;
                m_readholderp->CopyBits(bit_from, bit_to, guts, 0, (2*m_kmer_len+63)/64);
                
                return kmer;
            }

            // iterator advance
            kmer_iterator& operator++() {
                if(m_position == 2*(m_readholderp->m_total_seq-m_kmer_len)) {
                    m_position = 2*m_readholderp->m_total_seq;
                    return *this;
                }

                m_position += 2;
                if(++m_position_in_read == m_readholderp->m_read_length[m_read]-m_kmer_len+1) {
                    m_position += 2*(m_kmer_len-1);
                    ++m_read;
                    m_position_in_read = 0;
                    SkipShortReads();
                } 

                return *this;
            }
            // doesn't check read boundaries - should be used only if landing in the SAME read
            kmer_iterator& operator+=(int l) {
                m_position += 2*l;
                m_position_in_read += l;

                return *this;
            }

            friend bool operator==(kmer_iterator const& li, kmer_iterator const& ri) { return li.m_position == ri.m_position && li.m_readholderp == ri.m_readholderp; }
            friend bool operator!=(kmer_iterator const& li, kmer_iterator const& ri) { return li.m_position != ri.m_position || li.m_readholderp != ri.m_readholderp; }
            friend class CReadHolder;

        private:
            kmer_iterator(int kmer_len, const CReadHolder& rholder, size_t position = 0, size_t position_in_read = 0, size_t read = 0) : m_readholderp(&rholder), m_read(read), m_position(position), m_kmer_len(kmer_len), m_position_in_read(position_in_read) {
                SkipShortReads();
            }

            void SkipShortReads() {
                while(m_position < 2*m_readholderp->m_total_seq && m_read < m_readholderp->m_read_length.size() && m_readholderp->m_read_length[m_read] < m_kmer_len)
                    m_position += 2*m_readholderp->m_read_length[m_read++];                               
            }
            const CReadHolder* m_readholderp;
            size_t m_read;               // read number
            size_t m_position;           // BIT num in concatenated sequence
            uint32_t m_kmer_len;
            uint32_t m_position_in_read; // SYMBOL in read
        };

        // iterator-type clas to access reads
        class string_iterator;
        string_iterator send() const { return string_iterator(*this, 2*m_total_seq, m_read_length.size()); }
        string_iterator sbegin() const { return string_iterator(*this); }

        enum {eSingle = 0, eFirstMate = 1, eSecondMate = 2};
        class string_iterator {
        public:

            string operator*() const {
                int read_length = m_readholderp->m_read_length[m_read];
                string read;
                read.reserve(read_length);
                size_t position = m_position+m_readholderp->m_front_shift+2*(read_length-1);
                for(int i = 0; i < read_length; ++i) {
                    read.push_back(bin2NT[(m_readholderp->m_storage[position/64] >> position%64) & 3]);
                    position -= 2;
                }
                return read;
            }
            string_iterator& operator++() {
                if(m_read == m_readholderp->m_read_length.size())
                    return *this;
                m_position +=  2*m_readholderp->m_read_length[m_read++]; 
                return *this;
            }
            size_t ReadLen() const { return m_readholderp->m_read_length[m_read]; }
            kmer_iterator KmersForRead(int kmer_len) const {
                if(kmer_len <= (int)m_readholderp->m_read_length[m_read])
                    return kmer_iterator(kmer_len, *m_readholderp, m_position, 0, m_read); 
                else
                    return m_readholderp->kend();
            }
            size_t Hash() const { 
                hash<const CReadHolder*> h1;
                hash<size_t> h2;
                return h1(m_readholderp)^h2(m_position); 
            }
            bool HasMate() const { return m_readholderp->m_contains_paired; }
            int PairType() const {
                if(!m_readholderp->m_contains_paired)
                    return eSingle;
                else if(m_read%2) // odd
                    return eSecondMate;
                else              // even
                    return eFirstMate;
            }
            string_iterator GetMate() const { // undefined behavior if not paired container
                if(m_read%2) // odd
                    return string_iterator(*m_readholderp, m_position-2*m_readholderp->m_read_length[m_read-1], m_read-1);
                else         // even
                    return string_iterator(*m_readholderp, m_position+2*m_readholderp->m_read_length[m_read], m_read+1);
            }

            friend bool operator==(const string_iterator& li, const string_iterator& ri) { return li.m_read == ri.m_read && li.m_readholderp == ri.m_readholderp; }
            friend bool operator!=(const string_iterator& li, const string_iterator& ri) { return li.m_read != ri.m_read || li.m_readholderp != ri.m_readholderp; }
            friend class CReadHolder;
        

        private:
            string_iterator(const CReadHolder& rholder, size_t position = 0, size_t read = 0) : m_readholderp(&rholder), m_position(position), m_read(read) {}

            const CReadHolder* m_readholderp;
            size_t m_position;
            size_t m_read;
        };

    private:
        // efficiently copies sequence to destination without converting it to string
        // assumes that destination is extended properly and filled with 0; destination_size - number of 'used' 8-byte words in destination after copy
        template <typename Dest>
        void CopyBits(size_t bit_from, size_t bit_to, Dest& destination, size_t destination_bit_from, size_t destination_size) const {
            if(bit_to <= bit_from)
                return;

            size_t word = bit_from/64;
            size_t last_word = (bit_to-1)/64;
            unsigned shift = bit_from%64;
            size_t destination_word = destination_bit_from/64; 
            unsigned destination_shift = destination_bit_from%64;
            if(shift > 0) {                                                               // first word partial
                uint64_t chunk = (m_storage[word++] >> shift);
                if(destination_shift > 0) {                                               // destination word partial
                    destination[destination_word] += (chunk << destination_shift);
                    if(shift <= destination_shift)                                        // we used all remaining destination word
                        ++destination_word;
                    if(shift < destination_shift && destination_word < destination_size)  // first word spills out
                        destination[destination_word] += (chunk >> (64-destination_shift));
                } else {                                                                  // desination word is not partial - it is bigger than chunk
                    destination[destination_word] = chunk;
                }
                destination_shift = (destination_shift+64-shift)%64;
            }
            for( ; word <= last_word; ++word, ++destination_word) {
                if(destination_shift > 0) {
                    destination[destination_word] += (m_storage[word] << destination_shift);
                    if(destination_word+1 < destination_size)
                        destination[destination_word+1] += (m_storage[word] >> (64-destination_shift));
                } else {
                    destination[destination_word] = m_storage[word];
                }
            }
            int partial_bits = (destination_bit_from+bit_to-bit_from)%64;
            if(partial_bits > 0) {
                uint64_t mask = (uint64_t(1) << partial_bits) - 1;
                destination[destination_size-1] &= mask;
            }
        }


        deque<uint64_t> m_storage;
        deque<uint32_t> m_read_length;
        size_t m_total_seq;
        int m_front_shift;
        bool m_contains_paired;
    };

    typedef list<string> TStrList;
    typedef deque<CDBGraph::Successor> TBases;
    class SContig;
    typedef list<SContig> TContigList;

}; // namespace
#endif /* _DeBruijn_Graph_ */
