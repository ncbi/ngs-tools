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

#ifndef _Concurrent_Hash_
#define _Concurrent_Hash_


#include "Integer.hpp"

// This file contains classes which facilitate basic operation of storing reads, counting kmers,
// and creating and traversing a de Bruijn graph

using namespace std;
namespace DeBruijn {

    template<int BlockSize>
    class CConcurrentBlockedBloomFilter {
    public:
        enum EInsertResult {eNewKmer = 0, eAboveThresholdKmer = 1, eExistingKmer = 2};

        // table_size - number of elements in bloom filter
        // counter_bit_size - number of bith per counter (2, 4, 8...)
        // block_size - in bits (multiple of 64)
        // hash_num - number of has functions (generated from two) 
        CConcurrentBlockedBloomFilter(size_t table_size, int counter_bit_size, int hash_num, size_t min_count) {
            Reset(table_size, counter_bit_size, hash_num, min_count);
        } 
        void Reset(size_t table_size, int counter_bit_size, int hash_num, size_t min_count) {
            CStopWatch timer;
            timer.Restart();
            
            m_count_table.clear();
            m_counter_bit_size = counter_bit_size;
            m_block_size = BlockSize;
            m_hash_num = hash_num;           
            m_max_element = (size_t(1) << m_counter_bit_size) - 1;
            m_min_count = min(min_count, m_max_element);
            m_elements_in_block = BlockSize/m_counter_bit_size;
            m_blocks = ceil((double)table_size/m_elements_in_block);
            m_table_size = m_blocks*m_elements_in_block;
            size_t table_length = m_table_size*m_counter_bit_size/(8*sizeof(TCell));
            m_count_table.resize(table_length, 0);

            cerr << "Resize time: " << timer.Elapsed();
        }
        EInsertResult Insert(size_t hashp, size_t hashm, unsigned min_count) {
            unsigned mc = Test(hashp, hashm);
            if(mc >= min_count)
                return eExistingKmer;

            size_t blk_pos = hashp%m_blocks*m_block_size;
            size_t count = numeric_limits<size_t>::max();
            for(int h = 1; h < m_hash_num; ++h) {
                hashp += hashm;
                size_t pos = blk_pos+(hashp&(m_elements_in_block-1))*m_counter_bit_size;
                auto& cell = m_count_table[pos >> m_bits_in_cell_log].m_atomic;
                int shift = pos&(m_bits_in_cell-1);                
                size_t mask = m_max_element << shift;
                size_t one = size_t(1) << shift;                 
                auto existing_value = cell.load();
                if(((existing_value&mask) >> shift) > mc)
                    continue;
                count = min(count, (existing_value&mask) >> shift);
                while((existing_value&mask) < mask && !cell.compare_exchange_strong(existing_value, existing_value+one)); // check for overflow and try to increment                    
            } 
            if(count == 0)
                return eNewKmer;
            else if(count == m_min_count-1)
                return eAboveThresholdKmer;
            else
                return eExistingKmer;
        }
        int Test(size_t hashp, size_t hashm) const {
            int count = numeric_limits<int>::max();
            size_t blk_pos = hashp%m_blocks*m_block_size;
            for(int h = 1; h < m_hash_num; ++h) {
                hashp += hashm;
                size_t pos = blk_pos+(hashp&(m_elements_in_block-1))*m_counter_bit_size;
                auto& cell = m_count_table[pos >> m_bits_in_cell_log].m_atomic;
                int shift = pos&(m_bits_in_cell-1);
                size_t mask = m_max_element << shift;
                int cn = (cell.load()&mask) >> shift;
                /*!!!!!!!!!!!!!!!REMOVE
                if(cn <= 1)
                    return cn;
                */
                count = min(count, cn);
            }

            return count;
        }
        size_t MaxElement() const { return m_max_element; }
        int HashNum() const { return m_hash_num; }
        size_t TableSize() const { return m_table_size; } // number of counters
        size_t TableFootPrint() const { return sizeof(TCell)*m_count_table.size(); } // bytes

    private:
        size_t m_table_size;
        int m_counter_bit_size;
        int m_block_size;
        int m_hash_num;
        size_t m_blocks;
        size_t m_elements_in_block;
        size_t m_max_element;
        size_t m_min_count;

        size_t m_bits_in_cell = 8*sizeof(TCell);
        int m_bits_in_cell_log = log(m_bits_in_cell)/log(2);
        typedef SAtomic<uint64_t> TCell;
        vector<TCell> m_count_table;
    };
    typedef CConcurrentBlockedBloomFilter<1024> TConcurrentBlockedBloomFilter;

    // minimalistic multithread safe forward list
    // allows reading and inserts in the beginning
    // reading thread will not see new entries inserted after reading started
    template <typename E>
    class CForwardList {
    public:
        struct SNode {
            E m_data;
            SNode* m_next;
        };

        template<bool is_const>
        class iterator : public std::iterator<forward_iterator_tag, E> {
        public:
            template <bool flag, class IsTrue, class IsFalse> struct choose;
            template <class IsTrue, class IsFalse> struct choose<true, IsTrue, IsFalse> {  typedef IsTrue type; };
            template <class IsTrue, class IsFalse> struct choose<false, IsTrue, IsFalse> {  typedef IsFalse type; };
            typedef typename choose<is_const, const E&, E&>::type reference;
            typedef typename choose<is_const, const E*, E*>::type pointer;
            iterator(SNode* node = nullptr) : m_node(node) {};
            iterator& operator++() {
                m_node = m_node->m_next;
                return *this;
            }
            reference& operator*() { return m_node->m_data; }
            pointer operator->() { return &m_node->m_data; }
            bool operator!=(const iterator& i) const { return i.m_node != m_node; }

        private:
            friend class CForwardList;
            SNode* m_node;
        };
        iterator<false> begin() { return iterator<false>(m_head.load()); }
        iterator<false> end() { return iterator<false>(); }
        iterator<true> begin() const { return iterator<true>(m_head.load()); }
        iterator<true> end() const { return iterator<true>(); }

        CForwardList() { m_head.store(nullptr); }
        CForwardList(const CForwardList& other) { Delete(); m_head.store(other.m_head.load()); }
        ~CForwardList() { Delete(); }
        E& front() { return m_head.load()->m_data; }
        SNode* Head() { return m_head; }
        SNode* NewNode(const E& e) {
            SNode* p = new SNode;
            p->m_data = e;
            p->m_next = m_head;
            return p;
        }
        SNode* NewNode() {
            SNode* p = new SNode;
            p->m_next = m_head;
            return p;
        }
        bool TryPushFront(SNode* nodep) { return m_head.compare_exchange_strong(nodep->m_next, nodep); }
        E* Emplace() {
            SNode* p = NewNode();
            while (!m_head.compare_exchange_strong(p->m_next, p));
            return &(p->m_data);
        }
        void PushFront(const E& e) {
            SNode* p = NewNode(e);
            while (!m_head.compare_exchange_strong(p->m_next, p));
        }
        // not mutithread safe
        template <typename P> void remove_if(const P& prob) {
            while(m_head.load() != nullptr && prob(m_head.load()->m_data)) {
                SNode* p = m_head;
                m_head = p->m_next;
                delete p;
            }
            for(SNode* p = m_head; p != nullptr; ) {
                SNode* after = p->m_next;
                if(after != nullptr && prob(after->m_data)) {
                    p->m_next = after->m_next;
                    delete after;
                } else {
                    p = p->m_next;
                }
            }
        }

    private:
        void Delete() {
            for(SNode* p = m_head; p != nullptr; ) {
                auto tmp = p->m_next;
                delete p;
                p = tmp;
            }
            m_head.store(nullptr);
        }
        atomic<SNode*> m_head;
    };

    // minimalistic deque-type container allowing multithread initialization of a large hash_table
    template <typename E>
    class CDeque {
    public:
        typedef E value_type;
        CDeque(size_t size = 0) : m_chunks(1) {
            m_data.resize(m_chunks);
            Reset(size, m_chunks);
        }
        E& operator[](size_t index) { return m_data[index/m_chunk_size][index%m_chunk_size]; }        
        const E& operator[](size_t index) const { return m_data[index/m_chunk_size][index%m_chunk_size]; }        
        size_t size() const { return m_size; }
        void Reset(size_t size, size_t chunks) {
            m_chunks = chunks;
            m_data.resize(m_chunks);
            m_size = size;
            m_chunk_size = (size+m_chunks-1)/m_chunks;
            if(m_chunks == 1) {
                ResetChunk(0, m_chunk_size);
            } else {
                list<function<void()>> jobs;
                for(unsigned chunk = 0; chunk < m_chunks; ++chunk) {
                    size_t chunk_size = min(m_chunk_size, size);
                    jobs.push_back(bind(&CDeque::ResetChunk, this, chunk, chunk_size));
                    size -= min(m_chunk_size, size);                    
                }
                RunThreads(m_chunks, jobs);                            
            }
        }
    private:
        void ResetChunk(size_t chunk, size_t chunk_size) {
            m_data[chunk].clear();
            m_data[chunk].resize(chunk_size);
        }
        size_t m_chunks = 0;
        size_t m_size = 0;
        size_t m_chunk_size = 0;
        vector<vector<E>> m_data;
    };

    // Moderate value of BucketBlock will improve memory cache use
    // Larger values will reduce the number of entries in the spillover lists but eventually will increase the search time
    // 0 (all entries in the lists) is permitted and could be used for reduction of the table size for large sizeof(V)
    #define StatusNum 8
    template<int N, class V, int BucketBlock>
    struct SHashBlock {
        static_assert(StatusNum >= BucketBlock, "");
        typedef LargeInt<N> large_t;
        typedef V mapped_t;
        typedef pair<large_t,V> element_t;
        typedef CForwardList<element_t> list_t;
        typedef uint8_t status_t;
        enum {eAssigned = 1, eKeyExists = 2};
        SHashBlock() {
            for(auto& status : m_status)
                status = 0;
        }
        SHashBlock(const SHashBlock& other) : m_data(other.m_data) { // used for table initialisation only
            for(unsigned i = 0; i < m_status.size(); ++i)
                m_status[i] = other.m_status[i].load();
        }
        
        array<element_t, BucketBlock> m_data;
        list_t m_extra;

        bool Lock(int shift, const large_t& kmer) {
            status_t expected = 0;
            if(!m_status[shift].compare_exchange_strong(expected, eAssigned))
                return false;

            m_data[shift].first = kmer;
            m_status[shift] |= eKeyExists;
            return true;
        }
        void Wait(int shift) { while(!(m_status[shift].load()&eKeyExists)); }
        bool isEmpty(int shift) const { return m_status[shift].load() == 0; }

        void Move(element_t& cell, int to) {
            m_data[to] = cell;
            m_status[to] = (eAssigned|eKeyExists);
            cell.second = V();
        }
        void Move(int from, int to) {
            Move(m_data[from], to);
            m_status[from] = 0;
        }
        void Clear(int shift) {
            m_data[shift].second = V();
            m_status[shift] = 0;
        }

    private:
        array<atomic<status_t>,StatusNum>  m_status;
    };

    template <typename MappedV, int BucketBlock> 
    class CKmerHashMap {
    public:
        CKmerHashMap(int kmer_len, size_t size = 0) : m_kmer_len(kmer_len) {
            m_hash_table = CreateVariant<TKmerHashTable<MappedV>, THashBlockVec, MappedV>((m_kmer_len+31)/32);
            Reset(size, 1);
        }
        void Reset(size_t size, size_t chunks) {
            size_t blocks = size/max(1,BucketBlock);
            if(size%max(1,BucketBlock))
                ++blocks;
            m_table_size = max(1,BucketBlock)*blocks;
            apply_visitor(resize(blocks, chunks), m_hash_table);
        }
        //returns pointer to mapped value if exists, otherwise nullptr
        const MappedV* Find(const TKmer& kmer) const { return apply_visitor(find(kmer), m_hash_table); }        
        // if kmer already included returns pointer to mapped value
        // if not inserts a new entry and returns pointer to default value
        // caller MUST update the mapped value
        // assumes that any updates will be atomic
        MappedV* FindOrInsertInBucket(const TKmer& kmer, size_t index) {return apply_visitor(find_or_insert(kmer,index), m_hash_table); }
        int KmerLen() const { return m_kmer_len; }
        size_t TableSize() const { return m_table_size; }
        size_t TableFootPrint() const { return apply_visitor(foot_print(), m_hash_table); }
        size_t BucketsNum() const { return m_table_size/max(1,BucketBlock); }
        void Info() const { apply_visitor(info(), m_hash_table); }

    protected:
        template<int N, class V> using THashBlockVec = CDeque<SHashBlock<N,V,BucketBlock>>;
        template<class V> using TKmerHashTable = BoostVariant<THashBlockVec,V>;

        struct foot_print : public boost::static_visitor<size_t> {
            template <typename T> const size_t operator()(T& v) const { return sizeof(typename  T::value_type)*v.size(); }
        };
        
        //returns pointer to mapped value if exists, otherwise nullptr
        struct find : public boost::static_visitor<const MappedV*> {
            find(const TKmer& k) : kmer(k) {}
            template <typename T> const MappedV* operator()(T& v) const {
                typedef typename  T::value_type::large_t large_t;
                const large_t& k = kmer.get<large_t>();
                size_t index = k.oahash()%(v.size()*max(1,BucketBlock));
                auto& bucket = v[index/max(1,BucketBlock)];

                //try exact position first
                if(BucketBlock > 0) {
                    unsigned exact_pos = index%BucketBlock;
                    auto& cell = bucket.m_data[exact_pos];                
                    if(bucket.isEmpty(exact_pos)) 
                        return nullptr;
                    else if(cell.first == k)
                        return &cell.second;

                    //scan array    
                    for(unsigned shift = 0; shift < bucket.m_data.size(); ++shift) {
                        if(shift != exact_pos) {
                            auto& cell = bucket.m_data[shift]; 
                            if(bucket.isEmpty(shift)) 
                                return nullptr;
                            else if(cell.first == k)
                                return &cell.second;
                        }
                    }       
                }

                //scan spillover list
                for(auto& cell : bucket.m_extra) {
                    if(cell.first == k)
                        return &cell.second;
                }

                return nullptr;            
            } 
            const TKmer& kmer;
        };        

        template <typename T> 
        static MappedV* find_or_insert_in_bucket(T& v, size_t index, const typename T::value_type::large_t& k) {
            typedef typename T::value_type hash_block_t;
            typedef typename hash_block_t::mapped_t mapped_t;
            typedef typename hash_block_t::list_t list_t;
            typedef typename hash_block_t::large_t large_t;

            auto& bucket = v[index/max(1,BucketBlock)];

            auto TryCell = [&](int shift) {
                auto& cell = bucket.m_data[shift];

                //try to grab
                if(bucket.Lock(shift, k))
                    return &cell.second;
                
                //already assigned to some kmer
                //wait if kmer is not stored yet
                bucket.Wait(shift);
               
                if(cell.first == k) // kmer matches
                    return &cell.second; 
                else
                    return (mapped_t*)nullptr; // other kmer   
            };

            if(BucketBlock > 0) {
                //try exact position first 
                unsigned exact_pos = index%BucketBlock;
                auto rslt = TryCell(exact_pos);
                if(rslt != nullptr)
                    return rslt;            

                //scan remaining array  
                for(unsigned shift = 0; shift < bucket.m_data.size(); ++shift) {
                    if(shift != exact_pos) {
                        auto rslt = TryCell(shift);
                        if(rslt != nullptr)
                            return rslt;
                    }
                }
            }

            //scan spillover list   
            auto existing_head = bucket.m_extra.Head();
            for(auto p = existing_head; p != nullptr; p = p->m_next) {
                if(p->m_data.first == k)
                    return &(p->m_data.second); 
            }

            typename list_t::SNode* nodep = new typename list_t::SNode;
            nodep->m_data.first = k;
            nodep->m_next = existing_head;
            while(!bucket.m_extra.TryPushFront(nodep)) {
                //check if a new elemet matches 
                for(auto p = nodep->m_next; p != existing_head; p = p->m_next) {
                    if(p->m_data.first == k) {
                        delete nodep;
                        return &(p->m_data.second); 
                    }
                }
                existing_head = nodep->m_next;
            }                

            return &(nodep->m_data.second);
        }

        // if kmer already included returns pointer to mapped value
        // if not inserts a new entry and returns pointer to default value
        // caller MUST update the mapped value
        // assumes that any updated will be atomic
        struct find_or_insert : public boost::static_visitor<MappedV*> {
            find_or_insert(const TKmer& k, size_t i) : kmer(k), index(i) {}
            template <typename T> MappedV* operator()(T& v) const {
                typedef typename T::value_type::large_t large_t;
                const large_t& k = kmer.get<large_t>();
                return find_or_insert_in_bucket(v, index, k);
            }
            const TKmer& kmer;
            size_t index;
        };

        struct info : public boost::static_visitor<> {
            template <typename T> void operator()(T& v) const {
                
                map<int,int> numbers;
                for(size_t i = 0; i < v.size(); ++i) {
                    auto& bucket= v[i];
                    int num = distance(bucket.m_extra.begin(), bucket.m_extra.end());  
                    for(auto& cell : bucket.m_data) {
                        if(!(cell.second == MappedV()))
                            ++num;
                    }                                        
                    ++numbers[num];
                }
                for(auto& rslt : numbers)
                    cerr << "Bucket:\t" << rslt.first << "\t" << rslt.second << endl;
            }            
        };

        struct resize : public boost::static_visitor<> { 
            resize(size_t s, size_t c) : size(s), chunks(c) {}
            template <typename T> void operator()(T& v) const { v.Reset(size, chunks); }            
            size_t size;
            size_t chunks;
        };

        TKmerHashTable<MappedV> m_hash_table;
        size_t m_table_size;
        int m_kmer_len;
    };

    struct SKmerCounter {
        SKmerCounter() : m_count(0) {}
        bool operator==(const SKmerCounter& kc) const { return kc.m_count == m_count; }
        uint32_t Increment(bool is_plus) { return (m_count.m_atomic += (is_plus ? 0x100000001 : 1)); }
        uint32_t Count() const { return m_count; } // clips plus part
        SAtomic<uint64_t> m_count;
    };
    
    class CKmerHashCount : public CKmerHashMap<SKmerCounter, 8> {
    public:
        CKmerHashCount(int kmer_len, size_t size = 0) : CKmerHashMap(kmer_len, size) {}
        // returns true if kmer was new
        bool UpdateCount(const TKmer& kmer, bool is_plus) {
            size_t index = kmer.oahash()%m_table_size;
            return (FindOrInsertInBucket(kmer, index)->Increment(is_plus) == 1);
        } 
        size_t UpdateCounts(const CReadHolder::string_iterator& is, const TConcurrentBlockedBloomFilter& bloom, int min_count) {
            return apply_visitor(update_counts(is, bloom, TableSize(), min_count, m_kmer_len), m_hash_table);
        }
        // rehash bucket from other container
        void RehashOtherBuckets(CKmerHashCount& other, size_t bucket_from, size_t bucket_to) {
            apply_visitor(rehash_bucket(bucket_from, bucket_to, *this), m_hash_table, other.m_hash_table); 
        }
        //remove false positives
        size_t CleanBuckets(int min_count, size_t bucket_from, size_t bucket_to) {
            return apply_visitor(clean_buckets(min_count, bucket_from, bucket_to, TableSize()), m_hash_table);
        }
    private:
        struct update_counts : public boost::static_visitor<size_t> {
            update_counts(const CReadHolder::string_iterator& i, const TConcurrentBlockedBloomFilter& bl, size_t ts, int mc, unsigned kl) : is(i), bloom(bl), table_size(ts), min_count(mc), kmer_len(kl) {}
            template <typename T> size_t operator() (T& v) const {
                typedef typename T::value_type::large_t large_t;

                size_t read_len = is.ReadLen();
                if(read_len < kmer_len)
                    return 0;

                unsigned kmer_bytes = (2*kmer_len+7)/8;                 //number of whole bytes in kmer
                unsigned kmer_size = (2*kmer_len+63)/64;                //number of whole 8-byte words in kmer
                int partial_bits = (2*kmer_len)%64;                     //number of used bits in partial 8 byte word (if any)
                uint64_t mask = numeric_limits<uint64_t>::max();
                if(partial_bits > 0)
                    mask = (uint64_t(1) << partial_bits) - 1;
                size_t buf_size = (2*read_len+63)/64+1;
                uint64_t* read_buf = new uint64_t[buf_size]; //(enough + 1) 8-byte words for read (one extra because we'll copy kmers using whole bytes which can go beyond the sequence)

                size_t new_kmers = 0;
                large_t kmer(0);
                for(int shift = 0; shift < 4 && read_len-shift >= kmer_len; ++shift) {
                    memset(read_buf, 0, 8*buf_size);
                    is.BSeq(shift, read_buf);
                    for(unsigned k = 0; k <= read_len-shift-kmer_len; k += 4) { // every 4th kmer on the byte boundary
                        memcpy(kmer.getPointer(), (uint8_t*)read_buf+k/4, kmer_bytes);
                        kmer.getPointer()[kmer_size-1] &= mask;
                        large_t rkmer = revcomp(kmer, kmer_len);
                        large_t* min_kmerp = &rkmer; 
                        bool is_plus = false;
                        size_t hashp = rkmer.oahash();
                        size_t hashm = kmer.oahash();
                        if(kmer < rkmer) {
                            min_kmerp = &kmer;
                            is_plus = true;
                            swap(hashp, hashm);
                        }

                        if(min_count > 1 && bloom.Test(hashp, hashm) < min(min_count, (int)bloom.MaxElement()))
                            continue;  

                        size_t index = hashp%table_size;
                        if(find_or_insert_in_bucket(v, index, *min_kmerp)->Increment(is_plus) == 1)
                            ++new_kmers;
                    }
                }

                delete[] read_buf;
                return new_kmers;
            }
            const CReadHolder::string_iterator& is;
            const TConcurrentBlockedBloomFilter& bloom;
            size_t table_size;
            int min_count;
            unsigned kmer_len;
        };

        struct rehash_bucket : public boost::static_visitor<> {
            rehash_bucket(size_t bf, size_t bt, CKmerHashCount& h) : bucket_from(bf), bucket_to(bt), hash(h) {}
            template <typename T> void operator() (T& a, T& b) const {
                typedef typename T::value_type::element_t element_t;
                for(size_t indexb = bucket_from; indexb <= bucket_to; ++indexb) {
                    auto& bucket_b = b[indexb];

                    for(auto& cell : bucket_b.m_data) {
                        if(cell.second.Count() != 0) {
                            auto& kmer = cell.first;
                            size_t indexa = kmer.oahash()%hash.TableSize();
                            *hash.FindOrInsertInBucket(TKmer(kmer), indexa) = cell.second;
                        }
                    }
                    
                    for(element_t& cell : bucket_b.m_extra) {
                        auto& kmer = cell.first;
                        size_t indexa = kmer.oahash()%hash.TableSize();
                        *hash.FindOrInsertInBucket(TKmer(kmer), indexa) = cell.second;
                    }
                }
            }
            template <typename T, typename U> void operator() (T& a, U& b) const { throw runtime_error("Can't rehash from different type container"); }
            
            size_t bucket_from;
            size_t bucket_to;
            CKmerHashCount& hash;
        };

        struct clean_buckets : public boost::static_visitor<size_t> {
            clean_buckets(int mc, size_t bf, size_t bt, size_t tb) : min_count(mc), bucket_from(bf), bucket_to(bt), table_size(tb) {}            
            template <typename T> size_t operator()(T& v) const {
                typedef typename T::value_type::element_t element_t;
                typedef typename T::value_type::mapped_t mapped_t;

                size_t num = 0;
                for(size_t bind = bucket_from; bind <= bucket_to; ++bind) {
                    auto& bucket = v[bind];
                    int empty_cells = 0;

                    auto Reposition = [&bucket, this](element_t& cell, unsigned limit) {
                        size_t index = cell.first.oahash()%table_size;
                        size_t orig_shift = index%bucket.m_data.size();
                        if(orig_shift == limit)
                            return false;

                        if(bucket.m_data[orig_shift].second.Count() < min_count) {
                            if(limit < bucket.m_data.size())
                                bucket.Move(limit, orig_shift);
                            else
                                bucket.Move(cell, orig_shift);

                            return orig_shift > limit;
                        } 

                        for(unsigned shift = 0; shift < limit; ++shift) {
                            if(shift != orig_shift && bucket.m_data[shift].second.Count() < min_count) {
                                if(limit < bucket.m_data.size())
                                    bucket.Move(limit, shift);
                                else
                                    bucket.Move(cell, shift);

                                return false;
                            }                         
                        }

                        return false;
                    };

                                      
                    for(unsigned shift = 0; shift < bucket.m_data.size(); ++shift) {
                        auto& cell = bucket.m_data[shift];
                        auto count = cell.second.Count();
                        if(count < min_count) {
                            ++empty_cells;                            
                            if(count > 0)
                                bucket.Clear(shift);
                        } else {
                            if(Reposition(cell, shift))
                                ++empty_cells;          // moved down and created new empty cell (will be counted later)
                            else
                                ++num;                  // stayed or moved up
                        }
                    }
 
                    for(auto& cell : bucket.m_extra) {
                        if(cell.second.Count() >= min_count) {
                            ++num;
                            if(empty_cells > 0) {
                                Reposition(cell, bucket.m_data.size());
                                --empty_cells;
                            }
                        }
                    }
                                         
                    bucket.m_extra.remove_if([this](const element_t& elem) {return elem.second.Count() < min_count;});
                 }

                return num;
            }            
            unsigned min_count;
            size_t bucket_from;
            size_t bucket_to;
            size_t table_size;
        };
    };


    class CKmerHashCounter {
    public:
        CKmerHashCounter(const list<array<CReadHolder,2>>& reads, int kmer_len, int min_count, size_t estimated_kmer_num, bool is_stranded, int ncores, bool skip_bloom) : 
            m_kmer_len(kmer_len), m_min_count(min_count), m_is_stranded(is_stranded), m_ncores(ncores), m_skip_bloom(skip_bloom), m_hash_table(m_kmer_len), 
            m_estimated_table_size(0),  m_estimated_uniq_kmers(0), m_kmer_num(0), m_kmer_num_raw(0), m_kmer_count(0), m_rehash_status(false) {

            m_kmer_step = max(1., 0.1*m_hash_table.TableSize()/m_ncores);
            for(auto& rholder : reads) 
                m_start_position.push_back(make_pair(0, rholder[0].sbegin()));
            
            CStopWatch timer;

            TConcurrentBlockedBloomFilter bloom(0, 2, 1, m_min_count);
            if(!m_skip_bloom)  {
                m_estimated_uniq_kmers.store(estimated_kmer_num);
                while(true) {
                    timer.Restart();
                    int counter_bit_size = 2;
                    for( ; counter_bit_size <= 8 &&  (1 << counter_bit_size)-1 < m_min_count; counter_bit_size *= 2);
                    double false_positive_rate = 0.03;
                    //                    size_t bloom_table_size = -1.5*(double)m_estimated_uniq_kmers.load()*log(false_positive_rate)/log(2.)/log(2.); // 50% extra because blocked
                    size_t bloom_table_size = -(double)m_estimated_uniq_kmers.load()*log(false_positive_rate)/log(2.)/log(2.); // 50% extra because blocked
                    int hash_num = ceil(-log(false_positive_rate)/log(2.));
                    bloom.Reset(bloom_table_size, counter_bit_size, hash_num, m_min_count);

                    cerr << "Bloom table size: " << bloom.TableSize() << "(" << 0.1*(bloom.TableFootPrint()/100000) << "MB)" << " Counter bit size: " << counter_bit_size << " Hash num: " << hash_num << endl;  

                    m_estimated_table_size.store(0);
                    m_estimated_uniq_kmers.store(0);

                    list<function<void()>> jobs;
                    for(auto& job_input : reads) {
                        if(job_input[0].ReadNum() > 0 || job_input[1].ReadNum() > 0) {   // not empty                   
                            jobs.push_back(bind(&CKmerHashCounter::InsertInBloomJob, this, ref(job_input), ref(bloom)));
                        }
                    }
                    RunThreads(m_ncores, jobs);            
            
                    if(m_min_count == 1)
                        m_estimated_table_size.store(m_estimated_uniq_kmers.load());
                    double kmers = m_estimated_uniq_kmers.load();
                    false_positive_rate = pow(1.-exp(-hash_num*kmers/bloom_table_size), hash_num);
                    cerr << "Estimated kmers above threshold: " << m_estimated_table_size.load() << " Estimated uniq kmers: " << m_estimated_uniq_kmers.load() << " Estimated bloom false positive rate " << false_positive_rate << endl;
                    cerr << "Bloom filter in " << timer.Elapsed();
                    
                    if(false_positive_rate < 0.15)
                        break;

                    cerr << "\nBloom filter false positive rate is too high - increasing the bloom filter size and recalculating" << endl;
                }
            } else {
                m_estimated_table_size.store(estimated_kmer_num);
            }

            timer.Restart();
            m_hash_table.Reset(1.5*m_estimated_table_size.load(), m_ncores);            
            while(true) {
                {
                    CStopWatch timer;
                    timer.Restart();
                    list<function<void()>> jobs;
                    auto start_pos = m_start_position.begin();
                    for(auto& job_input : reads) {
                        if(job_input[0].ReadNum() > 0 || job_input[1].ReadNum() > 0) {   // not empty               
                            jobs.push_back(bind(&CKmerHashCounter::CountKmersJob, this, ref(job_input), ref(*start_pos), ref(bloom)));
                        }
                        ++start_pos;
                    }
                    RunThreads(m_ncores, jobs);
               }
                
                if(!m_rehash_status.load())
                    break;

                //Rehash
                {
                    CStopWatch timer;
                    timer.Restart();
                    size_t new_size = m_hash_table.TableSize()*m_increase_factor;
                    cerr << "Rehash new size: " << new_size << endl;
                    CKmerHashCount hash_table_tmp(m_kmer_len);
                    hash_table_tmp.Reset(new_size, m_ncores);
                    swap(m_hash_table, hash_table_tmp);
                    m_kmer_step = max(1., 0.1*m_hash_table.TableSize()/m_ncores);
                    m_rehash_status.store(false);
                
                    list<function<void()>> jobs;
                    size_t step = ceil((double)hash_table_tmp.BucketsNum()/m_ncores);
                    for(int thr = 0; thr < m_ncores; ++thr) {
                        size_t from = step*thr;
                        size_t to = min(hash_table_tmp.BucketsNum()-1,from+step-1);
                        if(to >= from)
                            jobs.push_back(bind(&CKmerHashCounter::RehashJob, this, ref(hash_table_tmp), from, to));
                    }
                    RunThreads(m_ncores, jobs);
                    cerr << "Rehashing in " << timer.Elapsed();
                }
            }

            cerr << "Create hash in " << timer.Elapsed();
            timer.Restart();            
                        
            //Remove false positives
            {
                list<function<void()>> jobs;
                size_t step = ceil((double)m_hash_table.BucketsNum()/m_ncores);
                for(int thr = 0; thr < m_ncores; ++thr) {
                    size_t from = step*thr;
                    size_t to = min(m_hash_table.BucketsNum()-1,from+step-1);
                    if(to >= from)
                        jobs.push_back(bind(&CKmerHashCounter::CleanJob, this, from, to));
                }
                RunThreads(m_ncores, jobs);
            }                       

            cerr << "Clean hash in " << timer.Elapsed();
            timer.Restart();            

            cerr << "Initial kmers: " << m_kmer_num_raw.load() << " Kmers above threshold: " << m_kmer_num.load() << " Total kmers: " << m_kmer_count.load() << " Hash table size: " <<  m_hash_table.TableSize() << "(" << 0.1*(m_hash_table.TableFootPrint()/100000) << "MB)" << endl;

        }
        void Info() const {
            m_hash_table.Info();
        }
        CKmerHashCount& HasCount() { return m_hash_table; }
        

    private:
        
        void CleanJob(size_t bucket_from, size_t bucket_to) {
             m_kmer_num += m_hash_table.CleanBuckets(m_min_count, bucket_from, bucket_to);
        }

        class CBloomInserter : public TKmer {
        public:
            CBloomInserter(int kmer_len) : TKmer(kmer_len, 0), m_kmer_len(kmer_len) {}
            pair<size_t,size_t> InsertInBloom(const CReadHolder::string_iterator& is, TConcurrentBlockedBloomFilter& bloom, unsigned min_count) {
                return apply_visitor(insert_in_bloom(is, bloom, m_kmer_len, min_count), v);
            }

        private:
            unsigned m_kmer_len;

            struct insert_in_bloom : public boost::static_visitor<pair<size_t,size_t>> {
                insert_in_bloom(const CReadHolder::string_iterator& i, TConcurrentBlockedBloomFilter& bl, unsigned kl, unsigned mc) : is(i), bloom(bl), kmer_len(kl), min_count(mc) {}
                template <typename large_t> pair<size_t,size_t> operator() (large_t& kmer) const {
                    size_t above_threshold_kmers = 0;
                    size_t uniq_kmers = 0;

                    size_t read_len = is.ReadLen();
                    if(read_len < kmer_len)
                        return make_pair(above_threshold_kmers, uniq_kmers);

                    unsigned kmer_bytes = (2*kmer_len+7)/8;                 //number of whole bytes in kmer
                    unsigned kmer_size = (2*kmer_len+63)/64;                //number of whole 8-byte words in kmer
                    int partial_bits = (2*kmer_len)%64;                     //number of used bits in partial 8 byte word (if any)
                    uint64_t mask = numeric_limits<uint64_t>::max();
                    if(partial_bits > 0)
                        mask = (uint64_t(1) << partial_bits) - 1;
                    size_t buf_size = (2*read_len+63)/64+1;
                    uint64_t* read_buf = new uint64_t[buf_size]; //(enough + 1) 8-byte words for read (one extra because we'll copy kmers using whole bytes which can go beyond the sequence)

                    for(int shift = 0; shift < 4 && read_len-shift >= kmer_len; ++shift) {
                        memset(read_buf, 0, 8*buf_size);
                        is.BSeq(shift, read_buf);
                        for(unsigned k = 0; k <= read_len-shift-kmer_len; k += 4) { // every 4th kmer on the byte boundary
                            memcpy(kmer.getPointer(), (uint8_t*)read_buf+k/4, kmer_bytes);
                            kmer.getPointer()[kmer_size-1] &= mask;
                            large_t rkmer = revcomp(kmer, kmer_len);

                            size_t hashp = rkmer.oahash();
                            size_t hashm = kmer.oahash();
                            if(kmer < rkmer)
                                swap(hashp, hashm);  

                            switch(bloom.Insert(hashp, hashm, min_count)) {
                            case TConcurrentBlockedBloomFilter::eNewKmer : ++uniq_kmers; continue;
                            case TConcurrentBlockedBloomFilter::eAboveThresholdKmer : ++above_threshold_kmers; continue;
                            default : continue;
                            }
                        }
                    }
                    delete[] read_buf;
                    return make_pair(above_threshold_kmers, uniq_kmers);
                }                
                const CReadHolder::string_iterator& is;
                TConcurrentBlockedBloomFilter& bloom;
                unsigned kmer_len;
                unsigned min_count;
            };
        };

        void InsertInBloomJob(const array<CReadHolder,2>& rholder, TConcurrentBlockedBloomFilter& bloom) {
            size_t above_threshold_kmers = 0;
            size_t uniq_kmers = 0;
 
            CBloomInserter bloom_inserter(m_kmer_len);
            for(int p = 0; p < 2; ++p) {
               for(CReadHolder::string_iterator is = rholder[p].sbegin(); is != rholder[p].send(); ++is) {
                   auto rslt = bloom_inserter.InsertInBloom(is, bloom, m_min_count);
                   above_threshold_kmers += rslt.first;
                   uniq_kmers += rslt.second;
               }
            }

            m_estimated_table_size += above_threshold_kmers;
            m_estimated_uniq_kmers += uniq_kmers;
        }
        void RehashJob(CKmerHashCount& other_hash_table, size_t bucket_from, size_t bucket_to) {
            m_hash_table.RehashOtherBuckets(other_hash_table, bucket_from, bucket_to);
        }

        void CountKmersJob(const array<CReadHolder,2>& rholder, pair<int,CReadHolder::string_iterator>& start_pos, const TConcurrentBlockedBloomFilter& bloom) {
            size_t kmer_num = 0;
            size_t kmer_count = 0;
            for(int p = start_pos.first; p < 2; ++p) {
                CReadHolder::string_iterator from = start_pos.second;
                if(p != start_pos.first)
                    from = rholder[p].sbegin();
                for(CReadHolder::string_iterator is = from; is != rholder[p].send(); ++is) {
                    size_t read_len = is.ReadLen();
                    if(read_len >= (unsigned)m_kmer_len)
                        kmer_count += read_len-m_kmer_len+1;
                    else
                        continue;

                    kmer_num += m_hash_table.UpdateCounts(is, bloom, m_skip_bloom ? 0 : m_min_count);

                    if(kmer_num >= m_kmer_step) {
                        m_kmer_num_raw += kmer_num;
                        m_kmer_count += kmer_count;
                        kmer_num = 0;
                        kmer_count = 0;
                        if(m_kmer_num_raw.load() > m_hash_table.TableSize()*m_max_load_factor)
                            m_rehash_status.store(true); 
                        if(m_rehash_status.load()) {
                            start_pos.first = p;
                            ++is;
                            start_pos.second = is;
                            return;
                        }
                    }                    
                }
            }
            m_kmer_num_raw += kmer_num;            
            m_kmer_count += kmer_count;
            start_pos.first = 2;
        }

        int m_kmer_len;
        int m_min_count;
        bool m_is_stranded;
        int m_ncores;
        bool m_skip_bloom;

        CKmerHashCount m_hash_table;
        atomic<size_t> m_estimated_table_size;
        atomic<size_t> m_estimated_uniq_kmers;
        atomic<size_t> m_kmer_num;
        atomic<size_t> m_kmer_num_raw;
        atomic<size_t> m_kmer_count;
        atomic<bool> m_rehash_status;
        size_t m_kmer_step;
        double m_max_load_factor = 1;
        int m_increase_factor = 2;
        list<pair<int,CReadHolder::string_iterator>> m_start_position;
    };


    /* TODO new graph based on hash
    struct SGraphData {
        SGraphData() : m_count(0), m_plus_fraction(0), m_branches(0), m_status(0) {}
        bool operator==(const SGraphData& kc) const { 
            return kc.m_count == m_count && 
                kc.m_plus_fraction == m_plus_fraction &&
                kc.m_branches == m_branches;
        }
        
        uint32_t m_count;
        uint16_t m_plus_fraction;
        uint8_t m_branches;
        SAtomic<uint8_t> m_status;
    };
    class CDBHashGraph : public CKmerHashMap<SGraphData> {
    public:
        CDBHashGraph(const CKmerHashCount& counter) : CKmerHashMap(counter.KmerLen(), counter.TableSize()) {

            //            const TKmerHashTable<SKmerCounter>& other_table = counter.m_hash_table;

            string max_kmer(KmerLen(), bin2NT[3]);
            m_max_kmer = TKmer(max_kmer);
        }

    private:
        TKmer m_max_kmer;             // contains 1 in all kmer_len bit positions  
        TBins m_bins;
        bool m_is_stranded;
    };
    */

}; // namespace
#endif /*_Concurrent_Hash_*/
