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

    class CCentinel {
    public:
        CCentinel(size_t size) { m_centinel.resize(size, 0); }
        void GrabBucket(size_t index) {            
            while(!m_centinel[index].Set(1, 0));
        }
        void ReleaseBucket(size_t index) { 
            m_centinel[index] = 0;
        }
    protected:
        vector<SAtomic<uint8_t>> m_centinel;
    };

    //hash CANNOT CONTAIN default value of MappedValue
    template <typename MappedValue, typename Inserter> class CKmerHashMapSimple : public CCentinel {
    public:
        CKmerHashMapSimple(int kmer_len, size_t size) : CCentinel(size), m_table_size(size), m_kmer_len(kmer_len) {
            m_hash_table = CreateVariant<TKmerAndV<MappedValue>, TLargeIntVecV, MappedValue>((m_kmer_len+31)/32);
            apply_visitor(resize(m_table_size), m_hash_table);
        }
        int KmerLen() const { return m_kmer_len; }
        size_t TableSize() const { return m_table_size; }
        // returns false if table is full
        bool Insert(const TKmer& kmer, const Inserter& inserter) { return apply_visitor(insert(kmer, inserter, this), m_hash_table); }
        MappedValue* Find(const TKmer& kmer) { return apply_visitor(find(kmer, this), m_hash_table); }
    
    private:
        //returns pointer to mapped value if exists, otherwise nullptr
        struct find : public boost::static_visitor<MappedValue*> {
            find(const TKmer& k, CKmerHashMapSimple* p) : kmer(k), callerp(p) {}
            template <typename T> MappedValue* operator()(T& v) const {
                typedef typename  T::value_type::first_type large_t;
                const large_t& k = kmer.get<large_t>();
                size_t index = kmer.oahash()%callerp->m_table_size;
                for(size_t step = 0; step < callerp->m_table_size; ++step, index = (index+1)%callerp->m_table_size) {
                    auto& bucket = v[index];
                    if(bucket.second == MappedValue())   // empty
                        return nullptr;
                    else if(bucket.first == k)
                        return &bucket.second;
                }
                return nullptr;
            } 
            const TKmer& kmer;
            CKmerHashMapSimple* callerp;
        };
        struct insert : public boost::static_visitor<bool> {
            insert(const TKmer& k, const Inserter& i, CKmerHashMapSimple* p) : kmer(k), inserter(i), callerp(p) {}
            template <typename T> bool operator()(T& v) const {
                typedef typename  T::value_type::first_type large_t;
                const large_t& k = kmer.get<large_t>();
                size_t index = kmer.oahash()%callerp->m_table_size;
                for(size_t step = 0; step < callerp->m_table_size; ++step, index = (index+1)%callerp->m_table_size) {
                    auto& bucket = v[index];
                    callerp->GrabBucket(index);
                    if(bucket.second == MappedValue()) {  // empty
                        bucket.first = k;
                        inserter(bucket.second);
                        callerp->ReleaseBucket(index);
                        return true;
                    } else if(bucket.first == k) {        // kmer match
                        inserter(bucket.second);
                        callerp->ReleaseBucket(index);
                        return true;
                    }                        
                    callerp->ReleaseBucket(index);
                }

                return false;
            }

            const TKmer& kmer;
            const Inserter& inserter;
            CKmerHashMapSimple* callerp;
        };
        struct resize : public boost::static_visitor<> { 
            resize(size_t s) : size(s) {}
            template <typename T> void operator()(T& v) const { v.resize(size); }            
            size_t size;
        };

        TKmerAndV<MappedValue> m_hash_table;
        size_t m_table_size;
        int m_kmer_len;
    };

    //hash CANNOT CONTAIN default value of V
    template <typename V> class CKmerHashMap {
    public:
        CKmerHashMap(int kmer_len, size_t size) : m_table_size(size), m_kmer_len(kmer_len) {
            m_hash_table = CreateVariant<TKmerHashTable<V>, TOneWayListVec, V>((m_kmer_len+31)/32);
            apply_visitor(resize(m_table_size), m_hash_table);
        }
        //OK to use in multiple threads if noone changes the hash
        const V* Find(const TKmer& kmer) const { return apply_visitor(find(kmer), m_hash_table); }
        // MUST be sandwiched between GrabBucket/ReleaseBucket
        V* FindOrInsertInBucket(const TKmer& kmer, size_t index) {return apply_visitor(find_or_insert(kmer,index), m_hash_table); }
        int KmerLen() const { return m_kmer_len; }
        size_t TableSize() const { return m_table_size; }

    protected:
        // if kmer already included returns pointer to mapped value
        // if not inserts a new entry and returns pointer to default value
        // caller MUST update the mapped value
        struct find_or_insert : public boost::static_visitor<V*> {
            find_or_insert(const TKmer& k, size_t i) : kmer(k), index(i) {}
            template <typename T> V* operator()(T& v) const {
                typedef typename  T::value_type::large_t large_t;
                const large_t& k = kmer.get<large_t>();
                auto& bucket = v[index];
                if(bucket.m_data.second == V()) {
                    bucket.m_data.first = k;
                    return &bucket.m_data.second;
                } else if(bucket.m_data.first == k) {
                    return &bucket.m_data.second;
                }

                for(auto& elem : bucket.m_extra) {
                    if(elem.first == k)
                        return &elem.second;
                }

                bucket.m_extra.emplace_front(k, V());
                return &bucket.m_extra.front().second;                
            }
            const TKmer& kmer;
            size_t index;
        };
        //returns pointer to mapped value if exists, otherwise nullptr
        struct find : public boost::static_visitor<const V*> {
            find(const TKmer& k) : kmer(k) {}
            template <typename T> const V* operator()(T& v) const {
                typedef typename  T::value_type::large_t large_t;
                auto& bucket = v[kmer.oahash()%v.size()];
                if(bucket.m_data.second == V())
                    return nullptr;
                const large_t& k = kmer.get<large_t>();
                if(bucket.m_data.first == k)
                    return &bucket.m_data.second;
                for(auto& elem : bucket.m_extra) {
                    if(elem.first == k)
                        return &elem.second;
                }
                return nullptr;            
            } 
            const TKmer& kmer;
        };


        // returns position of element if found
        template <typename Bucket, typename element_t = typename Bucket::element_t>
        static element_t* find_in_bucket(Bucket& bucket, const TKmer& kmer) {
            typedef typename Bucket::large_t large_t;

            if(bucket.m_data.second == V())
                return nullptr;
            const large_t& k = kmer.get<large_t>();
            if(bucket.m_data.first == k)
                return &bucket.m_data;
            for(auto& elem : bucket.m_extra) {
                if(elem.first == k)
                    return &elem;
            }
            return nullptr;            
        }
        // we KNOW it is a new entry (used for rehashing)
        template <typename Bucket, typename element_t = typename Bucket::element_t>
        static void insert_in_bucket(Bucket& bucket, const element_t& elem) {
            if(bucket.m_data.second == V())
                bucket.m_data = elem;
            else            
                bucket.m_extra.push_front(elem);
        }
        // if kmer already included returns position of the element
        // if not inserts a new entry and returns position
        // caller MUST update the mapped value
        // also returns the position before the entry (caller may decide to swap them to put the most used element first)
        template <typename Bucket, typename element_t = typename Bucket::element_t> 
        static pair<element_t*,element_t*> check_and_insert_in_bucket(Bucket& bucket, const TKmer& kmer) {
            typedef typename Bucket::large_t large_t;
            const large_t& k = kmer.get<large_t>();

            pair<element_t*,element_t*> rslt;

            if(bucket.m_data.second == V()) {
                bucket.m_data.first = k;
                rslt.first = &bucket.m_data;
                return rslt;
            } else if(bucket.m_data.first == k) {
                rslt.first = &bucket.m_data;
                return rslt;
            }

            rslt.second = &bucket.m_data;
            for(auto& elem : bucket.m_extra) {
                if(elem.first == k) {
                    rslt.first = &elem;
                    return rslt;
                }
                rslt.second = &elem;
            }

            rslt.second = &bucket.m_data;
            bucket.m_extra.emplace_front(k, V());
            rslt.first = &bucket.m_extra.front();                

            return rslt;
        }

        struct resize : public boost::static_visitor<> { 
            resize(size_t s) : size(s) {}
            template <typename T> void operator()(T& v) const { v.resize(size); }            
            size_t size;
        };

        TKmerHashTable<V> m_hash_table;
        size_t m_table_size;
        int m_kmer_len;
    };

    struct SKmerCounter {
        SKmerCounter() : m_count(0), m_plus_count(0) {}
        bool operator==(const SKmerCounter& kc) const { return kc.m_count == m_count && kc.m_plus_count == m_plus_count; }
        uint32_t m_count;
        uint32_t m_plus_count;
    };
    class CKmerHashCount : public CKmerHashMap<SKmerCounter>, CCentinel {
    public:
        CKmerHashCount(int kmer_len, size_t size) : CKmerHashMap(kmer_len, size), CCentinel(size) {}

        // returnts true if kmer was new
        bool UpdateCount(const TKmer& kmer, bool is_plus) {
            size_t index = kmer.oahash()%m_table_size;
            //grab the bucket
            GrabBucket(index);
            //update count or insert new element
            bool rslt = apply_visitor(update_count(kmer, is_plus, index), m_hash_table);
            //release the bucket
            ReleaseBucket(index); 
            return rslt;
        }        
        // rehash bucket from other container
        void RehashOtherBuckets(CKmerHashCount& other, size_t bucket_from, size_t bucket_to) {
            apply_visitor(rehash_bucket(bucket_from, bucket_to, *this), m_hash_table, other.m_hash_table); 
        }
        //remove false positives
        size_t CleanBuckets(int min_count, size_t bucket_from, size_t bucket_to) {
            return apply_visitor(clean_buckets(min_count, bucket_from, bucket_to), m_hash_table);
        }
        void Info() const { apply_visitor(info(), m_hash_table); }
        friend class CDBHashGraph;
    private:
        struct info : public boost::static_visitor<> {
            template <typename T> void operator()(T& v) const { 
                map<int,int> numbers;
                for(auto& bucket : v) {
                    int num = distance(bucket.m_extra.begin(), bucket.m_extra.end());
                    if(bucket.m_data.second.m_count != 0)
                        ++num;
                    ++numbers[num];
                }
                for(auto& rslt : numbers)
                    cerr << "Bucket:\t" << rslt.first << "\t" << rslt.second << endl;
            }            
        };
        struct clean_buckets : public boost::static_visitor<size_t> {
            clean_buckets(int mc, size_t bf, size_t bt) : min_count(mc), bucket_from(bf), bucket_to(bt) {}
            template <typename T> size_t operator()(T& v) const {
                typedef typename T::value_type::element_t element_t;

                size_t num = 0;
                unsigned mc = min_count;
                for(size_t index = bucket_from; index <= bucket_to; ++index) {
                    auto& bucket = v[index];
                    bucket.m_extra.remove_if([mc](const element_t& elem) {return elem.second.m_count < mc;});
                    num += distance(bucket.m_extra.begin(), bucket.m_extra.end());
                    if(bucket.m_data.second.m_count >= mc) {
                        ++num;
                    } else if(!bucket.m_extra.empty()) {
                        bucket.m_data = bucket.m_extra.front();
                        bucket.m_extra.erase_after(bucket.m_extra.before_begin());
                    } else {
                        bucket.m_data.second.m_count = 0;
                        bucket.m_data.second.m_plus_count = 0;
                    }
                }

                return num;
            }            
            unsigned min_count;
            size_t bucket_from;
            size_t bucket_to;
        };
        friend struct rehash_bucket;
        struct rehash_bucket : public boost::static_visitor<> {
            rehash_bucket(size_t bf, size_t bt, CKmerHashCount& h) : bucket_from(bf), bucket_to(bt), hash(h) {}
            template <typename T> void operator() (T& a, T& b) const {
                typedef typename T::value_type::element_t element_t;
                for(size_t indexb = bucket_from; indexb <= bucket_to; ++indexb) {
                    auto& bucket_b = b[indexb];
                    if(bucket_b.m_data.second.m_count == 0)
                        continue;

                    auto& kmer = bucket_b.m_data.first;
                    size_t indexa = kmer.oahash()%a.size();
                    //grab the bucket   
                    hash.GrabBucket(indexa);
                    insert_in_bucket(a[indexa], bucket_b.m_data);
                    //release the bucket    
                    hash.ReleaseBucket(indexa);
                    
                    for(element_t& elem : bucket_b.m_extra) {
                        auto& kmer = elem.first;
                        size_t indexa = kmer.oahash()%a.size();
                        //grab the bucket   
                        hash.GrabBucket(indexa);
                        insert_in_bucket(a[indexa],elem );
                        //release the bucket    
                        hash.ReleaseBucket(indexa);
                    }
                }
            }
            template <typename T, typename U> void operator() (T& a, U& b) const { throw runtime_error("Can't rehash from different type container"); }
            
            size_t bucket_from;
            size_t bucket_to;
            CKmerHashCount& hash;
        };
        struct update_count : public boost::static_visitor<bool> { 
            update_count(const TKmer& k, bool isp, size_t i) : kmer(k), index(i), is_plus(isp) {}
            template <typename T> bool operator()(T& v) const {
                typedef typename T::value_type::element_t element_t;

                pair<element_t*,element_t*> rslt = check_and_insert_in_bucket(v[index], kmer);
                element_t* elemp = rslt.first;
                bool new_entry = elemp->second.m_count == 0;
                ++elemp->second.m_count;
                if(is_plus)
                    ++elemp->second.m_plus_count;
                element_t* prev_elemp = rslt.second;
                if(prev_elemp != nullptr && elemp->second.m_count > prev_elemp->second.m_count)
                    swap(*elemp, *prev_elemp);

                return new_entry;
            }

            const TKmer& kmer;
            size_t index;
            bool is_plus;
        };
    };


    class CKmerHashCounter {

    public:
        CKmerHashCounter(const list<array<CReadHolder,2>>& reads, int kmer_len, int min_count, bool is_stranded, int ncores) : 
            m_kmer_len(kmer_len), m_min_count(min_count), m_is_stranded(is_stranded), m_ncores(ncores), m_hash_table(m_kmer_len, 10000000), 
            m_kmer_num(0), m_kmer_num_raw(0), m_kmer_count(0), m_rehash_status(false) {

            m_kmer_step = max(1., 0.1*m_hash_table.TableSize()/m_ncores);
            for(auto& rholder : reads) 
                m_start_position.push_back(make_pair(0, rholder[0].kbegin(m_kmer_len)));

            CStopWatch timer;
            timer.Restart();

            size_t estimated_kmer_num = 100000000;
            double false_positive_rate = 0.05;
            size_t bloom_table_size = -(double)estimated_kmer_num*log(false_positive_rate)/log(2.)/log(2.);
            int hash_num = ceil(log(2.)*bloom_table_size/estimated_kmer_num);
            cerr << "Bloom table size: " << bloom_table_size << " Hash num:" << hash_num << endl;

            //REMOVE            CConcurrentBlockedBloomFilter bloom(bloom_table_size, 2, 512, hash_num);
            CConcurrentBlockedBloomFilter bloom(bloom_table_size, 8, 512, hash_num);
            {
                list<function<void()>> jobs;
                for(auto& job_input : reads) {
                    if(job_input[0].ReadNum() > 0 || job_input[1].ReadNum() > 0) {   // not empty               
                        jobs.push_back(bind(&CKmerHashCounter::InsertInBloomJob, this, ref(job_input), ref(bloom)));
                    }
                }
                RunThreads(m_ncores, jobs);
            }

            cerr << "Bloom filter in " << timer.Elapsed();
            timer.Restart();
            
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
                    CKmerHashCount hash_table_tmp(m_kmer_len, new_size);
                    swap(m_hash_table, hash_table_tmp);
                    m_kmer_step = max(1., 0.1*m_hash_table.TableSize()/m_ncores);
                    m_rehash_status.store(false);
                
                    list<function<void()>> jobs;
                    size_t step = ceil((double)hash_table_tmp.TableSize()/m_ncores);
                    for(int thr = 0; thr < m_ncores; ++thr) {
                        size_t from = step*thr;
                        size_t to = min(hash_table_tmp.TableSize()-1,from+step-1);
                        if(to >= from)
                            jobs.push_back(bind(&CKmerHashCounter::RehashJob, this, ref(hash_table_tmp), from, to));
                    }
                    RunThreads(m_ncores, jobs);
                    //                    cerr << "Rehashing in " << timer.Elapsed();
                }
            }

            cerr << "Create hash in " << timer.Elapsed();
            timer.Restart();            

            //Remove false positives
            {
                list<function<void()>> jobs;
                size_t step = ceil((double)m_hash_table.TableSize()/m_ncores);
                for(int thr = 0; thr < m_ncores; ++thr) {
                    size_t from = step*thr;
                    size_t to = min(m_hash_table.TableSize()-1,from+step-1);
                    if(to >= from)
                        jobs.push_back(bind(&CKmerHashCounter::CleanJob, this, from, to));
                }
                RunThreads(m_ncores, jobs);
            }

            cerr << "Clean hash in " << timer.Elapsed();
            timer.Restart();            

            cerr << "Initial kmers: " << m_kmer_num_raw.load() << " Distinct kmers: " << m_kmer_num.load() << " Total kmers: " << m_kmer_count.load() << " Hash table: " <<  m_hash_table.TableSize() << endl;

        }
        void Info() const {
            m_hash_table.Info();
        }
        

    private:
        
        class CConcurrentBlockedBloomFilter {
        public:
            // table_size - number of elements in bloom filter
            // counter_bit_size - number of bith per counter (2, 4, 8...)
            // block_size - in bits (multiple of 64)
            // hash_num - number of has functions (generated from two)
            CConcurrentBlockedBloomFilter(size_t table_size, int counter_bit_size, int block_size, int hash_num) : m_counter_bit_size(counter_bit_size), m_block_size(block_size), m_hash_num(hash_num),
                                                                                                                   m_max_element((size_t(1) << m_counter_bit_size) - 1) {
                m_elements_in_block = block_size/m_counter_bit_size;
                m_blocks = ceil((double)table_size/m_elements_in_block);
                m_table_size = m_blocks*m_elements_in_block;
                size_t table_length = m_table_size*m_counter_bit_size/(8*sizeof(TElement));
                m_count_table.resize(table_length, 0);
            }
            void Insert(size_t hashp, size_t hashm) {
                size_t blk_pos = hashp%m_blocks*m_block_size;
                for(int h = 1; h < m_hash_num; ++h) {
                    size_t hash = hashp+h*hashm;
                    size_t pos = blk_pos+hash%m_elements_in_block*m_counter_bit_size;
                    auto& cell = m_count_table[pos/(8*sizeof(TElement))].m_atomic;

                    int shift = pos%(8*sizeof(TElement));
                    size_t mask = m_max_element << shift;
                    size_t one = size_t(1) << shift;                               
                    auto existing_value = cell.load();
                    while((existing_value&mask) < mask && !cell.compare_exchange_strong(existing_value, existing_value+one)); // check for overflow and try to increment                    
                }            
            }
            int Test(size_t hashp, size_t hashm) {
                int count = numeric_limits<int>::max();
                size_t blk_pos = hashp%m_blocks*m_block_size;
                for(int h = 1; h < m_hash_num; ++h) {
                    size_t hash = hashp+h*hashm;
                    size_t pos = blk_pos+hash%m_elements_in_block*m_counter_bit_size;
                    auto& cell = m_count_table[pos/(8*sizeof(TElement))].m_atomic;

                    int shift = pos%(8*sizeof(TElement));
                    size_t mask = m_max_element << shift;
                    int cn = (cell.load()&mask) >> shift;
                    if(cn <= 1)
                        return cn;
                    count = min(count, cn);
                }

                return count;
            }
            size_t MaxElement() const { return m_max_element; }
            int HashNum() const { return m_hash_num; }

        private:
            size_t m_table_size;
            int m_counter_bit_size;
            int m_block_size;
            int m_hash_num;
            size_t m_blocks;
            int m_elements_in_block;
            size_t m_max_element;
            typedef SAtomic<size_t> TElement;
            vector<TElement> m_count_table;
        };

        void CleanJob(size_t bucket_from, size_t bucket_to) {
             m_kmer_num += m_hash_table.CleanBuckets(m_min_count, bucket_from, bucket_to);
        }
        void InsertInBloomJob(const array<CReadHolder,2>& rholder, CConcurrentBlockedBloomFilter& bloom) {
            for(int p = 0; p < 2; ++p) {
                for(CReadHolder::kmer_iterator itk = rholder[p].kbegin(m_kmer_len); itk != rholder[p].kend(); ++itk) {
                    TKmer kmer = *itk;
                    TKmer rkmer = revcomp(kmer, m_kmer_len);
                    size_t hashp = kmer.oahash();
                    size_t hashm = rkmer.oahash();
                    if(rkmer < kmer)
                        swap(hashp, hashm);
                    bloom.Insert(hashp, hashm);                     
                }
            }
        }
        void RehashJob(CKmerHashCount& other_hash_table, size_t bucket_from, size_t bucket_to) {
            m_hash_table.RehashOtherBuckets(other_hash_table, bucket_from, bucket_to);
        }
        void CountKmersJob(const array<CReadHolder,2>& rholder, pair<int,CReadHolder::kmer_iterator>& start_pos, CConcurrentBlockedBloomFilter& bloom) {
            size_t kmer_num = 0;
            size_t kmer_count = 0;
            for(int p = start_pos.first; p < 2; ++p) {
                CReadHolder::kmer_iterator from = start_pos.second;
                if(p != start_pos.first)
                    from = rholder[p].kbegin(m_kmer_len);
                for(CReadHolder::kmer_iterator itk = from; itk != rholder[p].kend(); ++itk) {
                    TKmer kmer = *itk;
                    TKmer rkmer = revcomp(kmer, m_kmer_len);
                    TKmer* min_kmerp = &rkmer; 
                    bool is_plus = false;
                    size_t hashp = rkmer.oahash();
                    size_t hashm = kmer.oahash();
                    if(kmer < rkmer) {
                        min_kmerp = &kmer;
                        is_plus = true;
                        swap(hashp, hashm);
                    }
                    ++kmer_count;                    

                    if(bloom.Test(hashp, hashm) < min(m_min_count, (int)bloom.MaxElement()))
                        continue;                    

                    if(m_hash_table.UpdateCount(*min_kmerp, is_plus))
                        ++kmer_num;
                    if(kmer_num == m_kmer_step) {
                        m_kmer_num_raw += kmer_num;
                        m_kmer_count += kmer_count;
                        kmer_num = 0;
                        kmer_count = 0;
                        if(m_kmer_num_raw.load() > m_hash_table.TableSize()*m_max_load_factor)
                            m_rehash_status.store(true); 
                        if(m_rehash_status.load()) {
                            start_pos.first = p;
                            ++itk;
                            start_pos.second = itk;
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

        CKmerHashCount m_hash_table;
        atomic<size_t> m_kmer_num;
        atomic<size_t> m_kmer_num_raw;
        atomic<size_t> m_kmer_count;
        atomic<bool> m_rehash_status;
        size_t m_kmer_step;
        double m_max_load_factor = 1;
        int m_increase_factor = 2;
        list<pair<int,CReadHolder::kmer_iterator>> m_start_position;
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
