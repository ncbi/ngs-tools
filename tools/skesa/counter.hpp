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

#ifndef _KmerCounter_
#define _KmerCounter_

#include "DBGraph.hpp"

namespace DeBruijn {

    // CKmerCounter counts kmers in reads using multiple threads and stores them in TKmerCount
    // It also finds neighbors (in GetBranches) if a user wants to use this class to build a CDBGraph (de Bruijn graph)
    // As Kmer counting could be memory expensive, CKmerCounter accepts an upper limit for the memory available and will 
    // subdivide the task, if needed.
    // If the number of subtasks exceeds 10, it will throw an exception asking for more memory.

    class CKmerCounter {
    public:

        // reads - raw reads (ncores or more elements in the list)
        // kmer_len - size of kmer
        // min_count - minimal count for accepted kmers
        // is_stranded - flag indicating whether kmers are from input reads where strand is informative or from connected paired
        //               reads generated internally by the program where strand is not a meaningful observation
        // mem_available - allowed memory in bytes
        // ncores - number of cores
        CKmerCounter(const list<array<CReadHolder,2>>& reads, int kmer_len, int min_count, bool is_stranded, int64_t mem_available, int ncores) : 
            m_kmer_len(kmer_len), m_min_count(min_count), m_is_stranded(is_stranded), m_mem_available(mem_available), m_ncores(ncores), m_reads(reads) {

            cerr << endl << "Kmer len: " << m_kmer_len << endl;
            CStopWatch timer;
            timer.Restart();

            int64_t raw_kmer_num = 0;
            for(const auto& reads : m_reads)
                raw_kmer_num += reads[0].KmerNum(m_kmer_len)+reads[1].KmerNum(m_kmer_len);

            int64_t GB = 1000000000;
            int kmer_size = TKmerCount(m_kmer_len).ElementSize();
            int64_t mem_needed = 1.2*raw_kmer_num*kmer_size;

            int max_cycles = 10;  // maximum cycles allowed
            int64_t mbuf = 2*GB;  // memory buffer for allocation uncertainity
            if(mem_needed >= max_cycles*(mem_available-mbuf)) {
                throw runtime_error("Memory provided is insufficient to do runs in 10 cycles for the read coverage. We find that 16 Gb for 20x coverage of a 5 Mb genome is usually sufficient");
            }
            int cycles = ceil(double(mem_needed)/(mem_available-mbuf));

            cerr << "Raw kmers: " << raw_kmer_num  << " Memory needed (GB): " << double(mem_needed)/GB << " Memory available (GB): " << double(mem_available-mbuf)/GB << " " << cycles << " cycle(s) will be performed" << endl;
        
            int njobs = 8*m_reads.size();   // many buckets reduce short-lived memory overhead spike in SortAndMergeJob    
            int kmer_buckets = cycles*njobs; 
    
            for(int cycl = 0; cycl < cycles; ++cycl) {
                pair<int,int> bucket_range(cycl*njobs, (cycl+1)*njobs-1);
                list<vector<TKmerCount>> raw_kmers;

                list<function<void()>> jobs;
                for(auto& job_input : m_reads) {
                    if(job_input[0].ReadNum() > 0 || job_input[1].ReadNum() > 0) {   // not empty       
                        raw_kmers.push_back(vector<TKmerCount>());
                        jobs.push_back(bind(&CKmerCounter::SpawnKmersJob, this, ref(job_input), kmer_buckets, bucket_range, ref(raw_kmers.back())));
                    }
                }
                RunThreads(ncores, jobs);

                // size_t total = 0;
                // for(auto& v : raw_kmers) {
                //     for(auto& tc : v) 
                //         total += tc.MemoryFootprint();
                // }
                SortAndMergeKmers(raw_kmers);
            }
    
            size_t utotal = 0;
            for(auto& c : m_uniq_kmers)
                utotal += c.Size();

            cerr << "Distinct kmers: " << utotal << endl;    
            cerr << "Kmer count in " << timer.Elapsed();
            
            MergeSortedKmers();
            if(m_uniq_kmers.empty())
                m_uniq_kmers.push_back(TKmerCount(m_kmer_len));                        
        }
        virtual ~CKmerCounter() {}

        // reference to counted kmers
        TKmerCount& Kmers() { return m_uniq_kmers.front(); }
        const TKmerCount& Kmers() const { return m_uniq_kmers.front(); }

        // average count of kmers in the histogram with the main peak
        double AverageCount() const {
            map<int,size_t> bins;
            for(size_t index = 0; index < Kmers().Size(); ++index) {
                ++bins[Kmers().GetCount(index)];                  // count clipped to integer automatically
            }
            TBins hist(bins.begin(), bins.end());
            pair<int,int> grange =  HistogramRange(hist);
            if(grange.first < 0)
                grange.first = 0;

            int genome = 0;
            size_t kmers = 0;
            for(int i = grange.first; i <= grange.second; ++i) {
                genome += hist[i].second;
                kmers += hist[i].first*hist[i].second;
            }

            if(genome > 0)
                return double(kmers)/genome;
            else
                return 0.;
        }

        // prepares kmer counts to be used in  CDBGraph (de Bruijn graph)
        // runs multiple instances of GetBranchesJob
        void GetBranches() {
            CStopWatch timer;
            timer.Restart();
            if(Kmers().Size() > 0) {
                vector<uint8_t> branches(Kmers().Size());
                size_t bucket_size = Kmers().Size()/m_ncores+1;
                list<function<void()>> jobs;
                for(int i = 0; i < m_ncores; ++i) {
                    pair<size_t,size_t> range(bucket_size*i,min(bucket_size*(i+1)-1,Kmers().Size()-1));
                    if(range.second >= range.first)
                        jobs.push_back(bind(&CKmerCounter::GetBranchesJob, this, range, ref(branches)));
                }
                RunThreads(m_ncores, jobs);

                for(size_t index = 0; index < Kmers().Size(); ++index) {
                    size_t b = branches[index];
                    size_t count = Kmers().GetCount(index);
                    uint32_t total_count = count;
                    uint32_t plus_count = (count >> 32);
                    size_t plusf = uint16_t(double(plus_count)/total_count*numeric_limits<uint16_t>::max()+0.5);
                    Kmers().UpdateCount((plusf << 48)+(b << 32)+total_count, index);  // we put strand info and branching in the high half of the count!!!!!                   
                }
            }

            cerr << "Kmers branching in " << timer.Elapsed();
        }

        bool IsStranded() const { return m_is_stranded; }              // indicates if contains stranded information

    private:

        // one-thread worker producing kmers and putting them in multiple non-overlapping buckets
        // rholder - input reads 
        // buckets - total number of buckets
        // bucket_range - range of buckets used by this worker
        // kmers - output kmers
        void SpawnKmersJob(const array<CReadHolder,2>& rholder, int buckets, pair<int,int> bucket_range,  vector<TKmerCount>& kmers) {
            size_t total = rholder[0].KmerNum(m_kmer_len)+rholder[1].KmerNum(m_kmer_len);
            size_t reserve = 1.1*total/buckets;
            int active_buckets = bucket_range.second-bucket_range.first+1;
            kmers.resize(active_buckets, TKmerCount(m_kmer_len));
            for(auto& k : kmers)
                k.Reserve(reserve);

            for(int p = 0; p < 2; ++p) {
                for(CReadHolder::kmer_iterator itk = rholder[p].kbegin(m_kmer_len); itk != rholder[p].kend(); ++itk) {
                    TKmer kmer = *itk;
                    TKmer rkmer = revcomp(kmer, m_kmer_len);
                    size_t count = 1;
                    TKmer* min_kmerp = &rkmer;
                    if(kmer < rkmer) {
                        min_kmerp = &kmer;
                        count += (size_t(1) << 32);
                    }
                    int bucket = min_kmerp->oahash()%buckets;
                    if(bucket < bucket_range.first || bucket > bucket_range.second)
                        continue;
                    // good to go   
                    int ind = bucket - bucket_range.first;
                    if(kmers[ind].Size() == kmers[ind].Capacity()) { //expensive plan B for the case of failed hash uniformity          
                        //            cerr << "Warning: Hash uniformity problem" << endl;           
                        TKmerCount bigger(m_kmer_len);
                        bigger.Reserve(kmers[ind].Size()*1.2);
                        bigger.PushBackElementsFrom(kmers[ind]);
                        bigger.Swap(kmers[ind]);
                    }
                    kmers[ind].PushBack(*min_kmerp, count);            
                }
            }
        }

        //SortAndMergeJob briefly doubles the input memory - should be executed in small chunks!!!!!!   
        // one-thread worker which accepts all containers for a given bucket and merges, sorts and counts them
        // group - list of containers
        // ukmers - counted kmers
        typedef list<TKmerCount*> TContainerPList;
        void SortAndMergeJob(TContainerPList group, TKmerCount& ukmers) {
            TKmerCount all_kmers(group.front()->KmerLen());
            if(group.size() == 1) {
                all_kmers.Swap(*group.front());
            } else {
                size_t total = 0;
                for(auto p : group)
                    total += p->Size();
                   
                all_kmers.Reserve(total); // doubles the input memory!!!!       
            
                for(auto p : group) {
                    all_kmers.PushBackElementsFrom(*p);
                    TKmerCount(p->KmerLen()).Swap(*p);
                }
            }

            all_kmers.SortAndExtractUniq(m_min_count, ukmers);         
        }

        // runs multiple instances of SortAndMergeJob and stores results in m_uniq_kmers
        // raw_kmers - input kmers
        void SortAndMergeKmers(list<vector<TKmerCount>>& raw_kmers) {

            list<function<void()>> jobs;
            int bucken_num = raw_kmers.front().size();

            for(int bucket = 0; bucket < bucken_num; ++bucket) {
                TContainerPList job_input;
                for(auto& vec : raw_kmers)
                    job_input.push_back(&vec[bucket]);
                m_uniq_kmers.push_back(TKmerCount());
                jobs.push_back(bind(&CKmerCounter::SortAndMergeJob, this, job_input, ref(m_uniq_kmers.back())));
            }
            RunThreads(m_ncores, jobs);
        }

        // one-thread worker which merges two sorted buckets
        static void MergeSortedJob(TKmerCount& akmers, TKmerCount& bkmers) {
            akmers.MergeTwoSorted(bkmers);
            TKmerCount(bkmers.KmerLen()).Swap(bkmers); // release bkmers memory
        }

        // runs multiple instances of MergeSortedJob
        // at the end m_uniq_kmers has only one element with final kmers
        void MergeSortedKmers() {
            CStopWatch timer;
            timer.Restart();
            while(m_uniq_kmers.size() > 1) {
                list<function<void()>> jobs;
                for(list<TKmerCount>::iterator first = m_uniq_kmers.begin(); first != m_uniq_kmers.end(); ++first) {
                    list<TKmerCount>::iterator second = first;
                    if(++second != m_uniq_kmers.end()) {
                        jobs.push_back(bind(MergeSortedJob, ref(*first), ref(*second)));
                        first = second;
                    }
                }
                RunThreads(m_ncores, jobs);
                for(auto iloop = m_uniq_kmers.begin(); iloop != m_uniq_kmers.end(); ) {
                    auto it = iloop++;
                    if(it->Size() == 0)
                        m_uniq_kmers.erase(it);
                }
            }
            cerr << "Uniq kmers merging in " << timer.Elapsed();
        }

        // one-thread worker which calculates the branching information (neighbors) for a range of kmers
        // range - from,to indexes for kmers
        // branches - vector of branching information (one bit is used for each of the eight possible neighbors)  
        void GetBranchesJob(pair<size_t,size_t> range, vector<uint8_t>& branches) {
            TKmer max_kmer(string(m_kmer_len, bin2NT[3]));
            for(size_t index = range.first; index <= range.second; ++index) {
                pair<TKmer,size_t> kmer_count = Kmers().GetKmerCount(index);
                //direct        
                TKmer shifted_kmer = (kmer_count.first << 2) & max_kmer;
                //inverse       
                TKmer shifted_rkmer = (revcomp(kmer_count.first, m_kmer_len) << 2) & max_kmer;
                for(int nt = 0; nt < 4; ++nt) {
                    TKmer k = shifted_kmer + TKmer(m_kmer_len, nt);
                    size_t new_index = Kmers().Find(min(k, revcomp(k, m_kmer_len)));
                    // New kmer is a neighbor if it exists in reads and is not same as current kmer
                    if(new_index != Kmers().Size() && new_index != index) 
                        branches[index] |= (1 << nt);

                    k = shifted_rkmer + TKmer(m_kmer_len, nt);
                    new_index = Kmers().Find(min(k, revcomp(k, m_kmer_len)));
                    if(new_index != Kmers().Size() && new_index != index) 
                        branches[index] |= (1 << (nt+4));
                }
            }
        }

        int m_kmer_len;
        int m_min_count;
        bool m_is_stranded;
        size_t m_mem_available;
        int m_ncores;
        const list<array<CReadHolder,2>>& m_reads;
        list<TKmerCount> m_uniq_kmers;                       // storage for kmer buckets; at the end will have one element which is the result     
    };

}; // namespace
#endif /* _KmerCounter_ */
