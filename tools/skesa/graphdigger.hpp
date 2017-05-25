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

#ifndef _GraphDigger_
#define _GraphDigger_

#include <stack> 
#include "DBGraph.hpp"

namespace DeBruijn {

/************************
General description

Class CDBGraphDigger, defined in this file, performs most of the actual assembling work. 

struct SContig is used to hold assembled sequences. Members are 
        deque<char> m_seq;               // sequence representing contig
        deque<CDBGraph::Node> m_kmers;   // all kmers of this contig sequence
        int m_kmer_len;                  // size of kmers used for building the contig
        CDBGraph::Node m_next_left;      // denied left kmer (connection possible but it is already owned)
        CDBGraph::Node m_next_right;     // denied right kmer (connection possible but it is already owned)
        SContig* m_left_link;            // if set points to 'left' contig
        SContig* m_right_link;           // if set points to 'right' contig
        int m_left_shift;                // shift+1 for m_next_left in this contig (positive for the right end)
        int m_right_shift;               // shift+1 for m_next_right in this contig (positive for the right end)
        int m_left_extend;               // number of newly assembled bases which could be clipped
        int m_right_extend;              // number of newly assembled bases which could be clipped
        SAtomic<uint8_t> m_is_taken;

There are three main scenarios how this structure may be created

1) From a previously assembled contig represented by a c++ string. 
   In this case, no members other than m_seq, m_kmers, and m_kmer_len change their default zero value. 

2) Assembled starting from one of the kmers which has not been used so far. Because assembly is done in multiple threads,
   two or more threads could start assembling the same contig from different starting kmers. At some point they will
   collide with each other and will try to obtain a kmer which has already been used by the other contig in a different
   thread. When this happens, the thread stops extending the contig and assigns the denied kmer to m_next_left or m_next_right.
   These partially assembled contigs (which internally are called fragments) could be connected to each other using
   m_next_left/m_next_right. It is done in ConnectFragments().

3) When we increase the kmer size, some previously assembled contigs could be extended or connected because
   the longer kmer could resolve some repeats. To achieve this, we assemble new contigs starting from each of the
   flank kmers. When these contigs are started, m_left_link/m_right_link are assigned to point to the
   parent contig and m_left_shift/m_right_shift are assigned to indicate the start position. Because the work is done
   in mutiple threads, the contigs could come in two fragments if they started from different contigs and come together.
   The connected fragments will have links on both sides. These contigs are 'connectors'. The rest of contigs are
   'extenders'. They are used by ConnectAndExtendContigs() to form a new contig set. There is an important corner case
   for connectors. A thread could start from contig A and finish assembling a connector all the way to contig B before
   some other threads starts dealing with B. In this case the sequence starting from B will not contain any real bases
   but will only have m_next_left/m_next_right and a link. Those are mentioned in the code as 'empty linkers' and should
   be treated as special cases. 

A note about multiprocessing in 2) and 3): In both cases, the threads should be able to decide which kmers or contigs are
still available for work. For communicating this information between the threads, the code uses lock-free c++ atomic
variables. For the kmers, this is stored in m_visited vector in CDBGraph. For the contigs, this is stored in m_is_taken.

The length of newly assembled sequence is stored in  m_left_extend/m_right_extend.

************************/


    mutex out_mutex;


    struct SContig {
        typedef forward_list<CDBGraph::Node> TNodeList;
        SContig(CDBGraph& graph) :m_graph(graph), m_kmer_len(graph.KmerLen()) {}
        SContig(const CContigSequence& contig, CDBGraph& graph) : m_seq(contig), m_graph(graph), m_kmer_len(graph.KmerLen()) { GenerateKmersAndCleanSNPs(); }        

        void GenerateKmersAndCleanSNPs() {
            m_seq.RemoveShortUniqIntervals(m_kmer_len);
            for(int i = m_seq.size()-1; i >= 0; ) {
                if(i == (int)m_seq.size()-1) {
                    if((int)m_seq.ChunkLenMax(i) >= m_kmer_len) {  // last chunk size >= kmer_len
                        CReadHolder rh(false);
                        rh.PushBack(m_seq.back().front());
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(m_kmer_len) ; ik != rh.kend(); ++ik) {
                            CDBGraph::Node node = m_graph.GetNode(*ik);
                            if(node && !m_graph.SetVisited(node))
                                m_graph.SetMultContig(node);
                        }
                    }
                    
                    --i;
                } else { // all uniq chunks >= kmer_len-1
                    if((int)m_seq.ChunkLenMax(i-1) >= m_kmer_len) {
                        CReadHolder rh(false);
                        rh.PushBack(m_seq[i-1].front());
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(m_kmer_len); ik != rh.kend(); ++ik) {
                            CDBGraph::Node node = m_graph.GetNode(*ik);
                            if(node && !m_graph.SetVisited(node))
                                m_graph.SetMultContig(node);
                        }
                    }

                    unordered_set<CDBGraph::Node> kmers;
                    list<deque<CDBGraph::Node>> failed_nodes;
                    list<TLocalVariants::iterator> failed_variants;
                    for(auto prev = m_seq[i].before_begin(); ;++prev) {
                        auto current = prev;
                        if(++current == m_seq[i].end())
                            break;

                        auto& variant = *current;
                        int left = min(m_kmer_len-1, (int)m_seq.ChunkLenMax(i-1));
                        TVariation var_seq(m_seq[i-1].front().end()-left, m_seq[i-1].front().end());
                        var_seq.insert(var_seq.end(), variant.begin(), variant.end());
                        int right = min(m_kmer_len-1, (int)m_seq.ChunkLenMax(i+1));
                        var_seq.insert(var_seq.end(), m_seq[i+1].front().begin(), m_seq[i+1].front().begin()+right);
                        CReadHolder rh(false);
                        rh.PushBack(var_seq);
                        deque<CDBGraph::Node> var_nodes;
                        bool failed = false;
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(m_kmer_len) ; ik != rh.kend(); ++ik)  {
                            var_nodes.emplace_back(m_graph.GetNode(*ik));
                            if(!var_nodes.back())
                                failed = true;
                        }
                        if(failed) {
                            failed_nodes.push_back(move(var_nodes));
                            failed_variants.push_front(prev);           // reverse order for deleting
                        } else {
                            kmers.insert(var_nodes.begin(), var_nodes.end()); 
                        }
                    }

                    if((int)failed_variants.size() == m_seq.VariantsNumber(i)) { // all failed
                        for(auto& nodes : failed_nodes)
                            kmers.insert(nodes.begin(), nodes.end());
                    } else { // some are good
                        /*
                        if(!failed_variants.empty()) {
                            lock_guard<mutex> guard(out_mutex);
                            cerr << "Removed " << failed_variants.size() << " SNPs from " << m_seq.VariantsNumber(i) << endl;
                        }
                        */

                        for(auto prev : failed_variants)
                            m_seq[i].erase_after(prev);
                        if(m_seq.UniqueChunk(i)) { // only one left
                            m_seq[i-1].front().insert(m_seq[i-1].front().end(), m_seq[i].front().begin(), m_seq[i].front().end());
                            m_seq[i-1].front().insert(m_seq[i-1].front().end(), m_seq[i+1].front().begin(), m_seq[i+1].front().end());
                            m_seq.erase(m_seq.begin()+i, m_seq.begin()+i+2);
                        }
                    }

                    for(auto& node : kmers) {
                        if(node && !m_graph.SetVisited(node))
                            m_graph.SetMultContig(node);                        
                    }
                     
                    i -= 2;   
                }
            }
            m_seq.ContractVariableIntervals();
            m_seq.IncludeRepeatsInVariableIntervals();
            m_seq.RemoveShortUniqIntervals(m_kmer_len);
            m_seq.StabilizeVariantsOrder();
        }

        SContig(const SContig& to_left, const SContig& to_right, CDBGraph::Node initial_node, CDBGraph::Node lnode, CDBGraph::Node rnode, CDBGraph& graph) :  
        m_next_left(lnode), m_next_right(rnode), m_graph(graph),  m_kmer_len(graph.KmerLen()) {                                                                                                                                                   
//          initial_node - the starting kmer
//          to_left - left extension of the starting kmer
//          to_right - right extension of the starting kmer
//          lnode - left denied node
//          rnode - right denied node
//          graph - de Bruijn graph

                                                                       
            // take parts of the assembled sequence and put them together in SContig

            if(!to_left.m_seq.empty()) {
                m_seq = to_left.m_seq;
                ReverseComplement();
            }

            // could be changed by ReverseComplement
            m_next_left = lnode;
            m_next_right = rnode;

            string ikmer = graph.GetNodeSeq(initial_node);
            if(m_seq.empty() || m_seq.VariableChunk(m_seq.size()-1))  // empty or variant
                m_seq.InsertNewChunk(ikmer.begin(), ikmer.end());
            else
                m_seq.ExtendTopVariant(ikmer.begin(), ikmer.end());

            if(!to_right.m_seq.empty()) {
                if(to_right.m_seq.UniqueChunk(0)) {
                    m_seq.ExtendTopVariant(to_right.m_seq.front().front().begin(), to_right.m_seq.front().front().end());
                    m_seq.insert(m_seq.end(), to_right.m_seq.begin()+1, to_right.m_seq.end());
                } else {
                    m_seq.insert(m_seq.end(), to_right.m_seq.begin(), to_right.m_seq.end());
                }
            }
            m_seq.StabilizeVariantsOrder();
                
            m_left_extend = m_right_extend = LenMax();
        }
        SContig(SContig* link, int shift, CDBGraph::Node takeoff_node, const SContig& extension, CDBGraph::Node rnode, CDBGraph& graph) :
            m_next_left(takeoff_node), m_next_right(rnode), m_left_link(link), m_left_shift(shift),  m_graph(graph), m_kmer_len(graph.KmerLen()) {

            /*
            m_seq.m_left_repeat = 0;
            m_seq.m_right_repeat = m_kmer_len-1; 
            */

            string kmer = graph.GetNodeSeq(takeoff_node);
            m_seq.InsertNewChunk(kmer.begin()+1, kmer.end()); // don't include first base
            if(!extension.m_seq.empty()) {
                if(extension.m_seq.UniqueChunk(0)) {
                    m_seq.ExtendTopVariant(extension.m_seq.front().front().begin(), extension.m_seq.front().front().end());
                    m_seq.insert(m_seq.end(), extension.m_seq.begin()+1, extension.m_seq.end());
                } else {
                    m_seq.insert(m_seq.end(), extension.m_seq.begin(), extension.m_seq.end());
                }
            }
            m_seq.StabilizeVariantsOrder();

            m_left_extend = m_right_extend = LenMax();
        }

        CDBGraph::Node FrontKmer() const {
            if(m_seq.VariableChunk(0) || (int)m_seq.ChunkLenMax(0) < m_kmer_len)
                return 0;

            TKmer kmer(m_seq.front().front().begin(), m_seq.front().front().begin()+m_kmer_len); // front must be unambiguous
            return m_graph.GetNode(kmer);
        }
        CDBGraph::Node BackKmer() const { 
            int last = m_seq.size()-1;
            if(m_seq.VariableChunk(last) || (int)m_seq.ChunkLenMax(last) < m_kmer_len)
                return 0;

            TKmer kmer(m_seq.back().front().end()-m_kmer_len, m_seq.back().front().end());
            return m_graph.GetNode(kmer);
        }

        // don't 'own' any kmers
        bool EmptyLinker() const { return ((int)max(m_seq.ChunkLenMax(0), m_seq.ChunkLenMax(m_seq.size()-1)) < m_kmer_len && m_seq.size() <= 3); }

        bool RightSNP() const { return (m_seq.size() >= 3 && m_seq.UniqueChunk(m_seq.size()-1) && (int)m_seq.ChunkLenMax(m_seq.size()-1) < m_kmer_len); }
        bool LeftSNP() const { return (m_seq.size() >= 3 && m_seq.UniqueChunk(0) &&  (int)m_seq.ChunkLenMax(0) < m_kmer_len); }

        CDBGraph::Node RightConnectingNode() const {
            int last_index = m_seq.size()-1;
            if((int)m_seq.ChunkLenMax(last_index) >= m_kmer_len) {   // normal end
                return BackKmer(); 
            } else if(m_seq.size() >= 3) {                           // snp
                if((int)m_seq.ChunkLenMax(last_index-2) >= m_kmer_len) {
                    TKmer kmer(m_seq[last_index-2].front().end()-m_kmer_len, m_seq[last_index-2].front().end());
                    return m_graph.GetNode(kmer);
                }
            }

            return m_next_left;                        // empty linker
        }
        CDBGraph::Node LeftConnectingNode() const {
            if((int)m_seq.ChunkLenMax(0) >= m_kmer_len) {  // normal end
                return FrontKmer();
            } else if(m_seq.size() >= 3) {                 // snp
                if((int)m_seq.ChunkLenMax(2) >= m_kmer_len) {
                    TKmer kmer(m_seq[2].front().begin(), m_seq[2].front().begin()+m_kmer_len); // front must be unambiguous
                    return m_graph.GetNode(kmer);
                }
            }

            return m_next_right;                        // empty linker
        }
        
        void ReverseComplement() {  
            m_seq.ReverseComplement();
            swap(m_next_left, m_next_right);
            m_next_left = CDBGraph::ReverseComplement(m_next_left);
            m_next_right = CDBGraph::ReverseComplement(m_next_right);
            swap(m_left_link, m_right_link);
            swap(m_left_shift, m_right_shift);
            swap(m_left_extend, m_right_extend);
        }
        void AddToRight(const SContig& other) {
            m_seq.m_circular = false;

            m_next_right = other.m_next_right;
            m_right_link = other.m_right_link;
            m_right_shift = other.m_right_shift; 
            if(EmptyLinker() && other.EmptyLinker())
                return;            
 
            auto& last_chunk = m_seq.back().front();
            int last_chunk_len = last_chunk.size();
            int overlap = m_kmer_len-1;
            auto first_other_chunk_it = other.m_seq.begin(); 
            if(RightSNP() && other.LeftSNP()) {    // skip snp chunk
                overlap = last_chunk_len+other.m_seq.ChunkLenMax(1)+first_other_chunk_it->front().size();
                first_other_chunk_it += 2;
            }

            if(other.m_right_extend < (int)other.LenMax()) {
                m_right_extend = other.m_right_extend;
            } else {
                m_right_extend += other.m_right_extend-overlap;
                if(m_left_extend == (int)LenMax())
                    m_left_extend = m_right_extend;
            }
            
            /*
            m_seq.m_right_repeat = other.m_seq.m_right_repeat;
            */

            auto& first_other_chunk = first_other_chunk_it->front();
            last_chunk.insert(last_chunk.end(), first_other_chunk.begin()+min(m_kmer_len-1,last_chunk_len), first_other_chunk.end());  // combine overlapping chunks
            m_seq.insert(m_seq.end(), first_other_chunk_it+1, other.m_seq.end());                                                       // insert remaining chunks
        }
        void AddToLeft(const SContig& other) {
            m_seq.m_circular = false;

            m_next_left = other.m_next_left;
            m_left_link = other.m_left_link;
            m_left_shift = other.m_left_shift;
            if(EmptyLinker() && other.EmptyLinker())
                return;            
 
            auto& first_chunk = m_seq.front().front();
            int first_chunk_len = first_chunk.size();
            int overlap = m_kmer_len-1;
            auto last_other_chunk_it = other.m_seq.end()-1; 
            if(LeftSNP() && other.RightSNP()) {    // skip snp chunk
                overlap = first_chunk_len+other.m_seq.ChunkLenMax(other.m_seq.size()-2)+last_other_chunk_it->front().size();
                last_other_chunk_it -= 2;
            }
                
            if(other.m_left_extend < (int)other.LenMax()) {
                m_left_extend = other.m_left_extend;
            } else {
                m_left_extend += other.m_left_extend-overlap; 
                if(m_right_extend == (int)LenMax())
                    m_right_extend = m_left_extend;
            }

            /*
            m_seq.m_left_repeat = other.m_seq.m_left_repeat;
            */

            auto& last_other_chunk = last_other_chunk_it->front();
            first_chunk.insert(first_chunk.begin(),last_other_chunk.begin(), last_other_chunk.end()-min(m_kmer_len-1,first_chunk_len));  // combine overlapping chunks
            m_seq.insert(m_seq.begin(), other.m_seq.begin(), last_other_chunk_it);                                     // insert remaining chunks
        }
        
        void ClipRight(int clip) {
            if(clip <= 0)
                return;

            m_seq.m_circular = false;
            m_next_right = 0;
            m_right_link = nullptr;
            m_right_shift = 0;            

            while(!m_seq.empty() && (m_seq.VariableChunk(m_seq.size()-1) || (int)m_seq.ChunkLenMax(m_seq.size()-1) <= clip)) {
                int chunk_len = m_seq.ChunkLenMax(m_seq.size()-1);
                clip -= chunk_len;
                m_right_extend = max(0, m_right_extend-chunk_len);
                /*
                m_seq.m_right_repeat = max(0, m_seq.m_right_repeat-chunk_len);
                */
                m_seq.pop_back();
            }
            if(clip > 0 && !m_seq.empty()) {
                m_right_extend = max(0, m_right_extend-clip);
                /*
                m_seq.m_right_repeat = max(0, m_seq.m_right_repeat-clip);
                */
                m_seq.back().front().erase(m_seq.back().front().end()-clip, m_seq.back().front().end());
            }

            if((int)LenMin() < m_kmer_len-1)
                m_seq.clear();
        }
        void ClipLeft(int clip) {
            if(clip <= 0)
                return;

            m_seq.m_circular = false;
            m_next_left = 0;
            m_left_link = nullptr;
            m_left_shift = 0; 

            while(!m_seq.empty() && (m_seq.VariableChunk(0) || (int)m_seq.ChunkLenMax(0) <= clip)) {
                int chunk_len = m_seq.ChunkLenMax(0);
                clip -= chunk_len;
                m_left_extend = max(0, m_left_extend-chunk_len);
                /*
                m_seq.m_left_repeat = max(0, m_seq.m_left_repeat-chunk_len);
                */
                m_seq.pop_front();
            }
            if(clip > 0 && !m_seq.empty()) {
                m_left_extend = max(0, m_left_extend-clip);
                /*
                m_seq.m_left_repeat = max(0, m_seq.m_left_repeat-clip);
                */
                m_seq.front().front().erase(m_seq.front().front().begin(), m_seq.front().front().begin()+clip);
            }         

            if((int)LenMin() < m_kmer_len-1)
                m_seq.clear();
        }        

        size_t LenMax() const { return m_seq.LenMax(); }
        size_t LenMin() const { return m_seq.LenMin(); }

        // find position of the minimal non-zero kmer
        tuple<int, int, CDBGraph::Node> MinKmerPosition() const {  //chunk, position in chunk, kmer
            unordered_map<CDBGraph::Node, tuple<int, int, CDBGraph::Node>> kmers;

            for(int i = m_seq.size()-1; i >= 0; i -= 2) {
                deque<forward_list<CDBGraph::Node>> chunk_kmers;
                if(i == (int)m_seq.size()-1) {
                    if((int)m_seq.ChunkLenMax(i) >= m_kmer_len) { // last chunk could be short
                        chunk_kmers.resize(m_seq.ChunkLenMax(i)-m_kmer_len+1);
                        CReadHolder rh(false);
                        rh.PushBack(m_seq.back().front());
                        int pos = chunk_kmers.size();
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(m_kmer_len) ; ik != rh.kend(); ++ik) // iteration from last kmer to first     
                            chunk_kmers[--pos].push_front(m_graph.GetNode(*ik));
                    }
                } else { // all uniq chunks in the middle >= kmer_len-1; first/last could be short
                    chunk_kmers.resize(m_seq.ChunkLenMax(i)+m_seq.ChunkLenMax(i+1));
                    if((int)m_seq.ChunkLenMax(i) >= m_kmer_len) {
                        TVariation seq(m_seq[i].front().begin(), m_seq[i].front().end());
                        CReadHolder rh(false);
                        rh.PushBack(seq);
                        int pos = seq.size()-m_kmer_len+1;
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(m_kmer_len) ; ik != rh.kend(); ++ik)  // iteration from last kmer to first 
                            chunk_kmers[--pos].push_front(m_graph.GetNode(*ik));
                    }
                    for(auto& variant : m_seq[i+1]) {
                        TVariation seq;
                        if((int)m_seq.ChunkLenMax(i) >= m_kmer_len-1)
                            seq.insert(seq.end(), m_seq[i].front().end()-m_kmer_len+1, m_seq[i].front().end());
                        else
                            seq.insert(seq.end(), m_seq[i].front().begin(), m_seq[i].front().end());
                        seq.insert(seq.end(), variant.begin(), variant.end());
                        if((int)m_seq.ChunkLenMax(i+2) >= m_kmer_len-1)
                            seq.insert(seq.end(), m_seq[i+2].front().begin(), m_seq[i+2].front().begin()+m_kmer_len-1);
                        else
                            seq.insert(seq.end(), m_seq[i+2].front().begin(), m_seq[i+2].front().end());
                        CReadHolder rh(false);
                        rh.PushBack(seq);
                        int pos = seq.size()-m_kmer_len+1;
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(m_kmer_len) ; ik != rh.kend(); ++ik)  // iteration from last kmer to first 
                            chunk_kmers[--pos].push_front(m_graph.GetNode(*ik));
                    }
                }

                for(unsigned pos = 0; pos < chunk_kmers.size(); ++pos) {
                    int k = pos;
                    int chunk = i;
                    if(pos >= m_seq.ChunkLenMax(i)) {
                        k = pos-m_seq.ChunkLenMax(i);
                        chunk = i+1;
                    }
                    for(auto& kmer : chunk_kmers[pos]) {
                        if(kmer) {
                            auto mn = min(kmer, m_graph.ReverseComplement(kmer));
                            auto rslt = kmers.insert(make_pair(mn, make_tuple(chunk, k, kmer)));
                            if(!rslt.second)
                                get<2>(rslt.first->second) = 0;
                        }
                    }
                }
            }

            tuple<int, int, CDBGraph::Node> rslt(0, 0, 0);
            for(auto& elem : kmers) {
                if(!get<2>(elem.second))
                    continue;
                if(!get<2>(rslt) || elem.first < min(get<2>(rslt), m_graph.ReverseComplement(get<2>(rslt))))
                    rslt = elem.second;
            }
            return rslt;
        }

        // stabilize contig orientation using minimal kmer in the contig
        void SelectMinDirection() {
            CDBGraph::Node minkmer = get<2>(MinKmerPosition());
            if(minkmer && minkmer%2)
                ReverseComplement();            
        }

        // finds stable origin for circular contigs by placing minimal kmer at the beginning of the sequence
        void RotateCircularToMinKmer() { // assumes that the next extension of sequence would give the first kmer (m_next_right == m_kmers.front())
            m_seq.back().front().erase(m_seq.back().front().end()-m_kmer_len+1, m_seq.back().front().end());
            auto rslt = MinKmerPosition();
            if(!get<2>(rslt))
                return;

            size_t first_chunk = get<0>(rslt);
            size_t first_base = get<1>(rslt);

            if(get<2>(rslt)%2) {
                first_base += m_kmer_len;
                while(first_base >= m_seq.ChunkLenMax(first_chunk)) {
                    first_base -= m_seq.ChunkLenMax(first_chunk);
                    first_chunk = (first_chunk+1)%m_seq.size();
                }
                if(m_seq.VariableChunk(first_chunk)) {               // ambiguous interval - we don't want to cut it
                    ++first_chunk;  // variable chunk cant be last
                    first_base = 1;
                } else if(first_chunk > 0 && first_base == 0) {       // we want some uniq intervals on both ends
                    first_base = 1;
                }
            } else {
                if(m_seq.VariableChunk(first_chunk)) {                // ambiguous interval - we don't want to cut it
                    --first_chunk;  // variable chunk cant be first
                    first_base = m_seq.ChunkLenMax(first_chunk)-1;     // leave one base
                } else if(first_chunk > 0 && first_base == 0) {       // we want some uniq intervals on both ends
                    first_chunk -= 2;
                    first_base = m_seq.ChunkLenMax(first_chunk)-1;     // leave one base
                }
            }

            if(m_seq.size() == 1) {
                rotate(m_seq.front().front().begin(), m_seq.front().front().begin()+first_base, m_seq.front().front().end()); 
            } else {
                if(first_chunk > 0) {
                    auto& last_seq = m_seq.back().front();
                    last_seq.insert(last_seq.end(), m_seq.front().front().begin(), m_seq.front().front().end());
                    m_seq.pop_front();
                    rotate(m_seq.begin(), m_seq.begin()+first_chunk-1, m_seq.end());
                }
                if(first_base > 0) {
                    if(m_seq.VariableChunk(m_seq.size()-1)) {
                        m_seq.InsertNewChunk();
                        m_seq.InsertNewVariant();
                    }
                    auto& last_seq = m_seq.back().front();
                    last_seq.insert(last_seq.end(), m_seq.front().front().begin(), m_seq.front().front().begin()+first_base);
                    m_seq.front().front().erase(m_seq.front().front().begin(), m_seq.front().front().begin()+first_base);
                }
            }
            
            //clean edges   
            m_next_left = 0;
            m_next_right = 0;
            m_left_link = nullptr;
            m_left_shift = 0;
            m_right_link = nullptr;
            m_right_shift = 0;
            m_left_extend = 0;   // prevents any further clipping    
            m_right_extend = 0;  // prevents any further clipping    
            /*
            m_seq.m_left_repeat = 0; 
            m_seq.m_right_repeat = 0;
            */
            m_seq.m_circular = true;
        }

        bool operator<(const SContig& other) const { return m_seq < other.m_seq; }

        // connects fragments created in different threads and combines doubled 'empty' linkers
        static TContigList ConnectFragments(vector<TContigList>& fragments, const CDBGraph& graph) {

            int total = 0;
            int len = 0;
            for(auto& ns : fragments) {
                for(auto& seq : ns) {
                    ++total;
                    len += seq.LenMax()-graph.KmerLen()+1;
                }
            }
            
            cerr << "Fragments before: " << total << " " <<  len << endl;

            TContigList connected;        
            unordered_map<CDBGraph::Node, TContigList::iterator> denied_left_nodes;
            unordered_map<CDBGraph::Node, TContigList::iterator> denied_right_nodes;
            for(auto& ns : fragments) {
                for(auto iloop = ns.begin(); iloop != ns.end(); ) {
                    auto ic = iloop++;
                    connected.splice(connected.begin(), ns, ic);
                    SContig& contig = *connected.begin();                
                    if(contig.m_next_left > contig.m_next_right)  // need this to pair two identical empty links
                        contig.ReverseComplement();                

                    if(contig.m_next_left) {
                        auto rslt = denied_left_nodes.insert(make_pair(contig.m_next_left, connected.begin()));
                        if(!rslt.second) {
                            TContigList::iterator other = rslt.first->second;
                            if(contig.EmptyLinker() && contig.m_left_link && !other->m_left_link && contig.m_next_right == other->LeftConnectingNode()) {
                                other->AddToLeft(contig); // add left link to other
                                connected.pop_front();
                                continue;
                            }else if(other->EmptyLinker() && other->m_left_link && !contig.m_left_link && other->m_next_right == contig.LeftConnectingNode()) {
                                contig.AddToLeft(*other); // add left link to contig    
                                rslt.first->second = connected.begin();
                                denied_right_nodes.erase(other->m_next_right);
                                connected.erase(other);
                            } else {
                                cerr << "Unexpected left fork: " << graph.GetNodeSeq(contig.m_next_left) << " " << contig.m_next_left << endl;
                                /*
                                cerr << "Other link: " << other->m_left_link << " " << other->m_right_link << endl;
                                cerr << "Contig link: " << contig.m_left_link << " " << contig.m_right_link << endl;

                                cerr << "Contig: " << contig.m_seq.size() << " " << contig.m_seq.m_circular << " ";
                                if(contig.m_next_right)
                                    cerr << contig.m_graph.GetNodeSeq(contig.m_next_right);
                                cerr << endl;
                                for(auto& chunk : contig.m_seq) {
                                    cerr << "Chunk: " << distance(chunk.begin(), chunk.end()) << endl;
                                    for(auto& seq : chunk) {
                                        for(char c : seq)
                                            cerr << c;
                                        cerr << endl;
                                    }
                                }

                                cerr << "Other: " << other->m_seq.size() << " " << other->m_seq.m_circular << " ";
                                if(other->m_next_right)
                                    cerr << other->m_graph.GetNodeSeq(other->m_next_right);
                                cerr << endl;
                                for(auto& chunk : other->m_seq) {
                                    cerr << "Chunk: " << distance(chunk.begin(), chunk.end()) << endl;
                                    for(auto& seq : chunk) {
                                        for(char c : seq)
                                            cerr << c;
                                        cerr << endl;
                                    }
                                }
                                */

                            }
                        }
                    }
                    if(contig.m_next_right) {
                        auto rslt = denied_right_nodes.insert(make_pair(contig.m_next_right, connected.begin()));
                        if(!rslt.second) {
                            TContigList::iterator other = rslt.first->second;
                            if(contig.EmptyLinker() && contig.m_right_link && !other->m_right_link && contig.m_next_left == other->RightConnectingNode()) {
                                other->AddToRight(contig); // add right link to other
                                denied_left_nodes.erase(contig.m_next_left);
                                connected.pop_front();
                            } else if (other->EmptyLinker() && other->m_right_link &&  !contig.m_right_link && other->m_next_left == contig.RightConnectingNode()) {
                                contig.AddToRight(*other); // add right link to contig
                                rslt.first->second = connected.begin();
                                denied_left_nodes.erase(other->m_next_left);
                                connected.erase(other);
                            } else {
                                cerr << "Unexpected right fork: " << graph.GetNodeSeq(contig.m_next_right) << " " << contig.m_next_right << endl;                                
                            }
                        }
                    }
                }
            }

            for(SContig& contig : connected) {
                if(contig.EmptyLinker())
                    continue; 

                if(contig.m_next_right)
                    denied_right_nodes.erase(contig.m_next_right);
                if(contig.m_next_left)
                    denied_left_nodes.erase(contig.m_next_left);
                bool keep_doing = true;
                while(keep_doing) {
                    keep_doing = false;
                    if(contig.m_next_right) {
                        CDBGraph::Node rnode = contig.RightConnectingNode();
                        auto rslt = denied_left_nodes.find(rnode);
                        if(rslt != denied_left_nodes.end()) {
                            keep_doing = true;
                            SContig& rcontig = *rslt->second;
                            if(rcontig.m_next_right) 
                                denied_right_nodes.erase(rcontig.m_next_right);
                            contig.AddToRight(rcontig);
                            connected.erase(rslt->second);
                            denied_left_nodes.erase(rslt);
                        } else if((rslt = denied_right_nodes.find(CDBGraph::ReverseComplement(rnode))) != denied_right_nodes.end()) {
                            keep_doing = true;
                            SContig& rcontig = *rslt->second;
                            if(rcontig.m_next_left) 
                                denied_left_nodes.erase(rcontig.m_next_left);
                            rcontig.ReverseComplement();
                            contig.AddToRight(rcontig);
                            connected.erase(rslt->second);
                            denied_right_nodes.erase(rslt);
                        }
                    }               
                    if(contig.m_next_left) {
                        CDBGraph::Node lnode = contig.LeftConnectingNode();
                        auto rslt = denied_right_nodes.find(lnode);
                        if(rslt != denied_right_nodes.end()) {
                            keep_doing = true;
                            SContig& lcontig = *rslt->second;
                            if(lcontig.m_next_left) 
                                denied_left_nodes.erase(lcontig.m_next_left);
                            contig.AddToLeft(lcontig);
                            connected.erase(rslt->second);
                            denied_right_nodes.erase(rslt);
                        } else if((rslt = denied_left_nodes.find(CDBGraph::ReverseComplement(lnode))) != denied_left_nodes.end()) {
                            keep_doing = true;
                            SContig& lcontig = *rslt->second;
                            if(lcontig.m_next_right) 
                                denied_right_nodes.erase(lcontig.m_next_right);
                            lcontig.ReverseComplement();
                            contig.AddToLeft(lcontig);
                            connected.erase(rslt->second);
                            denied_left_nodes.erase(rslt);
                        }
                    }                                
                }
            
                if(contig.m_next_right == contig.LeftConnectingNode() && (int)contig.LenMax() >= 2*graph.KmerLen()-1)  // circular and not very short  
                    contig.RotateCircularToMinKmer();                     
            }
            
            total = 0;
            len = 0;
            for(auto& seq : connected) {
                ++total;
                len += seq.LenMax()-graph.KmerLen()+1;
            }            
            
            cerr << "Fragments after: " << total << " " <<  len << endl;

            return connected;                
        }

        // connects and extends contigs from previous iteration using a longer kmer
        // scontigs - previous contigs
        // extensions - connectors and extenders produced by longer kmer
        
        static void ConnectAndExtendContigs(TContigList& scontigs, TContigList& extensions) {
            if(scontigs.empty())
                return;

            int kmer_len = scontigs.front().m_kmer_len;
            typedef unordered_map<SContig*, SContig*> TExtensionsMap; // connections to left contig sides   (pointer to contig; pointer to connection)
            TExtensionsMap left_connections;
            TExtensionsMap right_connections;
            TExtensionsMap left_extensions;
            TExtensionsMap right_extensions;
            int connectors = 0;
            int extenders = 0;

            for(auto& ex : extensions) {
                if(ex.m_left_link && ex.m_right_link) {
                    ++connectors;

                    if(ex.m_left_shift < 0)
                        left_connections[ex.m_left_link] = &ex;
                    else
                        right_connections[ex.m_left_link] = &ex;
                    if(ex.m_right_shift < 0)
                        left_connections[ex.m_right_link] = &ex;
                    else
                        right_connections[ex.m_right_link] = &ex;
                } else if(ex.m_left_link) {
                    ++extenders;

                    if(ex.m_left_shift < 0)
                        left_extensions[ex.m_left_link] = &ex;
                    else
                        right_extensions[ex.m_left_link] = &ex;
                } else if(ex.m_right_link) {
                    ++extenders;

                    if(ex.m_right_shift < 0)
                        left_extensions[ex.m_right_link] = &ex;
                    else
                        right_extensions[ex.m_right_link] = &ex;
                }
            
            }
            cerr << "Connectors: " << connectors << " Extenders: " << extenders << endl;

            for(auto& contig : scontigs)
                contig.m_is_taken = 0;

            for(auto& contig : scontigs) {
                if(contig.m_is_taken)
                    continue;

                bool circular = false;
                for(int p = 0; p < 2; ++p) {
                    bool plus = (p == 0);
                    SContig* fragment = &contig;

                    while(true) {
                        TExtensionsMap::iterator rslt;
                        // check if connection to other contigs is possible
                        if((plus && (rslt = right_connections.find(fragment)) != right_connections.end()) ||
                           (!plus && (rslt = left_connections.find(fragment)) != left_connections.end())) {
                        
                            SContig* connector = rslt->second;
                            if(connector->m_right_link == fragment) { // either reversed or circular            
                                if(CDBGraph::ReverseComplement(connector->m_next_right) == contig.RightConnectingNode())
                                    connector->ReverseComplement();
                            }
                            if(connector->m_left_link != fragment || contig.RightConnectingNode() != connector->m_next_left)
                                cerr << "Corrupted connectionA" << endl;

                            /*
                            if(test_contig) {
                                cerr << "ConnectorA link: " << connector->m_left_link << " " << fragment << endl;

                                cerr << "Contig: " << contig.m_seq.size() << " " << contig.m_seq.m_circular << " ";
                                if(contig.m_next_right)
                                    cerr << contig.m_graph.GetNodeSeq(contig.m_next_right);
                                cerr << endl;
                                for(auto& chunk : contig.m_seq) {
                                    cerr << "Chunk: " << distance(chunk.begin(), chunk.end()) << endl;
                                    for(auto& seq : chunk) {
                                        for(char c : seq)
                                            cerr << c;
                                        cerr << endl;
                                    }
                                }

                                cerr << "Connector: " << connector->m_seq.size() << " " << connector->m_seq.m_circular << " ";
                                if(connector->m_next_left)
                                    cerr << connector->m_graph.GetNodeSeq(connector->m_next_left);
                                cerr << endl;
                                for(auto& chunk : connector->m_seq) {
                                    cerr << "Chunk: " << distance(chunk.begin(), chunk.end()) << endl;
                                    for(auto& seq : chunk) {
                                        for(char c : seq)
                                            cerr << c;
                                        cerr << endl;
                                    }
                                }

                            }
                            */

                            contig.AddToRight(*connector);

                            fragment = connector->m_right_link;
                            if(fragment->m_is_taken)      // don't connect already used contig (this is result of multiple connection)          
                                break;
                        
                            fragment->m_is_taken = 1;     // fragment will be removed
                            plus = connector->m_right_shift < 0;
                            if(!plus)
                                fragment->ReverseComplement();

                            if(fragment->LeftConnectingNode() != connector->m_next_right)
                                cerr << "Corrupted connectionB" << endl;
                            
                            /*
                            if(test_contig) {
                                cerr << "ConnectorB link: " << connector->m_right_link << " " << fragment << endl;

                                cerr << "Contig: " << contig.m_seq.size() << " " << contig.m_seq.m_circular << " ";
                                if(contig.m_next_right)
                                    cerr << contig.m_graph.GetNodeSeq(contig.m_next_right);
                                cerr << endl;
                                for(auto& chunk : contig.m_seq) {
                                    cerr << "Chunk: " << distance(chunk.begin(), chunk.end()) << endl;
                                    for(auto& seq : chunk) {
                                        for(char c : seq)
                                            cerr << c;
                                        cerr << endl;
                                    }
                                }

                                cerr << "Fragment: " << fragment->m_seq.size() << " " << fragment->m_seq.m_circular << " ";
                                if(fragment->m_next_left)
                                    cerr << fragment->m_graph.GetNodeSeq(fragment->m_next_left);
                                cerr << endl;
                                for(auto& chunk : fragment->m_seq) {
                                    cerr << "Chunk: " << distance(chunk.begin(), chunk.end()) << endl;
                                    for(auto& seq : chunk) {
                                        for(char c : seq)
                                            cerr << c;
                                        cerr << endl;
                                    }
                                }

                            }
                            */

                            circular = (fragment == &contig);

                            if(!circular) {  // not circular 
                                contig.AddToRight(*fragment);
                                continue;
                            } else if((int)contig.LenMax() >= 2*kmer_len-1) { //stabilize circular contig            
                                contig.RotateCircularToMinKmer();
                                break;                        
                            }
                        } 
                        if((plus && (rslt = right_extensions.find(fragment)) != right_extensions.end()) ||
                           (!plus && (rslt = left_extensions.find(fragment)) != left_extensions.end())) {

                            SContig* extender = rslt->second;
                            if((int)extender->LenMax() >= kmer_len) {
                                if(extender->m_right_link && extender->m_right_link == fragment)
                                    extender->ReverseComplement();
                                //                                if(extender->m_left_link != fragment || contig.BackKmer() != extender->m_next_left)
                                if(extender->m_left_link != fragment || contig.RightConnectingNode() != extender->m_next_left) 
                                    cerr << "Corrupted extension" << endl;

                                contig.AddToRight(*extender);

                                /*
                                if(test_contig) {
                                    cerr << "Ex links: " << extender->m_left_link << " " << extender->m_right_link << " " << fragment << endl;
                                    cerr << "Extender: " << extender->m_seq.size() << " " << extender->m_seq.m_circular << " ";
                                    cerr << (extender->m_next_left ? extender->m_graph.GetNodeSeq(extender->m_next_left) : "-");
                                    cerr << (extender->m_next_right ? extender->m_graph.GetNodeSeq(extender->m_next_right) : "-");
                                    cerr << endl;
                                    for(auto& chunk : extender->m_seq) {
                                        cerr << "Chunk:     " << distance(chunk.begin(), chunk.end()) << endl;
                                        for(auto& seq : chunk) {
                                            for(char c : seq)
                                                cerr << c;
                                            cerr << endl;
                                        }
                                    }

                                    cerr << "Contig links: " << contig.m_left_link << " " << contig.m_right_link << " " << fragment << endl;
                                    cerr << "Contig: " << contig.m_seq.size() << " " << contig.m_seq.m_circular << " ";
                                    cerr << (contig.m_next_left ? contig.m_graph.GetNodeSeq(contig.m_next_left) : "-") << " ";
                                    cerr << (contig.m_next_right ? contig.m_graph.GetNodeSeq(contig.m_next_right) : "-");
                                    cerr << endl;
                                    for(auto& chunk : contig.m_seq) {
                                        cerr << "Chunk:     " << distance(chunk.begin(), chunk.end()) << endl;
                                        for(auto& seq : chunk) {
                                            for(char c : seq)
                                                cerr << c;
                                            cerr << endl;
                                        }
                                    }
                                }
                                */

                            }
                        }                        

                        break;
                    }
                    contig.m_is_taken = 2;  // final contig will be kept
                    if(circular)
                        break;
                    contig.ReverseComplement();
                } 
                //clip flanks which are not 'double' checked 

                auto& graph = contig.m_graph;

                for(int low_abundance_clip = 10; low_abundance_clip > 0 && contig.m_left_extend > 0; --low_abundance_clip) {
                    auto kmer = contig.FrontKmer();
                    if(!kmer || graph.Abundance(kmer) > 5)
                        break;
                    
                    contig.ClipLeft(1);
                }
                int left_clip = min(kmer_len,contig.m_left_extend);
                contig.ClipLeft(left_clip);
                if(contig.m_left_extend > 0)
                    contig.m_seq.m_left_repeat = min(kmer_len-1, contig.m_left_extend+contig.m_seq.m_left_repeat);

                for(int low_abundance_clip = 10; low_abundance_clip > 0 && contig.m_right_extend > 0; --low_abundance_clip) {
                    auto kmer = contig.BackKmer();
                    if(!kmer || graph.Abundance(kmer) > 5)
                        break;
                    
                    contig.ClipRight(1);
                }
                int right_clip = min(kmer_len,contig.m_right_extend); 
                contig.ClipRight(right_clip);
                if(contig.m_right_extend > 0)
                    contig.m_seq.m_right_repeat = min(kmer_len-1, contig.m_right_extend+contig.m_seq.m_right_repeat);
            }

            //remove fragments; stabilize orientation and order which are random in multithreading          
            for(auto iloop = scontigs.begin(); iloop != scontigs.end(); ) {
                auto ic = iloop++;
                if(ic->m_is_taken != 2 || (int)ic->LenMin() < kmer_len)
                    scontigs.erase(ic);
                else
                    ic->SelectMinDirection();
            }
            scontigs.sort();
        }        

        CContigSequence m_seq;               // sequence

        CDBGraph::Node m_next_left = 0;      // denied left kmer (connection possible but it is already owned)
        CDBGraph::Node m_next_right = 0;     // denied right kmer (connection possible but it is already owned)

        SContig* m_left_link = nullptr;      // if set points to 'left' contig
        int m_left_shift = 0;                // shift+1 for m_next_left in this contig (positive for the right end)
        SContig* m_right_link = nullptr;     // if set points to 'right' contig
        int m_right_shift = 0;               // shift+1 for m_next_right in this contig (positive for the right end)

        int m_left_extend = 0;               // number of newly assembled bases which could be clipped
        int m_right_extend = 0;              // number of newly assembled bases which could be clipped

        CDBGraph& m_graph;
        int m_kmer_len;
        SAtomic<uint8_t> m_is_taken = 0;
    };


    // This is a very lightweight class holding a reference to de Bruijn graph and main assembling parameters
    // It provides function used in assembling
    class CDBGraphDigger {
    public:

        CDBGraphDigger(CDBGraph& graph, double fraction, int jump, int low_count, bool allow_snps = false) : m_graph(graph), m_fraction(fraction), m_jump(jump), m_hist_min(graph.HistogramMinimum()), m_low_count(low_count), m_allow_snps(allow_snps) { 
            m_max_branch = 200; // maximum number of paths explored before quitting
        }

    private:
        typedef tuple<TStrList::iterator,int> TContigEnd;

    public:

        // starting from a node, find an extension of len l with maximal abundance
        string MostLikelyExtension(CDBGraph::Node node, unsigned len) const { //don't do FilterNeighbors because it is called in it     
            string s;
            while(s.size() < len) {
                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(node);
                if(successors.empty())
                    return s;            
                sort(successors.begin(), successors.end(), [&](const CDBGraph::Successor& a, const CDBGraph::Successor& b) {return m_graph.Abundance(a.m_node) > m_graph.Abundance(b.m_node);}); 
                node = successors[0].m_node;
                s.push_back(successors[0].m_nt);
            }
            return s;
        }            

        string MostLikelySeq(CDBGraph::Successor base, unsigned len) const {
            string s(1, base.m_nt);
            return s+MostLikelyExtension(base.m_node, len-1);
        }

        // starting from a node, find an extension of len l without forks (simple path); returns true if hit dead end
        pair<string, bool> StringentExtension(CDBGraph::Node node, unsigned len) const {
            string s;
            while(s.size() < len) {
                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(node);
                FilterNeighbors(successors, false);
                if(successors.empty())
                    return make_pair(s, true);
                if(successors.size() != 1)
                    return make_pair(s, false);            
                node = successors[0].m_node;
                s.push_back(successors[0].m_nt);
            }
            return make_pair(s, false);
        }

        bool ExtendableSuccessor(const CDBGraph::Successor& initial_suc) const {
            int kmer_len = m_graph.KmerLen();
            int total_len = max(100, kmer_len);

            unordered_map<CDBGraph::Node,int> node_len;
            node_len.emplace(initial_suc.m_node,0);

            stack<pair<CDBGraph::Node,int>> active_nodes;
            active_nodes.emplace(initial_suc.m_node,0);

            while(!active_nodes.empty()) {

                CDBGraph::Node node = active_nodes.top().first;
                int len = active_nodes.top().second;
                active_nodes.pop();
                
                if(len == kmer_len) {
                    vector<CDBGraph::Successor> step_back = m_graph.GetNodeSuccessors(m_graph.ReverseComplement(node));
                    FilterLowAbundanceNeighbors(step_back);
                    bool found = false;
                    for(auto& back : step_back) {
                        if(back.m_nt == Complement(initial_suc.m_nt)) {
                            found = true;
                            break;
                        }
                    }
                    if(!found)
                        continue;
                }                

                if(len == total_len)
                    return true;

                if(len > kmer_len) {
                    int& l = node_len[node];
                    if(len > l)
                        l = len;
                    else
                        continue;
                }
        
                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(node);
                FilterLowAbundanceNeighbors(successors);

                if(!successors.empty()) {
                    for(int i = successors.size()-1; i >= 0; --i)
                        active_nodes.emplace(successors[i].m_node, len+1);
                } 
            }

            return false;
        }

        vector<CDBGraph::Successor> GetReversibleNodeSuccessors(CDBGraph::Node node) const {
            vector<CDBGraph::Successor> neighbors = m_graph.GetNodeSuccessors(node);
            FilterNeighbors(neighbors, true);
            for(auto& neighbor : neighbors) {
                vector<CDBGraph::Successor> step_back = m_graph.GetNodeSuccessors(m_graph.ReverseComplement(neighbor.m_node));
                FilterNeighbors(step_back, true);
                bool found = false;
                for(auto& back : step_back) {
                    if(back.m_node == m_graph.ReverseComplement(node)) {
                        found = true;
                        break;
                    }
                }
                if(!found) {
                    neighbors.clear();
                    return neighbors;
                }
            }
            return neighbors;
        }


        bool GoodNode(const CDBGraph::Node& node) const { return m_graph.Abundance(node) >= m_low_count; }
        int HistMin() const { return m_hist_min; }

        // removes noise forks
        void FilterLowAbundanceNeighbors(vector<CDBGraph::Successor>& successors) const {
            // low abundance forks
            if(successors.size() > 1) {
                int abundance = 0;
                for(auto& suc : successors) {
                    abundance += m_graph.Abundance(suc.m_node);
                }
                sort(successors.begin(), successors.end(), [&](const CDBGraph::Successor& a, const CDBGraph::Successor& b) {return m_graph.Abundance(a.m_node) > m_graph.Abundance(b.m_node);});
                for(int j = successors.size()-1; j > 0 && m_graph.Abundance(successors.back().m_node) <= m_fraction*abundance; --j) 
                    successors.pop_back();            
            }

            // strand specific noise reduction for Illumina issue of GGT->GG[ACG]
            if(m_graph.GraphIsStranded() && successors.size() > 1) {

                double fraction = 0.1*m_fraction;
            
                int target = -1;
                for(int j = 0; target < 0 && j < (int)successors.size(); ++j) {
                    if(m_graph.GetNodeSeq(successors[j].m_node).substr(m_graph.KmerLen()-3) == "GGT") 
                        target = j;
                }
                if(target >= 0) {
                    int abundance = m_graph.Abundance(successors[target].m_node);
                    if(abundance > 5) {
                        double am = abundance*(1-m_graph.PlusFraction(successors[target].m_node));
                        for(int j = 0; j < (int)successors.size(); ) {
                            if(m_graph.Abundance(successors[j].m_node)*(1-m_graph.PlusFraction(successors[j].m_node)) < fraction*am)
                                successors.erase(successors.begin()+j);
                            else
                                ++j;
                        }
                    }
                }
            }

        }
        
        void FilterNeighbors(vector<CDBGraph::Successor>& successors, bool check_extension) const {
            // low abundance forks
            FilterLowAbundanceNeighbors(successors);

            //not extendable forks
            if(check_extension && successors.size() > 1 && m_graph.Abundance(successors.front().m_node) > 5) {
                for(int i = 0; i < (int)successors.size(); ) {
                    if(ExtendableSuccessor(successors[i]))
                        ++i;
                    else
                        successors.erase(successors.begin()+i);
                }
            }

            // strand specific noise reduction for Illumina issue of GGT->GG[ACG] for negative strand and low coverage (the prev loop didn't work)
            if(m_graph.GraphIsStranded() && successors.size() > 1 && (!check_extension || m_graph.Abundance(successors.front().m_node) <= 5)) {

                double fraction = 0.1*m_fraction;
                int target = -1;

                for(int j = 0; target < 0 && j < (int)successors.size(); ++j) {
                    if(MostLikelySeq(successors[j], 3) == "ACC") 
                        target = j;
                }
                if(target >= 0) {
                    int abundance =  m_graph.Abundance(successors[target].m_node);
                    if(abundance > 5) {
                        double ap = abundance*m_graph.PlusFraction(successors[target].m_node);
                        for(int j = 0; j < (int)successors.size(); ) {
                            if(m_graph.Abundance(successors[j].m_node)*m_graph.PlusFraction(successors[j].m_node) < fraction*ap)
                                successors.erase(successors.begin()+j);
                            else
                                ++j;
                        }
                    }
                }
            }

            // strand balance issue
            if(m_graph.GraphIsStranded() && successors.size() > 1) {

                double fraction = 0.1*m_fraction;
                                      
                bool has_both = false;
                for(int j = 0; !has_both && j < (int)successors.size(); ++j) {
                    double plusf = m_graph.PlusFraction(successors[j].m_node);
                    double minusf = 1.- plusf;
                    has_both = GoodNode(successors[j].m_node) && (min(plusf,minusf) > 0.25);
                }

                if(has_both) {
                    for(int j = 0; j < (int)successors.size(); ) {
                        double plusf = m_graph.PlusFraction(successors[j].m_node);
                        double minusf = 1.- plusf;
                        if(min(plusf,minusf) < fraction*max(plusf,minusf))
                            successors.erase(successors.begin()+j);
                        else
                            ++j;
                    }
                } 
            }
        }

        CDBGraph& Graph() { return m_graph; }

        enum EConnectionStatus {eSuccess, eNoConnection, eAmbiguousConnection};

        struct SElement {
            SElement (CDBGraph::Successor suc, SElement* link) : m_link(link), m_suc(suc) {}
            struct SElement* m_link;   // previous element  
            CDBGraph::Successor m_suc;
        };
        // connects two nodes in a finite number of steps
        pair<TBases, EConnectionStatus> ConnectTwoNodes(const CDBGraph::Node& first_node, const CDBGraph::Node& last_node, int steps) const {
        
            pair<TBases, EConnectionStatus> bases(TBases(), eNoConnection);

            deque<SElement> storage; // will contain ALL extensions (nothing is deleted)    
            typedef unordered_map<CDBGraph::Node,SElement*> TElementMap;  //pointer to its own element OR zero if ambiguous path    
            TElementMap current_elements;

            vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(first_node);
            FilterNeighbors(successors, false);
            for(auto& suc : successors) {
                storage.push_back(SElement(suc, 0));
                current_elements[suc.m_node] = &storage.back();
            }

            list<SElement> connections;
            for(int step = 1; step < steps && !current_elements.empty(); ++step) {
                TElementMap new_elements;
                for(auto& el : current_elements) {
                    vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(el.first);
                    FilterNeighbors(successors, false);
                    if(el.second == 0) {  // ambiguous path 
                        for(auto& suc : successors) {
                            new_elements[suc.m_node] = 0;
                            if(suc.m_node == last_node) {
                                bases.second = eAmbiguousConnection;
                                return bases;
                            }
                        }
                    } else {
                        for(auto& suc : successors) {
                            storage.push_back(SElement(suc, el.second));
                            if(suc.m_node == last_node) {
                                if(!connections.empty()) {
                                    bases.second = eAmbiguousConnection;
                                    return bases;
                                } else {
                                    connections.push_back(storage.back()); 
                                }                           
                            }
                            pair<TElementMap::iterator, bool> rslt = new_elements.insert(make_pair(suc.m_node, &storage.back()));
                            if(!rslt.second || !GoodNode(suc.m_node))
                                rslt.first->second = 0;
                        }
                    }                    
                }
                swap(current_elements, new_elements);
                if(current_elements.size() > m_max_branch)
                    return bases;
            }

            if(connections.empty())
                return bases;

            SElement el = connections.front();
            while(el.m_link != 0) {
                bases.first.push_front(el.m_suc);
                el = *el.m_link;
            }
            bases.first.push_front(el.m_suc);
            bases.second = eSuccess;
            return bases;
        }

        typedef list<TBases> TBasesList;
        typedef unordered_map<CDBGraph::Node, forward_list<TBasesList::iterator>> TBranch;  // all 'leaves' will have the same length  

        // makes one step extension for DiscoverSNPs()
        // sequences - keep the actual assembled sequences
        // branch - a lightweight map of the last kmer to the sequence
        // indels - allowed diff in length
        void OneStepBranchExtend(TBranch& branch, TBasesList& sequences, unsigned indels) {
            TBranch new_branch;
            for(auto& leaf : branch) {
                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(leaf.first);
                FilterNeighbors(successors, true);
                if(successors.empty()) {
                    for(auto is : leaf.second)
                        sequences.erase(is);
                    continue;
                }
                for(int i = successors.size()-1; i >= 0; --i) {
                    CDBGraph::Node& node = successors[i].m_node;
                    auto& lst = new_branch[node];
                    for(auto is : leaf.second) {
                        if(i > 0) {  // copy sequence if it is a fork               
                            sequences.push_front(*is);
                            is = sequences.begin();
                        }
                        TBases& bases = *is;
                        bases.push_back(successors[i]);
                        lst.emplace_front(is);
                    }
                }
            }

            int kmer_len = m_graph.KmerLen();
            for(auto it_loop = new_branch.begin(); it_loop != new_branch.end(); ) {
                auto it = it_loop++;
                bool merged = false;
                for(auto jt = new_branch.begin(); !merged && jt != new_branch.end(); ++jt) {
                    if(it == jt)
                        continue;

                    for(unsigned k = 1; !merged && k <= indels; ++k) {
                        bool all_same = true;
                        auto& jlst = jt->second;
                        auto& jseq = *jlst.front();
                        for(auto& is : jlst) {
                            if(is->size() < k+kmer_len || *(is->end()-k-1) != *(jseq.end()-k-1)) {
                                all_same = false;
                                break;
                            }                            
                        }
                        if(!all_same)
                            break;
                        
                        if(it->first == (jseq.end()-k-1)->m_node) {
                            auto& ilst = it->second;
                            for(auto is : ilst) {
                                is->insert(is->end(), jseq.end()-k, jseq.end());
                                jlst.push_front(is);
                            }
                            new_branch.erase(it);
                            merged = true;
                        }
                    }
                }
            }

            swap(branch, new_branch);
        }

        // for a given fork (successors.size() > 1) perform no more than max_extent extension steps allowinf snps until only one path is left
        // returns the extension if successful, returns empty sequence otherwise
        TBasesList DiscoverSNPs(const vector<CDBGraph::Successor>& successors, int max_extent, int min_extent, int indels) {
            TBranch extensions;
            TBasesList sequences; // assembled seqs 

            if(max_extent == 0)
                return sequences;

            for(auto& suc : successors) {
                sequences.emplace_front(1,suc);
                extensions[suc.m_node].emplace_front(sequences.begin());         
            }

            while(!extensions.empty() && extensions.size() < m_max_branch) {
                int len = 0;
                for(auto& seq : sequences)
                    len = max(len, (int)seq.size());
                if(len >= max_extent)
                    break;
            
                OneStepBranchExtend(extensions, sequences, indels);
                if(extensions.empty()) { // can't extend 
                    sequences.clear();
                    return sequences;
                }

                if(extensions.size() == 1 && len+1 >= min_extent)
                    break;
            }

            if(extensions.size() == 1 && sequences.size() > 1) { // found snp
                // clip extra matches from the end
                int matches = 0;
                while(true) {
                    bool all_same = true;
                    for(auto& seq : sequences) {
                        if(matches == (int)seq.size() || (seq.end()-matches-1)->m_nt != (sequences.front().end()-matches-1)->m_nt) {
                            all_same = false;
                            break;
                        }
                    }
                    if(all_same)
                        ++matches;
                    else
                        break;
                }
                int kmer_len = m_graph.KmerLen();
                if(matches > kmer_len) {
                    for(auto& seq : sequences)
                        seq.erase(seq.end()-(matches-kmer_len), seq.end());
                }

                return sequences;            
            }

            sequences.clear();
            return sequences;
        }


        // starting from initial_node assembles the right extension        
        tuple<SContig, CDBGraph::Node, int> ExtendToRight(const CDBGraph::Node& initial_node, int allowed_intrusion) { // initial_node may be not owned 
            CDBGraph::Node node = initial_node;
            SContig extension(m_graph);
            int max_extent = m_jump;
            int kmer_len = m_graph.KmerLen();

            int initial_node_intrusion = 0;

            while(true) {
                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(node);                    
                FilterNeighbors(successors, true);
                if(successors.empty())                    // no extensions 
                    break;                                             

                int allowed_indels = 15;

                TBasesList step;
                if(successors.size() == 1) {  // simple extension  
                    step.emplace_back(1, successors.front());
                } else if(m_allow_snps) {     // try snps
                    step = DiscoverSNPs(successors, max_extent, 0, allowed_indels);  //TODO make sure that rhe branch is at the same point
                    if(step.empty())        // no snp
                        break;
                    
                } else {
                    break;                    
                }

                bool all_good = true;
                for(auto& stp : step) {
                    for(auto& s : stp) {
                        if(!GoodNode(s.m_node)) {
                            all_good = false;
                            break;
                        }
                    }
                }
                if(!all_good)
                    break;

                int step_size = 0;  // not all step seqs have same length
                for(auto& var : step)
                    step_size = max(step_size, (int)var.size());                

                CDBGraph::Node rev_node = CDBGraph::ReverseComplement(step.front().back().m_node);  // all step seqs have same last kmer
                vector<CDBGraph::Successor> predecessors = m_graph.GetNodeSuccessors(rev_node);
                FilterNeighbors(predecessors, true);
                if(predecessors.empty())                            // no extensions   
                    break;

                if(step_size == 1 && predecessors.size() == 1) {
                    if(CDBGraph::ReverseComplement(predecessors[0].m_node) != node) // no return
                        break;

                    if(!m_graph.SetVisited(successors.front().m_node))  // node is taken
                        return make_tuple(extension, successors.front().m_node, initial_node_intrusion);                    

                    // normal extension   
                    if(extension.m_seq.empty()) {
                        extension.m_seq.InsertNewChunk();
                        extension.m_seq.InsertNewVariant();
                    }
                    extension.m_seq.ExtendTopVariant(successors.front().m_nt);

                    node = successors.front().m_node;
                    continue;
                }
                
                if(step_size == 1 && predecessors.size() > 1) {           // end of unique seq before repeat

                    /* bad idea - randomly releases nodes after clip
                    bool returned = false;
                    for(auto& pred : predecessors) {                        // check if branch doesn't come back
                        if(CDBGraph::ReverseComplement(pred.m_node) == node)
                            returned = true;
                    }
                    if(!returned) 
                        break;                    

                    if(extension.m_seq.empty()) {
                        extension.m_seq.InsertNewChunk();
                        extension.m_seq.InsertNewVariant();
                    }
                    extension.m_seq.ExtendTopVariant(successors.front().m_nt);     // keep one extra base
                    */
                    
                    break; 
                }

                if(!m_allow_snps)
                    break;
                
                // only snps after this

                TBasesList step_back = DiscoverSNPs(predecessors, max_extent, step_size, allowed_indels);

                if(step.size() != step_back.size())
                    break;

                int step_back_size = 0;  // not all step seqs have same length
                for(auto& var : step_back)
                    step_back_size = max(step_back_size, (int)var.size());

                if(step_back_size != step_size)
                    break;

                for(auto& seq : step_back) {
                    reverse(seq.begin(), seq.end());
                    for(auto& suc : seq) {
                        suc.m_nt = Complement(suc.m_nt); // not correct but it doesn't matter because we compare nodes  
                        suc.m_node = CDBGraph::ReverseComplement(suc.m_node);
                    }
                }
                if(step_back.front().front().m_node != node)
                    break;                

                for(auto& seq : step_back) {
                    seq.pop_front();
                    seq.push_back(step.front().back());
                }

                step.sort();
                step_back.sort();
                if(step != step_back)
                    break; 

                bool has_empty_variant = false;
                for(auto& seq : step) {
                    if((int)seq.size() == kmer_len) {
                        has_empty_variant = true;
                        break;
                    }
                }
                TLocalVariants extra(1);
                if(has_empty_variant) {
                    int last_chunk = extension.m_seq.empty() ? 0 : extension.m_seq.back().front().size();
                    string initial_node_seq;
                    if(extension.m_seq.size() <= 1 && last_chunk < kmer_len)
                        initial_node_seq = m_graph.GetNodeSeq(initial_node);
                    TLocalVariants seqs;
                    for(auto& seq : step) {
                        seqs.push_front(TVariation());
                        seqs.front().insert(seqs.front().end(), initial_node_seq.begin(), initial_node_seq.end());
                        if(!extension.m_seq.empty()) {
                            if(last_chunk <= kmer_len)
                                seqs.front().insert(seqs.front().end(), extension.m_seq.back().front().begin(), extension.m_seq.back().front().end());
                            else
                                seqs.front().insert(seqs.front().end(), extension.m_seq.back().front().end()-kmer_len, extension.m_seq.back().front().end());
                        }
                        for(unsigned i = 0; i < seq.size()-kmer_len; ++i)
                            seqs.front().push_back(seq[i].m_nt);
                    }

                    int shift = 0;
                    bool all_same = true;
                    while(all_same) {
                        for(auto& seq : seqs) {
                            if(shift == (int)seq.size() || *(seq.end()-shift-1) != *(seqs.front().end()-shift-1)) {
                                all_same = false;
                                break;
                            }
                        } 
                        if(all_same)
                            ++shift;
                    }                                    
                    if(shift > 0) {
                        if(shift >= kmer_len)
                            break;
                        if(shift > last_chunk) {
                            initial_node_intrusion = shift-last_chunk;
                            if(initial_node_intrusion > allowed_intrusion)
                                break;
                            extra.front().insert(extra.front().end(), initial_node_seq.end()-initial_node_intrusion, initial_node_seq.end());
                            shift = last_chunk;
                            /*
                            lock_guard<mutex> guard(out_mutex);
                            cerr << "ShiftedSNP intrusion: " << initial_node_intrusion << endl;
                            */
                        }
                        if(shift > 0) {
                            if(extension.m_seq.size() > 1 && last_chunk-shift < kmer_len) { // interfere with existing snp
                                int existing_snp = extension.m_seq.size()-2;
                                if(kmer_len+(int)extension.m_seq.ChunkLenMax(existing_snp)+step_size-shift > max_extent) // merged snp is too long
                                    break;
                                extra.clear();
                                for(auto& var : extension.m_seq[existing_snp]) {
                                    extra.push_front(var);
                                    extra.front().insert(extra.front().end(), extension.m_seq.back().front().begin(), extension.m_seq.back().front().end());
                                }
                                extension.m_seq.pop_back();
                                extension.m_seq.pop_back();
                                /*
                                lock_guard<mutex> guard(out_mutex);
                                cerr << "ShiftedSNP merge" << endl;
                                */
                            } else {
                                extra.front().insert(extra.front().end(), extension.m_seq.back().front().end()-shift, extension.m_seq.back().front().end());
                                extension.m_seq.back().front().resize(extension.m_seq.back().front().size()-shift);
                                if(extension.m_seq.back().front().empty())
                                    extension.m_seq.pop_back();
                                /*
                                lock_guard<mutex> guard(out_mutex);
                                cerr << "ShiftedSNP standalone" << endl;
                                */
                            }
                        }
                    }
                }
                

                auto last_node = step.front().back().m_node;
                char last_base = step.front().back().m_nt;
                bool my_snp = m_graph.SetVisited(last_node);                                

                extension.m_seq.InsertNewChunk();       // empty chunk for variable part
                for(auto& extra_seq : extra) {
                    for(auto& seq : step) {
                        extension.m_seq.InsertNewVariant(); // empty seq for new variant    
                        extension.m_seq.ExtendTopVariant(extra_seq.begin(), extra_seq.end());
                        for(unsigned i = 0; i < seq.size()-kmer_len; ++i)
                            extension.m_seq.ExtendTopVariant(seq[i].m_nt);
                    }
                }

                extension.m_seq.InsertNewChunk();   // empty chunk for matching kmer-1 bases    
                extension.m_seq.InsertNewVariant(); // empty seq for matching kmer-1 bases 
                for(auto i = step.front().end()-kmer_len; i != step.front().end()-1; ++i)
                    extension.m_seq.ExtendTopVariant(i->m_nt);

                if(my_snp) {
                    extension.m_seq.ExtendTopVariant(last_base);
                    for(auto& seq : step) {
                        for(unsigned i = 0; i < seq.size()-1; ++i)
                            m_graph.SetVisited(seq[i].m_node);
                    }

                    /*
                    if(!extra.front().empty()) {
                        lock_guard<mutex> guard(out_mutex);
                        cerr << "ShiftedSNP finished" << endl;
                        cerr << "Initial node: " << m_graph.GetNodeSeq(initial_node) << endl;
                        int num = 0;
                        for(auto& chunk : extension.m_seq) {
                            cerr << "Chunk" << ++num << endl;
                            for(auto& seq : chunk) {
                                for(char c : seq)
                                    cerr << c;
                                cerr << endl;
                            }
                        }
                    }
                    */
                } else {
                    /*
                    if(!extra.front().empty()) {
                        lock_guard<mutex> guard(out_mutex);
                        cerr << "ShiftedSNP for connection" << endl;
                        cerr << "Initial node: " << m_graph.GetNodeSeq(initial_node) << endl;
                        int num = 0;
                        for(auto& chunk : extension.m_seq) {
                            cerr << "Chunk" << ++num << endl;
                            for(auto& seq : chunk) {
                                for(char c : seq)
                                    cerr << c;
                                cerr << endl;
                            }
                        }
                    }
                    */

                    return make_tuple(extension, last_node, initial_node_intrusion);
                }                

                node = last_node;
            }

            return make_tuple(extension, 0, initial_node_intrusion);
        }        

        // assembles a contig starting from initial_node 
        // min_len - minimal length for accepted contigs
        // changes the state of all used nodes to 'visited' or 'temporary holding'   
        SContig GetContigForKmer(CDBGraph::Node initial_node, int min_len) {
            if(m_graph.Abundance(initial_node) < m_hist_min || !GoodNode(initial_node) || !m_graph.SetVisited(initial_node))
                return SContig(m_graph);

            //node is good and this thread owns it  

            // don't allow intrusion of snps in the initial kmer
            tuple<SContig, CDBGraph::Node, int> to_right = ExtendToRight(initial_node, 0);
            tuple<SContig, CDBGraph::Node, int> to_left = ExtendToRight(CDBGraph::ReverseComplement(initial_node), 0);

            SContig scontig(get<0>(to_left), get<0>(to_right), initial_node, CDBGraph::ReverseComplement(get<1>(to_left)), get<1>(to_right), m_graph);
            
            if(!scontig.m_next_left && !scontig.m_next_right && (int)scontig.LenMin() < min_len) {
                int kmer_len = m_graph.KmerLen();
                for(int i = scontig.m_seq.size()-1; i >= 0; i -= 2) {
                    if(i == (int)scontig.m_seq.size()-1) { // last chunk size >= kmer_len
                        CReadHolder rh(false);
                        rh.PushBack(scontig.m_seq.back().front());
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)
                            m_graph.SetTempHolding(m_graph.GetNode(*ik));
                    } else {
                        if((int)scontig.m_seq.ChunkLenMax(i) >= kmer_len) {
                            TVariation seq(scontig.m_seq[i].front().begin(), scontig.m_seq[i].front().end());
                            CReadHolder rh(false);
                            rh.PushBack(seq);
                            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)
                                m_graph.SetTempHolding(m_graph.GetNode(*ik));
                        }
                        for(auto& variant : scontig.m_seq[i+1]) {
                            TVariation seq(scontig.m_seq[i].front().end()-kmer_len+1, scontig.m_seq[i].front().end()); // all uniq chunks >= kmer_len-1
                            seq.insert(seq.end(), variant.begin(), variant.end());
                            seq.insert(seq.end(), scontig.m_seq[i+2].front().begin(), scontig.m_seq[i+2].front().begin()+kmer_len-1);
                            CReadHolder rh(false);
                            rh.PushBack(seq);
                            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)
                                m_graph.SetTempHolding(m_graph.GetNode(*ik));
                        }
                    }
                }

                return SContig(m_graph);
            } else {
                return scontig;
            }
        }

        void CheckRepeats(TContigList& scontigs) {
            int kmer_len = m_graph.KmerLen();

            for(auto it = scontigs.begin(); it != scontigs.end(); ++it) {
                auto& contig = it->m_seq;

                if(contig.m_left_repeat >= kmer_len && contig.m_left_repeat < (int)contig.LenMin()) {
                    int last_chunk = 0;
                    for(int len = contig.ChunkLenMin(last_chunk); len < contig.m_left_repeat+1; len += contig.ChunkLenMin(++last_chunk));

                    vector<forward_list<CDBGraph::Node>> kmers(contig.m_left_repeat+1-kmer_len+1);

                    stack<pair<TVariation*, int>> active_chunks;
                    active_chunks.emplace(&contig[0].front(), 0);
                    deque<TVariation*> current_seqs;
                    while(!active_chunks.empty()) {
                        TVariation* seqp = active_chunks.top().first;
                        int chunk_num = active_chunks.top().second;
                        active_chunks.pop();

                        current_seqs.resize(chunk_num);
                        current_seqs.push_back(seqp);
                        for(int chunk = chunk_num+1; chunk <= last_chunk; ++chunk) {
                            auto it = contig[chunk].begin();
                            current_seqs.push_back(&(*it));
                            for(++it; it != contig[chunk].end(); ++it) 
                                active_chunks.emplace(&(*it), chunk);
                        }
                        TVariation seq;
                        for(unsigned i = 0; i < current_seqs.size()-1; ++i)
                            seq.insert(seq.end(), current_seqs[i]->begin(), current_seqs[i]->end());
                        seq.insert(seq.end(), current_seqs.back()->begin(), current_seqs.back()->begin()+contig.m_left_repeat+1-seq.size());
                        CReadHolder rh(false);
                        rh.PushBack(seq);
                        int pos = kmers.size()-1;
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) {
                            CDBGraph::Node node = m_graph.GetNode(*ik);
                            if(find(kmers[pos].begin(), kmers[pos].end(), node) == kmers[pos].end())
                                kmers[pos].push_front(node);
                        }
                    }
                    
                    for( ; contig.m_left_repeat >= kmer_len; --contig.m_left_repeat) {
                        int p = contig.m_left_repeat-kmer_len;

                        bool bad_node = false;
                        for(auto& kmer : kmers[p]) {
                            if(!kmer || !GoodNode(kmer)) {
                                bad_node = true;
                                break;
                            }
                        }
                        if(bad_node)
                            break;
                                                 
                        bool no_step = false;
                        for(auto& kmer : kmers[p]) {
                            if(kmer) {
                                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(kmer);
                                FilterNeighbors(successors, true);
                                if(successors.empty()) {
                                    no_step = true;
                                    break;
                                }
                                auto& next_lst = kmers[p+1];
                                for(auto& suc : successors) {
                                    if(find_if(next_lst.begin(), next_lst.end(), [suc](const CDBGraph::Node& node) {return node == suc.m_node; }) == next_lst.end()) {
                                        no_step = true;
                                        break;
                                    }
                                }
                            }
                        }
                        if(no_step)
                            break;
                    }                    
                }
            

                if(contig.m_right_repeat >= kmer_len && contig.m_right_repeat < (int)contig.LenMin()) {
                    int first_chunk = contig.size()-1;
                    for(int len = contig.ChunkLenMin(first_chunk); len < contig.m_right_repeat+1; len += contig.ChunkLenMin(--first_chunk));

                    vector<forward_list<CDBGraph::Node>> kmers(contig.m_right_repeat+1-kmer_len+1);

                    stack<pair<TVariation*, int>> active_chunks;
                    active_chunks.emplace(&contig[contig.size()-1].front(), contig.size()-1);
                    deque<TVariation*> current_seqs;
                    while(!active_chunks.empty()) {
                        TVariation* seqp = active_chunks.top().first;
                        int chunk_num = active_chunks.top().second;
                        active_chunks.pop();

                        if(!current_seqs.empty())
                            current_seqs.erase(current_seqs.begin(), current_seqs.begin()+chunk_num-first_chunk+1);
                        current_seqs.push_front(seqp);
                        for(int chunk = chunk_num-1; chunk >= first_chunk; --chunk) {
                            auto it = contig[chunk].begin();
                            current_seqs.push_front(&(*it));
                            for(++it; it != contig[chunk].end(); ++it) 
                                active_chunks.emplace(&(*it), chunk);
                        }
                        TVariation seq;
                        for(unsigned i = current_seqs.size()-1; i > 0; --i)
                            seq.insert(seq.begin(), current_seqs[i]->begin(), current_seqs[i]->end());
                        seq.insert(seq.begin(), current_seqs.front()->end()-(contig.m_right_repeat+1-seq.size()), current_seqs.front()->end());
                        CReadHolder rh(false);
                        rh.PushBack(seq);
                        int pos = kmers.size()-1;
                        for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) {
                            CDBGraph::Node node = m_graph.GetNode(*ik);
                            if(find(kmers[pos].begin(), kmers[pos].end(), node) == kmers[pos].end())
                                kmers[pos].push_front(node);
                        }
                    }
                    
                    for( ; contig.m_right_repeat >= kmer_len; --contig.m_right_repeat) {
                        int p = kmers.size()-(contig.m_right_repeat-kmer_len+1);

                        bool bad_node = false;
                        for(auto& kmer : kmers[p]) {
                            if(!kmer || !GoodNode(kmer)) {
                                bad_node = true;
                                break;
                            }
                        }
                        if(bad_node)
                            break;

                        bool no_step = false;
                        for(auto& kmer : kmers[p]) {
                            if(kmer) {
                                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(CDBGraph::ReverseComplement(kmer));
                                FilterNeighbors(successors, true);
                                if(successors.empty()) {
                                    no_step = true;
                                    break;
                                }
                                auto& prev_lst = kmers[p-1];
                                for(auto& suc : successors) {
                                    if(find_if(prev_lst.begin(), prev_lst.end(), [suc](const CDBGraph::Node& node) {return node == CDBGraph::ReverseComplement(suc.m_node); }) == prev_lst.end()) {
                                        no_step = true;
                                        break;
                                    }
                                }
                            }
                        }                        
                        if(no_step)
                            break;
                    }                    
                }
            }


            /*
            for(auto it = scontigs.begin(); it != scontigs.end(); ++it) {
                SContig& contig = *it;

                if(contig.m_seq.m_left_repeat >= kmer_len && contig.m_seq.m_left_repeat < (int)contig.m_seq.ChunkLenMax(0)) {
                    deque<CDBGraph::Node> left_kmers;
                    TVariation left(contig.m_seq[0].front().begin(), contig.m_seq[0].front().begin()+contig.m_seq.m_left_repeat+1);
                    CReadHolder rh(false);
                    rh.PushBack(left);
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik) {
                        CDBGraph::Node node = m_graph.GetNode(*ik);
                        left_kmers.push_front(node);
                    }

                    for( ; contig.m_seq.m_left_repeat >= kmer_len; --contig.m_seq.m_left_repeat) {
                        int p = contig.m_seq.m_left_repeat-kmer_len;
                        auto& kmer = left_kmers[p];
                        if(!kmer || !GoodNode(kmer))
                            break;
                        vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(kmer);
                        FilterNeighbors(successors, true);
                        if(successors.size() != 1 || successors[0].m_node != left_kmers[p+1])
                           break;
                    }

                }

                int last_chunk = contig.m_seq.size()-1;
                if(contig.m_seq.m_right_repeat >= kmer_len && contig.m_seq.m_right_repeat < (int)contig.m_seq.ChunkLenMax(last_chunk)) {
                    deque<CDBGraph::Node> right_kmers;
                    TVariation right(contig.m_seq[last_chunk].front().end()-contig.m_seq.m_right_repeat-1, contig.m_seq[last_chunk].front().end());
                    CReadHolder rh(false);
                    rh.PushBack(right);
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik) {
                        CDBGraph::Node node = m_graph.GetNode(*ik);
                        right_kmers.push_front(node);
                    }

                    for( ; contig.m_seq.m_right_repeat >= kmer_len; --contig.m_seq.m_right_repeat) {
                        int p = right_kmers.size()+kmer_len-1-contig.m_seq.m_right_repeat;
                        auto& kmer = right_kmers[p];
                        if(!kmer || !GoodNode(kmer))
                            break;
                        vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(CDBGraph::ReverseComplement(kmer));
                        if(successors.size() != 1 || CDBGraph::ReverseComplement(successors[0].m_node) != right_kmers[p-1])
                           break;
                    }
                }
            }
            */
        }

        void ConnectOverlappingContigs(TContigList& scontigs) {
            int kmer_len = m_graph.KmerLen();
            unordered_map<CDBGraph::Node, forward_list<pair<TContigList::iterator, int>>> kmers;
            for(auto it = scontigs.begin(); it != scontigs.end(); ++it) {
                SContig& contig = *it;
                if((int)contig.m_seq.ChunkLenMax(0) > kmer_len) {
                    CReadHolder rh(false);
                    rh.PushBack(contig.m_seq[0].front());
                    int pos = contig.m_seq.ChunkLenMax(0)-kmer_len;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) {
                        if(pos < (int)contig.m_seq.ChunkLenMax(0)-kmer_len) {
                            CDBGraph::Node node = m_graph.GetNode(*ik);
                            if(node)
                                kmers[node].emplace_front(it, pos);
                        }
                    }
                }
                if(contig.m_seq.size() > 1 && (int)contig.m_seq.ChunkLenMax(contig.m_seq.size()-1) > kmer_len) {
                    CReadHolder rh(false);
                    rh.PushBack(contig.m_seq[contig.m_seq.size()-1].front());
                    int pos = contig.m_seq.LenMax()-kmer_len;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) {
                        if(pos > (int)contig.m_seq.LenMax()-(int)contig.m_seq.ChunkLenMax(contig.m_seq.size()-1)) {
                            CDBGraph::Node node = m_graph.GetNode(*ik);
                            if(node)
                                kmers[node].emplace_front(it, pos);
                        }
                    }
                }
            }

            list<tuple<TContigList::iterator, TContigList::iterator, int, int, int>> overlaps; // first contig, second contig, start/end, start/end, len    

            for(auto it = scontigs.begin(); it != scontigs.end(); ++it) {
                SContig& icontig = *it;

                // right overlap    
                {
                    list<tuple<TContigList::iterator, TContigList::iterator, int, int, int>> contig_overlaps;

                    auto& irchunk = icontig.m_seq.back().front();
                    auto rslt = kmers.find(icontig.BackKmer()); // rightend to left end 
                    if(rslt != kmers.end()) {
                        for(auto& hit : rslt->second) {
                            auto jt = hit.first;
                            if(jt == it)
                                continue;
                            int overlap_len = hit.second+kmer_len;
                            auto& jlchunk = jt->m_seq.front().front();
                            if(overlap_len > (int)irchunk.size() || overlap_len > (int)jlchunk.size())
                                continue;
                            if(!equal(jlchunk.begin(), jlchunk.begin()+hit.second, irchunk.end()-overlap_len))
                                continue;

                            contig_overlaps.emplace_back(it, jt, 1, -1, overlap_len);
                        }
                    }
                    rslt = kmers.find(CDBGraph::ReverseComplement(icontig.BackKmer())); // right end to right end   
                    if(rslt != kmers.end()) {
                        for(auto& hit : rslt->second) {
                            auto jt = hit.first;
                            if(jt == it)
                                continue;
                            int overlap_len = jt->m_seq.LenMax()-hit.second;
                            auto& jrchunk = jt->m_seq.back().front();
                            if(overlap_len > (int)irchunk.size() || overlap_len > (int)jrchunk.size())
                                continue;
                            TVariation seq(irchunk.end()-overlap_len, irchunk.end()-kmer_len);
                            ReverseComplementSeq(seq.begin(), seq.end());
                            if(!equal(seq.begin(), seq.end(), jrchunk.end()-overlap_len+kmer_len))
                                continue;

                            contig_overlaps.emplace_back(it, jt, 1, 1, overlap_len);
                        }
                    }
                    if(contig_overlaps.size() == 1)
                        overlaps.splice(overlaps.end(), contig_overlaps);
                }

                //left overlap  
                {
                    list<tuple<TContigList::iterator, TContigList::iterator, int, int, int>> contig_overlaps;

                    auto& ilchunk = it->m_seq.front().front();
                    auto rslt = kmers.find(icontig.FrontKmer()); // left end to right end   
                    if(rslt != kmers.end()) {
                        for(auto& hit : rslt->second) {
                            auto jt = hit.first;
                            if(jt == it)
                                continue;
                            int overlap_len = jt->m_seq.LenMax()-hit.second;
                            auto& jrchunk = jt->m_seq.back().front();
                            if(overlap_len > (int)ilchunk.size() || overlap_len > (int)jrchunk.size())
                                continue;
                            if(!equal(jrchunk.end()-overlap_len+kmer_len, jrchunk.end(), ilchunk.begin()+kmer_len))
                                continue;

                            contig_overlaps.emplace_back(it, jt, -1, 1, overlap_len);
                        }
                    }
                    rslt = kmers.find(CDBGraph::ReverseComplement(icontig.FrontKmer())); // left end to left end    
                    if(rslt != kmers.end()) {
                        for(auto& hit : rslt->second) {
                            auto jt = hit.first;
                            if(jt == it)
                                continue;
                            int overlap_len = hit.second+kmer_len;
                            auto& jlchunk = jt->m_seq.front().front();
                            if(overlap_len > (int)ilchunk.size() || overlap_len > (int)jlchunk.size())
                                continue;
                            TVariation seq(ilchunk.begin()+kmer_len, ilchunk.begin()+overlap_len);
                            ReverseComplementSeq(seq.begin(), seq.end());
                            if(!equal(jlchunk.begin(), jlchunk.begin()+hit.second, seq.begin()))
                                continue;

                            contig_overlaps.emplace_back(it, jt, -1, -1, overlap_len);
                        }
                    }
                    if(contig_overlaps.size() == 1)
                        overlaps.splice(overlaps.end(), contig_overlaps);
                }
            }
                
            for(auto it = overlaps.begin(); it != overlaps.end(); ) {
                auto overlap = *it;
                swap(get<0>(overlap), get<1>(overlap));
                swap(get<2>(overlap), get<3>(overlap));
                auto jt = find(it, overlaps.end(), overlap);
                if(jt == overlaps.end()) {
                    auto tmp = it++;
                    overlaps.erase(tmp);
                } else {
                    overlaps.erase(jt);
                    ++it;
                }
            }                


            for(auto it_loop = overlaps.begin(); it_loop != overlaps.end(); ) {
                auto it = it_loop++;
                auto& overlap = *it;
                int overlap_len = get<4>(overlap);
                auto icontigp = get<0>(overlap);
                auto jcontigp = get<1>(overlap);
                int diri = get<2>(overlap);
                int dirj = get<3>(overlap);
                
                auto NextIBase = [&]() {
                    CDBGraph::Node node = diri > 0 ? icontigp->BackKmer(): CDBGraph::ReverseComplement(icontigp->FrontKmer());
                    auto forward = m_graph.GetNodeSuccessors(node);
                    FilterNeighbors(forward, true);
                    if(forward.size() == 1) {
                        auto backward = m_graph.GetNodeSuccessors(CDBGraph::ReverseComplement(forward.front().m_node));
                        FilterNeighbors(backward, true);
                        if(backward.size() == 1 && CDBGraph::ReverseComplement(backward.front().m_node) == node)
                            return forward.front().m_nt;                            
                    }
                    return 'N';
                };
                auto NextJBase = [&]() {
                    return dirj < 0 ? *(jcontigp->m_seq.front().front().begin()+overlap_len) : Complement(*(jcontigp->m_seq.back().front().end()-overlap_len-1));
                };

                
                bool connected;
                if(diri > 0)
                    connected = (icontigp->m_seq.m_right_repeat < kmer_len);
                else
                    connected = (icontigp->m_seq.m_left_repeat < kmer_len);
                if(dirj > 0)
                    connected = connected && (jcontigp->m_seq.m_right_repeat < kmer_len);
                else
                    connected = connected && (jcontigp->m_seq.m_left_repeat < kmer_len);
                connected = connected && (NextIBase() == NextJBase());
                if(connected) {
                    swap(icontigp, jcontigp);
                    swap(diri, dirj);
                    connected = (NextIBase() == NextJBase());
                }
                if(!connected)
                    overlaps.erase(it);
            }
            
            cerr << "Overlap connections: " << overlaps.size() << " " << kmer_len << endl;
                
            while(!overlaps.empty()) {
                auto& overlap = overlaps.front();
                auto icontigp = get<0>(overlap);
                auto jcontigp = get<1>(overlap);
                int diri = get<2>(overlap);
                int dirj = get<3>(overlap);
                int overlap_len = get<4>(overlap);
                
                if(diri > 0) {
                    if(dirj > 0)
                        jcontigp->ReverseComplement();
                    jcontigp->ClipLeft(overlap_len-kmer_len+1);  // AddToRight assumes kmer-1 overlap   
                    icontigp->AddToRight(*jcontigp);
                } else {
                    if(dirj < 0)
                        jcontigp->ReverseComplement();
                    jcontigp->ClipRight(overlap_len-kmer_len+1);  // AddToLeft assumes kmer-1 overlap   
                    icontigp->AddToLeft(*jcontigp);
                }
                overlaps.pop_front();

                for(auto& overlap : overlaps) {
                    if(get<0>(overlap) == jcontigp) {
                        get<0>(overlap) = icontigp;
                        get<2>(overlap) = diri;
                    } else if(get<1>(overlap) == jcontigp) {
                        get<1>(overlap) = icontigp;
                        get<3>(overlap) = diri;
                    }
                }

                scontigs.erase(jcontigp);
            }                
        }


        // Starting from available graph nodes, generates all contigs >= min_len_for_new_seeds. Uses ncores threads.
        TContigList GenerateNewSeeds(int min_len_for_new_seeds, int ncores, CDBGraphDigger* test_graphdiggerp) {
            //assemble new seeds
            vector<TContigList> new_seeds_for_threads(ncores);
            list<function<void()>> jobs;
            for(auto& ns : new_seeds_for_threads) {
                jobs.push_back(bind(&CDBGraphDigger::NewSeedsJob, this, ref(ns), min_len_for_new_seeds));
            }
            RunThreads(ncores, jobs);

            //connect fragments 
            Graph().ClearHoldings();
            TContigList new_seeds = SContig::ConnectFragments(new_seeds_for_threads, Graph());
            
            int kmer_len = Graph().KmerLen();
            CReadHolder removed_seq(false);
            for(auto iloop = new_seeds.begin(); iloop != new_seeds.end(); ) {
                auto ic = iloop++;
                if((int)ic->LenMin() < min_len_for_new_seeds+2*kmer_len) {
                    removed_seq.PushBack(ic->m_seq[0].front());
                    new_seeds.erase(ic);
                    continue;
                } 

                if(test_graphdiggerp != nullptr) {
                    CReadHolder rh(false);
                    rh.PushBack(ic->m_seq[0].front());
                    double abundance = 0;
                    int knum = 0;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(test_graphdiggerp->Graph().KmerLen()) ; ik != rh.kend(); ++ik, ++knum) {
                        auto node = test_graphdiggerp->Graph().GetNode(*ik);
                        abundance += test_graphdiggerp->Graph().Abundance(node);
                    }
                    if(abundance < knum*test_graphdiggerp->m_hist_min) {
                        removed_seq.PushBack(ic->m_seq[0].front());
                        new_seeds.erase(ic);
                        continue;
                    }
                }

                string left(ic->m_seq[0].front().begin(), ic->m_seq[0].front().begin()+2*kmer_len-1);
                removed_seq.PushBack(left);
                string right(ic->m_seq[0].front().end()-2*kmer_len+1, ic->m_seq[0].front().end());
                removed_seq.PushBack(right);
                ic->ClipLeft(kmer_len);
                ic->ClipRight(kmer_len); 
                ic->m_seq.m_left_repeat = kmer_len-1;
                ic->m_seq.m_right_repeat = kmer_len-1;
            }
            for(CReadHolder::kmer_iterator ik = removed_seq.kbegin(kmer_len) ; ik != removed_seq.kend(); ++ik)            
                Graph().ClearVisited(Graph().GetNode(*ik));            
        
            return new_seeds;
        }
        // Using a longer kmer generates connectors and extenders and improves previously assembled contigs
        // scontigs - contigs (input/output)
        // ncores - number of threads
        void ConnectAndExtendContigs(TContigList& scontigs, int ncores) {

            vector<TContigList> extensions_for_jobs(ncores);
            list<function<void()>> jobs;
            for(auto& ex : extensions_for_jobs) {
                jobs.push_back(bind(&CDBGraphDigger::ExtendContigsJob, this, ref(scontigs), ref(ex)));
            }
            RunThreads(ncores, jobs);

            TContigList extensions = SContig::ConnectFragments(extensions_for_jobs, Graph()); 
            SContig::ConnectAndExtendContigs(scontigs, extensions);  
        }

        /*
        void ExtendContigsInRepeats(TContigList& scontigs, int ncores) {
            list<function<void()>> jobs;
            for(int thr = ncores; thr > 0; --thr) 
                jobs.push_back(bind(&CDBGraphDigger::ExtendContigsInRepeatsJob, this, ref(scontigs)));
            RunThreads(ncores, jobs);
        }
        */
        
        list<array<CReadHolder,2>> ConnectPairs(const list<array<CReadHolder,2>>& mate_pairs, int insert_size, int ncores) {
            CStopWatch timer;
            timer.Restart();

            list<array<CReadHolder,2>> paired_reads;
            list<function<void()>> jobs;
            for(auto& reads : mate_pairs) {
                auto& job_input = reads[0];
                paired_reads.push_back(array<CReadHolder,2>({CReadHolder(false), CReadHolder(true)}));
                if(job_input.ReadNum() > 0)  // not empty       
                    jobs.push_back(bind(&CDBGraphDigger::ConnectPairsJob, this, insert_size, ref(job_input), ref(paired_reads.back())));            
            }
            RunThreads(ncores, jobs);
                       
            size_t connected = 0;
            size_t not_connected = 0;
            for(auto& rh : paired_reads) {
                connected += rh[0].ReadNum();
                not_connected += rh[1].ReadNum();
            }
            size_t mates = 0;
            for(auto& rh : mate_pairs)
                mates += rh[0].ReadNum();
            cerr << "Connected: " << connected << " ambiguously connected: " << not_connected/2 << " from " << mates/2 << " mate pairs" << endl;        
            cerr << "Connect pairs in " << timer.Elapsed(); 

            return paired_reads;
        }
        
    private:
        // Prepares one of the mates of a read pair for connection
        // Finds the longest stretch of the read which could be assembled from both ends and clips the rest
        // read - read input/output
        // nodes - kmers for the remaining part
        void CheckAndClipRead(string& read, deque<CDBGraph::Node>& nodes) {
            int kmer_len = m_graph.KmerLen();

            string lextend = MostLikelyExtension(CDBGraph::ReverseComplement(m_graph.GetNode(read.substr(0, kmer_len))), kmer_len);        
            ReverseComplementSeq(lextend.begin(), lextend.end());
            string rextend = MostLikelyExtension(m_graph.GetNode(read.substr(read.size()-kmer_len)), kmer_len);

            deque<CDBGraph::Node> extended_nodes;
            CReadHolder rh(false);
            rh.PushBack(lextend+read+rextend);
            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik)  // iteration from last kmer to first  
                extended_nodes.push_front(m_graph.GetNode(*ik));

            vector<int> bases(read.size(), 0);
            unsigned read_pos = kmer_len-lextend.size();
            for(int kk = 0; lextend.size()+read_pos+1 < extended_nodes.size() && read_pos < read.size(); ++kk, ++read_pos) {
                CDBGraph::Node left = extended_nodes[kk];
                CDBGraph::Node node = extended_nodes[kk+1];
                if(!left || !GoodNode(left) || !node || !GoodNode(node))
                    continue;
                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(left);
                FilterNeighbors(successors, false);
                if(find_if(successors.begin(),successors.end(),[node](const CDBGraph::Successor& s){return s.m_node == node;}) == successors.end())                
                    continue;
            
                CDBGraph::Node right = m_graph.ReverseComplement(extended_nodes[lextend.size()+read_pos+1]);
                node = m_graph.ReverseComplement(extended_nodes[read_pos+lextend.size()]);
                if(!right || !GoodNode(right) || !node || !GoodNode(node))
                    continue;
                successors = m_graph.GetNodeSuccessors(right);
                FilterNeighbors(successors, false);
                if(find_if(successors.begin(),successors.end(),[node](const CDBGraph::Successor& s){return s.m_node == node;}) == successors.end())
                    continue;                

                bases[read_pos] = 1;
            }        

            int left = 0;             // first kmer position    
            int len = 0;              // number of consecutive good kmers (the sequence is longer by kmer_len-1)    
            for(unsigned k = 0; k < read.size(); ++k) {
                for( ; k < read.size() && !bases[k]; ++k);         // skip bad bases    
                int current_left = k;
                int current_len = 0;
                for( ; k < read.size() && bases[k]; ++k, ++current_len);   // count adjacent good bases 
                if(current_len > len) {
                    left = current_left;
                    len = current_len;
                }
            }

            if(len < kmer_len) {
                read.clear();
                nodes.clear();
            } else {
                read = read.substr(left, len);
                nodes.resize(len-kmer_len+1);
                copy(extended_nodes.begin()+lextend.size()+left, extended_nodes.begin()+lextend.size()+left+len-kmer_len+1, nodes.begin());
            }
        }    

        // one-thread worker for paired reads connection
        // saves reads which were unambiguously connected; extends the ends of ambiguously connected reads and
        // keeps them for future; discards reads which don't have connection
        // insert_size - the maximal limit of the insert length
        // mate_pairs - pairs for connection (one mate after another)
        // paired_reads - [0] connected reads, [1] reads for future connection    
        void ConnectPairsJob(int insert_size, const CReadHolder& mate_pairs, array<CReadHolder,2>& paired_reads) {
            if(mate_pairs.ReadNum() < 2)
                return;

            int kmer_len = m_graph.KmerLen();

            for(CReadHolder::string_iterator is = mate_pairs.sbegin(); is != mate_pairs.send(); ++is) {
                string read1 = *is;
                string read2 = *(++is);
                if((int)min(read1.size(),read2.size()) < kmer_len)
                    continue;

                deque<CDBGraph::Node> nodes1;
                CheckAndClipRead(read1, nodes1);
                if(read1.empty())
                    continue;
                CDBGraph::Node last_node1 = nodes1.back();                        
                
                ReverseComplementSeq(read2.begin(), read2.end());
                deque<CDBGraph::Node> nodes2;
                CheckAndClipRead(read2, nodes2);
                if(read2.empty())
                    continue;
                CDBGraph::Node first_node2 = nodes2.front();

                int steps = insert_size;

                pair<TBases, EConnectionStatus> rslt = ConnectTwoNodes(last_node1, first_node2, steps);

                bool ambiguous = false;

                if(rslt.second == CDBGraphDigger::eAmbiguousConnection) {
                    ambiguous = true; 
                } else {
                    string read;
                    unordered_set<CDBGraph::Node> read_nodes;
                    if(rslt.second == eSuccess) {
                        string r1 = read1;
                        for(auto& suc : rslt.first) {
                            r1.push_back(suc.m_nt);
                            read_nodes.insert(suc.m_node);
                        }
                        r1 += read2.substr(kmer_len);
                        read_nodes.insert(nodes2.begin()+1, nodes2.end());
                        rslt = ConnectTwoNodes(CDBGraph::ReverseComplement(first_node2), CDBGraph::ReverseComplement(last_node1), steps);
                        if(rslt.second == eSuccess) {
                            string seq;
                            for(auto& suc : rslt.first)
                                seq.push_back(suc.m_nt);
                            ReverseComplementSeq(seq.begin(), seq.end());
                            string r2 = read1.substr(0, read1.size()-kmer_len)+seq+read2;
                            if(r1 == r2)
                                read = r1;
                        }
                    
                        if(read.empty())
                            ambiguous = true;
                    } else {
                        //check for long overlap with extension     
                        int hit = find(nodes2.begin(), nodes2.end(), last_node1) - nodes2.begin(); // first kmer position of the hit        
                        if(hit < (int)min(nodes1.size(),nodes2.size()) && equal(nodes2.begin(), nodes2.begin()+hit, nodes1.end()-hit-1)) {
                            read = read1+read2.substr(hit+kmer_len);
                            read_nodes.insert(nodes2.begin()+hit+1, nodes2.end());
                        }
                    }

                    if(!read.empty()) { 
                        read_nodes.insert(nodes1.begin(), nodes1.end());
                        if(read_nodes.size() == read.size()-kmer_len+1) {
                            paired_reads[0].PushBack(read);
                            continue;
                        } else {
                            ambiguous = true;
                        }
                    }
                }

                if(ambiguous) {
                    string lextend =  StringentExtension(CDBGraph::ReverseComplement(nodes1.front()), kmer_len).first;
                    ReverseComplementSeq(lextend.begin(), lextend.end());
                    paired_reads[1].PushBack(lextend+read1);
                    read2 += StringentExtension(nodes2.back(), kmer_len).first;
                    ReverseComplementSeq(read2.begin(), read2.end());
                    paired_reads[1].PushBack(read2);
                }
            }                                
        }

        // one-thread worker for generating new seeds
        // returns contigs sequences which are either >= min_len or are known fragments
        // contigs - generated contigs
        // min_len - minimal length for acceptable contigs
        void NewSeedsJob(TContigList& contigs, int min_len) {
            for(size_t index = 0; index < Graph().GraphSize(); ++index) {
                CDBGraph::Node initial_node = 2*(index+1);
                SContig contig = GetContigForKmer(initial_node, min_len);
                if(!contig.m_seq.empty())
                    contigs.push_back(contig);
            }
        }

        // one-thread worker for generating connectors and extenders for previously assembled contigs
        // scontigs - contigs (input/output)
        // extensions - generated sequences
        //               (K-1)mers from a window are used as contig sequence near end of contig may not be correct
        
        void ExtendContigsJob(TContigList& scontigs, TContigList& extensions) {
            for(auto& contig : scontigs) {
                if(!contig.m_is_taken.Set(1))  // grab contig
                    continue;

                int kmer_len = contig.m_kmer_len;
                int chunks = contig.m_seq.size();

                
                if(contig.m_seq.m_right_repeat < kmer_len) {
                    if((int)contig.m_seq[chunks-1].front().size() >= kmer_len) { // normal end 
                        CDBGraph::Node takeoff_node = contig.BackKmer();

                        if(takeoff_node && GoodNode(takeoff_node) && !Graph().IsMultContig(takeoff_node)) {         // valid uniq kmer  
                            int allowed_intrusion = max(0, (int)contig.m_seq.ChunkLenMax(chunks-1)-kmer_len);
                            tuple<SContig, CDBGraph::Node, int> extension = ExtendToRight(takeoff_node, allowed_intrusion);
                            if(!get<0>(extension).m_seq.empty() || get<1>(extension)) { // extension could be empty - starting kmer + landing kmer  
                                contig.ClipRight(get<2>(extension));
                                SContig sc(&contig, 1, contig.BackKmer(), get<0>(extension), get<1>(extension), Graph());
                                extensions.push_back(sc); 
                            }                   
                        }                
                    } else if(chunks > 2 && (int)contig.m_seq[chunks-3].front().size() >= kmer_len) {  // snp close to the end  
                        CDBGraph::Node takeoff_node = contig.RightConnectingNode();
                        if(takeoff_node && GoodNode(takeoff_node) && !Graph().IsMultContig(takeoff_node)) {         // valid uniq kmer  
                            int shift = 0;
                            for( ; contig.m_seq.AllSameL(chunks-2, shift); ++shift);
                            if(shift < kmer_len) {
                                if(shift > 0) {
                                    auto& takeoff_seq = contig.m_seq[chunks-3].front();
                                    auto& snp_seq = contig.m_seq[chunks-2].front();
                                    string kmer(takeoff_seq.end()-kmer_len+shift, takeoff_seq.end());
                                    kmer.insert(kmer.end(), snp_seq.begin(), snp_seq.begin()+shift);
                                    takeoff_node = Graph().GetNode(kmer);
                                }
                                if(takeoff_node && !Graph().IsMultContig(takeoff_node)) {      
                                    tuple<SContig, CDBGraph::Node, int> extension = ExtendToRight(takeoff_node, shift);
                                    int ext_chunks = get<0>(extension).m_seq.size();
                                    if(ext_chunks > 1 && get<2>(extension) == shift) {
                                        SContig sc(&contig, 1, contig.RightConnectingNode(), get<0>(extension), get<1>(extension), Graph());
                                        sc.m_seq.StabilizeVariantsOrder();
                                        if(sc.m_seq[1] == contig.m_seq[chunks-2])
                                            extensions.push_back(sc); 

                                        /*
                                        lock_guard<mutex> guard(out_mutex);
                                        cerr << "ShiftR: " << shift << endl;
                                        
                                        cerr << "Right ext: " << sc.m_seq.size() << " ";
                                        if(sc.m_next_right)
                                            cerr << m_graph.GetNodeSeq(sc.m_next_right);
                                        cerr << endl;
                                        for(auto& chunk: sc.m_seq) {
                                            cerr << "Chunk" << endl;
                                            for(auto& seq : chunk) {
                                                for(char c : seq)
                                                    cerr << c;
                                                cerr << endl;
                                            }
                                        }                                

                                        cerr << "Contig: " << contig.m_seq.size() << " " << contig.m_seq.m_circular << " ";
                                        cerr << (contig.m_next_left ? contig.m_graph.GetNodeSeq(contig.m_next_left) : "-") << " ";
                                        cerr << (contig.m_next_right ? contig.m_graph.GetNodeSeq(contig.m_next_right) : "-");
                                        cerr << endl;
                                        for(auto& chunk : contig.m_seq) {
                                            cerr << "Chunk:     " << distance(chunk.begin(), chunk.end()) << endl;
                                            for(auto& seq : chunk) {
                                                for(char c : seq)
                                                    cerr << c;
                                                cerr << endl;
                                            }
                                        }
                                        */

                                    }
                                }
                            }
                        }
                    }
                }
                                 

                
                if(contig.m_seq.m_left_repeat < kmer_len) {
                    if((int)contig.m_seq[0].front().size() >= kmer_len) { // normal end 
                        CDBGraph::Node takeoff_node = CDBGraph::ReverseComplement(contig.FrontKmer());
                        if(takeoff_node && GoodNode(takeoff_node) && !Graph().IsMultContig(takeoff_node)) {         // valid uniq kmer     
                            int allowed_intrusion = max(0, (int)contig.m_seq.ChunkLenMax(0)-kmer_len);
                            tuple<SContig, CDBGraph::Node, int> extension = ExtendToRight(takeoff_node, allowed_intrusion);
                            if(!get<0>(extension).m_seq.empty() || get<1>(extension)) { // extension could be empty - starting kmer + landing kmer   
                                contig.ClipLeft(get<2>(extension));
                                SContig sc(&contig, -1, CDBGraph::ReverseComplement(contig.FrontKmer()), get<0>(extension), get<1>(extension), Graph());
                                sc.ReverseComplement();
                                extensions.push_back(sc);                        
                            }
                        }
                    } else if(chunks > 2 && (int)contig.m_seq[2].front().size() >= kmer_len) {  // snp close to the end 
                        CDBGraph::Node takeoff_node = CDBGraph::ReverseComplement(contig.LeftConnectingNode());
                        
                        if(takeoff_node && GoodNode(takeoff_node) && !Graph().IsMultContig(takeoff_node)) {         // valid uniq kmer  
                            int shift = 0;
                            for( ; contig.m_seq.AllSameR(1, shift); ++shift);

                            if(shift < kmer_len) {
                                if(shift > 0) {
                                    auto& takeoff_seq = contig.m_seq[2].front();
                                    auto& snp_seq = contig.m_seq[1].front();
                                    string kmer(snp_seq.end()-shift, snp_seq.end());
                                    kmer.insert(kmer.end(), takeoff_seq.begin(), takeoff_seq.begin()+kmer_len-shift);
                                    takeoff_node = CDBGraph::ReverseComplement(Graph().GetNode(kmer));
                                }
                                if(takeoff_node && !Graph().IsMultContig(takeoff_node)) {                                          
                                    tuple<SContig, CDBGraph::Node, int> extension = ExtendToRight(takeoff_node, shift);
                                    int ext_chunks = get<0>(extension).m_seq.size();

                                    if(ext_chunks > 1 && get<2>(extension) == shift) {
                                        SContig sc(&contig, -1, CDBGraph::ReverseComplement(contig.LeftConnectingNode()), get<0>(extension), get<1>(extension), Graph());
                                        sc.ReverseComplement();
                                        sc.m_seq.StabilizeVariantsOrder();

                                        if(sc.m_seq[ext_chunks-1] == contig.m_seq[1]) {
                                            extensions.push_back(sc);  

                                            /*
                                            lock_guard<mutex> guard(out_mutex);
                                            cerr << "ShiftL: " << shift << endl;
                                            cerr << "Left ext: " << sc.m_seq.size() << " ";
                                            if(sc.m_next_left)
                                                cerr << m_graph.GetNodeSeq(sc.m_next_left);
                                            cerr << endl;
                                            for(auto& chunk: sc.m_seq) {
                                                cerr << "Chunk" << endl;
                                                for(auto& seq : chunk) {
                                                    for(char c : seq)
                                                        cerr << c;
                                                    cerr << endl;
                                                }
                                            }  
                                            
                                            cerr << "Contig: " << contig.m_seq.size() << " " << contig.m_seq.m_circular << " ";
                                            cerr << (contig.m_next_left ? contig.m_graph.GetNodeSeq(contig.m_next_left) : "-") << " ";
                                            cerr << (contig.m_next_right ? contig.m_graph.GetNodeSeq(contig.m_next_right) : "-");
                                            cerr << endl;
                                            for(auto& chunk : contig.m_seq) {
                                                cerr << "Chunk:         " << distance(chunk.begin(), chunk.end()) << endl;
                                                for(auto& seq : chunk) {
                                                    for(char c : seq)
                                                        cerr << c;
                                                    cerr << endl;
                                                }
                                            }
                                            */
                                      

                                        }                                                              
                                    }
                                }
                            }
                        }
                    }
                }
                  
            }
        }        


        CDBGraph& m_graph;
        double m_fraction;
        int m_jump;
        int m_hist_min;
        int m_low_count;
        size_t m_max_branch;
        bool m_allow_snps;
    };

}; // namespace
#endif /* _GraphDigger_ */
