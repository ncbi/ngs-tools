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
   scan_window flank kmers. When these contigs are started, m_left_link/m_right_link are assigned to point to the
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
    
    struct SContig {
        SContig() : m_next_left(0), m_next_right(0), m_left_extend(0), m_right_extend(0) {}
        SContig(const string& contig, CDBGraph& graph) : 
            m_seq(contig.begin(), contig.end()), m_next_left(0), m_next_right(0), m_left_link(nullptr), m_left_shift(0), m_right_link(nullptr), m_right_shift(0), 
            m_left_extend(0), m_right_extend(0), m_kmer_len(graph.KmerLen()), m_is_taken(0) {
            CReadHolder rh(false);
            rh.PushBack(contig);
            for(CReadHolder::kmer_iterator itk = rh.kbegin(m_kmer_len); itk != rh.kend(); ++itk) {  // gives kmers in reverse order!
                CDBGraph::Node node = graph.GetNode(*itk);
                m_kmers.push_front(node);  // may be 0
                if(node)
                    graph.SetVisited(node);
            }          
        }
        SContig(const TBases& to_left, const TBases& to_right, CDBGraph::Node initial_node, CDBGraph::Node lnode, CDBGraph::Node rnode, const CDBGraph& graph) :  
        m_next_left(lnode), m_next_right(rnode), m_left_link(nullptr), m_left_shift(0), m_right_link(nullptr), m_right_shift(0), m_kmer_len(graph.KmerLen()), m_is_taken(0) {                                                                                                                                                   
//          initial_node - the starting kmer
//          to_left - left extension of the starting kmer
//          to_right - right extension of the starting kmer
//          lnode - left denied node
//          rnode - right denied node
//          graph - de Bruijn graph
                                                                       
            // take parts of the assembled sequence and put them together in SContig
            for(const auto& base : to_left) {
                m_seq.push_front(Complement(base.m_nt));
                m_kmers.push_front(CDBGraph::ReverseComplement(base.m_node));
            }
            m_kmers.push_back(initial_node);
            string ikmer = graph.GetNodeSeq(initial_node);
            m_seq.insert(m_seq.end(), ikmer.begin(), ikmer.end());
            for(const auto& base : to_right) {
                m_seq.push_back(base.m_nt);
                m_kmers.push_back(base.m_node);
            }

            m_left_extend = m_right_extend = m_seq.size();
        }
        SContig(SContig* link, int shift, CDBGraph::Node takeoff_node, const TBases& extension, CDBGraph::Node rnode, const CDBGraph& graph) :
            m_next_left(takeoff_node), m_next_right(rnode), m_left_link(link), m_left_shift(shift), m_right_link(nullptr), m_right_shift(0), 
            m_kmer_len(graph.KmerLen()), m_is_taken(0) {
            string kmer = graph.GetNodeSeq(takeoff_node);
            m_seq.insert(m_seq.end(), kmer.begin()+1, kmer.end());
            for(const auto& base : extension) {
                m_seq.push_back(base.m_nt);
                m_kmers.push_back(base.m_node);
            }

            m_left_extend = m_right_extend = m_seq.size();
        }

        void ReverseComplement() {
            ReverseComplementSeq(m_seq.begin(), m_seq.end());
            reverse(m_kmers.begin(), m_kmers.end());
            for(auto& kmer : m_kmers)
                kmer = CDBGraph::ReverseComplement(kmer);
            swap(m_next_left, m_next_right);
            m_next_left = CDBGraph::ReverseComplement(m_next_left);
            m_next_right = CDBGraph::ReverseComplement(m_next_right);
            swap(m_left_link, m_right_link);
            swap(m_left_shift, m_right_shift);
            swap(m_left_extend, m_right_extend);
        }
        void AddToRight(const SContig& other) { 
            if(other.m_right_extend < (int)other.m_seq.size()) {
                m_right_extend = other.m_right_extend;
            } else {
                m_right_extend += other.m_right_extend-m_kmer_len+1;
                if(m_left_extend == (int)m_seq.size())
                    m_left_extend = m_right_extend;
            }
            m_seq.insert(m_seq.end(), other.m_seq.begin()+m_kmer_len-1, other.m_seq.end());
            m_kmers.insert(m_kmers.end(), other.m_kmers.begin(), other.m_kmers.end());
            m_next_right = other.m_next_right;
            m_right_link = other.m_right_link;
            m_right_shift = other.m_right_shift;            
        }
        void AddToLeft(const SContig& other) {
            if(other.m_left_extend < (int)other.m_seq.size()) {
                m_left_extend = other.m_left_extend;
            } else {
                m_left_extend += other.m_left_extend-m_kmer_len+1; 
                if(m_right_extend == (int)m_seq.size())
                    m_right_extend = m_left_extend;
            }               
            m_seq.insert(m_seq.begin(), other.m_seq.begin(), other.m_seq.end()-m_kmer_len+1);
            m_kmers.insert(m_kmers.begin(), other.m_kmers.begin(), other.m_kmers.end());
            m_next_left = other.m_next_left;
            m_left_link = other.m_left_link;
            m_left_shift = other.m_left_shift;
        }
        void ClipRight(int clip) {
            if(clip > 0) {
                m_right_extend = max(0, m_right_extend-clip);
                m_seq.erase(m_seq.end()-clip, m_seq.end());
                m_kmers.erase(m_kmers.end()-clip, m_kmers.end());
                m_next_right = 0;
                m_right_link = nullptr;
                m_right_shift = 0;            
            }
        }
        void ClipLeft(int clip) {
            if(clip > 0) {
                m_left_extend = max(0, m_left_extend-clip);
                m_seq.erase(m_seq.begin(), m_seq.begin()+clip);
                m_kmers.erase(m_kmers.begin(), m_kmers.begin()+clip);
                m_next_left = 0;
                m_left_link = nullptr;
                m_left_shift = 0;   
            }         
        }

        size_t Len() const { return m_seq.size(); }

        // find position of the minimal non-zero kmer
        size_t MinKmerPosition() const {
            size_t mkp = 0;
            for(size_t i = 0; i < m_kmers.size(); ++i) {
                if(m_kmers[i] != 0 && (m_kmers[mkp] == 0 || m_kmers[i] < m_kmers[mkp]))
                    mkp = i;
            }
            return mkp;
        }

        // stabilize contig orientation using minimal kmer in the contig
        void SelectMinDirection() {
            CDBGraph::Node minkmer = m_kmers[MinKmerPosition()];
            if(minkmer && minkmer%2)
                ReverseComplement();            
        }

        // finds stable origin for circular contigs by placing minimal kmer at the beginning of the sequence
        void RotateCircularToMinKmer() { // assumes that the next extension of sequence would give the first kmer (m_next_right == m_kmers.front())
            m_seq.erase(m_seq.end()-m_kmer_len+1, m_seq.end());
            size_t first_base = MinKmerPosition();
            // cut after the minimal kmer so rotate brings it to the beginning
            if(m_kmers[first_base]%2)
                first_base = (first_base+m_kmer_len)%m_kmers.size();
            rotate(m_seq.begin(), m_seq.begin()+first_base, m_seq.end());
            rotate(m_kmers.begin(), m_kmers.begin()+first_base, m_kmers.end());
            m_next_left = m_kmers.back();
            m_next_right = *(m_kmers.end()-m_kmer_len+1);

            //remove extra kmers
            m_kmers.erase(m_kmers.end()-m_kmer_len+1, m_kmers.end());
            
            //clean edges   
            m_left_link = nullptr;
            m_left_shift = 0;
            m_right_link = nullptr;
            m_right_shift = 0;
            m_left_extend = 0;   // prevents any further clipping    
            m_right_extend = 0;  // prevents any further clipping    
        }

        bool operator<(const SContig& other) const { return m_seq < other.m_seq; }

        // connects fragments created in different threads and combines doubled 'empty' linkers
        static TContigList ConnectFragments(vector<TContigList>& fragments, const CDBGraph& graph) {
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
                            if(contig.m_kmers.empty() && contig.m_left_link && !other->m_left_link && contig.m_next_right == (other->m_kmers.empty() ? other->m_next_right:other->m_kmers.front())) {
                                other->AddToLeft(contig); // add left link to other
                                connected.pop_front();
                                continue;
                            } else if(other->m_kmers.empty() && other->m_left_link && !contig.m_left_link && other->m_next_right == (contig.m_kmers.empty() ? contig.m_next_right:contig.m_kmers.front())) {
                                contig.AddToLeft(*other); // add left link to contig    
                                rslt.first->second = connected.begin();
                                denied_right_nodes.erase(other->m_next_right);
                                connected.erase(other);
                            } else {
                                cerr << "Unexpected left fork: " << graph.GetNodeSeq(contig.m_next_left) << " " << contig.m_next_left << endl;
                            }
                        }
                    }
                    if(contig.m_next_right) {
                        auto rslt = denied_right_nodes.insert(make_pair(contig.m_next_right, connected.begin()));
                        if(!rslt.second) {
                            TContigList::iterator other = rslt.first->second;
                            if(contig.m_kmers.empty() && contig.m_right_link && !other->m_right_link && contig.m_next_left == (other->m_kmers.empty() ? other->m_next_left:other->m_kmers.back())) {
                                other->AddToRight(contig); // add right link to other
                                denied_left_nodes.erase(contig.m_next_left);
                                connected.pop_front();
                            } else if(other->m_kmers.empty() && other->m_right_link &&  !contig.m_right_link && other->m_next_left == (contig.m_kmers.empty() ? contig.m_next_left:contig.m_kmers.back())) {
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
                if(contig.m_kmers.empty())
                    continue; 

                if(contig.m_next_right)
                    denied_right_nodes.erase(contig.m_next_right);
                if(contig.m_next_left)
                    denied_left_nodes.erase(contig.m_next_left);
                bool keep_doing = true;
                while(keep_doing) {
                    keep_doing = false;
                    if(contig.m_next_right) {
                        CDBGraph::Node rnode = contig.m_kmers.back();
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
                        CDBGraph::Node lnode = contig.m_kmers.front();
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
            
                if(contig.m_next_right == contig.m_kmers.front() && (int)contig.m_kmers.size() >= graph.KmerLen())  // circular and not very short  
                    contig.RotateCircularToMinKmer();                     
            }

            return connected;                
        }

        // connects and extends contigs from previous iteration using a longer kmer
        // scontigs - previous contigs
        // extensions - connectors and extenders produced by longer kmer
        static void ConnectAndExtendContigs(TContigList& scontigs, TContigList& extensions) {
            if(scontigs.empty())
                return;

            int kmer_len = scontigs.front().m_kmer_len;
            typedef unordered_map<SContig*, map<int, SContig*>> TExtensionsDoubleMap; // connections to left contig sides   (pointer to contig; shift on contig's left side (sorted); pointer to connection)
            TExtensionsDoubleMap left_connections;
            TExtensionsDoubleMap right_connections;
            TExtensionsDoubleMap left_extensions;
            TExtensionsDoubleMap right_extensions;
            int connectors = 0;
            int extenders = 0;
            for(auto& ex : extensions) {
                if(ex.m_left_link && ex.m_right_link) {
                    ++connectors;
                    if(ex.m_left_shift < 0)
                        left_connections[ex.m_left_link][-(ex.m_left_shift+1)] = &ex;
                    else
                        right_connections[ex.m_left_link][ex.m_left_shift-1] = &ex;
                    if(ex.m_right_shift < 0)
                        left_connections[ex.m_right_link][-(ex.m_right_shift+1)] = &ex;
                    else
                        right_connections[ex.m_right_link][ex.m_right_shift-1] = &ex;
                } else if(ex.m_left_link) {
                    ++extenders;
                    if(ex.m_left_shift < 0)
                        left_extensions[ex.m_left_link][-(ex.m_left_shift+1)] = &ex;
                    else
                        right_extensions[ex.m_left_link][ex.m_left_shift-1] = &ex;
                } else if(ex.m_right_link) {
                    ++extenders;
                    if(ex.m_right_shift < 0)
                        left_extensions[ex.m_right_link][-(ex.m_right_shift+1)] = &ex;
                    else
                        right_extensions[ex.m_right_link][ex.m_right_shift-1] = &ex;
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
                        TExtensionsDoubleMap::iterator rslt;
                        // check if connection to other contigs is possible
                        if((plus && (rslt = right_connections.find(fragment)) != right_connections.end()) ||
                           (!plus && (rslt = left_connections.find(fragment)) != left_connections.end())) {

                            if(rslt->second.size() > 1) 
                                cerr << "Multiple connections" << endl;
                        
                            SContig* connector = rslt->second.begin()->second;
                            if(connector->m_right_link == fragment) { // either reversed or circular            
                                int rshift = connector->m_right_shift;
                                rshift = rshift > 0 ? rshift-1 : -(rshift+1);
                                if(rshift < (int)contig.m_kmers.size() && CDBGraph::ReverseComplement(connector->m_next_right) == *(contig.m_kmers.end()-rshift-1))
                                    connector->ReverseComplement();
                            }
                            int lshift = connector->m_left_shift;
                            contig.ClipRight(lshift > 0 ? lshift-1 : -(lshift+1));
                            if(connector->m_left_link != fragment || contig.m_kmers.back() != connector->m_next_left)
                                cerr << "Corrupted connectionA" << endl;
                            contig.AddToRight(*connector);

                            fragment = connector->m_right_link;
                            if(fragment->m_is_taken)      // don't connect already used contig (this is result of multiple connection)          
                                break;
                        
                            fragment->m_is_taken = 1;     // fragment will be removed
                            int rshift = connector->m_right_shift;
                            plus = rshift < 0;
                            rshift = plus ? -(rshift+1) : rshift-1;
                            if(!plus)
                                fragment->ReverseComplement();
                            fragment->ClipLeft(rshift);
                            if(fragment->m_kmers.front() != connector->m_next_right)
                                cerr << "Corrupted connectionB" << endl;
                            circular = (fragment == &contig);

                            if(!circular) {  // not circular            
                                contig.AddToRight(*fragment);
                                continue;
                            } else { //stabilize circular contig            
                                contig.RotateCircularToMinKmer();
                                break;                        
                            }
                        } 
                        if((plus && (rslt = right_extensions.find(fragment)) != right_extensions.end()) ||
                           (!plus && (rslt = left_extensions.find(fragment)) != left_extensions.end())) {
                            int extra_len = 0;
                            for(auto& ext : rslt->second) {
                                int shift = ext.first;
                                SContig* extender = ext.second;
                                if((int)extender->m_kmers.size()-shift > extra_len) {
                                    if(extender->m_right_link && extender->m_right_link == fragment)
                                        extender->ReverseComplement();
                                    contig.ClipRight(extra_len+shift);
                                    if(extender->m_left_link != fragment || contig.m_kmers.back() != extender->m_next_left)
                                        cerr << "Corrupted extension" << endl;
                                    contig.AddToRight(*extender);
                                    extra_len = extender->m_kmers.size()-shift;
                                }
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
                contig.ClipLeft(min(kmer_len,contig.m_left_extend));
                contig.ClipRight(min(kmer_len,contig.m_right_extend)); 
            }

            //remove fragments; stabilize orientation and order which are random in multithreading          
            for(auto iloop = scontigs.begin(); iloop != scontigs.end(); ) {
                auto ic = iloop++;
                if(ic->m_is_taken != 2)
                    scontigs.erase(ic);
                else
                    ic->SelectMinDirection();
            }
            scontigs.sort();
        }

        // m_seq.size() == m_kmer.size()+kmer_len-1
        // Extreme case: m_kmer.size() == 0; m_seq.size == kmer_len-1 (represents two connected 'next' kmers)
        // SContig is not very convenient for 'circular' sequences which should have m_seq.size() == m_kmer.size()
        deque<char> m_seq;               // sequence
        deque<CDBGraph::Node> m_kmers;   // kmers

        CDBGraph::Node m_next_left;      // denied left kmer (connection possible but it is already owned)
        CDBGraph::Node m_next_right;     // denied right kmer (connection possible but it is already owned)

        SContig* m_left_link;  // if set points to 'left' contig
        int m_left_shift;      // shift+1 for m_next_left in this contig (positive for the right end)
        SContig* m_right_link; // if set points to 'right' contig
        int m_right_shift;     // shift+1 for m_next_right in this contig (positive for the right end)

        int m_left_extend;     // number of newly assembled bases which could be clipped
        int m_right_extend;    // number of newly assembled bases which could be clipped

        int m_kmer_len;
        SAtomic<uint8_t> m_is_taken;
    };


    // This is a very lightweight class holding a reference to de Bruijn graph and main assembling parameters
    // It provides function used in assembling
    class CDBGraphDigger {
    public:

        CDBGraphDigger(CDBGraph& graph, double fraction, int jump, int low_count) : m_graph(graph), m_fraction(fraction), m_jump(jump), m_hist_min(graph.HistogramMinimum()), m_low_count(low_count) { 
            m_max_branch = 200; // maximum number of paths explored before quitting
            cerr << "Valley: " << m_hist_min << endl; 
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
                FilterNeighbors(successors);
                if(successors.empty())
                    return make_pair(s, true);
                if(successors.size() != 1)
                    return make_pair(s, false);            
                node = successors[0].m_node;
                s.push_back(successors[0].m_nt);
            }
            return make_pair(s, false);
        }

        bool IsSimpleNode(CDBGraph::Node initial_node, int len) const {
            stack<pair<CDBGraph::Node,int>> q;
            q.emplace(initial_node,len);
            while(!q.empty()) {
                CDBGraph::Node node = q.top().first;
                int remaining_len = q.top().second;
                q.pop();
                vector<CDBGraph::Successor> successors = GetReversibleNodeSuccessors(node);
                if(successors.empty())
                    return false;

                if(remaining_len > 1) {
                    for(auto& suc : successors)
                      q.emplace(suc.m_node,remaining_len-1);
                }
            }
            return true;
        }


        bool IsSimpleNodeOld(CDBGraph::Node initial_node, int len) const {
            stack<pair<CDBGraph::Node,int>> q;
            q.emplace(initial_node,len);
            while(!q.empty()) {
                CDBGraph::Node node = q.top().first;
                int remaining_len = q.top().second;
                q.pop();
                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(node);
                FilterNeighbors(successors);
                if(successors.empty())
                    return false;
                for(auto& suc : successors) {
                    vector<CDBGraph::Successor> step_back = m_graph.GetNodeSuccessors(m_graph.ReverseComplement(suc.m_node));
                    FilterNeighbors(step_back);
                    if(step_back.size() != 1 || m_graph.ReverseComplement(step_back.front().m_node) != node) 
                        return false;
                    else if(remaining_len > 1)
                        q.emplace(suc.m_node,remaining_len-1);
                }
            }
            return true;
        }

        vector<CDBGraph::Successor> GetReversibleNodeSuccessors(CDBGraph::Node node) const {
            vector<CDBGraph::Successor> neighbors = m_graph.GetNodeSuccessors(node);
            FilterNeighbors(neighbors);
            for(auto& neighbor : neighbors) {
                vector<CDBGraph::Successor> step_back = m_graph.GetNodeSuccessors(m_graph.ReverseComplement(neighbor.m_node));
                FilterNeighbors(step_back);
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
        void FilterNeighbors(vector<CDBGraph::Successor>& successors) const {
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
        
            // strand specific noise reduction for Illumina issue of GGT->GG[ACG] or strand balance issue
            if(m_graph.GraphIsStranded() && successors.size() > 1) {

                double fraction = 0.1*m_fraction;
            
                int target = -1;
                for(int j = 0; target < 0 && j < (int)successors.size(); ++j) {
                    if(m_graph.GetNodeSeq(successors[j].m_node).substr(m_graph.KmerLen()-3) == "GGT") 
                        target = j;
                }
                if(target >= 0 && GoodNode(successors[target].m_node)) {
                    double am = m_graph.Abundance(successors[target].m_node)*(1-m_graph.PlusFraction(successors[target].m_node));
                    for(int j = 0; j < (int)successors.size(); ) {
                        if(m_graph.Abundance(successors[j].m_node)*(1-m_graph.PlusFraction(successors[j].m_node)) < fraction*am)
                            successors.erase(successors.begin()+j);
                        else
                            ++j;
                    }
                    return;
                }

                for(int j = 0; target < 0 && j < (int)successors.size(); ++j) {
                    if(MostLikelySeq(successors[j], 3) == "ACC") 
                        target = j;
                }
                if(target >= 0 && GoodNode(successors[target].m_node)) {
                    double ap = m_graph.Abundance(successors[target].m_node)*m_graph.PlusFraction(successors[target].m_node);
                    for(int j = 0; j < (int)successors.size(); ) {
                        if(m_graph.Abundance(successors[j].m_node)*m_graph.PlusFraction(successors[j].m_node) < fraction*ap)
                            successors.erase(successors.begin()+j);
                        else
                            ++j;
                    }
                    return;
                }
                       
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
            FilterNeighbors(successors);
            for(auto& suc : successors) {
                storage.push_back(SElement(suc, 0));
                current_elements[suc.m_node] = &storage.back();
            }

            list<SElement> connections;
            for(int step = 1; step < steps && !current_elements.empty(); ++step) {
                TElementMap new_elements;
                for(auto& el : current_elements) {
                    vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(el.first);
                    FilterNeighbors(successors);
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

        typedef tuple<TBases,size_t> TSequence;
        typedef list<TSequence> TSeqList;
        typedef unordered_map<CDBGraph::Node, tuple<TSeqList::iterator, bool>> TBranch;  // all 'leaves' will have the same length  

        // makes one step extension for JumpOver()
        // sequences - keep the actual assembled sequences
        // branch - a lightweight map of the last kmer to the sequence
        void OneStepBranchExtend(TBranch& branch, TSeqList& sequences) {
            TBranch new_branch;
            for(auto& leaf : branch) {
                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(leaf.first);
                FilterNeighbors(successors);
                if(successors.empty()) {
                    sequences.erase(get<0>(leaf.second));
                    continue;
                }
                for(int i = successors.size()-1; i >= 0; --i) {
                    TSeqList::iterator is = get<0>(leaf.second);
                    if(i > 0) {  // copy sequence if it is a fork               
                        sequences.push_front(*is);
                        is = sequences.begin();
                    }
                    TBases& bases = get<0>(*is);
                    size_t& abundance = get<1>(*is);
                    bases.push_back(successors[i]);
                    CDBGraph::Node& node = successors[i].m_node;
                    abundance += m_graph.Abundance(node);

                    pair<TBranch::iterator, bool> rslt = new_branch.insert(make_pair(node, make_tuple(is, get<1>(leaf.second))));
                    if(!rslt.second) {  // we already have this kmer - select larger abundance and mark as ambiguous 
                        get<1>(rslt.first->second) = true;                      // mark as ambiguous    
                        TSeqList::iterator& js = get<0>(rslt.first->second);    // existing one     
                        if(abundance > get<1>(*js)) {                           // new is better        
                            sequences.erase(js);
                            js = is;
                        } else {                                                // existing is better or equal      
                            sequences.erase(is);                                        
                        }
                    }
                }
            }

            swap(branch, new_branch);
        }

        // for a given fork (successors.size() > 1) perform at least min_extent but no more than max_extent extension steps until only one unambiguous path is left
        // returns the extension if successful, returns empty sequence otherwise
        TBases JumpOver(vector<CDBGraph::Successor>& successors, int max_extent, int min_extent) {
            //if(max_extent == 0 || m_hist_min == 0)    
            if(max_extent == 0)
                return TBases();

            TBranch extensions;
            TSeqList sequences; // storage  
            for(auto& suc : successors) {
                sequences.push_front(TSequence(TBases(1,suc), m_graph.Abundance(suc.m_node)));
                extensions[suc.m_node] = make_tuple(sequences.begin(), false);         
            }

            while(!extensions.empty() && extensions.size() < m_max_branch) {
                TSeqList::iterator is = get<0>(extensions.begin()->second);
                int len = get<0>(*is).size();
                if(len == max_extent)
                    break;
            
                OneStepBranchExtend(extensions, sequences);
                if(extensions.empty())  // can't extend 
                    return TBases();

                if(extensions.size() == 1 && len+1 >= min_extent)
                    break;
            }

            if(extensions.size() == 1 && !get<1>(extensions.begin()->second)) {
                TSeqList::iterator is = get<0>(extensions.begin()->second);
                TBases& bases = get<0>(*is);
                bool all_good = true;
                for(auto& base : bases) {
                    if(!GoodNode(base.m_node)) {
                        all_good = false;
                        break;
                    }
                }
                if(all_good)
                    return bases;            
            }

            return TBases();        
        }

        // starting from initial_node assembles the right extension
        pair<TBases, CDBGraph::Node> ExtendToRight(const CDBGraph::Node& initial_node) { // initial_node may be not owned 
            CDBGraph::Node node = initial_node;
            TBases extension;
            int max_extent = m_jump;

            while(true) {
                vector<CDBGraph::Successor> successors = m_graph.GetNodeSuccessors(node);
                FilterNeighbors(successors);
                if(successors.empty())                    // no extensions  
                    break;             

                TBases step;
                if(successors.size() > 1) {                // test for dead end     
                    step = JumpOver(successors, max_extent, 0);
                } else {                                   // simple extension      
                    step.push_back(successors.front());
                }
                if(step.empty())                    // multiple extensions  
                    break;

                bool all_good = true;
                for(auto& s : step) {
                    if(!GoodNode(s.m_node)) {
                        all_good = false;
                        break;
                    }
                }
                if(!all_good)
                    break;

                int step_size = step.size();
        
                CDBGraph::Node rev_node = CDBGraph::ReverseComplement(step.back().m_node);
                vector<CDBGraph::Successor> predecessors = m_graph.GetNodeSuccessors(rev_node);
                FilterNeighbors(predecessors);
                if(predecessors.empty())                     // no extensions   
                    break; 

                TBases step_back;
                if(predecessors.size() > 1 || step_size > 1)
                    step_back = JumpOver(predecessors, max_extent, step_size);
                else
                    step_back.push_back(predecessors.front());

                int step_back_size = step_back.size();
                if(step_back_size < step_size)
                    break;

                bool good = true;
                for(int i = 0; i <= step_size-2 && good; ++i)
                    good = (CDBGraph::ReverseComplement(step_back[i].m_node) == step[step_size-2-i].m_node);
                if(!good)
                    break;

                int overshoot = step_back_size-step_size;  // >= 0  
                if(CDBGraph::ReverseComplement(step_back[step_back_size-1-overshoot].m_node) != node) 
                    break;
            
                if(overshoot > 0) { // overshoot    
                    CDBGraph::Node over_node = CDBGraph::ReverseComplement(step_back.back().m_node);
                    vector<CDBGraph::Successor> oversuc = m_graph.GetNodeSuccessors(over_node);
                    FilterNeighbors(oversuc);
                    if(oversuc.empty())
                        break;
                    TBases step_over;
                    if(oversuc.size() > 1 || overshoot > 1)
                        step_over = JumpOver(oversuc, max_extent, overshoot);
                    else
                        step_over.push_back(oversuc.front());

                    if((int)step_over.size() < overshoot)
                        break;
                    for(int i = 0; i < overshoot && good; ++i)
                        good = (CDBGraph::ReverseComplement(step_over[i].m_node) == step_back[step_back_size-2-i].m_node);
                    if(!good)
                        break;
                }            

                for(auto& s : step) {
                    if(!m_graph.SetVisited(s.m_node))
                        return make_pair(extension, s.m_node);

                    extension.push_back(s);
                }

                node = extension.back().m_node;
            }

            return make_pair(extension, 0);
        }

        // assembles a contig starting from initial_node 
        // min_len - minimal length for accepted contigs
        // changes the state of all used nodes to 'visited' or 'temporary holding'   
        SContig GetContigForKmer(CDBGraph::Node initial_node, int min_len) {
            if(m_graph.Abundance(initial_node) < m_hist_min || !GoodNode(initial_node) || !m_graph.SetVisited(initial_node))
                return SContig();

            //node is good and this thread owns it  

            pair<TBases, CDBGraph::Node> to_right = ExtendToRight(initial_node);
            pair<TBases, CDBGraph::Node> to_left = ExtendToRight(CDBGraph::ReverseComplement(initial_node));
            
            if(!to_left.second && !to_right.second && (int)(to_left.first.size()+m_graph.KmerLen()+to_right.first.size()) < min_len) {
                m_graph.SetVisited(initial_node, 2, 1);
                for(auto& base : to_right.first)
                    m_graph.SetVisited(base.m_node, 2, 1);
                for(auto& base : to_left.first)
                    m_graph.SetVisited(base.m_node, 2, 1);

                return SContig();
            } else {
                return SContig(to_left.first, to_right.first, initial_node, CDBGraph::ReverseComplement(to_left.second), to_right.second, m_graph);
            }
        }
        // Starting from available graph nodes, generates all contigs >= min_len_for_new_seeds. Uses ncores threads.
        TContigList GenerateNewSeeds(int min_len_for_new_seeds, int ncores) {
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

            for(auto iloop = new_seeds.begin(); iloop != new_seeds.end(); ) {
                auto ic = iloop++;
                if((int)ic->Len() < min_len_for_new_seeds) {
                    for(auto& kmer : ic->m_kmers)
                        Graph().ClearVisited(kmer);
                    new_seeds.erase(ic);
                }
            }
        
            return new_seeds;
        }
        // Using a longer kmer generates connectors and extenders and improves previously assembled contigs
        // scontigs - contigs (input/output)
        // scan_window - the size-1 of the contig's flank area used for extensions and connections
        // ncores - number of threads
        void ConnectAndExtenContigs(TContigList& scontigs, int scan_window, int ncores) {
            vector<TContigList> extensions_for_jobs(ncores);
            list<function<void()>> jobs;
            for(auto& ex : extensions_for_jobs) {
                jobs.push_back(bind(&CDBGraphDigger::ExtendContigsJob, this, ref(scontigs), ref(ex), scan_window));
            }
            RunThreads(ncores, jobs);
            TContigList extensions = SContig::ConnectFragments(extensions_for_jobs, Graph()); 
            SContig::ConnectAndExtendContigs(scontigs, extensions);  
        }
        
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
                FilterNeighbors(successors);
                if(find_if(successors.begin(),successors.end(),[node](const CDBGraph::Successor& s){return s.m_node == node;}) == successors.end())                
                    continue;
            
                CDBGraph::Node right = m_graph.ReverseComplement(extended_nodes[lextend.size()+read_pos+1]);
                node = m_graph.ReverseComplement(extended_nodes[read_pos+lextend.size()]);
                if(!right || !GoodNode(right) || !node || !GoodNode(node))
                    continue;
                successors = m_graph.GetNodeSuccessors(right);
                FilterNeighbors(successors);
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
        // scan_window - the size-1 of the contig's flank area used for extensions and connections
        //               (K-1)mers from a window are used as contig sequence near end of contig may not be correct
        void ExtendContigsJob(TContigList& scontigs, TContigList& extensions, int scan_window) {
            int kmer_len = Graph().KmerLen();
            for(auto& contig : scontigs) {
                if(!contig.m_is_taken.Set(1))  // grab contig
                    continue;

                int len = contig.m_seq.size();
                int klen = contig.m_kmers.size();
                int sw = min(scan_window,(len-kmer_len)/2);

                for(int shift = 0; shift <= sw; ++shift) {
                    CDBGraph::Node takeoff_node = contig.m_kmers[klen-shift-1];
                    if(takeoff_node && GoodNode(takeoff_node)) {         // valid kmer 
                        pair<TBases, CDBGraph::Node> extension = ExtendToRight(takeoff_node);
                        if(!extension.first.empty() || extension.second) { // extension could be empty - starting kmer + landing kmer
                            if(shift == 0 || !extension.second || extension.second != contig.m_kmers[klen-shift]) {
                                SContig sc(&contig, shift+1, takeoff_node, extension.first, extension.second, Graph());
                                extensions.push_back(sc);
                            }
                        }                   
                    }
                }
            
                for(int shift = 0; shift <= sw; ++shift) {
                    CDBGraph::Node takeoff_node = CDBGraph::ReverseComplement(contig.m_kmers[shift]);
                    if(takeoff_node && GoodNode(takeoff_node)) {         // valid kmer     
                        pair<TBases, CDBGraph::Node> extension = ExtendToRight(takeoff_node);
                        if(!extension.first.empty() || extension.second) { // extension could be empty - starting kmer + landing kmer   
                            if(shift == 0 || !extension.second  || extension.second != CDBGraph::ReverseComplement(contig.m_kmers[shift-1])) {
                                SContig sc(&contig, -(shift+1), takeoff_node, extension.first, extension.second, Graph());
                                sc.ReverseComplement();
                                extensions.push_back(sc);
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
    };

}; // namespace
#endif /* _GraphDigger_ */
