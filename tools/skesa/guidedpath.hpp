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

#ifndef _GuidedPath_
#define _GuidedPath_

#include "graphdigger.hpp"
#include "glb_align.hpp"
#include <stack> 

using namespace std;
namespace DeBruijn {

    class CGudedPath {
    public:

        struct SPathChunk {
            SPathChunk() : m_score(numeric_limits<int>::min()), m_tlen(0) {}
            int m_score;
            int m_tlen;
            string m_seq;
        };

        CGudedPath(CDBGraph::Node initial_node, const string& target_extension, CDBGraphDigger& graph_digger, SMatrix& delta, int gapopen, int gapextend, int dropoff, int max_len = numeric_limits<int>::max()) : 
            m_snp_detected(false), m_b(target_extension), m_fork_count(0), m_path_end(false), m_num(0), m_graph_digger(graph_digger), m_delta(delta), m_rho(gapopen), m_sigma(gapextend), m_dropoff(dropoff), m_max_len(max_len) {

            int bignegative = numeric_limits<int>::min()/2;
            int nb = m_b.size();
            m_branch.m_sm.resize(nb+1,bignegative);
            m_branch.m_gapb.resize(nb+1,bignegative);

            m_branch.m_jmin = 0;
            m_branch.m_jmax = nb-1;
            m_branch.m_sm[0] = 0;
            m_branch.m_sm[1] = -m_rho-m_sigma;                                        // scores for --------------      
            for(int i = 2; i <= nb && m_branch.m_sm[i-1]-m_sigma > -m_dropoff; ++i) { //            BBBBBBBBBBBBBB      
                m_branch.m_sm[i] = m_branch.m_sm[i-1]-m_sigma;
                m_branch.m_jmax = min(nb-1,i);
            }    
                
            m_branch.m_maxscore = 0;
            m_branch.m_maxposa = -1;
            m_branch.m_maxposb = -1;
            m_branch.m_na = 0;

            m_branch.m_node = initial_node;
            m_last_step_nodes.push_back(m_branch.m_node);
            if(m_branch.m_node) {
                vector<CDBGraph::Successor> neighbors = m_graph_digger.GetReversibleNodeSuccessors(m_branch.m_node);
                if(neighbors.size() == 2 && m_max_len < numeric_limits<int>::max())   // detect SNP
                    DetectSNP(neighbors);
                if(!m_path_end) {
                    for(auto& neighbor : neighbors)
                        m_edges.push(neighbor);
                    if(neighbors.size() > 1) {
                        int fcount = neighbors.size()-1;
                        m_forks.push(make_pair(m_branch, fcount));
                        m_fork_count += fcount;
                    }
                }
            }        
            m_path_end = m_edges.empty();
        }
        bool ProcessNextEdge() {
            if(!m_snp_detected || m_edges.empty() || m_path_end)
                m_last_step_nodes.clear();

            m_snp_detected = false;

            if(m_edges.empty())
                return false;            
        
            ++m_num;

            if(m_path_end) 
                DeleteLastBranch();

            m_branch.m_node = m_edges.top().m_node;
            m_last_step_nodes.push_back(m_branch.m_node);
            m_path_end = !AddOneBase(m_edges.top().m_nt) || (m_branch.m_na == m_max_len);
            m_edges.pop();
            if(m_path_end)
                return true;

            vector<CDBGraph::Successor> neighbors = m_graph_digger.GetReversibleNodeSuccessors(m_branch.m_node);
            m_path_end = neighbors.empty();
            if(m_path_end)
                return true;

            if(neighbors.size() == 2 && m_max_len < numeric_limits<int>::max()) {  // detect SNP
                DetectSNP(neighbors);
                if(m_path_end)
                    return true;
            }
            
            for(auto& neighbor : neighbors)
                m_edges.push(neighbor);
            if(neighbors.size() > 1) {
                int fcount = neighbors.size()-1;
                m_forks.push(make_pair(m_branch, fcount));
                m_fork_count += fcount;
            }            

            return true;
        }
        SPathChunk GetBestPart() const {
            SPathChunk rslt;
            rslt.m_score = m_branch.m_maxscore;
            rslt.m_tlen = m_branch.m_maxposb+1;
            rslt.m_seq = m_a.substr(0, m_branch.m_maxposa+1);
            return rslt;
        }
        string GetAssembledSeq() const { return m_a; }

        bool PathEnd() const { return m_path_end; }
        int Num() const { return m_num; }
        int ForkCount() const { return m_fork_count; }
        list<CDBGraph::Node> LastStepNodes() const { return m_last_step_nodes; }

    private:

        struct SBranch {
            vector<int> m_sm;    // best scores in previous a-row
            vector<int> m_gapb;  // best score with b-gap
            CDBGraph::Node m_node;         // last node
            int m_na;
            int m_maxscore;
            int m_maxposa;
            int m_maxposb;
            int m_jmin;          // b-interval evaluted using prevous row results
            int m_jmax;
        };

        bool AddOneBase(char c) {
            m_a.push_back(c);
            ++m_branch.m_na;
            return UpdateScore();
        }
        bool UpdateScore() {
            int bignegative = numeric_limits<int>::min()/2;
            int nb = m_b.size();
            int rs = m_rho+m_sigma;
            int next_jmax = -1;
            int next_jmin = nb;

            vector<int> s(nb+1,bignegative);             // best scores in current a-raw    
            if(-m_rho-m_branch.m_na*m_sigma > m_branch.m_maxscore-m_dropoff) {
                next_jmin = 0;
                s[0] = -m_rho-m_branch.m_na*m_sigma;      // score for AAAAAAAAAAA  
            }                                             //           -----------  
            
            int gapa = bignegative;
            int ai = m_a[m_branch.m_na-1];
            const char* matrix = (m_delta.matrix)[ai];
            int* sp = &s[m_branch.m_jmin];

            for(int j = m_branch.m_jmin; j <= m_branch.m_jmax; ) { // here j is 'real' position in b
                int ss = m_branch.m_sm[j]+matrix[(int)m_b[j]];

                gapa -= m_sigma;
                if(*sp-rs > gapa)
                    gapa = *sp-rs;
			
                int& gapbj = m_branch.m_gapb[++j]; // here j is one-shifted to account for extra element in vectors (nb+1)
                gapbj -= m_sigma;
                if(m_branch.m_sm[j]-rs > gapbj)
                    gapbj = m_branch.m_sm[j]-rs;

                if(gapa > gapbj) {
                    if(ss > gapa) {
                        *(++sp) = ss;
                        if(ss > m_branch.m_maxscore) {
                            m_branch.m_maxscore = ss;
                            m_branch.m_maxposa = m_branch.m_na-1;
                            m_branch.m_maxposb = j-1;
                        }
                    } else {
                        *(++sp) = gapa;
                    }
                } else {
                    if(ss > gapbj) {
                        *(++sp) = ss;
                        if(ss > m_branch.m_maxscore) {
                            m_branch.m_maxscore = ss;
                            m_branch.m_maxposa = m_branch.m_na-1;
                            m_branch.m_maxposb = j-1;
                        }
                    } else {
                        *(++sp) = gapbj;
                    }
                }

                if(max(*sp,gapbj) > m_branch.m_maxscore-m_dropoff) {
                    next_jmin = min(next_jmin, j-1);
                    next_jmax = min(nb-1,j);
                }
            }
            swap(m_branch.m_sm,s);
            // jmin never decreases
            m_branch.m_jmin = next_jmin;
            //right may decrease
            for(int l = next_jmax+1; l <= m_branch.m_jmax; ++l) {
                m_branch.m_gapb[l+1] = bignegative;
                m_branch.m_sm[l+1] = bignegative;
            }
            m_branch.m_jmax = next_jmax;

            return m_branch.m_jmax >= m_branch.m_jmin;
        }
        void DeleteLastBranch() {
            if(m_edges.empty())
                return;    

            if(!m_path_end) {
                int neigbors_num = 1;
                if(!m_forks.empty() && m_forks.top().first.m_na == m_branch.m_na) {  // last base is a fork 
                    neigbors_num += m_forks.top().second;
                    m_fork_count -= m_forks.top().second;
                    m_forks.pop();
                }
                while(neigbors_num-- > 0)
                    m_edges.pop();    
            }

            if(!m_forks.empty()) {
                m_path_end = false;
                m_branch = m_forks.top().first;
                --m_fork_count;
                m_a.resize(m_branch.m_na);
                if(--m_forks.top().second == 0)
                    m_forks.pop(); 
            } else {
                m_a.clear();
            }
        }

        bool CheckIfSNP(CDBGraph::Node nodea, CDBGraph::Node nodeb, vector<CDBGraph::Successor>& patha, vector<CDBGraph::Successor>& pathb) {
            int kmer = m_graph_digger.Graph().KmerLen();
            patha.reserve(kmer);
            pathb.reserve(kmer);
            for(int s = 0; s < kmer; ++s) {
                vector<CDBGraph::Successor> nbra = m_graph_digger.GetReversibleNodeSuccessors(nodea);
                if(nbra.size() != 1)
                    break;
                vector<CDBGraph::Successor> nbrb = m_graph_digger.GetReversibleNodeSuccessors(nodeb);
                if(nbrb.size() != 1)
                    break;
                patha.push_back(nbra.front());
                nodea = nbra.front().m_node;
                pathb.push_back(nbrb.front());
                nodeb = nbrb.front().m_node;                             
            }

            return ((int)patha.size() == kmer && nodea == nodeb);
        }

        void DetectSNP(vector<CDBGraph::Successor>& neighbors) {
            int kmer = m_graph_digger.Graph().KmerLen();
            vector<CDBGraph::Successor> patha;
            vector<CDBGraph::Successor> pathb;
            CDBGraph::Node nodea = neighbors[0].m_node;
            CDBGraph::Node nodeb = neighbors[1].m_node;
            if(!CheckIfSNP(nodea, nodeb, patha, pathb))
                return;
                        
            { // check in reverse
                nodea = m_graph_digger.Graph().ReverseComplement(patha[kmer-2].m_node);
                nodeb = m_graph_digger.Graph().ReverseComplement(pathb[kmer-2].m_node);
                vector<CDBGraph::Successor> pathaa;
                vector<CDBGraph::Successor> pathbb;
                if(!CheckIfSNP(nodea, nodeb, pathaa, pathbb))
                    return;
            }
            
            m_snp_detected = true;
            m_branch.m_node = neighbors[0].m_node;
            m_last_step_nodes.push_back(neighbors[0].m_node);
            m_last_step_nodes.push_back(neighbors[1].m_node);
            string ambig{neighbors[0].m_nt, neighbors[1].m_nt};
            m_path_end = !AddOneBase(ToAmbiguousIUPAC[AmbiguousString(ambig)]);
            if(m_path_end)
                return;
            for(int s = 0; s < kmer-1; ++s) {
                m_branch.m_node = patha[s].m_node;
                m_last_step_nodes.push_back(patha[s].m_node);
                m_last_step_nodes.push_back(pathb[s].m_node);
                m_path_end = !AddOneBase(patha[s].m_nt);
                if(m_path_end)
                    return;
            }
            neighbors.clear();
            neighbors.push_back(patha.back());
        }

        SBranch m_branch;
        list<CDBGraph::Node> m_last_step_nodes;
        bool m_snp_detected;
        string m_a;
        string m_b;
        stack<pair<SBranch,int> > m_forks;
        stack<CDBGraph::Successor> m_edges;
        int m_fork_count;
        bool m_path_end;
        int m_num;

        CDBGraphDigger m_graph_digger; 
        SMatrix& m_delta;
        int m_rho;
        int m_sigma;
        int m_dropoff;
        int m_max_len;
    };


}; // namespace
#endif /* _GuidedPath_ */
