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
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/seek.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <iomanip>  

#include "DBGraph.hpp"
#include "glb_align.hpp"

using namespace boost::program_options;
using namespace DeBruijn;

class CwgMLST {
public:
    typedef list<tuple<string,string>> TLocus;
    typedef map<string,TLocus> TAlles;       // [locus], list<allele,sequence>    
    typedef vector<tuple<string,string>> TGenome; // <accession,contig>
    typedef CKmerMap<list<tuple<int, int, int>>> TKmerToGenome;  // tuple<contig position, strand, contig's vector index>; list for repeated kmers

    CwgMLST (boost::iostreams::filtering_istream& alleles_file, ofstream& output_mappings, ofstream& output_loci, int ncores) : 
        m_output_mappings(output_mappings), m_output_loci(output_loci), m_ncores(ncores), m_delta(m_match, m_mismatch) {  
        // read alleles
        string locus;
        while(alleles_file >> locus) {
            string allele;
            string seq;
            if(!(alleles_file >> allele >> seq))
                throw runtime_error("Invalid alleles file");

            //convert to upper case
            for(char& c : seq) c = toupper(c);

            // check valid sequence
            if(seq.find_first_not_of("ACGTYRWSKMDVHBXN-") != string::npos)
                throw runtime_error("Invalid sequence in the alleles file");

            m_alleles[locus].push_back(make_tuple(allele,seq));
        } 
    }

    void CheckAndCleanAlleles() {
        vector<pair<TLocus*, SAtomic<uint8_t>>> loci;
        for(auto& locus : m_alleles)
            loci.push_back(make_pair(&locus.second,SAtomic<uint8_t>(0)));       

        size_t before = 0;
        for(auto& loc : m_alleles)
            before += loc.second.size();

        list<function<void()>> jobs;
        for(int thr = m_ncores; thr > 0; --thr) 
            jobs.push_back(bind(&CwgMLST::CheckAndCleanAllelesJob, this, ref(loci)));
        RunThreads(m_ncores, jobs);

        vector<string> removed_loci;
        size_t after = 0;        
        for(auto& loc : m_alleles) {
            if(loc.second.empty())
                removed_loci.push_back(loc.first);
            else
                after += loc.second.size();
        }
        if(before > after) {
            cerr << "Removed " << before-after << " alleles" << endl;
            cerr << "Completely removed " << removed_loci.size() << " loci:"  << endl;
            for(auto& locus : removed_loci) {
                cerr << locus << endl;
                m_alleles.erase(locus);
            }
        }
    }

    void ReadGenome(boost::iostreams::filtering_istream& genome_file) {
        m_genome.clear();

        // read genome
        char c;
        if(!(genome_file >> c) || c != '>')
            throw runtime_error("Invalid fasta file format for genome");

        string record;
        while(getline(genome_file, record, '>')) {
            size_t first_ret = min(record.size(), record.find('\n'));
            if(first_ret == string::npos)
                throw runtime_error("Invalid fasta file format for genome");

            size_t first_space = min(first_ret, record.find_first_of(" \t"));
            string acc = record.substr(0, first_space);
            string contig = record.substr(first_ret+1);
            contig.erase(remove(contig.begin(),contig.end(),'\n'),contig.end());
            //convert to upper case
            for(char& c : contig) c = toupper(c);
            if(contig.find_first_not_of("ACGTYRWSKMDVHBXN-") != string::npos)
                throw runtime_error("Invalid fasta file format for genome");
            m_genome.push_back(make_tuple(acc, contig));            
            m_genome_acc_to_index[acc] = m_genome.size();
        }
    }

    void ReadBadBases(ifstream& bad_bases_file) {
        m_bad_bases.clear();
        string line;
        int line_num = 0;
        while(getline(bad_bases_file, line)) {
            ++line_num;
            istringstream ss(line);
            vector<string> tokens;
            string tk;
            while(ss >> tk)
                tokens.push_back(tk);
            if(tokens.size() < 2)
                throw runtime_error("Invalid format for bad bases in line "+to_string(line_num));
            string::size_type idx;
            int base = stoi(tokens[1], &idx);
            if(idx != tokens[1].size())
                throw runtime_error("Invalid format for bad bases in line "+to_string(line_num));
            string& contig_acc = tokens[0];
            int index = m_genome_acc_to_index[contig_acc]-1; // it starts from 1
            if(index < 0)
                throw runtime_error("Unknown genome accession in bad bases file");
            if(!m_bad_bases[index].insert(base-1).second)
                throw runtime_error("Do not expect to have same position present multiple times in bad bases input");
        }
    }

    bool AllGoodBases(int contig_index, int from, int to) {
        if(!m_bad_bases.count(contig_index))
            return true;

        auto& bad_bases = m_bad_bases[contig_index];
        auto lb = bad_bases.lower_bound(from);
        return (lb == bad_bases.end() || *lb > to);
    }

    void ReadBlastHits(boost::iostreams::filtering_istream& blast_file) {
        m_blast_hits.clear();

        // read blast hits
        string qid, sid;
        double ident;
        int length, mism, gap, qs, qe, ss, se;
        double score, evalue;
        string btop;

        // allele - query
        // contig - subject 
        // btop is ALREADY REVERSED for minus strand
        while(blast_file >> qid >> sid >> ident >> length >> mism >> gap >> qs >> qe >> ss >> se >> evalue >> score >> btop) {
            string remain;
            getline(blast_file, remain);

            if(!m_alleles.count(qid))
                continue;

            int strand = 1;
            if(ss > se) {
                swap(ss, se);
                swap(qs, qe);
                int qlen = get<1>(m_alleles[qid].front()).size();
                qs = qlen-qs+1;  // coordinates on REVERSED query
                qe = qlen-qe+1;
                strand = -1;
            }
            int matches = ident*length/100+0.5;
            int contig_key = m_genome_acc_to_index[sid]*strand;
            if(contig_key == 0)
                throw runtime_error("Unknown genome accession");
            m_blast_hits[qid][contig_key].push_back(SLinkedHit(qs-1, qe-1, ss-1, se-1, matches, contig_key, btop)); 

            /* debug output
            const string& contig = get<1>(m_genome[abs(contig_key)-1]);
            string query = get<1>(m_alleles[qid].front());
            if(strand < 0)
                ReverseComplementSeq(query.begin(), query.end());
            CCigar cigar(qs-2, ss-2);
            for(auto& elem : m_blast_hits[qid][contig_key].back().m_btop) {
                cigar.PushBack(elem);
                cerr << elem.m_len << elem.m_type;
            }
            cerr << endl;
            cerr << qid << " " << sid << " " << qs << " " << qe << " " << ss << " " << se << " " << strand << " " << query.size() << " " << m_blast_hits[qid][contig_key].back().m_btop.size() << " " << btop << endl;
            cerr << cigar.CigarString(0, query.size()) << endl;

            TCharAlign align = cigar.ToAlign(query.c_str(), contig.c_str());
            for(unsigned i = 0; i < align.first.size(); ++i)
                cerr << align.first[i];
            cerr << endl;
            for(unsigned i = 0; i < align.first.size(); ++i) {
                if(align.first[i] == align.second[i])
                    cerr << "|";
                else
                    cerr << " ";
            }
            cerr << endl;
            for(unsigned i = 0; i < align.first.size(); ++i)
                cerr << align.second[i];
            cerr << endl << endl;                            
            */
        }
        if(m_blast_hits.empty())
            throw runtime_error("Blast hits don't correspond to alleles");
    }

    void PrepareKmerMap(int kmer_len) {
        m_genome_kmers = TKmerToGenome(kmer_len);

        size_t total = 0;
        for(auto& genome : m_genome) {
            const string& contig = get<1>(genome);
            if((int)contig.size() >= kmer_len)
                total += contig.size()-kmer_len+1;
        }
        m_genome_kmers.Reserve(total);

        for(int k = 0; k < (int)m_genome.size(); ++k) {
            const string& contig = get<1>(m_genome[k]);
            size_t start = 0;
            while(start < contig.size()) {
                size_t stop = min(contig.size(),contig.find_first_not_of("ACGT", start));
                int len = stop-start;
                if(len >= kmer_len) {
                    CReadHolder rh(false);
                    rh.PushBack(contig.substr(start, len));
                    int pos = stop-kmer_len;
                    for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --pos) { // iteration from last kmer to first      
                        TKmer kmer = *ik;
                        TKmer rkmer = revcomp(kmer, kmer_len);
                        if(kmer < rkmer)
                            m_genome_kmers[kmer].push_back(make_tuple(pos, +1, k));
                        else
                            m_genome_kmers[rkmer].push_back(make_tuple(pos, -1, k));
                    }
                }
                start = contig.find_first_of("ACGT", stop);
            }
        }
    }

    void AnalyzeAlleles(int min_kmer_bases, double min_fraction_of_matches, int match, int mismatch, int gap_open, int gap_extend) {
        m_min_bases_to_consider = min_kmer_bases;
        m_min_deep_search_fraction = min_fraction_of_matches;
        m_match = match;
        m_mismatch = mismatch;
        m_gap_open = gap_open;
        m_gap_extend = gap_extend;

        vector<tuple<string,list<SFragment>,SAtomic<uint8_t>>> loci;
        for(auto& loc : m_alleles) {
            if(m_blast_hits.empty() || m_blast_hits.count(loc.first))
                loci.push_back(make_tuple(loc.first,list<SFragment>(), 0));
        }
        
        list<function<void()>> jobs;
        for(int thr = m_ncores; thr > 0; --thr) 
            jobs.push_back(bind(&CwgMLST::AnalyzeAllelesJob, this, ref(loci)));
        RunThreads(m_ncores, jobs);

        if(m_output_loci.is_open()) {
            // remove partials and overlapping with partials
            vector<pair<int,int>> partial_margins(m_genome.size());
            for(unsigned index = 0; index < m_genome.size(); ++index)
                partial_margins[index] = make_pair(-1, get<1>(m_genome[index]).size()); 
            for(auto& loc : loci) {
                list<SFragment>& lst = get<1>(loc);
                for(auto it_loop = lst.begin(); it_loop != lst.end(); ) {
                    auto it = it_loop++;
                    bool left_partial = (it->m_strand > 0 && it->m_not_aligned_fivep > 0) || (it->m_strand < 0 && it->m_not_aligned_threep > 0);
                    bool right_partial = (it->m_strand < 0 && it->m_not_aligned_fivep > 0) || (it->m_strand > 0 && it->m_not_aligned_threep > 0);
                    if(left_partial || right_partial) {
                        if(left_partial) {
                            partial_margins[it->m_contig].first = max(partial_margins[it->m_contig].first, it->m_to);
                            if(partial_margins[it->m_contig].second <= it->m_to)
                                partial_margins[it->m_contig].second = 0;
                        }
                        if(right_partial) {
                            partial_margins[it->m_contig].second = min(partial_margins[it->m_contig].second, it->m_from);
                            if(partial_margins[it->m_contig].first >= it->m_from)
                                partial_margins[it->m_contig].first = get<1>(m_genome[it->m_contig]).size(); 
                        }
                        lst.erase(it);
                    }
                }
            }
            for(auto& loc : loci) {
                list<SFragment>& lst = get<1>(loc);
                for(auto it_loop = lst.begin(); it_loop != lst.end(); ) {
                    auto it = it_loop++;
                    pair<int,int>& margins = partial_margins[it->m_contig];
                    int from = it->m_from;
                    int to = it->m_to;
                    double critical_overlap = min((double)m_consistent_bases, m_percent_contained*(to-from+1));
                    if(margins.first-from+1 >= critical_overlap) {
                        cerr << "Removing " << get<0>(m_genome[it->m_contig]) << " " << from+1 << " " << to+1 << " because of " << 1 << " " << margins.first << endl;
                        lst.erase(it);
                    } else if(to-margins.second+1 >= critical_overlap) {
                        cerr << "Removing " << get<0>(m_genome[it->m_contig]) << " " << from+1 << " " << to+1 << " because of " << margins.second+1 << " " << get<1>(m_genome[it->m_contig]).size() << endl;
                        lst.erase(it);
                    }                                                                    
                }
            }
            
            list<tuple<string,string,SFragment>> overlap_data; // 'allele', 'locus', fragment - apparently now allele and locus are swapped in the schema???
            for(auto& loc : loci) {
                for(auto& fragment : get<1>(loc)) {
                    string& locus = get<0>(loc);
                    string& allele = get<0>(m_alleles[locus].front());
                    overlap_data.push_back(make_tuple(allele, locus, fragment)); // we swapped allele and locus here
                }
                get<1>(loc).clear();
            }
            overlap_data.sort(SCompareOverlapsByLocation);            

            /* remove shorter overlapping ones */
            for(auto ind = overlap_data.begin(); ind != overlap_data.end(); ++ind) {
                string& ilocus = get<0>(*ind);
                SFragment& ifrag = get<2>(*ind);                
                for(auto jnd = ind; ++jnd != overlap_data.end(); ) {
                    string& jlocus = get<0>(*jnd);
                    SFragment& jfrag = get<2>(*jnd);
                    if(ifrag.m_contig != jfrag.m_contig || jfrag.m_from > ifrag.m_to)
                        break;

                    /* long enough overlap or one strictly contained and shorter than other */
                    if((ifrag.m_to-jfrag.m_from+1 > m_consistent_bases && jfrag.m_from > ifrag.m_from) ||
                       (jfrag.m_from > ifrag.m_from && jfrag.m_to <= ifrag.m_to) ||
                       (jfrag.m_from == ifrag.m_from && jfrag.m_to < ifrag.m_to)) {
                        int ilength = ifrag.m_to-ifrag.m_from+1;
                        int jlength = jfrag.m_to-jfrag.m_from+1;
                        if(ilength < jlength)
                            ilocus.clear();
                        if(jlength < ilength)
                            jlocus.clear();                        
                    }
                }
            }
            for(auto ind = overlap_data.begin(); ind != overlap_data.end(); ) {
                auto i = ind++;
                if(get<0>(*i).empty())
                    overlap_data.erase(i);
            }

            /* All overlapping ones left now are of same length but may still be only partially overlapping */
            
            /* don't care about allele name any more -- keep one of the best identity one for same location and locus */
            for(auto ind = overlap_data.begin(); ind != overlap_data.end(); ++ind) {
                SFragment& ifrag = get<2>(*ind);
                auto jnd = ind;
                for(++jnd; jnd != overlap_data.end() && get<0>(*ind) == get<0>(*jnd); ) {
                    SFragment& jfrag = get<2>(*jnd);
                    auto j = jnd++;
                    if(ifrag.m_contig == jfrag.m_contig && ifrag.m_from == jfrag.m_from && ifrag.m_to == jfrag.m_to)
                        overlap_data.erase(j);                    
                    else
                        break;                    
                }
            }
            
            /* find overlapping ones more than expected */
            int same_count = 0;
            int different_count = 0;
            for(auto ind = overlap_data.begin(); ind != overlap_data.end(); ++ind) {
                string& ilocus = get<0>(*ind);
                string& iallele = get<1>(*ind);
                SFragment& ifrag = get<2>(*ind);
                for(auto jnd = ind; ++jnd != overlap_data.end(); ) {
                    string& jlocus = get<0>(*jnd);
                    string& jallele = get<1>(*jnd);
                    SFragment& jfrag = get<2>(*jnd);
                    if(ifrag.m_contig != jfrag.m_contig || jfrag.m_from > ifrag.m_to)
                        break;

                    int length;
                    if(ifrag.m_to <= jfrag.m_to)
                        length = ifrag.m_to-jfrag.m_from+1;
                    else
                        length = jfrag.m_to-jfrag.m_from+1;
                    if(length > m_consistent_bases ||(ifrag.m_from == jfrag.m_from && ifrag.m_to == jfrag.m_to)) {
                        ifrag.m_print = false;
                        jfrag.m_print = false;

                        double first_percent = 100.*length/(ifrag.m_to-ifrag.m_from+1);
                        double second_percent = 100.*length/(jfrag.m_to-jfrag.m_from+1);
                        if(first_percent > second_percent)
                            swap(first_percent, second_percent);
                        
                        if(ifrag.m_from == jfrag.m_from && ifrag.m_to == jfrag.m_to)
                            cerr << "BOTH_SAME\t";
                        else if(ifrag.m_from == jfrag.m_from || ifrag.m_to == jfrag.m_to)
                            cerr << "ONE_SAME\t";
                        else if((ifrag.m_from <= jfrag.m_from && ifrag.m_to >= jfrag.m_to) ||
                                (ifrag.m_from >= jfrag.m_from && ifrag.m_to <= jfrag.m_to))
                            cerr << "CONTAINED\t";
                        else
                            cerr << "PARTIAL\t";

                        cerr << "OVERLAP\t";
                        cerr << length << '\t';
                        cerr << fixed << setprecision(2) << first_percent << '\t';
                        cerr << fixed << setprecision(2) << second_percent << '\t';
                        cerr << get<0>(m_genome[ifrag.m_contig]) << '\t';
                        cerr << ifrag.m_from+1 << '\t';
                        cerr << ifrag.m_to+1 << '\t';
                        cerr << fixed << setprecision(2) << ifrag.m_identity << '\t';
                        cerr << iallele << '\t';                 
                        cerr << ilocus << '\t';
                        cerr << jfrag.m_from+1 << '\t';
                        cerr << jfrag.m_to+1 << '\t';
                        cerr << fixed << setprecision(2) << jfrag.m_identity << '\t';
                        cerr << jallele << '\t';             
                        cerr << jlocus << endl;

                        if(ilocus == jlocus)
                            ++same_count;
                        else
                            ++different_count;
                    }
                }                
            }
            cerr << same_count+different_count << " overlapping regions with " << same_count << " from same locus and " << different_count << " from different loci" << endl;
            
            overlap_data.sort(SCompareOverlapsByLocus);
            /* print loci with multiple placements */
            same_count = 0;
            for(auto ind = overlap_data.begin(); ind != overlap_data.end(); ) {
                string& ilocus = get<0>(*ind);
                SFragment& ifrag = get<2>(*ind);
                auto jnd = ind;
                bool multiple = false;
                while(++jnd != overlap_data.end()) {
                    string& jlocus = get<0>(*jnd);
                    SFragment& jfrag = get<2>(*jnd);
                    if(ilocus == jlocus) {
                        ifrag.m_print = false;
                        jfrag.m_print = false;
                        multiple = true;
                    } else {
                        break;
                    }
                }
                if(multiple) {
                    for(auto i = ind; i != jnd; ++i) {
                        SFragment& fragment = get<2>(*i); 
                        cerr << "MULTIPLE\t";
                        cerr << get<0>(m_genome[fragment.m_contig]) << '\t';
                        cerr << fragment.m_from+1 << '\t';
                        cerr << fragment.m_to+1 << '\t';
                        cerr << get<0>(*i) << '\t';
                        cerr << get<1>(*i) << '\t';
                        cerr << fixed << setprecision(2) << fragment.m_identity << '\t';
                        cerr << fragment.m_to-fragment.m_from+1 << endl;
                    }
                    ++same_count;
                }
                ind = jnd;
            }
            cerr << same_count << " loci present in multiple places on the genome" << endl;

            /* check if region overlaps bad bases */
            int bad_count = 0;
            for(auto i = overlap_data.begin(); i != overlap_data.end(); ++i) {
                SFragment& fragment = get<2>(*i); 
                if(fragment.m_print && !AllGoodBases(fragment.m_contig, fragment.m_from, fragment.m_to)) {
                    ++bad_count;
                    fragment.m_print = false;
                    cerr << "BAD_BASE\t";
                    cerr << get<0>(m_genome[fragment.m_contig]) << '\t';
                    cerr << fragment.m_from+1 << '\t';
                    cerr << fragment.m_to+1 << '\t';
                    cerr << get<0>(*i) << '\t';
                    cerr << get<1>(*i) << '\t';
                    cerr << fixed << setprecision(2) << fragment.m_identity << '\t';
                    cerr << fragment.m_to-fragment.m_from+1 << endl;
                }
            }
            cerr << bad_count << " loci overlap bad bases on the genome" << endl;

            overlap_data.sort(SCompareOverlapsByLocation);            
            /* print kept locations to m_output_loci : results */
            /* print removed equally good locations to stderr */
            for(auto i = overlap_data.begin(); i != overlap_data.end(); ++i) {
                SFragment& fragment = get<2>(*i); 
                ostream& out_loci = fragment.m_print ? m_output_loci : cerr;
                out_loci << get<0>(*i) << '\t';
                out_loci << fragment.m_seq << '\t';
                out_loci << get<0>(m_genome[fragment.m_contig]) << '\t';
                out_loci << fragment.m_from+1 << "\t";                                                                         
                out_loci << fragment.m_to+1 << "\t";
                out_loci << fixed << setprecision(2) << fragment.m_identity << '\t';
                out_loci << fragment.m_to-fragment.m_from+1 << "\t";
                out_loci << get<1>(*i) << endl;
            }

            /* get total coverage on genome in loci kept;
            * we may still have pairs of hits with overlaps that have are at most CONSISTENT_BASES
            * doing it destructively so save if original ones needed
            */
            
            for(auto ind = overlap_data.begin(); ind != overlap_data.end(); ++ind) {
                SFragment& ifrag = get<2>(*ind);
                if(ifrag.m_print) {
                    for(auto jnd = ind; ++jnd != overlap_data.end(); ) {
                        SFragment& jfrag = get<2>(*jnd);
                        if(ifrag.m_contig != jfrag.m_contig || jfrag.m_from > ifrag.m_to)
                            break;
                        else if(jfrag.m_print)
                            jfrag.m_from = ifrag.m_to+1;
                    }
                }
            }
            int length = 0;
            for(auto ind = overlap_data.begin(); ind != overlap_data.end(); ++ind) {
                SFragment& ifrag = get<2>(*ind);
                if(ifrag.m_print && ifrag.m_from <= ifrag.m_to)
                    length += ifrag.m_to-ifrag.m_from+1;
            }
            int total = 0;
            for(auto& contig : m_genome)
                total += get<1>(contig).size();
            cerr << fixed << setprecision(2) << 100.*length/total << " coverage at " << length << " out of total " << total << " genome length" << endl;
        }        
    }

private:

    struct SHit {
        SHit(int qfrom, int qto, int sfrom, int sto, int score, int contig_key, const string& btop = "") : m_qfrom(qfrom), m_qto(qto), m_sfrom(sfrom), m_sto(sto), m_score(score), m_contig_key(contig_key) {
            istringstream istr_btop(btop);
            char c;
            string digits("0123456789");
            while(istr_btop >> c) {
                istr_btop.putback(c);
                if(digits.find(c) != string::npos) {
                    int l = 0;
                    istr_btop >> l;
                    if(m_btop.empty() || m_btop.back().m_type != 'M')
                        m_btop.push_back(CCigar::SElement(l, 'M'));
                    else
                        m_btop.back().m_len += l; 
                } 
                char a, b;
                if(istr_btop >> a >> b) {
                    char type = 'M';
                    if(a == '-')
                        type = 'D';
                    else if(b == '-')
                        type = 'I';
                    if(m_btop.empty() ||  m_btop.back().m_type != type)
                        m_btop.push_back(CCigar::SElement(1, type));
                    else
                        ++m_btop.back().m_len;
                }                
            }
        }
        bool operator<(const SHit& other) const {
            if(m_contig_key != other.m_contig_key)
                return m_contig_key < other.m_contig_key;
            else if(m_sfrom != other.m_sfrom)
                return m_sfrom < other.m_sfrom;
            else
                return m_qfrom < other.m_qfrom;
        }
        
        int m_qfrom;
        int m_qto;
        int m_sfrom;
        int m_sto;
        int m_score;
        int m_contig_key;
        list<CCigar::SElement> m_btop;        
    };
    struct SLinkedHit : public SHit {
        SLinkedHit(int qfrom, int qto, int sfrom, int sto, int score, int contig_key, const string& btop = "") : SHit(qfrom, qto, sfrom, sto, score, contig_key, btop), m_left(nullptr) {}
        SLinkedHit* m_left;
    };

    typedef vector<SHit> THits;
    typedef list<vector<SHit>> THitTandems;
    typedef vector<SLinkedHit> TLinkedHits;

    struct SFragment {
        string m_seq;
        int m_contig = numeric_limits<int>::max();
        int m_strand = 0;
        int m_from = numeric_limits<int>::max();
        int m_to = numeric_limits<int>::min();
        int m_not_aligned_fivep = 0;
        int m_not_aligned_threep = 0;
        bool m_consistent = false;
        bool m_has_start = false;
        bool m_has_stop = false;
        int m_matches = 0;
        string m_allele_id;
        double m_identity = 0;
        bool m_print = true;
    };

    struct SLocusRslt {
        SLocusRslt(const string& locus) : m_locus_id(locus) {}
            
        list<SFragment> m_fragments;
        const string& m_locus_id;
    };

    static bool SCompareOverlapsByLocation (const tuple<string,string,SFragment>& a, const tuple<string,string,SFragment>& b) 
    {
        const string& loca = get<0>(a);
        const string& locb = get<0>(b);
        const string& allelea = get<1>(a);
        const string& alleleb = get<1>(b);
        const SFragment& afrg = get<2>(a);
        const SFragment& bfrg = get<2>(b);
        
        /* genomic location and then locus and then identity */

        if(afrg.m_contig != bfrg.m_contig)
            return afrg.m_contig < bfrg.m_contig;
        if(afrg.m_from != bfrg.m_from)
            return afrg.m_from < bfrg.m_from;
        if(afrg.m_to != bfrg.m_to)
            return afrg.m_to > bfrg.m_to;

        if(loca != locb)
            return loca < locb;

        if(afrg.m_identity != bfrg.m_identity)
            return afrg.m_identity > bfrg.m_identity;

        return allelea < alleleb;
    }

    static bool SCompareOverlapsByLocus (const tuple<string,string,SFragment>& a, const tuple<string,string,SFragment>& b) 
    {
        const string& loca = get<0>(a);
        const string& locb = get<0>(b);
        const string& allelea = get<1>(a);
        const string& alleleb = get<1>(b);
        const SFragment& afrg = get<2>(a);
        const SFragment& bfrg = get<2>(b);

        if(loca != locb)
            return loca < locb;
        if(allelea != alleleb)
            return allelea < alleleb;

        if(afrg.m_identity != bfrg.m_identity)
            return afrg.m_identity > bfrg.m_identity;
        if(afrg.m_contig != bfrg.m_contig)
            return afrg.m_contig < bfrg.m_contig;
        if(afrg.m_from != bfrg.m_from)
            return afrg.m_from < bfrg.m_from;

        return afrg.m_to > bfrg.m_to;
    }

    void PrintAlignments(SLocusRslt& locus_rslt) const {
        const string& locus_id = locus_rslt.m_locus_id;
        ostream& out = m_output_mappings.is_open() ? m_output_mappings : cout;
        for(auto& fragment : locus_rslt.m_fragments) {
            out << locus_id << "\t";                                                                                   // locus
            out << (fragment.m_allele_id.empty() ? "-" : fragment.m_allele_id) << "\t";                                // allele
            if(fragment.m_not_aligned_fivep > 0 || fragment.m_not_aligned_threep > 0)
                out << "Partial\t";                                                                                    // partial
            else
                out << (fragment.m_consistent ? "Consistent" : "Notconsistent") << "\t";                               // consistency 
            if(fragment.m_has_start)
                out << "Start\t";                                                                                      // start 
            else
                out << fragment.m_not_aligned_fivep << "\t";                                                           // not aligned 5'
            if(fragment.m_has_stop)
                out << "Stop\t";                                                                                       // stop 
            else
                out << fragment.m_not_aligned_threep << "\t";                                                          // not aligned 3'
            out << get<0>(m_genome[fragment.m_contig]) << "\t";                                                        // contig
            out << fragment.m_from+1 << "\t";                                                                          // from on contig
            out << fragment.m_to+1 << "\t";                                                                            // to on contig
            out << (fragment.m_strand > 0 ? "+" : "-")  << "\t";                                                       // strand
            out << fragment.m_seq  << "\t";                                                                            // sequence
            out << fragment.m_matches << endl;
        }    
    }

    void AnalyzeAllelesJob(vector<tuple<string, list<SFragment>, SAtomic<uint8_t>>>& loci) {
        for(auto& loc : loci) {
            if(!get<2>(loc).Set(1))
                continue;

            const string& locus_id = get<0>(loc);
            const TLocus& locus = m_alleles[locus_id];
            auto& alleleid_seq = locus.front();
            const string& seq = get<1>(alleleid_seq);
            int qlen = seq.size();

            THitTandems tandems;
            if(m_blast_hits.empty())   // find hits from scratch
                tandems = GetMatchesForLocus(seq);
            else                       // use blast hits
                tandems = GetMatchesForLocusFromBlast(locus_id, qlen);            

            SLocusRslt locus_rslt(locus_id);
            list<SFragment>& fragments = locus_rslt.m_fragments; 
            set<tuple<int,int,int>> found_intervals;

            for(auto& tandem : tandems) {
                sort(tandem.begin(), tandem.end());
                list<SFragment> rslt = DeepSearch(tandem, seq, found_intervals); // hits are recalculated!
                for(auto& fragment : rslt) {
                    for(auto& alleleid_seq : locus) {
                        const string& allele_id = get<0>(alleleid_seq);
                        const string& seq = get<1>(alleleid_seq);
                        if(seq == fragment.m_seq) {
                            fragment.m_allele_id = allele_id;
                            break;
                        }
                    }                                
                    /* use longer length for identity */
                    if (fragment.m_to-fragment.m_from+1 < qlen)
                        fragment.m_identity = 100.*fragment.m_matches/qlen;
                    else
                        fragment.m_identity = 100.*fragment.m_matches/(fragment.m_to-fragment.m_from+1);
                    fragments.push_back(fragment);
                }
            }
            if(!fragments.empty()) {
                lock_guard<mutex> guard(m_out_mutex);
                PrintAlignments(locus_rslt);
            }
            if(m_output_loci.is_open()) {
                list<SFragment>& new_loci = get<1>(loc);
                for(auto& fragment : fragments) {
                    if((fragment.m_consistent && fragment.m_has_start && fragment.m_has_stop && fragment.m_identity >= m_minimum_identity) || // complete
                       (fragment.m_not_aligned_fivep > 0 || fragment.m_not_aligned_threep > 0)) {                                             // partial for filtering (will be deleted)
                        new_loci.push_back(fragment);
                    }
                }
            }
        }
    }

    list<SFragment> DeepSearch(vector<SHit>& locus_hits, const string& seq, set<tuple<int,int,int>>& found_intervals) {
        list<SFragment> rslt;

        int qlen = seq.size();
        if(qlen < m_flank)
            return rslt;

        int contig_key = locus_hits.front().m_contig_key;
        int strand = contig_key/abs(contig_key);
        int index = abs(contig_key)-1;
        const string& genome = get<1>(m_genome[index]);

        int left_shift = numeric_limits<int>::max();
        int max_right = 0;
        bool one_alignment = true;
        for(int i = 0; i < (int)locus_hits.size(); ++i) {
            auto& hit = locus_hits[i];
            int sfrom = hit.m_sfrom;
            int sto = hit.m_sto;
            int qfrom = hit.m_qfrom;
            int qto = hit.m_qto;
            left_shift = min(left_shift, sfrom-qfrom);
            max_right = max(max_right, sto+qlen-qto-1);
            if(i > 0 && !(locus_hits[i-1].m_qfrom < qfrom && locus_hits[i-1].m_qto < qto && locus_hits[i-1].m_sfrom < sfrom && locus_hits[i-1].m_sto < sto))
                one_alignment = false ;
        }
        left_shift = max(0, left_shift-qlen/2);
        max_right = min(max_right+qlen/2, (int)genome.size()-1);
        bool left_restricted = (left_shift == 0);
        bool right_restricted = (max_right == (int)genome.size()-1); 

        string contig = genome.substr(left_shift, max_right-left_shift+1);
        if(strand < 0) {
            ReverseComplementSeq(contig.begin(), contig.end());
            swap(left_restricted, right_restricted);
        }
        for(auto& hit : locus_hits) { // convert hits to clipped and possibly reversed contig
            hit.m_sfrom -= left_shift;
            hit.m_sto -= left_shift;
            if(strand < 0) {
                swap(hit.m_sfrom, hit.m_sto);
                hit.m_sfrom = contig.size()-hit.m_sfrom-1;
                hit.m_sto = contig.size()-hit.m_sto-1;
                swap(hit.m_qfrom, hit.m_qto);
                hit.m_qfrom = qlen-hit.m_qfrom-1;
                hit.m_qto = qlen-hit.m_qto-1;
                reverse(hit.m_btop.begin(), hit.m_btop.end());
            }
        }
        if(strand < 0)
            reverse(locus_hits.begin(), locus_hits.end());
        
        list<SHit> aligns;

        if(m_blast_hits.empty() || locus_hits.size() > 1) {
            { // main alignment 
                CCigar cigar;
                int left;
                int right;
                int length;
                int qleft = 0;
                int qright = qlen-1;
                int qlength = qlen;

                if(!m_blast_hits.empty() && locus_hits.size() > 1 && one_alignment) {  // single alignment with multiple blast hits
                    left = locus_hits.front().m_sfrom;
                    right = locus_hits.back().m_sto;
                    length = right-left+1; 
                    qleft = locus_hits.front().m_qfrom;
                    qright = locus_hits.back().m_qto;
                    qlength = qright-qleft+1;
                    vector<TRange> band(qlength, TRange(length-1,0));
                    int delta = 20;
                    //prepare variable band                    
                    for(int i = 0; i < (int)locus_hits.size(); ++i) {
                        auto hit = locus_hits[i];
                        hit.m_qfrom -= qleft;
                        hit.m_qto -= qleft;
                        hit.m_sfrom -= left;
                        hit.m_sto -= left;                        
                        if(i > 0) { //make a bulge
                            auto prev_hit = locus_hits[i-1];
                            prev_hit.m_qfrom -= qleft;
                            prev_hit.m_qto -= qleft;
                            prev_hit.m_sfrom -= left;
                            prev_hit.m_sto -= left;
                            int deltaq = hit.m_qfrom-prev_hit.m_qto;
                            int deltas = hit.m_sfrom-prev_hit.m_sto;
                            if(deltaq >= 0 && deltas >= 0) {
                                for(int qp = max(0,prev_hit.m_qto-delta); qp <= min(qlength-1,hit.m_qfrom+delta); ++qp) {
                                    band[qp].first = max(0, min(band[qp].first,prev_hit.m_sto-delta));
                                    band[qp].second = min(length-1, max(band[qp].second,hit.m_sfrom+delta));
                                }
                            } else if(deltaq < 0 && deltas >= 0) {
                                for(int qp = max(0,hit.m_qfrom-delta); qp <= min(qlength-1,prev_hit.m_qto+delta); ++qp) {
                                    band[qp].first = max(0, min(band[qp].first,prev_hit.m_sto-delta));
                                    band[qp].second = min(length-1, max(band[qp].second,hit.m_sfrom+delta));
                                }
                            } else if(deltaq >= 0 && deltas < 0) {
                                for(int qp = 0; qp <= hit.m_qfrom; ++qp)
                                    band[qp].first = max(0, min(band[qp].first,hit.m_sfrom-delta));
                                for(int qp = prev_hit.m_qto; qp < qlength; ++qp)
                                    band[qp].second = min(length-1, max(band[qp].second,prev_hit.m_sto+delta));
                            } else if(abs(deltaq) >= abs(deltas)) {
                                for(int qp = prev_hit.m_qto; qp <= min(qlength-1, prev_hit.m_qto+delta); ++qp)
                                    band[qp].first = max(0, min(band[qp].first,prev_hit.m_qto-delta));
                                for(int qp = max(0, hit.m_qfrom-delta); qp <= hit.m_qfrom; ++qp)
                                    band[qp].second = min(length-1, max(band[qp].second,hit.m_sfrom+delta));
                            } else {
                                for(int qp = 0; qp <= hit.m_qfrom; ++qp)
                                    band[qp].first = max(0, min(band[qp].first,hit.m_sfrom-delta));
                                for(int qp = prev_hit.m_qto; qp < qlength; ++qp)
                                    band[qp].second = min(length-1, max(band[qp].second,prev_hit.m_sto+delta));                                
                            }
                        }                        
                        
                        int qp = hit.m_qfrom-1;
                        int sp = hit.m_sfrom-1;                        
                        for(auto& elem : hit.m_btop) {
                            for(int k = 0; k < elem.m_len; ++k) {
                                if(elem.m_type == 'M') {
                                    ++qp;
                                    ++sp;
                                } else if(elem.m_type == 'I') {
                                    ++qp;
                                } else {
                                    ++sp;
                                }
                                for(int iq = max(0,qp-delta); iq <= min(qlength-1,qp+delta); ++iq) {
                                    band[iq].first = max(0, min(band[iq].first,sp-delta));
                                    band[iq].second = min(right-left,max(band[iq].second,sp+delta));
                                }
                            }
                        }
                    }
                    for(int i = 2; i < (int)band.size(); ++i)
                        band[i].second = max(band[i].second, band[i-1].second);
                    for(int i = (int)band.size()-2; i >= 0; --i)
                        band[i].first = min(band[i].first, band[i+1].first);
                    cigar = VariBandAlign(seq.c_str()+qleft, qlength, contig.c_str()+left, length, m_gap_open, m_gap_extend, m_delta.matrix, &band[0]); 
                } 
               
                if(cigar.QueryRange().first > cigar.QueryRange().second) { // no alignmnment from blast - do banded
                    //left test             
                    left = 0;
                    right = contig.size()-1;
                    length = right-left+1;                
                    qleft = 0;
                    qright = qlen-1;
                    qlength = qlen;
                    
                    CCigar lcigar = LclAlign(seq.c_str(), m_flank, contig.c_str()+left, length, m_gap_open, m_gap_extend, m_delta.matrix);
                    if(lcigar.Matches(seq.c_str(), contig.c_str()+left) >= m_high_frac*m_flank) { // good left              
                        left = left+lcigar.SubjectRange().first;
                        //right test                
                        length = right-left+1;
                        if(length >= m_flank) {
                            CCigar rcigar = LclAlign(seq.c_str()+qlen-m_flank, m_flank, contig.c_str()+left, length, m_gap_open, m_gap_extend, m_delta.matrix);
                            if(rcigar.Matches(seq.c_str()+qlen-m_flank, contig.c_str()+left) >= m_high_frac*m_flank) { // good right            
                                right = left+rcigar.SubjectRange().second;
                                //main a    lignment                            
                                length = right-left+1;
                                if(length > m_length_tolerance*qlen && qlen > m_length_tolerance*length) { // band alignment (we beleive start/stop)                
                                    vector<TRange> band;
                                    int width = (1-m_length_tolerance)*qlen+m_flank;
                                    for(int i = 0; i < qlen; ++i)
                                        band.push_back(TRange(max(0,i-width), min(right-left,i+width)));
                                    cigar = VariBandAlign(seq.c_str(), qlen, contig.c_str()+left, length, m_gap_open, m_gap_extend, m_delta.matrix, &band[0]);
                                }
                            }
                        }
                    }                                       
                }
                
                if(cigar.QueryRange().first > cigar.QueryRange().second) { // no banded alignment - do full alignmnment         
                    left = 0;
                    right = contig.size()-1;
                    length = right-left+1;
                    qleft = 0;
                    qright = qlen-1;
                    qlength = qlen;

                    cigar = LclAlign(seq.c_str(), qlen, contig.c_str()+left, length, m_gap_open, m_gap_extend, m_delta.matrix);
                } 
          
                TRange qrange = cigar.QueryRange();
                qrange.first += qleft;
                qrange.second += qleft;
                TRange srange = cigar.SubjectRange();
                srange.first += left;
                srange.second += left;

                /*
                cerr << qrange.first+1 << " " << qrange.second+1 << " " << srange.first+1 << " " << srange.second+1 << endl;
                TCharAlign align = cigar.ToAlign(seq.c_str()+qleft, contig.c_str()+left);
                for(unsigned i = 0; i < align.first.size(); ++i)
                    cerr << align.first[i];
                cerr << endl;
                for(unsigned i = 0; i < align.first.size(); ++i) {
                    if(align.first[i] == align.second[i])
                        cerr << "|";
                    else
                        cerr << " ";
                }
                cerr << endl;
                for(unsigned i = 0; i < align.first.size(); ++i)
                    cerr << align.second[i];
                cerr << endl << endl;
                */                                              
 	
                int aligned_len = qrange.second-qrange.first+1;
                int dist_to_edge = srange.first;
                //                  contig end            query not inside                         close to edge    
                bool left_partial = left_restricted && (qrange.first > dist_to_edge) && (dist_to_edge < (1-m_high_frac)*aligned_len);
                dist_to_edge = contig.size()-1-srange.second;
                bool right_partial = right_restricted && (qlen-qrange.second > dist_to_edge) && (dist_to_edge < (1-m_high_frac)*aligned_len);
                bool partial = (left_partial && right_partial) || (left_partial && (aligned_len > m_high_frac*(qlen-qrange.first))) || (right_partial && aligned_len > m_high_frac*qrange.second);
                if(aligned_len > m_high_frac*qlen || (partial && aligned_len > m_high_frac_for_partial*qlen)) {
                    aligns.push_back(SHit(qrange.first, qrange.second, srange.first, srange.second, cigar.Matches(seq.c_str()+qleft, contig.c_str()+left), contig_key));
                }
            }
            if(aligns.empty())
                return rslt;
                
            if(!one_alignment) {
                int left_end = aligns.front().m_sfrom;
                if(left_end >= (left_restricted ? m_high_frac_for_partial*qlen : m_high_frac*qlen)) {
                    int left = 0;
                    int right = left_end-1;
                    int length = right-left+1;
                    CCigar cigar = LclAlign(seq.c_str(), qlen, contig.c_str()+left, length, m_gap_open, m_gap_extend, m_delta.matrix);
                    TRange qrange = cigar.QueryRange();
                    TRange srange = cigar.SubjectRange();
                    int aligned_len = qrange.second-qrange.first+1;
                    int dist_to_edge = srange.first+left;
                    bool partial = left_restricted && (qrange.first > dist_to_edge) && (dist_to_edge < (1-m_high_frac)*aligned_len) && (aligned_len > m_high_frac*(qlen-qrange.first));
                    if(aligned_len > m_high_frac*qlen || (partial && aligned_len > m_high_frac_for_partial*qlen)) {
                        aligns.push_back(SHit(qrange.first, qrange.second, srange.first+left, srange.second+left, cigar.Matches(seq.c_str(), contig.c_str()+left), contig_key));
                    }
                }
                int right_end = aligns.front().m_sto;
                if((int)contig.size()-1-right_end >= (right_restricted ? m_high_frac_for_partial*qlen : m_high_frac*qlen)) {
                    int left = right_end+1;
                    int right = contig.size()-1;
                    int length = right-left+1;
                    CCigar cigar = LclAlign(seq.c_str(), qlen, contig.c_str()+left, length, m_gap_open, m_gap_extend, m_delta.matrix);
                    TRange qrange = cigar.QueryRange();
                    TRange srange = cigar.SubjectRange();
                    int aligned_len = qrange.second-qrange.first+1;
                    int dist_to_edge = contig.size()-1-srange.second-left;
                    bool partial = right_restricted && (qlen-qrange.second > dist_to_edge) && (dist_to_edge < (1-m_high_frac)*aligned_len) && (aligned_len > m_high_frac*qrange.second);
                    if(aligned_len > m_high_frac*qlen || (partial && aligned_len > m_high_frac_for_partial*qlen)) {
                        aligns.push_back(SHit(qrange.first, qrange.second, srange.first+left, srange.second+left, cigar.Matches(seq.c_str(), contig.c_str()+left), contig_key));
                    }
                }
            }
        } else {
            const SHit& hit = locus_hits.front();
            int aligned_len = hit.m_qto-hit.m_qfrom+1;
            int dist_to_edge = hit.m_sfrom;
            bool left_partial = left_restricted && (hit.m_qfrom > dist_to_edge) && (dist_to_edge < (1-m_high_frac)*aligned_len);
            dist_to_edge = contig.size()-1-hit.m_sto;
            bool right_partial = right_restricted && (qlen-hit.m_qto > dist_to_edge) && (dist_to_edge < (1-m_high_frac)*aligned_len);
            bool partial = (left_partial && right_partial) || (left_partial && (aligned_len > m_high_frac*(qlen-hit.m_qfrom))) || (right_partial && aligned_len > m_high_frac*hit.m_qto);
            if(aligned_len > m_high_frac*qlen || (partial && aligned_len > m_high_frac_for_partial*qlen)) {
                aligns.push_back(hit);
            }
        }

        for(auto& al : aligns) {
            int not_aligned_left = al.m_qfrom;
            int from = al.m_sfrom-not_aligned_left;
            int left_extend = not_aligned_left;
            if(from < 0 && left_restricted) {
                not_aligned_left = -from;
                left_extend -= not_aligned_left;
                from = 0;
            } else {
                not_aligned_left = 0;
            }

            int not_aligned_right = qlen-al.m_qto-1;
            int to = al.m_sto+not_aligned_right;
            int right_extend = not_aligned_right;
            if(to > (int)contig.size()-1 && right_restricted) {                
                not_aligned_right = to-contig.size()+1;
                right_extend -= not_aligned_right;
                to = contig.size()-1;
            } else {
                not_aligned_right = 0;
            }

            if(from < 0 || to > (int)contig.size()-1)
                continue;

            int matches = al.m_score;
            for(int i = 0; i < left_extend; ++i) {
                if(seq[i+not_aligned_left] == contig[i+from])
                    ++matches;
            }
            for(int i = 0; i < right_extend; ++i) {
                if(seq[qlen-1-i-not_aligned_right] == contig[to-i])
                    ++matches;
            }
            string fseq = contig.substr(from, to-from+1);

            if(strand < 0) {
                swap(from, to);
                from = contig.size()-1-from;
                to = contig.size()-1-to;
            }
            from += left_shift;
            to += left_shift;

            tuple<int,int,int> interval(from, to, locus_hits.front().m_contig_key);
            if(!found_intervals.insert(interval).second)
                continue;        

            SFragment fragment;
            fragment.m_contig = index;
            fragment.m_strand = strand;
            fragment.m_from = from;            
            fragment.m_to = to;
            fragment.m_not_aligned_fivep = not_aligned_left;
            fragment.m_not_aligned_threep = not_aligned_right;
            fragment.m_has_start = (not_aligned_left == 0 && isStart(fseq.begin()));
            fragment.m_has_stop = (not_aligned_right == 0 && isStop(fseq.end()-3));
            bool consistent = (to-from+1)%3 == 0;
            for(auto is = fseq.begin(); is <= fseq.end()-6 && consistent; is += 3)
                consistent = !isStop(is);
            fragment.m_seq = move(fseq);
            fragment.m_consistent = consistent;
            fragment.m_matches = matches;  
            rslt.push_back(fragment);
        }

        return rslt;
    }

    THitTandems GetMatchesForLocusFromBlast(const string& locus_id, int qlen) {
        THitTandems tandems;
        for(auto& blst : m_blast_hits[locus_id]) {
            int contig_key = blst.first;
            int contig_len = get<1>(m_genome[abs(contig_key)-1]).size();
            TLinkedHits& candidates = blst.second;
            sort(candidates.begin(), candidates.end());
            for(int i = 1; i < (int)candidates.size(); ++i) {
                int best_score = candidates[i].m_score;
                SLinkedHit* left = nullptr;
                for(int j = i-1; j >= 0; --j) {
                    if(candidates[i].m_sfrom-candidates[j].m_sto > qlen/2) // too far
                        continue;
                    if(candidates[i].m_sfrom > candidates[j].m_sfrom && candidates[i].m_sto > candidates[j].m_sto &&
                       candidates[i].m_qfrom > candidates[j].m_qfrom && candidates[i].m_qto > candidates[j].m_qto) {  // extension
                        int overlap = max(0, candidates[j].m_sto-candidates[i].m_sfrom+1);
                        int new_score = candidates[j].m_score+candidates[i].m_score-overlap;  // overlaps are overpenalized 
                        if(new_score > best_score) {
                            best_score = new_score;
                            left = &candidates[j];
                        }
                    }
                }
                if(left != nullptr) {
                    candidates[i].m_score = best_score;
                    candidates[i].m_left = left;
                }
            }

            /*
            map<pair<int,int>,list<SHit>> singles;
            for(int i = (int)candidates.size()-1; i >= 0; --i) {
                auto& candidate = candidates[i];
                if(candidate.m_score < 0) // already included
                    continue;
                SLinkedHit* left = &candidate;
                while(left->m_left != nullptr) {
                    left = left->m_left;
                    left->m_score = -1;  // prevent starting from already included
                }
                bool possibly_partial = (left->m_sfrom < left->m_qfrom || candidate.m_sto > contig_len-(qlen-candidate.m_qto));
                int len = candidate.m_qto-left->m_qfrom+1;
                if(candidate.m_score > 0.5*len && (possibly_partial || len > 0.75*qlen)) {
                    int a = left->m_sfrom;
                    int b = candidate.m_sto;
                    list<SHit> hits;
                    for(SLinkedHit* l = &candidate; l != nullptr; l = l->m_left)
                        hits.push_front(*l);
                    if(singles.empty()) {
                        singles[make_pair(a,b)] = hits;
                    } else {
                        for(auto i = singles.begin(); i != singles.end(); ++i) {
                            bool overlap = i->first.first <= b && i->first.second >= a;
                            if(overlap) {
                                if(candidate.m_score > i->second.back().m_score) {
                                    singles.erase(i);
                                    singles[make_pair(a,b)] = hits;
                                    break;
                                }
                            }
                            singles[make_pair(a,b)] = hits;
                        }
                    }
                }
            }
            for(auto& single : singles) {
                const SHit& last = single.second.back(); 
                if(tandems.empty() || tandems.back().back().m_contig_key != contig_key 
                   || last.m_sto+(qlen-last.m_qto)+qlen/2 < tandems.front().front().m_sfrom-tandems.front().front().m_qfrom
                   || (single.second.size() == 1 && last.m_qfrom==0 && last.m_qto == qlen-1)) {
                        tandems.push_back(vector<SHit>());
                }
                tandems.back().insert(tandems.back().end(), single.second.begin(), single.second.end());
            }
            */

            
            for(int i = (int)candidates.size()-1; i >= 0; --i) {
                auto& candidate = candidates[i];
                if(candidate.m_score < 0) // already included
                    continue;
                SLinkedHit* left = &candidate;
                while(left->m_left != nullptr) {
                    left = left->m_left;
                    left->m_score = -1;  // prevent starting from already included
                }
                bool possibly_partial = (left->m_sfrom < left->m_qfrom || candidate.m_sto > contig_len-(qlen-candidate.m_qto));
                int len = candidate.m_qto-left->m_qfrom+1;
                if(candidate.m_score > 0.5*len && (possibly_partial || len > 0.75*qlen)) {
                    if(tandems.empty() || tandems.front().back().m_contig_key != contig_key 
                       || candidate.m_sto+(qlen-candidate.m_qto)+qlen/2 < tandems.front().front().m_sfrom-tandems.front().front().m_qfrom
                       || (&candidate == left && candidate.m_qfrom==0 && candidate.m_qto == qlen-1)) {
                        tandems.push_front(vector<SHit>());
                    }
                    list<SHit> hits;
                    for(SLinkedHit* l = &candidate; l != nullptr; l = l->m_left)
                        hits.push_front(*l);
                    tandems.front().insert(tandems.front().begin(), hits.begin(), hits.end());
                }
            }
            
        }

        return tandems;
    }

    

    THitTandems GetMatchesForLocus(const string& seq) {
        map<int, TLinkedHits> kmers_for_locus;  // number of kmers for [contig_key]<hits>; contig_key represents strand and index in m_geneomes; hits are exact kmer matches
        int kmer_len = m_genome_kmers.KmerLen();
        int qlen = seq.size();
        if(qlen >= kmer_len) {
            CReadHolder rh(false);
            rh.PushBack(seq);
            int query_pos = qlen-kmer_len;
            for(CReadHolder::kmer_iterator ik = rh.kbegin(kmer_len) ; ik != rh.kend(); ++ik, --query_pos) { // iteration from last kmer to first  
                TKmer kmer = *ik;
                auto rslt = FindKmerInContigs(kmer);
                if(rslt.first != nullptr) {
                    // connect overlaping matches
                    int plus = rslt.second;
                    for(auto& hit : *rslt.first) {
                        int strand = plus*get<1>(hit);
                        int subject_pos = get<0>(hit);            // left kmer position regardless of strand
                        int contig_key = strand*(get<2>(hit)+1);  // represents strand and index in geneomes
                        TLinkedHits& matches = kmers_for_locus[contig_key];
                        int qp = (strand > 0) ? query_pos : qlen-kmer_len-query_pos;  // position for reversed query if strand < 0                        
                        bool extended = false;
                        if(!matches.empty()) {
                            SLinkedHit& last_match = matches.back();
                            int sdist = subject_pos-last_match.m_sfrom;
                            int qdist = qp-last_match.m_qfrom;
                            if(sdist == qdist) {
                                int len = last_match.m_sto-last_match.m_sfrom+1;
                                if(sdist > 0 && sdist <= len && sdist+kmer_len > len) {
                                    int extra_len = sdist+kmer_len-len;
                                    last_match.m_sto += extra_len;
                                    last_match.m_qto += extra_len;
                                    last_match.m_score += extra_len;
                                    extended = true;
                                } else if(sdist < 0 && -sdist <= kmer_len) {
                                    int extra_len = -sdist;
                                    last_match.m_sfrom -= extra_len;
                                    last_match.m_qfrom -= extra_len;
                                    last_match.m_score += extra_len;
                                    extended = true;
                                }
                            }
                        }
                        if(!extended) 
                            matches.push_back(SLinkedHit(qp, qp+kmer_len-1, subject_pos, subject_pos+kmer_len-1, kmer_len, contig_key));
                    }
                }
            }
        }

        THitTandems tandems;
        string revseq;  // uses lazy init
        for(auto& contig_matches : kmers_for_locus) {
            int contig_key = contig_matches.first;
            int index = abs(contig_key)-1;
            string& contig = get<1>(m_genome[index]);
            TLinkedHits& matches = contig_matches.second;
            sort(matches.begin(), matches.end());
            for(int i = 1; i < (int)matches.size(); ++i) {
                int best_score = matches[i].m_score;
                SLinkedHit* left = nullptr;
                for(int j = i-1; j >= 0; --j) {
                    if(matches[j].m_sfrom >= matches[i].m_sfrom || matches[j].m_sto >= matches[i].m_sto)     // not extension
                        continue;
                    int shift = matches[i].m_sfrom-matches[j].m_sfrom;
                    if(shift >  matches[i].m_qfrom) // too far
                        break;                                                           
                                                           
                    if(shift == matches[i].m_qfrom-matches[j].m_qfrom) {   // diagonal
                        int overlap = max(0, matches[j].m_sto-matches[i].m_sfrom+1);
                        int new_score = matches[j].m_score+matches[i].m_score-overlap;
                        if(new_score > best_score) {
                            best_score = new_score;
                            left = &matches[j];
                        }
                    }
                }
                if(left != nullptr) {
                    matches[i].m_score = best_score;
                    matches[i].m_left = left;
                }
            }

            TLinkedHits diag_candidates;
            for(int i = (int)matches.size()-1; i >= 0; --i) {
                if(matches[i].m_score >= m_min_bases_to_consider) {
                    SLinkedHit* left = &matches[i];
                    while(left->m_left != nullptr) {
                        left = left->m_left;
                        left->m_score = -1;  // prevent starting from already included
                    }
                    int qfrom = left->m_qfrom;
                    int sfrom = left->m_sfrom;
                    int len = matches[i].m_sto-sfrom+1;

                    int identicals = 0;
                    const string* sp = &seq;
                    if(contig_key < 0) {
                        if(revseq.empty()) {   // lazy init
                            revseq = seq;
                            ReverseComplementSeq(revseq.begin(), revseq.end());
                        }
                        sp = &revseq;                        
                    }
                    for(int p = 0; p < len; ++p) {
                        if((*sp)[qfrom+p] == contig[sfrom+p])
                            ++identicals;
                    }
                    if(identicals > 0.5*len)
                        diag_candidates.push_back(SLinkedHit(qfrom, qfrom+len-1, sfrom, sfrom+len-1, identicals, contig_key));
                }
            }

            THits candidates;
            sort(diag_candidates.begin(), diag_candidates.end());
            for(int i = 1; i < (int)diag_candidates.size(); ++i) {
                int best_score = diag_candidates[i].m_score;
                SLinkedHit* left = nullptr;
                for(int j = i-1; j >= 0; --j) {
                    if(diag_candidates[j].m_sfrom >= diag_candidates[i].m_sfrom || diag_candidates[j].m_sto >= diag_candidates[i].m_sto)     // not extension
                        continue;
                    int diag = diag_candidates[i].m_sfrom-diag_candidates[j].m_sfrom-diag_candidates[i].m_qfrom+diag_candidates[j].m_qfrom;
                    if(abs(diag) > m_max_off_diag)
                        continue;
                    int overlap = max(0, diag_candidates[j].m_sto-diag_candidates[i].m_sfrom+1);
                    int new_score = diag_candidates[j].m_score+diag_candidates[i].m_score-overlap;  // overpenilize for overlaps > kmer_len
                    if(new_score > best_score) {
                        best_score = new_score;
                        left = &diag_candidates[j];
                    }
                }
                if(left != nullptr) {
                    diag_candidates[i].m_score = best_score;
                    diag_candidates[i].m_left = left;
                }
            }
 
            for(int i = (int)diag_candidates.size()-1; i >= 0; --i) {
                if(diag_candidates[i].m_score >= m_min_deep_search_fraction*qlen) {
                    SLinkedHit* left = &diag_candidates[i];
                    while(left->m_left != nullptr) {
                        left = left->m_left;
                        left->m_score = -1;  // prevent starting from already included
                    }
                    candidates.push_back(SLinkedHit(left->m_qfrom, diag_candidates[i].m_qto, left->m_sfrom, diag_candidates[i].m_sto, diag_candidates[i].m_score, contig_key));
                }
            }

            sort(candidates.begin(), candidates.end());
            for(auto& candidate : candidates) {
                if(tandems.empty() || tandems.back().back().m_contig_key != contig_key || tandems.back().back().m_sto+(qlen-tandems.back().back().m_qto)+qlen/2 < candidate.m_sfrom-candidate.m_qfrom)
                    tandems.push_back(vector<SHit>(1, candidate));
                else if(tandems.back().back().m_sfrom < candidate.m_sfrom && tandems.back().back().m_sto < candidate.m_sto)
                    tandems.back().push_back(candidate);
            }
        }

        return tandems;
    }

    pair<TKmerToGenome::MappedType*, int> FindKmerInContigs(const TKmer& kmer) {
        TKmer rkmer = revcomp(kmer, m_genome_kmers.KmerLen());
        const TKmer* kmerp = &kmer;
        int plus = 1;
        if(rkmer < kmer) {
            kmerp = &rkmer;
            plus = -1;
        }
        return make_pair(m_genome_kmers.Find(*kmerp), plus);
    }

    bool isStart(string::const_iterator i) const {
        return find_if(m_start_codons.begin(), m_start_codons.end(), [i](const string& codon){return *i == codon[0] && *(i+1) == codon[1] && *(i+2) == codon[2];}) != m_start_codons.end();
    }

    bool isStop(string::const_iterator i) const {
        return find_if(m_stop_codons.begin(), m_stop_codons.end(), [i](const string& codon){return *i == codon[0] && *(i+1) == codon[1] && *(i+2) == codon[2];}) != m_stop_codons.end();
    }

    void CheckAndCleanAllelesJob(vector<pair<TLocus*, SAtomic<uint8_t>>>& loci) {
        for(auto& lpair : loci) {        
            if(!lpair.second.Set(1))
                continue;
            TLocus& locus = *lpair.first;
            for(auto iloop = locus.begin(); iloop != locus.end(); ) {
                auto it = iloop++;
                string& seq = get<1>(*it);
                if((int)seq.size() < m_flank) {  // too short for anything 
                    locus.erase(it);
                    continue;
                }

                if(seq.find_first_not_of("ACGT") != string::npos) {  // Ns are present
                    locus.erase(it);
                    continue;
                }
            }
        }
    }
    
    TAlles  m_alleles;
    TGenome m_genome;
    map<string, int> m_genome_acc_to_index;
    map<int, set<int>> m_bad_bases;
    TKmerToGenome m_genome_kmers;
    map<string, map<int,TLinkedHits>> m_blast_hits;

    double m_min_deep_search_fraction = 0.1;  //0.2;
    int m_min_bases_to_consider = 15;
    int m_max_off_diag = 15;

    ofstream& m_output_mappings;
    ofstream& m_output_loci;

    int m_ncores;
    //TODO: Other genetic code(s)?
    const array<string,7> m_start_codons = { {"ATG", "TTG", "CTG", "ATT", "ATC", "ATA", "GTG"} };
    const array<string,3> m_stop_codons = { {"TAA", "TAG", "TGA"} };

    int m_flank = 42;                   // pattern to align 5'/3' ends
    double m_length_tolerance = 0.9;    // length diff to apply banded alignment
    double m_high_frac = 0.75;
    double m_high_frac_for_partial = 0.1;
    double m_minimum_identity = 85.;    // identity at which we will call an allele as found
    int m_consistent_bases = 115;       // overlapping gene not allowed at this length or more
    double m_percent_contained = 0.75;

    int m_match = 1;
    int m_mismatch = 1;
    int m_gap_open = 8;
    int m_gap_extend = 2;
    SMatrix m_delta;

    mutex m_out_mutex;
};


int main(int argc, const char* argv[])
{
    options_description arguments("Program arguments");
    arguments.add_options()
        ("help,h", "Produce help message")
        ("genome", value<string>()->required(), "Assembled genome (required)")
        ("bad_bases", value<string>(), "Positions of low quality genome bases (optional)")
        ("alleles", value<string>()->required(), "Alleles (required)")
        ("output_mappings", value<string>(), "Output allele mappings (optional, default cout)")
        ("output_loci", value<string>(), "Output new loci (optional)")
        ("blast_hits", value<string>(), "Blast hits (optional)")
        ("kmer", value<int>()->default_value(15), "Kmer length for allele search")
        ("min_kmer_bases", value<int>()->default_value(15), "Minimal bases in exact diagonal kmer matches to consider for an alignment")
        ("min_fraction_of_matches", value<double>()->default_value(0.1, "0.1"), "Minimal fraction of the query found as matches in diagonal hits to consider for an alignment")
        ("match", value<int>()->default_value(1), "Bonus for match")
        ("mismatch", value<int>()->default_value(1), "Penalty for mismatch")
        ("gap_open", value<int>()->default_value(8), "Penalty for gap opening")
        ("gap_extend", value<int>()->default_value(2), "Penalty for gap extension")
        ("cores", value<int>()->default_value(0), "Number of cores to use (default all) [integer]");



    int kmer_len;
    int min_kmer_bases;
    double min_fraction_of_matches;
    boost::iostreams::filtering_istream genome_file;
    ifstream bad_bases_file;
    boost::iostreams::filtering_istream alleles_file;
    boost::iostreams::filtering_istream blast_file;
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
    bool blast_hits_present = false;
    ofstream output_mappings;
    ofstream output_loci;


    try {
        variables_map argmap;                                // boost arguments
        store(parse_command_line(argc, argv, arguments), argmap);

        if(argmap.count("help")) {
#ifdef SVN_REV
            cerr << "SVN revision:" << SVN_REV << endl << endl;
#endif
            cerr << arguments << "\n";
            return 1;
        }

        // must be after "help" if thre are 'required' options
        notify(argmap);    

        kmer_len = argmap["kmer"].as<int>();
        min_kmer_bases = argmap["min_kmer_bases"].as<int>();
        min_fraction_of_matches = argmap["min_fraction_of_matches"].as<double>();

        {
            string file = argmap["genome"].as<string>();
            boost::iostreams::file_source f{file};
            if(!f.is_open()) {
                cerr << "Can't open file " << file << endl;
                return 1; 
            }            
            if(file.size() > 3 &&  file.substr(file.size()-3) == ".gz")
                genome_file.push(boost::iostreams::gzip_decompressor());
            genome_file.push(f);
        }

        if(argmap.count("bad_bases")) {
            string file = argmap["bad_bases"].as<string>();
            bad_bases_file.open(file);
            if(!bad_bases_file.is_open()) {
                cerr << "Can't open file " << file << endl;
                return 1; 
            }    
        }

        {
            string file = argmap["alleles"].as<string>();
            boost::iostreams::file_source f{file};
            if(!f.is_open()) {
                cerr << "Can't open file " << file << endl;
                return 1;
            }
            if(file.size() > 3 &&  file.substr(file.size()-3) == ".gz")
                alleles_file.push(boost::iostreams::gzip_decompressor());
            alleles_file.push(f);
        }

        if(argmap.count("blast_hits")) {
            string file = argmap["blast_hits"].as<string>();
            boost::iostreams::file_source f{file};
            if(!f.is_open()) {
                cerr << "Can't open file " << file << endl;
                return 1;
            }
            if(file.size() > 3 &&  file.substr(file.size()-3) == ".gz")
                blast_file.push(boost::iostreams::gzip_decompressor());
            blast_file.push(f);
            blast_hits_present = true;
        }

        if(argmap.count("output_mappings")) {
            output_mappings.open(argmap["output_mappings"].as<string>());
            if(!output_mappings.is_open()) {
                cerr << "Can't open file " << argmap["output_mappings"].as<string>() << endl;
                exit(1);
            }
        }

        if(argmap.count("output_loci")) {
            output_loci.open(argmap["output_loci"].as<string>());
            if(!output_loci.is_open()) {
                cerr << "Can't open file " << argmap["output_loci"].as<string>() << endl;
                exit(1);
            }
        }

        int ncores = thread::hardware_concurrency();
        if(argmap["cores"].as<int>()) {
            int nc = argmap["cores"].as<int>();
            if(nc < 0) {
                cerr << "Value of --cores must be >= 0" << endl;
                exit(1);
            } else if(nc > ncores) {
                cerr << "WARNING: number of cores was reduced to the hardware limit of " << ncores << " cores" << endl;
            } else if(nc > 0) {
                ncores = nc;
            }
        }

        match = argmap["match"].as<int>();
        mismatch = argmap["mismatch"].as<int>();
        gap_open = argmap["gap_open"].as<int>();
        gap_extend = argmap["gap_extend"].as<int>();

        CStopWatch timer;
        timer.Restart();
        CwgMLST wg_mlst(alleles_file, output_mappings, output_loci, ncores); 
        //        cerr << "Alleles input in " << timer.Elapsed();
        timer.Restart();
        wg_mlst.CheckAndCleanAlleles();
        //        cerr << "Alleles cleaning in " << timer.Elapsed(); 

        timer.Restart();
        wg_mlst.ReadGenome(genome_file);
        if(bad_bases_file.is_open())
            wg_mlst.ReadBadBases(bad_bases_file);

        if(blast_hits_present) {
            wg_mlst.ReadBlastHits(blast_file);
            //            cerr << "Genome and blast hits input in " << timer.Elapsed();
        } else {
            wg_mlst.PrepareKmerMap(kmer_len);
            //            cerr << "Genome input and kmermap in " << timer.Elapsed(); 
        }
        
        timer.Restart();
        wg_mlst.AnalyzeAlleles(min_kmer_bases, min_fraction_of_matches, match, mismatch, gap_open, gap_extend);  
        //        cerr << "Alleles found in " << timer.Elapsed(); 
    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        return 1;
    }
}

