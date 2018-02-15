#ifndef __WGMLST_API_HPP
#define __WGMLST_API_HPP
/* 
 * ===========================================================================
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
 * Author:  Alex Kotliarov
 *
 * File Description:
 */
#include <fstream>
#include <functional>
#include <vector>
#include <string>
#include <tuple>

#include <boost/iostreams/filtering_stream.hpp>


#include "contigs.hpp"
#include "blast_hits.hpp"

namespace libwgmlst
{

class CAlleleAlignment
{
    std::string  m_Locus;
    std::pair<std::string, std::string>  m_Alleles;
    std::string  m_Consistency;
    int32_t m_Start;
    int32_t m_Stop;
    int32_t m_From;
    int32_t m_To;
    char    m_Strand;
    int32_t m_Matches;
public:
    CAlleleAlignment() = delete;
    CAlleleAlignment(std::string locus,
                     std::string allele_first,
                     std::string allele_second,
                     std::string consistency,
                     int32_t start,
                     int32_t stop,
                     int32_t from,
                     int32_t to,
                     int32_t matches,
                     char strand)
    : m_Locus(move(locus)),
      m_Alleles(allele_first, allele_second),
      m_Consistency(move(consistency)),
      m_Start(start), 
      m_Stop(stop),
      m_From(from),
      m_To(to), 
      m_Strand(strand),
      m_Matches(matches)
    {
    }

    std::string const& Locus() const { return m_Locus; }
    std::pair<std::string, std::string> const& Alleles() const { return m_Alleles; }
    int Start() const { return m_Start; }
    int Stop() const { return m_Stop; }
    int From() const { return m_From; }
    int To() const { return m_To; }
    int Matches() const { return m_Matches; }
    char Strand() const { return m_Strand; }
};


// Alignment(locus,allele_first, allele_second, consistency, start, stop, from, to, matches, strand, seq)
using TEmitAlignmentCallback = std::function<void(std::string const&, std::string const&, std::string const&, std::string const&, std::string const&, std::string const&, int32_t, int32_t, int32_t, char, std::string const&)>;
// AlleleCallback(locus, allele, contig-id, from, to, identity_pct, length, ref-allele-name); 
using TEmitAlleleCallback = std::function<void(std::string const&, std::string const&, std::string const&, int32_t, int32_t, double, int32_t, std::string const&)>;


// wgmlst_align_alleles
// //  Arguments:
// //   @locus: locus name
// //   @queries:   collection of (allele-name, allele-contig)
// //   @subjects:  collection of (allele-name, allele-contig)
// //   @kmer_length:                  default_value: 15,  Kmer length for allele search.
// //   @min_kmer_bases:              default_value: 15,  Minimal bases in exact diagonal kmer matches to consider for an alignment.
// //   @min_fraction_of_matches:     default_value: 0.1, Minimal fraction of the query found as matches in diagonal hits to consider for an alignment.
// //   @match:                       default_value: 1,   Bonus for match.
// //   @mismatch:                    default_value: 1,   Penalty for mismatch.
// //   @gap_open:                    default_value: 8,   Penalty for gap opening.
// //   @gap_extend:                  default_value: 2,   Penalty for gap extension.
//  
std::vector<CAlleleAlignment> wgmlst_allele_alignments(std::string locus,
                                                  std::vector<std::tuple<std::string, std::string> > queries,  // collection of (allele-id, contig) std::pairs
                                                  std::vector<std::tuple<std::string, std::string> > subjects, // collection of (allele-id, contig) std::pairs
                                                  int kmer_len,
                                                  int min_kmer_bases,
                                                  double min_fraction_of_matches,
                                                  int match,
                                                  int mismatch,
                                                  int gap_open,
                                                  int gap_extend);

void wgmlst(boost::iostreams::filtering_istream& genome_file, 
            boost::iostreams::filtering_istream& alleles_file,
            std::ifstream& bad_bases, 
            bool blast_hits_present,
            boost::iostreams::filtering_istream& blast_file,
            int kmer_len,
            int min_kmer_bases,
            double min_fraction_of_matches,
            int match,
            int mismatch,
            int gap_open,
            int gap_extend,
            int ncores,
            TEmitAlignmentCallback alignment_cb,
            TEmitAlleleCallback allele_cb);

void wgmlst(std::string const& genome_path,
            std::string const& alleles_path,
            std::string const& bad_bases_path,
            std::string const& blast_hits_path,
            std::string const& output_mappings_path,
            std::string const& output_loci_path,
            int kmer_len,
            int min_kmer_bases,
            double min_fraction_of_matches,
            int match,
            int mismatch,
            int gap_open,
            int gap_extend,
            int ncores);

using allele_seq_t = std::pair<std::string, std::string>;

std::vector<allele_seq_t> wgmlst(CContigsSource& genome,
                                 std::string const& alleles_path,
                                 std::string const& bad_bases_path, 
                                 CBlastHitsSource& blast_hits,
                                 int min_kmer_bases,
                                 double min_fraction_of_matches,
                                 int match,
                                 int mismatch,
                                 int gap_open,
                                 int gap_extend,
                                 int ncores);


//
} // namespace libwgmlst

#endif
