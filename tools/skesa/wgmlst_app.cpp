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
 * Author:  Alexander Souvorov
 *          
 * Edits:   Alex Kotliarov
 *
 * File Description:
 */
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <thread>

#include <boost/program_options.hpp>

#include "wgmlst.hpp"

using namespace std;

namespace po = boost::program_options;

int main(int argc, const char* argv[])
{
    po::options_description arguments("Program arguments");
    arguments.add_options()
        ("help,h", "Produce help message")
        ("version,v", "Print version")
        ("genome", po::value<string>()->required(), "Assembled genome (required)")
        ("bad_bases", po::value<string>(), "Positions of low quality genome bases (optional)")
        ("alleles", po::value<string>()->required(), "Alleles (required)")
        ("output_mappings", po::value<string>(), "Output allele mappings (optional, default cout)")
        ("output_loci", po::value<string>(), "Output new loci (optional)")
        ("blast_hits", po::value<string>(), "Blast hits (optional)")
        ("kmer", po::value<int>()->default_value(15), "Kmer length for allele search")
        ("min_kmer_bases", po::value<int>()->default_value(15), "Minimal bases in exact diagonal kmer matches to consider for an alignment")
        ("min_fraction_of_matches", po::value<double>()->default_value(0.1, "0.1"), "Minimal fraction of the query found as matches in diagonal hits to consider for an alignment")
        ("match", po::value<int>()->default_value(1), "Bonus for match")
        ("mismatch", po::value<int>()->default_value(1), "Penalty for mismatch")
        ("gap_open", po::value<int>()->default_value(8), "Penalty for gap opening")
        ("gap_extend", po::value<int>()->default_value(2), "Penalty for gap extension")
        ("cores", po::value<int>()->default_value(0), "Number of cores to use (default all) [integer]");



    int kmer_len;
    int min_kmer_bases;
    double min_fraction_of_matches;
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;

    try {
        po::variables_map argmap;                                // boost arguments
        po::store(po::parse_command_line(argc, argv, arguments), argmap);

        if(argmap.count("help")) {
#ifdef SVN_REV
            cerr << "SVN revision:" << SVN_REV << endl << endl;
#endif
            cerr << arguments << "\n";
            return 1;
        }

        if(argmap.count("version")) {
            cerr << "WGMLST v.1.1";
#ifdef SVN_REV
            cerr << "-SVN_" << SVN_REV;
#endif
            cerr << endl;
            return 0;
        }

        // must be after "help" if thre are 'required' options
        po::notify(argmap);    

        string genome_file_path = argmap["genome"].as<string>();
        string alleles_file_path = argmap["alleles"].as<string>();
        string bad_bases_path { (argmap.count("bad_bases") ? argmap["bad_bases"].as<string>() : "") };
        string blast_hits_path { (argmap.count("blast_hits") ? argmap["blast_hits"].as<string>() : "") };
        string output_mappings_path { (argmap.count("output_mappings") ? argmap["output_mappings"].as<string>() : "") };
        string output_loci_path { (argmap.count("output_loci") ? argmap["output_loci"].as<string>() : "") };

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

        kmer_len = argmap["kmer"].as<int>();
        min_kmer_bases = argmap["min_kmer_bases"].as<int>();
        min_fraction_of_matches = argmap["min_fraction_of_matches"].as<double>();

        match = argmap["match"].as<int>();
        mismatch = argmap["mismatch"].as<int>();
        gap_open = argmap["gap_open"].as<int>();
        gap_extend = argmap["gap_extend"].as<int>();

        libwgmlst::wgmlst(genome_file_path,
                          alleles_file_path,
                          bad_bases_path,
                          blast_hits_path,
                          output_mappings_path,
                          output_loci_path,
                          kmer_len,
                          min_kmer_bases,
                          min_fraction_of_matches,
                          match,
                          mismatch,
                          gap_open,
                          gap_extend,
                          ncores);

    } catch (exception &e) {
        cerr << endl << e.what() << endl;
        return 1;
    }
}
