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

#include <ncbi_pch.hpp>
#include <corelib/ncbienv.hpp>
#include <corelib/ncbiargs.hpp>

#include "kmercounts.hpp"

class CDeBruijnApplication : public CNcbiApplication
{
protected:
    virtual void Init();
    virtual int  Run(void);
};

void CDeBruijnApplication::Init()
{
    SetDiagPostLevel(eDiag_Info);

    auto_ptr<CArgDescriptions> arg_desc(new CArgDescriptions);
    arg_desc->SetUsageContext(GetArguments().GetProgramBasename(),"bact_kmer_counter");

    arg_desc->AddKey("fasta", "fasta", "List of reads fasta files", CArgDescriptions::eInputFile);
    arg_desc->AddKey("onefasta", "onefasta", "Fasta file of reads", CArgDescriptions::eString);
    arg_desc->AddKey("accessions", "accessions", "List of SRA accessions", CArgDescriptions::eInputFile);
    arg_desc->AddKey("oneaccession", "oneaccession", "SRA accession", CArgDescriptions::eString);
    arg_desc->SetDependency("fasta", CArgDescriptions::eExcludes, "onefasta");
    arg_desc->SetDependency("fasta", CArgDescriptions::eExcludes, "accessions");
    arg_desc->SetDependency("fasta", CArgDescriptions::eExcludes, "oneaccession");
    arg_desc->SetDependency("onefasta", CArgDescriptions::eExcludes, "accessions");
    arg_desc->SetDependency("onefasta", CArgDescriptions::eExcludes, "oneaccession");
    arg_desc->SetDependency("accessions", CArgDescriptions::eExcludes, "oneaccession");

    arg_desc->AddKey("kmer", "kmer", "Kmer length (minimal kmer length if steps > 1)", CArgDescriptions::eInteger);
    arg_desc->AddDefaultKey("steps", "steps", "Desired number of kmer steps", CArgDescriptions::eInteger, "1");
    arg_desc->AddKey("min_count", "min_count", "Minimal count for retained kmers", CArgDescriptions::eInteger);
    arg_desc->AddKey("out", "out", "Output file", CArgDescriptions::eOutputFile);
    arg_desc->AddOptionalKey("hist", "hist", "File for histogram", CArgDescriptions::eOutputFile);
    arg_desc->AddDefaultKey("cores", "cores", "Number of cores to use (default all)", CArgDescriptions::eInteger, "0");
    arg_desc->AddDefaultKey("memory", "memory", "Memory to use (GB)", CArgDescriptions::eInteger, "100");    

    SetupArgDescriptions(arg_desc.release());
}

void GetBranchesAndSaveKmers(TKmerCount& sorted_kmers, int ncores, const CArgs& args) {
    GetBranches(ncores, sorted_kmers);

    map<int,size_t> bins;
    for(size_t index = 0; index < sorted_kmers.Size(); ++index) {
        ++bins[sorted_kmers.GetCount(index)];                  // count clipped to integer automatically
    }

    CStopWatch timer;
    timer.Restart();
    ostream& out = args["out"].AsOutputFile();
    sorted_kmers.Save(out);

    int bin_num = bins.size();
    out.write(reinterpret_cast<const char*>(&bin_num), sizeof(bin_num));         
    for(auto& bin : bins)
        out.write(reinterpret_cast<const char*>(&bin), sizeof(bin)); 

    bool is_stranded = true;            
    out.write(reinterpret_cast<const char*>(&is_stranded), sizeof is_stranded);

    if(args["hist"]) {
        ostream& hist = args["hist"].AsOutputFile();
        for(auto& bin : bins)
            hist << sorted_kmers.KmerLen() << '\t' << bin.first << '\t' << bin.second << endl;
    }

    cerr << "Kmers saving in " << timer.Elapsed() << endl;
}
    
int CDeBruijnApplication::Run(void)
{
    const CArgs& args(GetArgs());

    int ncores = thread::hardware_concurrency();
    if(args["cores"].AsInteger())
        ncores = min(ncores, args["cores"].AsInteger());
    int steps = args["steps"].AsInteger();
    int min_count = args["min_count"].AsInteger();

    map<int,CDBGraph*> graphs;
    list<array<CReadHolder,2>> raw_reads;
    GetReads(ncores, args, raw_reads);

    list<TKmerCount> uniq_kmers;  // storage        
    
    int min_kmer = args["kmer"].AsInteger();
    TKmerCount* sorted_kmers = GetSortedKmers(min_kmer, min_count, raw_reads, uniq_kmers, ncores, args, graphs);
    GetBranchesAndSaveKmers(*sorted_kmers, ncores, args);        

    /*
    if(steps > 1) {
        int read_len = 0;
        for(auto& reads : raw_reads) 
            read_len = max(read_len, (int)reads.MaxLength());
        cerr << "Maximal read length: " << read_len << endl;
        size_t genome_size = GenomeSize(*sorted_kmers);
        cerr << "Genome size: " << genome_size << endl;

        int max_kmer = min(TKmer::MaxKmer(), read_len);
        max_kmer -= 1-max_kmer%2;
        while(max_kmer > min_kmer) {
            sorted_kmers = GetSortedKmers(max_kmer, raw_reads, uniq_kmers, ncores, args);
            size_t gsize = sorted_kmers->Size();
            if(gsize > 0.75*genome_size) {
                GetBranchesAndSaveKmers(*sorted_kmers, ncores, args);
                break;
            } else {
                max_kmer -= read_len/25;
                max_kmer -= 1-max_kmer%2;
            }
        }
        cerr << endl << "Max kmer: " << max_kmer << endl;

        double alpha = double(max_kmer-min_kmer)/(steps-1);
        for(int step = 1; step < steps-1; ++step) {
            int kmer_len = min_kmer+step*alpha+0.5;
            kmer_len -= 1-kmer_len%2;
            sorted_kmers = GetSortedKmers(kmer_len, raw_reads, uniq_kmers, ncores, args);
            GetBranchesAndSaveKmers(*sorted_kmers, ncores, args);
        }
        for(int j = 0; j < 2; ++j) {
            int kmer_len = max_kmer-(j+0.5)*alpha+0.5;
            kmer_len -= 1-kmer_len%2;
            sorted_kmers = GetSortedKmers(kmer_len, raw_reads, uniq_kmers, ncores, args);
            GetBranchesAndSaveKmers(*sorted_kmers, ncores, args);
        }
    }
    */

    return 0;
}

int main(int argc, const char* argv[])
{
    // Execute main application function
    return CDeBruijnApplication().AppMain(argc, argv, 0, eDS_ToStderr);
}
