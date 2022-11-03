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

#include "config_split_dbss.h"
#include "dbss.h"
#include "io.h"

using namespace std;

void write_header(string split_folder, int kmer_len)
{
    ofstream f(split_folder + "/header");
    f << kmer_len;
}

int main(int argc, char const *argv[])
{
	Config config(argc, argv);
    DBSS::DBSSFileReader dbss_reader(config.dbss);

    DBSS::DBSAnnotation annotation;
    auto annotation_filename = DBSS::DBSAnnot::annotation_filename(config.dbss);
    auto sum_offset = DBSS::load_dbs_annotation(annotation_filename, annotation);
    dbss_reader.check_consistency(sum_offset);

    string split_folder = config.dbss + ".split";
    if (!IO::create_folder(split_folder))
        throw std::runtime_error("failed to create " + split_folder);

    write_header(split_folder, dbss_reader.header.kmer_len);

    std::vector<hash_t> hashes;
    for (auto& annot : annotation) 
    {
        dbss_reader.load_kmers(hashes, annot.tax_id, annot);
        DBSIO::save_dbs(DBSS::DBSSFolderReader::tax_id_to_filename(split_folder, annot.tax_id), hashes, dbss_reader.header.kmer_len);
    }
}
