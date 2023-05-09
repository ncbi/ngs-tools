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

#pragma once

#include "aligns_to_dbs_job.h"
#include "dbss.h"

struct DBSSJob : public DBSJob
{
    DBSSJob(const std::string &dbss, const std::string &dbss_tax_list, int num_threads)
    {
        auto dbss_reader = DBSS::make_reader(dbss);
        kmer_len = dbss_reader->header.kmer_len;

        DBSS::DBSAnnotation annotation;
        auto sum_offset = DBSS::load_dbs_annotation(DBSS::DBSAnnot::annotation_filename(dbss), annotation);
        dbss_reader->check_consistency(sum_offset);

        auto tax_list = DBSS::load_tax_list(dbss_tax_list);
        DBSS::load_dbss(hash_array, dbss_reader, tax_list, annotation, num_threads);
    }
};
