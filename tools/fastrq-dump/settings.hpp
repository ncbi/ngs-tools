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
#ifndef _hpp_ngs_settings_
#define _hpp_ngs_settings_

#include <stdlib.h>
#include <stdint.h>


namespace fastrq
{
    enum FmtType
    {
        fmt_fastq,
        fmt_fasta
    };
    
    enum OpType
    {
        op_reads,
        op_refseq
    };
    
    struct FastRQSettings
    {
        void validate ();

        FastRQSettings ()
        : fmt ( fmt_fastq ), op ( op_reads ), min_length ( 0 ), n_count ( 0 )
        , num_accessions ( 0 )
        {
        }
                
        FmtType fmt; // an enum of output format, currently fastq or fasta
        OpType op;

        //Filter Settings
        size_t min_length;
        uint32_t n_count;
        uint32_t num_accessions;
    };
}

#endif // _hpp_ngs_settings_
