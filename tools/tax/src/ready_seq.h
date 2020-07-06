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

#include "fasta.h"
#include "seq_cleaner.h"

struct ReadySeq
{
    std::string seq;
    SeqCleaner::p_strings clean_strings;
};

static void swap(ReadySeq &a, ReadySeq &b)
{
    swap(a.seq, b.seq);
    swap(a.clean_strings, b.clean_strings);
}

static void load_sequence(Fasta *_fasta, ReadySeq *_seq)
{
    Fasta &fasta = *_fasta;
    ReadySeq &seq = *_seq;

    seq.seq.clear(); // for better performance clear instead of constructor
    seq.clean_strings.clear();

    if (!fasta.get_next_sequence(seq.seq))
        return;

    SeqCleaner cleaner(seq.seq);
    seq.clean_strings = move(cleaner.clean_strings);
}

