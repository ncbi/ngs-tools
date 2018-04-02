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

#include <iostream>
#include <fstream>
#include <stdexcept>

struct FilterDB
{
    static int predicted(const std::string &kmer)
    {
	    int pred = 0;

	    for (int i = 4; i < int(kmer.length()); i++)
		    if (kmer[i] == kmer[i - 1] && kmer[i - 1] == kmer[i - 2] && kmer[i - 2] == kmer[i - 3] && kmer[i - 3] == kmer[i - 4])
			    pred++;
		    else if (i >= 8 && kmer[i] == kmer[i - 2] && kmer[i - 2] == kmer[i - 4] && kmer[i - 4] == kmer[i - 6] && kmer[i - 6] == kmer[i - 8])
			    pred++;

	    return pred;
    }

    static int min_score(int kmer_len) { return 13 * kmer_len/32; }// const 13 was designed for 32 bp kmers
};
