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

#include <cstdint>

typedef uint64_t hash_t;

struct KmerHash
{
    typedef uint64_t hash_of_hash_t;

    static hash_of_hash_t hash_of(hash_t hash)
    {
	    return fnv1_hash(&hash, sizeof(hash));
    }

private:
    static uint64_t fnv1_hash(void *key, int n_bytes)
    {
        unsigned char *p = (unsigned char *)key;
        uint64_t h = 14695981039346656037UL;
    
        for (int i = 0; i < n_bytes; i++)
            h = (h * 1099511628211) ^ p[i];

        return h;
    }
};

