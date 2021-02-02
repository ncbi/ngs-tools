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

#include "reader.h"
#include "fasta_reader.h"
#ifndef NO_NGS_SUPPORT
#include "vdb_reader.h"
#endif
#include "mt_reader.h"
#include "aux_reader.h"
#include "omp_adapter.h"

template <typename ReaderImpl, typename... ReaderArgs>
static ReaderPtr create_wrapped(bool split_non_atgc, ReaderArgs... args) {
    if (split_non_atgc) {
        return ReaderPtr(new SplittingReader<ReaderImpl>(args...));
    } else {
        return ReaderPtr(new CuttingReader<ReaderImpl>(args...));
    }
}

template <typename ReaderImpl, typename... ReaderArgs>
static ReaderPtr create_filtered(const std::string& filter_file, bool exclude_filter, bool split_non_atgc, ReaderArgs... args) {
    if (!filter_file.empty()) {
        if (exclude_filter) {
            return create_wrapped<FilteringReader<ReaderImpl, ExcludeFileSpotFilter>>(split_non_atgc, filter_file, args...);
        } else {
            return create_wrapped<FilteringReader<ReaderImpl, IncludeFileSpotFilter>>(split_non_atgc, filter_file, args...);
        }
    } else {
        return create_wrapped<ReaderImpl>(split_non_atgc, args...);
    }
}

template <typename ReaderImpl, typename... ReaderArgs>
static ReaderPtr create_threaded(const std::string& filter_file, bool exclude_filter, bool split_non_atgc, int thread_count, size_t chunk_size, ReaderArgs... args) {
    if (thread_count < 0) {
//        thread_count = std::max(omp_get_max_threads() / 2, 1);
//        thread_count = std::max(std::thread::hardware_concurrency() / 2, 1u);
        thread_count = 1; //bug SRA-6371 - does not work with stdin std::max(std::thread::hardware_concurrency() / 2, 1u); workaround over ngs high overhead
    }
    if (chunk_size == 0) {
        chunk_size = Reader::DEFAULT_CHUNK_SIZE;
    }
    if (thread_count > 1) {
        return create_filtered<MTReader<ReaderImpl>>(filter_file, exclude_filter, split_non_atgc, thread_count, chunk_size, args...);
    } else {
        return create_filtered<ReaderImpl>(filter_file, exclude_filter, split_non_atgc, args...);
    }
}

ReaderPtr Reader::create(const std::string& path, const Reader::Params& params) {
    if (FastaReader::is_fasta(path)) {
        return create_threaded<FastaReader>(params.filter_file, params.exclude_filter, params.split_non_atgc, params.thread_count, params.chunk_size, path);
    } else {
#ifdef NO_NGS_SUPPORT
		throw std::runtime_error("Please rebuild the software with NGS library support to read runs directly or use pipes to stdin instead.");
#else
        if (!params.unaligned_only && AlignedVdbReader::is_aligned(path)) {
            return create_threaded<AlignedVdbReader>(params.filter_file, params.exclude_filter, params.split_non_atgc, params.thread_count, params.chunk_size, path, params.read_qualities);
        } else {
            return create_threaded<VdbReader>(params.filter_file, params.exclude_filter, params.split_non_atgc, params.thread_count, params.chunk_size, path, params.read_qualities, params.unaligned_only);
        }
#endif
    }
}
