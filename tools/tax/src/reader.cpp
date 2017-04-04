#include "reader.h"
#include "fasta_reader.h"
#include "vdb_reader.h"
#include "mt_reader.h"
#include "aux_reader.h"

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
        thread_count = std::max(std::thread::hardware_concurrency() / 2, 1u);
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
        if (!params.unaligned_only && AlignedVdbReader::is_aligned(path)) {
            return create_threaded<AlignedVdbReader>(params.filter_file, params.exclude_filter, params.split_non_atgc, params.thread_count, params.chunk_size, path, params.read_qualities);
        } else {
            return create_threaded<VdbReader>(params.filter_file, params.exclude_filter, params.split_non_atgc, params.thread_count, params.chunk_size, path, params.read_qualities, params.unaligned_only);
        }
    }
}
