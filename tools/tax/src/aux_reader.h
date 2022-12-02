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

#include "reader.h"
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include "checksum.h"
#include "log.h"

//test tc

static bool is_actg(char const ch) { return (ch == 'A') | (ch == 'C') | (ch == 'T') | (ch == 'G'); }
static bool non_actg(char const ch) { return !is_actg(ch); }

class SpotFilter {
public:
    virtual ~SpotFilter() {}
    virtual size_t expected_spot_count() const { return -1; }
    virtual bool is_good(const std::string& spotid) const = 0;
};

template <typename Callable>
class CallableSpotFilter final: public SpotFilter {
private:
    const Callable callable;
public:
    template <typename... Args>
    CallableSpotFilter(Args... args) : callable(args...) {}

    bool is_good(const std::string& spotid) const override { return callable(spotid); }
};

class BaseFileSpotFilter: public SpotFilter {
protected:
    std::unordered_set<std::string> file_spots;
public:
    BaseFileSpotFilter(const std::string& path) {
        std::ifstream is(path);
        std::string line;
        while (std::getline(is, line)) {
            auto dot_it = std::find(line.begin(), line.end(), '.');
            line.resize(dot_it - line.begin()); // strip dot and everything after
            file_spots.insert(line);
        }
    }
};

class IncludeFileSpotFilter final: public BaseFileSpotFilter {
public:
    IncludeFileSpotFilter(const std::string& path) : BaseFileSpotFilter(path) {}
    size_t expected_spot_count() const override { return file_spots.size(); }
    bool is_good(const std::string& spotid) const override { return file_spots.count(spotid) > 0; }
};

class ExcludeFileSpotFilter final: public BaseFileSpotFilter {
public:
    ExcludeFileSpotFilter(const std::string& path) : BaseFileSpotFilter(path) {}
    bool is_good(const std::string& spotid) const override { return file_spots.count(spotid) == 0; }
};

// only keeps spots with spotids passing the filter
template <typename ReaderType, typename FilterType>
class FilteringReader final: public Reader {
private:
    ReaderType reader;
    FilterType filter;
    
public:
    template <typename... ReaderArgs>
    FilteringReader(FilterType&& filter, ReaderArgs... reader_args)
        : reader(reader_args...)
        , filter(filter)
    {}
    template <typename... ReaderArgs>
    FilteringReader(const FilterType& filter, ReaderArgs... reader_args)
        : reader(reader_args...)
        , filter(filter)
    {}

    SourceStats stats() const override {
        auto stats = reader.stats();
        stats.expected_spot_count = std::min(stats.expected_spot_count, filter.expected_spot_count());
        return stats;
    }
    float progress() const override { return reader.progress(); }

    virtual bool read(Fragment* output) override {
        Fragment temp;
        if (!output) {
            output = &temp;
        }
        while (reader.read(output)) {
            if (filter.is_good(output->spotid)) {
                return true;
            }
        }
        return false;
    }
};

// splits reads by non-atgs values
template <typename ReaderType>
class SplittingReader final: public Reader {
private:
    ReaderType reader;
    Fragment last;
    size_t offset;

#if AUX_READER_LOG_CHECKSUM
    CheckSum checksum;
#endif

public:
    template <typename... ReaderArgs>
    SplittingReader(ReaderArgs... reader_args) : reader(reader_args...), offset(0){}

    SourceStats stats() const override { return reader.stats(); }
    float progress() const override { return reader.progress(); }

    virtual bool read(Fragment* output) override {
        while (true) {
            if (offset >= last.bases.size()) {
                if (!reader.read(&last)) {
#if AUX_READER_LOG_CHECKSUM
                    checksum.log();
#endif
                    last.bases.clear(); // just in case if read messed with last.bases
                    return false;
                }
#if AUX_READER_LOG_CHECKSUM
                checksum.update(last.bases);
#endif
                offset = 0;
            }
            auto from = std::find_if(last.bases.begin() + offset, last.bases.end(), is_actg);
            auto to = std::find_if(from, last.bases.end(), non_actg);
            if (from == last.bases.begin() && to == last.bases.end()) {
                if (output) {
#if AUX_READER_NO_BASES
                    last.bases = std::string();
#endif
                    std::swap(*output, last);
                }
                offset = last.bases.size();
                return true;
            } else {
                offset = to - last.bases.begin();
                if (from != to) {
                    if (output) {
                        output->spotid = last.spotid;
                        output->bases.assign(from, to);
                    }
                    return true;
                }
            }
        }
    }
};

// cuts reads at first non-atgc value, consumes empty reads
template <typename ReaderType>
class CuttingReader final: public Reader {
private:
    ReaderType reader;
    
public:
    template <typename... ReaderArgs>
    CuttingReader(ReaderArgs... reader_args) : reader(reader_args...) {}

    SourceStats stats() const override { return reader.stats(); }
    float progress() const override { return reader.progress(); }

    bool read(Fragment* output) override {
        while (true) {
            Fragment temp;
            if (!output) {
                output = &temp;
            }
            const bool res = reader.read(output);
            if (res) {
                auto it = find_if(output->bases.begin(), output->bases.end(), non_actg);
                if (it == output->bases.begin()) {
                    continue;
                } else if (it != output->bases.end()) {
                    output->bases.resize(it - output->bases.begin());
                }
            }
            return res;
        }
    }
};

// cuts reads at first non-atgc value, consumes empty reads
template <typename ReaderType>
class UltraFastSkipReader final: public Reader {
private:
    ReaderType reader;
    int skip_step = 1;

public:
    template <typename... ReaderArgs>
    UltraFastSkipReader(int skip_step, ReaderArgs... reader_args) : reader(reader_args...), skip_step(skip_step) 
    {
        LOG("UltraFastSkipReader " << skip_step);
    }

    SourceStats stats() const override { return reader.stats(); }
    float progress() const override { return reader.progress(); }

    bool read(Fragment* output) override {
        for (int i = 0; i < skip_step; i++)
            if (!reader.read(output))
                return false;
            
        return reader.read(output);
    }
};
