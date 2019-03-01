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

#include <stdint.h>
#include <string>
#include <vector>
#include <assert.h>
#include <memory>
#include <algorithm>

class Reader;
typedef std::unique_ptr<Reader> ReaderPtr;

class Reader {
public:
    static const size_t DEFAULT_CHUNK_SIZE = 1024; // todo: tune

    struct Fragment
    {
        std::string spotid; // unique spot identifier, multiple fragments can have the same spotid
        std::string bases; // must not be empty
        bool operator == (const Fragment& other) const { return spotid == other.spotid && bases == other.bases; }
    };

    struct SourceStats {
        size_t spot_count;
        size_t expected_spot_count; // hint, tries to account for filtering
        int frags_per_spot;

        SourceStats(size_t spot_count, int frags_per_spot = 1)
            : spot_count(spot_count)
            , expected_spot_count(spot_count)
            , frags_per_spot(frags_per_spot)
        {}
        SourceStats() : SourceStats(0) {}

        size_t frag_count() const { return spot_count * frags_per_spot; }
        size_t expected_frag_count() const { return expected_spot_count * frags_per_spot; }

        bool operator== (const SourceStats& other) const {
            return spot_count == other.spot_count
                && expected_spot_count == other.expected_spot_count
                && frags_per_spot == other.frags_per_spot;
        };
    };

    virtual ~Reader() {}

    // returns stats of original file
    virtual SourceStats stats() const = 0;

    // returns [0-1] value
    virtual float progress() const = 0;

    // if output is not null, reads one fragment into it
    // if output is null simply skips one fragment
    // returns true if fragment was read
    // returns false when at eof
    virtual bool read(Fragment* output) = 0;

    // read the most efficient count of fragments into output
    // replaces output content
    // returns true if anything was read (i.e. output is not empty)
    virtual bool read_many(std::vector<Fragment>& output) {
        output.resize(std::max(output.capacity(), size_t(DEFAULT_CHUNK_SIZE)));
        for (size_t i = 0; i < output.size(); ++i) {
            if (!read(&output[i])) {
                output.resize(i);
                break;
            }
        }
        return !output.empty();
    }

    // factory params aux struct
    struct Params {
        std::string filter_file;
        bool exclude_filter; // if true, inverse filtering, i.e. exclude spots listed in filter file
        bool read_qualities; // if true, low quality bases replaced with N
        bool split_non_atgc; // if true, splits reads by non-atgc characters, otherwise cuts reads at first non-atgc character
        bool unaligned_only; // if true, skips aligned reads
        int thread_count; // default means auto
        size_t chunk_size; ; // default means auto
        Params() : exclude_filter(false), read_qualities(false), split_non_atgc(false), unaligned_only(false), thread_count(-1), chunk_size(0) {}
    };
    // factory method, creates corresponding reader depending on file type
    static ReaderPtr create(const std::string& path, const Params& params = Params());
};
