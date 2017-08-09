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

#include <ngs/ncbi/NGS.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>
#include <klib/printf.h>

class BaseVdbReader: public Reader {
protected:
    const bool read_qualities;
	ngs::ReadCollection run;
    size_t spot_idx;

    BaseVdbReader(const std::string& acc, bool read_qualities)
        : read_qualities(read_qualities)
        , run(ncbi::NGS::openReadCollection( acc ))
        , spot_idx(0)
    {}

    static void apply_qualities(std::string& bases, const ngs::Fragment& fragment) {
        auto qualities = fragment.getFragmentQualities();
        assert(bases.size() == qualities.size());
        auto qualities_data = qualities.data();
        for (size_t i = 0; i < bases.size(); ++i) {
            static const int SAM_QUALITY_BASE = 33;
            static const int MIN_GOOD_QUALITY = SAM_QUALITY_BASE + 3;
            const unsigned char quality = qualities_data[i];
            if (quality < MIN_GOOD_QUALITY) {
                bases[i] = '!'; // not N for easy testing
            }
        }
    }

    template <typename NGSFragment>
    void copy(NGSFragment& fragment, Fragment* output, size_t numeric_spot_id = 0) {
        if (output) {
            auto bases = fragment.getFragmentBases();
            output->bases.assign(bases.data(), bases.size());
            if (read_qualities) {
                apply_qualities(output->bases, fragment);
            }

            if (numeric_spot_id) {
                // real spotid is way too slow, using index hack
                //auto spotid = it.getReadId();
                //output->spotid.assign(spotid.data(), spotid.size());
                char buffer[32];
                size_t written;
                string_printf(buffer, sizeof(buffer), &written, "%zd", numeric_spot_id);
                assert(written < sizeof(buffer));
                output->spotid = buffer;
            } else {
                auto spotid = fragment.getReadId();
                // leaving only last part of dot-separated spot it
                // to match non-aligned vdb reader spotid
                size_t start_pos = 0;
                for (size_t i = spotid.size() - 1; i > 0; --i) {
                    if (spotid.data()[i] == '.') {
                        start_pos = i + 1;
                        break;
                    }
                }
                output->spotid.assign(spotid.data() + start_pos, spotid.size() - start_pos);
            }
        }
    }

    SourceStats stats_for_category(ngs::Read::ReadCategory category) const {
        SourceStats res;
        res.spot_count = run.getReadCount(category);
        res.expected_spot_count = res.spot_count;

        auto it = run.getReads(ngs::Read::all);
        if (it.nextRead()) {
            res.frags_per_spot = it.getNumFragments();
        } else {
            res.frags_per_spot = 0;
        }

        return res;
    }

    bool next_fragment(ngs::ReadIterator& rit) {
        while (!rit.nextFragment()) {
            ++spot_idx;
            if (!rit.nextRead()) {
                return false;
            }
        }
        return true;
    }

public:

    static bool is_aligned(const std::string& acc) {
        auto run = ncbi::NGS::openReadCollection(acc);
        auto alit = run.getAlignments(ngs::Alignment::primaryAlignment);
        return alit.nextAlignment();
    }
};

class VdbReader final: public BaseVdbReader {
private:
    const ngs::Read::ReadCategory category;
    ngs::ReadIterator it;
    bool eof;
    size_t spot_count;

public:
	VdbReader(const std::string& acc, bool read_qualities = false, bool unaligned_only = false)
        : BaseVdbReader(acc, read_qualities)
        , category(unaligned_only ? ngs::Read::unaligned : ngs::Read::all)
        , it(run.getReads(category))
    {
        spot_count = run.getReadCount(category);
        eof = !it.nextRead();
    }

    SourceStats stats() const override { return stats_for_category(category); }

    float progress() const override {
        return spot_count ? float(spot_idx) / spot_count : 1;
    }

    bool read(Fragment* output) override {
        if (eof) {
            return false;
        }

        if (!next_fragment(it)) {
            eof = true;
            return false;
        }

        if (category != ngs::Read::all) { // unable to use spot_idx hack for unalgined-only reads
            copy(it, output);
        } else {
            copy(it, output, spot_idx + 1);
        }
        return true;
    }
};

class AlignedVdbReader final: public BaseVdbReader {
private:
    enum State {
        READING_ALIGNMNETS,
        READING_PARTIAL_READS,
        READING_UNALIGNED_READS,
        READING_EOF,
    };

    ngs::AlignmentIterator alit;
    ngs::ReadIterator pit;
	ngs::ReadIterator uit;
    State state;
    size_t alignment_idx;
    size_t alignment_count;
    size_t spot_count;

public:
	AlignedVdbReader(const std::string& acc, bool read_qualities = false)
        : BaseVdbReader(acc, read_qualities)
        , alit(run.getAlignments(ngs::Alignment::primaryAlignment))
        , pit(run.getReads(ngs::Read::partiallyAligned))
        , uit(run.getReads(ngs::Read::unaligned))
        , state(READING_ALIGNMNETS)
        , alignment_idx(0)
    {
        spot_count = run.getReadCount() - run.getReadCount(ngs::Read::fullyAligned);
        alignment_count = run.getAlignmentCount(ngs::Alignment::primaryAlignment);
    }

    SourceStats stats() const override { return stats_for_category(ngs::Read::ReadCategory::all); }

    float progress() const override {
        const size_t current = alignment_idx + spot_idx;
        const size_t total = alignment_count + spot_count;
        return total ? float(current) / total : 1;
    }

    bool read(Fragment* output) override {
        switch (state) {
        case READING_ALIGNMNETS:
            if (alit.nextAlignment()) {
                copy(alit, output);
                ++alignment_idx;
                return true;
            } else {
                if (pit.nextRead()) {
                    state = READING_PARTIAL_READS;
                } else if (uit.nextRead()) {
                    state = READING_UNALIGNED_READS;
                } else {
                    state = READING_EOF;
                }
                return read(output);
            }
        case READING_PARTIAL_READS:
            while (next_fragment(pit)) {
                if (!pit.isAligned()) {
                    copy(pit, output);
                    return true;
                }
            }
            if (uit.nextRead()) {
                state = READING_UNALIGNED_READS;
            } else {
                state = READING_EOF;
            }
            return read(output);
        case READING_UNALIGNED_READS:
            if (next_fragment(uit)) {
                copy(uit, output);
                return true;
            } else {
                state = READING_EOF;
                return false;
            }
        case READING_EOF:
            return false;
        default:
            assert(false);
            return false;
        }
    }
};
