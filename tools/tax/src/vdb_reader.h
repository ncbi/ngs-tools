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

    static bool is_aligned(const ngs::ReadCollection& run) {
        auto alit = run.getAlignments(ngs::Alignment::primaryAlignment);
        return alit.nextAlignment();
    }

public:
    static bool is_aligned(const std::string& acc) {
        return is_aligned(ncbi::NGS::openReadCollection(acc));
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
        , category((unaligned_only && is_aligned(run)) ? ngs::Read::unaligned : ngs::Read::all)
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

#include "vdb/cursor.h"
#include "vdb/table.h"
#include "vdb/database.h"
#include "vdb/manager.h"

class FastVdbReader final: public Reader {
    class VdbBaseReader {
    protected:
        VCursor const *curs;
        uint64_t row;
        uint64_t rowCount;
        int64_t firstRow;

        /* boilerplate code for accessing VDB */
        template <typename T>
        std::pair<T const *, uint32_t> columnDataNoThrow(uint32_t const cid)
        {
            void const *data = nullptr;
            uint32_t bits = 0;
            uint32_t boff = 0;
            uint32_t count = 0;
            if (0 != VCursorCellDataDirect(curs, row + firstRow, cid, &bits, &data, &boff, &count))
                return {nullptr, 0};
            assert(boff == 0);
            assert(bits == 8 * sizeof(T));
            return {reinterpret_cast<T const *>(data), count};
        }

        template <typename T>
        std::pair<T const *, uint32_t> columnData(uint32_t const cid)
        {
            void const *data = nullptr;
            uint32_t bits = 0;
            uint32_t boff = 0;
            uint32_t count = 0;
            if (0 != VCursorCellDataDirect(curs, row + firstRow, cid, &bits, &data, &boff, &count))
                throw std::runtime_error("Can't read row");
            assert(boff == 0);
            assert(bits == 8 * sizeof(T));
            return {reinterpret_cast<T const *>(data), count};
        }


        template <typename T>
        T const *columnData(uint32_t const cid, uint32_t const requiredCount)
        {
            void const *data = nullptr;
            uint32_t bits = 0;
            uint32_t boff = 0;
            uint32_t count = 0;
            if (0 != VCursorCellDataDirect(curs, row + firstRow, cid, &bits, &data, &boff, &count))
                throw std::runtime_error("Can't read row");
            assert(boff == 0);
            assert(bits == 8 * sizeof(T));
            assert(requiredCount == count);
            return reinterpret_cast<T const *>(data);
        }

        template <typename T>
        T readValue(uint32_t const cid) {
            return *columnData<T>(cid, 1);
        }

        uint32_t addColumn(char const *const col_expr) {
            auto cid = ~(uint32_t(0));
            if (0 != VCursorAddColumn(curs, &cid, col_expr))
                throw std::runtime_error("Can't add columns to cursor!");
            return cid;
        }

        void open() {
            if (0 != VCursorOpen(curs))
                throw std::runtime_error("Can't open cursor!");
        }

        void getRowRange(uint32_t const cid) {
            if (0 != VCursorIdRange(curs, cid, &firstRow, &rowCount))
                throw std::runtime_error("Can't cursor row range!");
        }

        VdbBaseReader(VTable const *tbl)
        : curs(nullptr)
        , row(0)
        , rowCount(0)
        , firstRow(0)
        {
            if (0 != VTableCreateCursorRead(tbl, &curs))
                throw std::runtime_error("Can't make a cursor!");
        }
        virtual ~VdbBaseReader() {
            VCursorRelease(curs);
        }
    public:
        virtual bool read(Fragment *) = 0;
    };
    class AlignReader final : public VdbBaseReader {
        uint32_t cid_read;
        uint32_t cid_spotId;
    public:
        AlignReader(VTable const *tbl)
        : VdbBaseReader(tbl)
        , cid_read(addColumn("RAW_READ"))
        , cid_spotId(addColumn("SEQ_SPOT_ID"))
        {
            open();
            getRowRange(cid_read);
        }
        bool read(Fragment *output) {
            auto const ps = columnDataNoThrow<char>(cid_read);
            if (ps.first) {
                auto const spotId = readValue<int64_t>(cid_spotId);

                ++row;
                output->bases = std::string(ps.first, ps.second);
                output->spotid = std::to_string(spotId);
                return true;
            }
            return false;
        }
    };

    /* This is for getting the reads of
     * *either* an unaligned run
     * *or* the unaligned reads of an
     * aligned run.
     */
    class SeqReader final : public VdbBaseReader {
        uint8_t const *readType;
        uint8_t const *readFilter;
        uint32_t const *readLen;
        char const *bases;
        uint32_t const *readStart;
        int64_t const *alignId;
        uint32_t sReadStart[2];
        uint32_t *hReadStart;

        uint64_t spotCount;
        uint64_t readCount;
        int64_t firstUnaligned;

        uint32_t cid_alignId;
        uint32_t cid_read;
        uint32_t cid_start;
        uint32_t cid_length;
        uint32_t cid_type;
        uint32_t cid_filter;

        unsigned readNo;
        unsigned numreads;
        unsigned max_numreads;

        bool const aligned;
        bool const useReadFilter;

        bool isFiltered(int const read) const {
            return useReadFilter && (readFilter[read] != 0);
        }
        bool isAligned(int const read) const {
            return aligned && (alignId[read] != 0);
        }

        bool loadReadType() {
            /* this is the first column we will load in any row,
             * if we can't read it, we won't be able to read the
             * rest of the columns.
             */
            auto const &ptr_count = columnDataNoThrow<uint8_t>(cid_type);
            readType = ptr_count.first;
            numreads = ptr_count.second;
            return readType != nullptr;
        }
        void loadBases() {
            auto const ps = columnData<char>(cid_read);
            bases = ps.first;
        }
        void loadFilters(bool loadAlignId) {
            if (useReadFilter)
                readFilter = columnData<uint8_t>(cid_filter, numreads);
            if (loadAlignId)
                alignId = columnData<int64_t>(cid_alignId, numreads);
        }
        void loadAligned() {
            loadFilters(true);
            readLen = columnData<uint32_t>(cid_start, numreads);
            loadBases();

            // Using CMP_READ instead of READ
            // READ_START is about READ
            // Need to compute own version of READ_START
            // that works for CMP_READ, which contains
            // only unaligned bases
            if (max_numreads < numreads) {
                max_numreads = numreads;
                if (hReadStart)
                    delete [] hReadStart;
                hReadStart = new uint32_t[max_numreads];
            }
            uint32_t *const rs = (numreads <= 2)
                               ? (&sReadStart[0])
                               : (hReadStart);
            uint32_t start = 0;
            for (unsigned i = 0; i < numreads; ++i) {
                rs[i] = start;
                if (alignId[i] == 0)
                    start += readLen[i];
            }
            readStart = rs;
        }
        void loadUnaligned() {
            loadFilters(false);
            readStart = columnData<uint32_t>(cid_start, numreads);
            readLen = columnData<uint32_t>(cid_start, numreads);
            loadBases();
        }
        void addAlignedColumns() {
            cid_type = addColumn("READ_TYPE");
            if (useReadFilter)
                cid_filter = addColumn("READ_FILTER");
            cid_alignId = addColumn("PRIMARY_ALIGNMENT_ID");
            cid_length = addColumn("READ_LEN");
            cid_read = addColumn("CMP_READ");
        }
        void addUnalignedColumns() {
            cid_type = addColumn("READ_TYPE");
            if (useReadFilter)
                cid_filter = addColumn("READ_FILTER");
            cid_start = addColumn("READ_START");
            cid_length = addColumn("READ_LEN");
            cid_read = addColumn("READ");
        }
    public:
        SeqReader(VTable const *tbl, bool is_Aligned, bool filtered = false)
        : VdbBaseReader(tbl)
        , hReadStart(nullptr)
        , spotCount(0)
        , readCount(0)
        , firstUnaligned(0)
        , max_numreads(2)
        , aligned(is_Aligned)
        , useReadFilter(filtered)
        {
            if (aligned)
                addAlignedColumns();
            else
                addUnalignedColumns();
            open();
            getRowRange(aligned ? cid_alignId : cid_read);

            auto foundFirstUnaligned = false;
            for (row = 0; row < rowCount; ++row) {
                if (!loadReadType())
                    continue;
                loadFilters(aligned);

                unsigned reads = 0;
                unsigned alignedReads = 0;
                for (unsigned readNo = 0; readNo < numreads; ++readNo) {
                    if ((readType[readNo] & 1) != 1)
                        continue;
                    if (isFiltered(readNo))
                        continue;
                    if (isAligned(readNo))
                        ++alignedReads;
                    ++reads;
                }
                readCount += reads;
                if (reads > 0) {
                    if (!foundFirstUnaligned && alignedReads < reads) {
                        firstUnaligned = row;
                        foundFirstUnaligned = true;
                    }
                    spotCount += 1;
                }
            }
            row = firstUnaligned;
        }
        ~SeqReader() {
            if (hReadStart)
                delete [] hReadStart;
        }
        SourceStats stats() const {
            return SourceStats(spotCount, readCount);
        }

        float progress(uint64_t const processed) const
        {
            auto const total = readCount;
            return total ? double(processed) / total : 1.0;
        }

        bool read(Fragment *output) {
            if (readNo >= numreads) {
                if (row >= rowCount)
                    return false;
                if (!loadReadType())
                    return false;
                if (aligned)
                    loadAligned();
                else
                    loadUnaligned();
                readNo = 0;
                ++row;
            }
            auto const thisRead = readNo++;

            if ((readType[thisRead] & 1) != 1)
                return read(output); // it's not biological

            if (isFiltered(thisRead) || isAligned(thisRead))
                return read(output);

            output->bases = std::string(bases + readStart[thisRead], readLen[thisRead]);
            output->spotid = std::to_string(row + firstRow - 1);
            return true;
        }
    };
    AlignReader *alg;
    SeqReader *seq;
    VdbBaseReader *current;
    uint64_t readsProcessed;
public:
    /* Gets aligned and unaligned reads
     * Never accesses quality scores
     */
    FastVdbReader(std::string const &acc, bool unaligned_only)
    {
        VDBManager const *mgr = nullptr;
        VDatabase const *db = nullptr;
        VTable const *seqtbl = nullptr;
        VTable const *algtbl = nullptr;
        rc_t rc = 0;

        rc = VDBManagerMakeRead(&mgr, nullptr);
        if (rc != 0)
            throw std::runtime_error("Can't make a VDB Manager!");

        /* it is okay for this to fail, e.g. it is not aligned */
        rc = VDBManagerOpenDBRead(mgr, &db, nullptr, "%s", acc.c_str());
        if (rc == 0) {
            VDBManagerRelease(mgr);

            if (!unaligned_only) {
                /* it is okay for this to fail */
                (void)VDatabaseOpenTableRead(db, &algtbl, "PRIMARY_ALIGNMENT");
            }

            /* it is not okay for this to fail */
            rc = VDatabaseOpenTableRead(db, &seqtbl, "SEQUENCE");
            VDatabaseRelease(db);
            if (rc != 0) {
                VTableRelease(algtbl);

                throw std::runtime_error(std::string("Can't read sequences from ") + acc);
            }
        }
        else {
            rc = VDBManagerOpenTableRead(mgr, &seqtbl, nullptr, "%s", acc.c_str());
            VDBManagerRelease(mgr);
            if (rc != 0)
                throw std::runtime_error(std::string("Can't open ") + acc);
        }
        seq = new SeqReader(seqtbl, algtbl != NULL);
        if (algtbl) {
            alg = new AlignReader(algtbl);
            current = alg;
        }
        else {
            alg = nullptr;
            current = seq;
        }
    }
    ~FastVdbReader() {
        delete seq;
        if (alg)
            delete alg;
    }

    SourceStats stats() const override
    {
        return seq->stats();
    }

    float progress() const override
    {
        return seq->progress(readsProcessed);
    }

    bool read(Fragment* output) override
    {
        if (current) {
            if (current->read(output)) {
                readsProcessed += 1;
                return true;
            }
            if (current == alg)
                current = seq;
            else
                current = nullptr;
            return read(output);
        }
        return false;
    }
};
