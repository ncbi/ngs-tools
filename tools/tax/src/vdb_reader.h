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
#include "log.h"

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

#include "kdb/meta.h"
#include "vdb/cursor.h"
#include "vdb/table.h"
#include "vdb/database.h"
#include "vdb/manager.h"

#include <map>
#include <algorithm>
#include <cstdlib>
#include <cinttypes>

struct FastVdbReader_Counts {
    uint64_t spots;
    uint64_t reads;
    uint64_t aligned;

    std::ostream &debugPrint(std::ostream &os) const {
        return os << "Spots: " << spots
                << ", Reads: " << reads
                << ", Aligned Reads: " << aligned;
    }
    friend std::ostream &operator <<(std::ostream &os, FastVdbReader_Counts const &counts) {
        return counts.debugPrint(os);
    }
};

class FastVdbReader final: public Reader {
    /**
     * Contains mosly boilerplate code for accessing VDB
     */
    class VdbBaseReader {
    protected:
        using RowId = int64_t;
        using RowCount = uint64_t;

        struct RowRange : public std::pair<RowId, RowCount> {
            using Base = std::pair<RowId, RowCount>;
            bool contains(RowId query) const {
                if (query < first || first + second <= query)
                    return false;
                return true;
            }
            explicit RowRange(Base const &b)
            : Base(b)
            {}
            explicit RowRange(Base &&b)
            : Base(std::move(b))
            {}
        };

        /**
         * Useful for binding data type to a column
         */
        template <typename T>
        class Column {
            static constexpr uint32_t nilColumnId() { return ~(uint32_t(0)); }
            static uint32_t add(VCursor const *curs, char const *const expr, std::string const &sourceName)
            {
                auto cid = nilColumnId();
                auto const rc = VCursorAddColumn(curs, &cid, expr);
                if (rc == 0)
                    return cid;

                throw std::runtime_error(  std::string("Can't add column '")
                                         + std::string(expr)
                                         + std::string("' to cursor on ")
                                         + sourceName);
            }
            std::string expr;
            uint32_t cid;

            [[noreturn]] void throwReadError(int64_t const row, std::string const &sourceName) const noexcept(false)
            {
                throw std::runtime_error(  std::string("Can't read row ")
                                         + std::to_string(row)
                                         + " of " + sourceName + "." + expr);
            }
        public:
            Column(VCursor const *const curs, char const *const expr, std::string const &sourceName)
            : expr(expr)
            , cid(add(curs, expr, sourceName))
            {}
            Column()
            : expr("no value")
            , cid(nilColumnId())
            {}

            using Value = T;
            using Ptr = Value const *;

            /** a single cell can not contain more than 2^32 elements */
            using Count = uint32_t;

            /** pointer + length */
            using Data = std::pair<Ptr, Count>;

            /**
             * Get the data for a row, does not throw if there is a problem.
             *
             * @return pointer + length or nullptr + 0
             */
            Data dataNoThrow(VCursor const *const curs, RowId const row) const noexcept(true)
            {
                void const *data = nullptr;
                uint32_t bits = 0;
                uint32_t boff = 0;
                uint32_t count = 0;
                assert(cid != nilColumnId());
                if (0 == VCursorCellDataDirect(curs, row, cid, &bits, &data, &boff, &count)) {
                    assert(boff == 0);
                    assert(bits == 8 * sizeof(T));
                    return {reinterpret_cast<T const *>(data), count};
                }
                return {nullptr, 0};
            }

            /**
             * Get the data for a row.
             *
             * @return pointer + length
             *
             * @throw std::runtime_error
             */
            Data data(VCursor const *const curs, RowId const row, std::string const &sourceName) const noexcept(false)
            {
                void const *data = nullptr;
                uint32_t bits = 0;
                uint32_t boff = 0;
                uint32_t count = 0;
                assert(cid != nilColumnId());
                if (0 == VCursorCellDataDirect(curs, row, cid, &bits, &data, &boff, &count)) {
                    assert(boff == 0);
                    assert(bits == 8 * sizeof(T));
                    return {reinterpret_cast<T const *>(data), count};
                }
                throwReadError(row, sourceName);
            }

            /**
             * Get the data for a row.
             *
             * @return pointer to data
             *
             * @throw std::runtime_error
             */
            Ptr data(Count const requiredCount, VCursor const *curs, RowId const row, std::string const &sourceName) const noexcept(false)
            {
                auto const result = data(curs, row, sourceName);
                if (result.second == requiredCount)
                    return result.first;
                throwReadError(row, sourceName);
            }

            /**
             * Get the data for a row.
             *
             * @return value for the row.
             *
             * @throw std::runtime_error
             */
            Value value(VCursor const *const curs, RowId const row, std::string const &sourceName) const noexcept(false)
            {
                if (sizeof(Value) > 1) {
                    union {
                        Value value;
                        uint8_t raw[sizeof(T)];
                    } u;
                    auto const p = reinterpret_cast<uint8_t const *>(data(1, curs, row, sourceName));
                    std::copy(p, p + sizeof(T), u.raw);
                    return u.value;
                }
                return *(data(1, curs, row, sourceName));
            }

            /**
             * Get the row range for the column.
             *
             * @Note VDB tables are not necessarily square, i.e. columns are not all required to be the same number of rows.
             *
             * @return row range for the column (start + count)
             *
             * @throw std::runtime_error (but very unlikely, probably programmer error)
             */
            RowRange rowRange(VCursor const *const curs, std::string const &sourceName) const {
                auto first = RowId(0);
                auto count = RowCount(0);
                assert(cid != nilColumnId());
                if (0 == VCursorIdRange(curs, cid, &first, &count))
                    return RowRange({first, count});
                throw std::runtime_error(std::string("Can't get row range of ") + sourceName);
            }
        };
        std::string const sourceName; ///< Table name
        VCursor const *curs;
        RowRange rowRange;
        RowId row; ///< current row number

        /**
         * Get the data (from the current row) for a column, does not throw if there is a problem.
         *
         * @return pointer + length or nullptr + 0
         */
        template <typename Col, typename Data = typename Col::Data>
        Data columnDataNoThrow(Col const &col) const noexcept(true)
        {
            return col.dataNoThrow(curs, row + rowRange.first);
        }

        /**
         * Get the data (from the current row) for a column.
         *
         * @return pointer + length
         *
         * @throw std::runtime_error
         */
        template <typename Col, typename Data = typename Col::Data>
        Data dataOf(Col const &col) const noexcept(false)
        {
            return col.data(curs, row + rowRange.first, sourceName);
        }

        /**
         * Get the data (from the current row) for a column.
         *
         * @return pointer to data
         *
         * @throw std::runtime_error
         */
        template <typename Col, typename Ptr = typename Col::Ptr, typename Count = typename Col::Count>
        Ptr dataOf(Col const &col, Count const requiredCount) const noexcept(false)
        {
            return col.data(requiredCount, curs, row + rowRange.first, sourceName);
        }

        /**
         * Get the value (from the current row) for a column.
         *
         * @return the value
         *
         * @throw std::runtime_error
         */
        template <typename Col, typename Value = typename Col::Value>
        Value valueOf(Col const &col)
        {
            return col.value(curs, row + rowRange.first, sourceName);
        }

        std::string stringValueOf(Column<char> const &col)
        {
            auto const &&data = dataOf(col);
            return std::string(data.first, data.second);
        }

        /**
         * Add a column (expression) to the cursor.
         *
         * @return object representing the new column.
         *
         * @throw std::runtime_error (probably programmer error)
         */
        template <typename T>
        T addColumn(char const *const col_expr) {
            return T(curs, col_expr, sourceName);
        }

        /**
         * Open the cursor (after adding all the columns)
         *
         * @throw std::runtime_error (probably programmer error)
         */
        void open() {
            if (0 != VCursorOpen(curs))
                throw std::runtime_error(std::string("Can't open cursor on ") + sourceName);
        }

        /**
         * Set the row range from a column.
         *
         * @Note VDB tables are not necessarily square, i.e. columns are not all required to be the same number of rows.
         * Pick a column that logically must have data for every row.
         *
         * @throw std::runtime_error (but very unlikely, probably programmer error)
         */
        template <typename T>
        void getRowRangeFrom(T const &col) {
            rowRange = col.rowRange(curs, sourceName);
        }

        /**
         * @param tbl the table to create the cursor on.
         * @param name informative.
         */
        VdbBaseReader(VTable const *tbl, char const *name)
        : sourceName(name)
        , curs(nullptr)
        , row(0)
        , rowRange({0, 0})
        {
            if (0 != VTableCreateCursorRead(tbl, &curs))
                throw std::runtime_error(std::string("Can't make a cursor on ") + sourceName);
        }

        /**
         * @param db the database to create the cursor in.
         * @param name the table name to create the cursor on.
         */
        VdbBaseReader(VDatabase const *const db, char const *const name)
        : sourceName(name)
        , curs(nullptr)
        , row(0)
        , rowRange({0, 0})
        {
            VTable const *tbl = nullptr;

            if (0 != VDatabaseOpenTableRead(db, &tbl, "%s", name))
                throw std::runtime_error(std::string("Can't open table ") + sourceName);

            if (0 != VTableCreateCursorRead(tbl, &curs))
                throw std::runtime_error(std::string("Can't make a cursor on ") + sourceName);

            VTableRelease(tbl);
        }

        /**
         * close the cursor
         */
        virtual ~VdbBaseReader() {
            VCursorRelease(curs);
        }

        /**
         * What table was the cursor created on.
         *
         * @Note this pointer needs to be released by the caller.
         */
        VTable const *table() const {
            VTable const *tbl = nullptr;

            VCursorOpenParentRead(curs, &tbl);

            return tbl;
        }

        /**
         * The metadata object for the table this cursor was created on.
         *
         * @Note this pointer needs to be released by the caller.
         */
        KMetadata const *metadata() const {
            KMetadata const *md = nullptr;
            VTable const *const tbl = table();

            VTableOpenMetadataRead(tbl, &md);
            VTableRelease(tbl);

            return md;
        }
    public:
        double progress() const {
            return rowRange.second ? double(row) / rowRange.second : -1.0;
        }
    };

    /**
     * A cursor on the database's Reference table.
     *
     * @see ncbi-vdb/libs/axf/align-restore-read.c
     */
    class ReferenceReader final : public VdbBaseReader {
        struct ReferenceInfo {
            using RowRange = VdbBaseReader::RowRange;
            using RowId = VdbBaseReader::RowId;
            std::string id;
            RowRange rowRange;
            uint32_t length;
            bool circular;

            bool operator <(ReferenceInfo const &rhs) const {
                return rowRange < rhs.rowRange;
            }

            ReferenceInfo(std::string id, RowId first, uint32_t length, bool circular)
            : id(id)
            , rowRange({first, 1})
            , length(length)
            , circular(circular)
            {}

            void addRow(uint32_t inc) {
                length += inc;
                rowRange.second += 1;
            }
        };
        using References = std::vector<ReferenceInfo>;
        using ReadCol = VdbBaseReader::Column<char>;

        ReadCol readCol;
        ReadCol::Data read;
        References references;
        References::const_iterator current;
        unsigned maxSeqLen; ///< this is the chunking size of the reference table; each row has logically this many bases (actual length may be less).

        /**
         * Use the absolute position to set the current row.
         *
         * @param start the absolute position (from GLOBAL_REF_START column)
         *
         * @return the offset of the absolute position in the current row.
         */
        unsigned loadAbsolute(uint64_t start) {
            auto const wantRow = start / maxSeqLen;
            auto const offset = start % maxSeqLen;

            if (row != wantRow) {
                row = wantRow;
                read = dataOf(readCol);

                if (!current->rowRange.contains(row)) {
                    auto next = current + 1;
                    if (next == references.end() || !next->rowRange.contains(row)) {
                        next = std::lower_bound(references.begin(), references.end(), *current, [&](ReferenceInfo const &a, ReferenceInfo const &b) { return a.rowRange.first < row; });
                        assert(next != references.end());
                        assert(next != current);
                        assert(next->rowRange.contains(row));
                    }
                    current = next;
                }
            }
            return (unsigned)offset;
        }
    public:
        using HasRefOffsetColData = VdbBaseReader::Column<uint8_t>::Data;
        using RefOffsetColData = VdbBaseReader::Column<int32_t>::Data;
        using HasMismatchColData = VdbBaseReader::Column<uint8_t>::Data;
        using MismatchColData = VdbBaseReader::Column<uint8_t>::Data;
        using GlobalRefStartColData = VdbBaseReader::Column<uint64_t>::Data;

        ReferenceReader(VDatabase const *const db)
        : VdbBaseReader(db, "REFERENCE")
        , readCol(addColumn<ReadCol>("(INSDC:dna:text)READ"))
        {
            auto const mslCol = addColumn<VdbBaseReader::Column<uint32_t>>("MAX_SEQ_LEN");
            auto const idCol = addColumn<VdbBaseReader::Column<char>>("SEQ_ID");
            auto const lenCol = addColumn<VdbBaseReader::Column<uint32_t>>("READ_LEN");
            auto const circCol = addColumn<VdbBaseReader::Column<uint8_t>>("CIRCULAR");

            open();
            getRowRangeFrom(idCol);
            maxSeqLen = valueOf(mslCol);
            for (row = 0; row < rowRange.second; ++row) {
                auto const id = stringValueOf(idCol);
                auto const len = valueOf(lenCol);

                if (references.empty() || id != references.back().id) {
                    auto const &&info = ReferenceInfo(id, row, len, valueOf(circCol));

                    references.emplace_back(info);
                }
                else {
                    references.back().addRow(len);
                }
            }
            current = references.begin();
            row = 0;
        }

        std::string restoreRead(  HasMismatchColData const &hasMismatchColData
                                , MismatchColData const &mismatchColData
                                , HasRefOffsetColData const &hasRefOffsetColData
                                , RefOffsetColData const &refOffsetColData
                                , uint64_t const globalStart
                                , bool const reversed
                                )
        {
            auto const readLen = hasMismatchColData.second;
            auto mi = decltype(mismatchColData.second)(0);
            auto roi = decltype(refOffsetColData.second)(0);
            auto const hasMismatch = hasMismatchColData.first;
            auto const hasRefOffset = hasRefOffsetColData.first;
            auto bi = decltype(refOffsetColData.first[0])(0);
            auto ri = decltype(bi)(0);
            auto result = std::string(readLen, 'N');
            auto f = result.begin();
            auto r = result.rbegin();

            for (auto i = decltype(readLen)(0); i < readLen; ++i, ++ri, ++bi) {
                if (hasRefOffset[i] && bi >= 0) {
                    assert(roi < refOffsetColData.second);
                    bi = refOffsetColData.first[roi++];
                    ri += bi;
                }
                auto base = 'N';
                if (!hasMismatch[i]) {
                    auto const offset = loadAbsolute(globalStart + ri);
                    if (offset < read.second)
                        base = read.first[offset];
                    else if (current->circular) {
                        auto const zero = current->rowRange.first * maxSeqLen;
                        auto const wrapped = (globalStart + ri - zero) % current->length;
                        auto const offset = loadAbsolute(wrapped + zero);
                        assert(offset < read.second);
                        base = read.first[offset];
                    }
                }
                else {
                    assert(mi < mismatchColData.second);
                    base = mismatchColData.first[mi++];
                }

                // Reference is always 5', but read may be either strand
                if (reversed) {
                    assert(r != result.rend());
                    switch (base) {
                    case 'A': *r++ = 'T'; break;
                    case 'C': *r++ = 'G'; break;
                    case 'G': *r++ = 'C'; break;
                    case 'T': *r++ = 'A'; break;
                    default:
                        *r++ = 'N';
                        break;
                    }
                }
                else {
                    assert(f != result.end());
                    *f++ = base;
                }
            }
            return result;
        }
    };

    /**
     * A cursor on the database's Primary Alignment table.
     *
     * @see ncbi-vdb/libs/axf/align-restore-read.c
     *
     * @Note The data here is also available from the table's
     * `RAW_READ` column, but that was causing heavy memory use.
     * So I switched to reconstructing the value here from the base columns.
     */
    class AlignReader final : public VdbBaseReader {
#define FULL_IMPLEMENTATION 1
#if FULL_IMPLEMENTATION
        using HasRefOffsetCol = VdbBaseReader::Column<uint8_t>;
        using RefOffsetCol = VdbBaseReader::Column<int32_t>;
        using HasMismatchCol = VdbBaseReader::Column<uint8_t>;
        using MismatchCol = VdbBaseReader::Column<uint8_t>;
        using GlobalRefStartCol = VdbBaseReader::Column<uint64_t>;
        using ReverseStrandCol = VdbBaseReader::Column<bool>;
#else
        using ReadCol = VdbBaseReader::Column<char>;
#endif
        using SpotIdCol = VdbBaseReader::Column<int64_t>;

#if FULL_IMPLEMENTATION
        ReferenceReader ref;
        HasRefOffsetCol hasRefOffsetCol;
        RefOffsetCol refOffsetCol;
        HasMismatchCol hasMismatchCol;
        MismatchCol mismatchCol;
        GlobalRefStartCol globalRefStartCol;
        ReverseStrandCol reverseStrandCol;
#else
        ReadCol readCol;
#endif
        SpotIdCol spotIdCol;
    public:
        AlignReader(VDatabase const *const db)
        : VdbBaseReader(db, "PRIMARY_ALIGNMENT")
#if FULL_IMPLEMENTATION
        , ref(db)
        , hasRefOffsetCol(addColumn<HasRefOffsetCol>("(bool)HAS_REF_OFFSET"))
        , refOffsetCol(addColumn<RefOffsetCol>("REF_OFFSET"))
        , hasMismatchCol(addColumn<HasMismatchCol>("(bool)HAS_MISMATCH"))
        , mismatchCol(addColumn<MismatchCol>("MISMATCH"))
        , globalRefStartCol(addColumn<GlobalRefStartCol>("GLOBAL_REF_START"))
        , reverseStrandCol(addColumn<ReverseStrandCol>("REF_ORIENTATION"))
#else
        , readCol(addColumn<ReadCol>("RAW_READ"))
#endif
        , spotIdCol(addColumn<SpotIdCol>("SEQ_SPOT_ID"))
        {
            open();
#if FULL_IMPLEMENTATION
            getRowRangeFrom(hasMismatchCol);
#else
            getRowRangeFrom(readCol);
#endif
        }
        bool read(Fragment *output) {
#if FULL_IMPLEMENTATION
            auto const hasMismatch = columnDataNoThrow(hasMismatchCol);
            if (hasMismatch.first) {
                auto const spotId = valueOf(spotIdCol);
                auto const hasRefOffset = dataOf(hasRefOffsetCol);
                auto const mismatch = dataOf(mismatchCol);
                auto const refOffset = dataOf(refOffsetCol);
                auto const globalRefStart = valueOf(globalRefStartCol);
                auto const reversed = valueOf(reverseStrandCol);

                output->bases = ref.restoreRead(hasMismatch, mismatch
                                                , hasRefOffset, refOffset
                                                , globalRefStart, reversed);
                output->spotid = std::to_string(spotId);
                ++row;
                return true;
            }
#else
            auto const read = columnDataNoThrow(readCol);
            if (read.first) {
                auto const spotId = valueOf(spotIdCol);

                output->bases = std::string(read.first, read.second);
                output->spotid = std::to_string(spotId);
                ++row;
                return true;
            }
#endif
            return false;
        }
    };

    /**
     * A cursor on the Sequence table.
     *
     * This is for getting the reads of
     * either an unaligned run
     * or the unaligned reads of an
     * aligned run.
     */
    class SeqReader final : public VdbBaseReader {
        using AlignIdCol = VdbBaseReader::Column<int64_t>;
        using ReadCol = VdbBaseReader::Column<char>;
        using ReadStartCol = VdbBaseReader::Column<uint32_t>;
        using ReadLenCol = VdbBaseReader::Column<uint32_t>;
        using ReadTypeCol = VdbBaseReader::Column<uint8_t>;
        using ReadFilterCol = VdbBaseReader::Column<uint8_t>;

        AlignIdCol alignIdCol;
        ReadCol readCol;
        ReadStartCol readStartCol;
        ReadLenCol readLenCol;
        ReadTypeCol readTypeCol;
        ReadFilterCol readFilterCol;

        ReadTypeCol::Ptr readType;
        ReadFilterCol::Ptr readFilter;
        ReadLenCol::Ptr readLen;
        ReadCol::Ptr bases;
        AlignIdCol::Ptr alignId;
        ReadStartCol::Ptr readStart;
        ReadStartCol::Value sReadStart[2];
        ReadStartCol::Value *hReadStart;

        uint64_t spotCount;
        uint64_t readCount;
        uint64_t alignedReadCount;

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

        /**
         * This is the first column we will load in any row,
         * If we can't read it, it is probably because we are
         * at the end of the table.
         */
        bool loadReadType() {
            auto const pc = columnDataNoThrow(readTypeCol);
            readType = pc.first;
            numreads = pc.second;
            return readType != nullptr;
        }

        /**
         * Create a usable `READ_START`
         *
         * Since we are using `CMP_READ` instead of `READ`,
         * and `READ_START` is about `READ`,
         * we need to create a `READ_START`
         * that divides up `CMP_READ`,
         * which contains the unaligned bases.
         */
        void computeReadStart() {
            if (max_numreads < numreads) {
                /* Generally there are no more than 2 reads and
                 * this code will never get triggered
                 */
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
        void loadBases() {
            auto const pc = dataOf(readCol);
            bases = pc.first;
        }
        void loadFilters() {
            if (useReadFilter)
                readFilter = dataOf(readFilterCol, numreads);
            if (aligned)
                alignId = dataOf(alignIdCol, numreads);
        }
        void addCommonColumns() {
            readTypeCol = addColumn<ReadTypeCol>("READ_TYPE");
            readLenCol = addColumn<ReadLenCol>("READ_LEN");
            if (useReadFilter)
                readFilterCol = addColumn<ReadFilterCol>("READ_FILTER");
        }
        void addAlignedColumns() {
            readCol = addColumn<ReadCol>("CMP_READ"); // notice, it is not the same name as the unaligned case
            alignIdCol = addColumn<AlignIdCol>("PRIMARY_ALIGNMENT_ID");
        }
        void addUnalignedColumns() {
            readCol = addColumn<ReadCol>("READ"); // notice, it is not the same name as the aligned case
            readStartCol = addColumn<ReadStartCol>("READ_START");
        }
        void addColumns() {
            addCommonColumns();
            if (aligned)
                addAlignedColumns();
            else
                addUnalignedColumns();
        }
        int64_t firstRowWithUnalignedRead() const {
            // NB. best effort using metadata, i.e. no table scan
            KMetadata const *meta = metadata();
            if (meta) {
                int64_t h = 0;
                int64_t u = 0;

                KMDataNode const *nodeH = nullptr;
                KMetadataOpenNodeRead(meta, &nodeH, "unaligned/first-half-aligned");
                KMDataNodeReadAsI64(nodeH, &h);
                KMDataNodeRelease(nodeH);

                KMDataNode const *nodeU = nullptr;
                KMetadataOpenNodeRead(meta, &nodeH, "unaligned/first-unaligned");
                KMDataNodeReadAsI64(nodeU, &u);
                KMDataNodeRelease(nodeU);

                KMetadataRelease(meta);

                if (nodeH && nodeU)
                    return std::min(h, u) - rowRange.first;
                if (nodeH)
                    return h - rowRange.first;
                if (nodeU)
                    return u - rowRange.first;
            }
            return 0;
        }
    public:
        FastVdbReader_Counts getCounts() const {
            auto result = FastVdbReader_Counts();
            result.spots = spotCount;
            result.reads = readCount;
            result.aligned = alignedReadCount;
            return result;
        }

        /**
         * Opens a cursor on the Sequence table of an aligned database.
         */
        SeqReader(VDatabase const *const db, bool is_Aligned, bool unalignedOnly, bool filtered = false)
        : VdbBaseReader(db, "SEQUENCE")
        , hReadStart(nullptr)
        , spotCount(0)
        , readCount(0)
        , alignedReadCount(0)
        , max_numreads(2)
        , aligned(is_Aligned)
        , useReadFilter(filtered)
        {
            addColumns();
            open();
            getRowRangeFrom(readCol);
            if (aligned && unalignedOnly)
                row = firstRowWithUnalignedRead();

            readNo = 0;
        }

        /**
         * Opens a cursor on an unaligned table.
         */
        SeqReader(VTable const *const tbl, bool filtered = false)
        : VdbBaseReader(tbl, "<implied>")
        , hReadStart(nullptr)
        , spotCount(0)
        , readCount(0)
        , alignedReadCount(0)
        , max_numreads(2)
        , aligned(false)
        , useReadFilter(filtered)
        {
            addColumns();
            open();
            getRowRangeFrom(readCol);
            readNo = 0;
        }
        ~SeqReader() {
            if (hReadStart)
                delete [] hReadStart;
        }
        float progress(AlignReader *other) const {
            auto const pctRows = VdbBaseReader::progress();
            if (pctRows < 0)
                return 1.0;

            if (!aligned || readCount == 0 || alignedReadCount == 0)
                return pctRows;

            auto const pctAligned = double(alignedReadCount) / readCount;
            return pctRows * (1.0 - pctAligned) + other->progress() * pctAligned;
        }

        bool read(Fragment *output, AlignReader *other) {
            while (row < rowRange.second) {
                if (readNo == 0) {
                    if (!loadReadType()) {
                        LOG("SeqReader stopped early at row " << row);
                        return false;
                    }

                    loadFilters();
                    readLen = dataOf(readLenCol, numreads);
                    if (aligned)
                        computeReadStart();
                    else
                        readStart = dataOf(readStartCol, numreads);
                    loadBases();
                    ++spotCount;
                }
                auto const thisRead = readNo;
                auto const thisRow = row + rowRange.first;

                if (++readNo == numreads) {
                    readNo = 0;
                    ++row;
                }
                ++readCount;
                if (isAligned(thisRead)) {
                    ++alignedReadCount;
                    continue;
                }
                if ((readType[thisRead] & 1) != 1 || isFiltered(thisRead))
                    continue;

                output->bases = std::string(bases + readStart[thisRead], readLen[thisRead]);
                output->spotid = std::to_string(thisRow);
                return true;
            }
            if (other)
                return other->read(output);
            return false;
        }
    };
    using Cache = std::map<std::string, FastVdbReader_Counts>;
    static std::pair<bool, Cache::const_iterator> statsFor(std::string const &accession, SeqReader const &seq) {
        static auto cache = Cache();
        auto const fnd = cache.find(accession);
        if (fnd != cache.end())
            return {false, fnd};

        auto const inserted = cache.insert({accession, seq.getCounts()});
        return {true, inserted.first};
    }
    std::string accession;
    AlignReader *alg;
    SeqReader *seq;
    bool unalignedOnly;
public:
    /**
     * Gets aligned and unaligned reads;
     * never accesses quality scores
     */
    FastVdbReader(std::string const &acc, bool unaligned_only)
    : accession(acc)
    , alg(nullptr)
    , seq(nullptr)
    , unalignedOnly(unaligned_only)
    {
        VDBManager const *mgr = nullptr;
        VDatabase const *db = nullptr;
        rc_t rc = 0;

        rc = VDBManagerMakeRead(&mgr, nullptr);
        if (rc != 0)
            throw std::runtime_error("Can't make a VDB Manager!");

        /* it is okay for this to fail, e.g. it is not aligned */
        rc = VDBManagerOpenDBRead(mgr, &db, nullptr, "%s", acc.c_str());
        if (rc == 0) {
            VDBManagerRelease(mgr);

            try {
                /* it is okay for this to fail */
                alg = new AlignReader(db);
            }
            catch (std::exception const &e) {
                alg = nullptr;
            }

            try {
                /* it is not okay for this to fail */
                seq = new SeqReader(db, alg != nullptr, unalignedOnly);
            }
            catch (...) {
                throw std::runtime_error(std::string("Can't read sequences from ") + acc);
            }
            VDatabaseRelease(db);
        }
        else {
            VTable const *seqtbl = nullptr;

            rc = VDBManagerOpenTableRead(mgr, &seqtbl, nullptr, "%s", acc.c_str());
            VDBManagerRelease(mgr);
            if (rc != 0)
                throw std::runtime_error(std::string("Can't open ") + acc);
            seq = new SeqReader(seqtbl);
            VTableRelease(seqtbl);
        }
    }
    ~FastVdbReader() {
        delete seq;
        if (alg)
            delete alg;
    }

    SourceStats stats() const override
    {
        auto const inserted = statsFor(accession, *seq);
        auto const &counts = inserted.second->second;
        auto const rps = counts.reads / counts.spots;

        if (unalignedOnly) {
            auto unaligned = counts.reads - counts.aligned;
            return SourceStats(size_t(0.5 + unaligned / rps), int(0.5 + rps));
        }
        else {
            return SourceStats(counts.spots, int(0.5 + rps));
        }
    }

    float progress() const override
    {
        return seq->progress(alg);
    }

    bool read(Fragment* output) override
    {
        return seq->read(output, alg);
    }
};
