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

#define NO_BASES 0

class FastVdbReader final: public Reader {
    class VdbBaseReader {
    protected:
        using RowId = int64_t;
        using RowCount = uint64_t;
        using RowRange = std::pair<RowId, RowCount>;

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
            using Count = uint32_t;
            using Data = std::pair<Ptr, Count>;

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
            Ptr data(Count const requiredCount, VCursor const *curs, RowId const row, std::string const &sourceName) const noexcept(false)
            {
                auto const result = data(curs, row, sourceName);
                if (result.second == requiredCount)
                    return result.first;
                throwReadError(row, sourceName);
            }
            Value value(VCursor const *const curs, RowId const row, std::string const &sourceName) const noexcept(false)
            {
                return *(data(1, curs, row, sourceName));
            }

            RowRange rowRange(VCursor const *const curs, std::string const &sourceName) const {
                auto first = RowId(0);
                auto count = RowCount(0);
                assert(cid != nilColumnId());
                if (0 == VCursorIdRange(curs, cid, &first, &count))
                    return {first, count};
                throw std::runtime_error(std::string("Can't get row range of ") + sourceName);
            }
        };
        std::string const sourceName;
        VCursor const *curs;
        RowId row;
        RowRange rowRange;

        /* boilerplate code for accessing VDB */
        template <typename Col, typename Data = typename Col::Data>
        Data columnDataNoThrow(Col const &col) const noexcept(true)
        {
            return col.dataNoThrow(curs, row + rowRange.first);
        }

        template <typename Col, typename Data = typename Col::Data>
        Data dataOf(Col const &col) const noexcept(false)
        {
            return col.data(curs, row + rowRange.first, sourceName);
        }

        template <typename Col, typename Ptr = typename Col::Ptr, typename Count = typename Col::Count>
        Ptr dataOf(Col const &col, Count const requiredCount) const noexcept(false)
        {
            return col.data(requiredCount, curs, row + rowRange.first, sourceName);
        }

        template <typename Col, typename Value = typename Col::Value>
        Value valueOf(Col const &col)
        {
            return col.value(curs, row + rowRange.first, sourceName);
        }

        template <typename T>
        T addColumn(char const *const col_expr) {
            return T(curs, col_expr, sourceName);
        }

        void open() {
            if (0 != VCursorOpen(curs))
                throw std::runtime_error(std::string("Can't open cursor on ") + sourceName);
        }

        template <typename T>
        void getRowRangeFrom(T const &col) {
            rowRange = col.rowRange(curs, sourceName);
        }

        VdbBaseReader(VTable const *tbl, char const *name)
        : sourceName(name)
        , curs(nullptr)
        , row(0)
        , rowRange({0, 0})
        {
            if (0 != VTableCreateCursorRead(tbl, &curs))
                throw std::runtime_error(std::string("Can't make a cursor on ") + sourceName);
        }
        virtual ~VdbBaseReader() {
            VCursorRelease(curs);
        }
    public:
        virtual bool read(Fragment *, VdbBaseReader *other) = 0;
    };
    class AlignReader final : public VdbBaseReader {
        using ReadCol = VdbBaseReader::Column<char>;
        using SpotIdCol = VdbBaseReader::Column<int64_t>;

        ReadCol readCol;
        SpotIdCol spotIdCol;
    public:
        AlignReader(VTable const *tbl)
        : VdbBaseReader(tbl, "PRIMARY_ALIGNMENT")
        , readCol(addColumn<ReadCol>("RAW_READ"))
        , spotIdCol(addColumn<SpotIdCol>("SEQ_SPOT_ID"))
        {
            open();
            getRowRangeFrom(readCol);
        }
        bool read(Fragment *output, VdbBaseReader *other) {
            auto const pc = columnDataNoThrow(readCol);
            if (pc.first) {
                auto const spotId = valueOf(spotIdCol);

                ++row;
                output->bases = std::string(pc.first, pc.second);
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

        bool loadReadType() {
            /* this is the first column we will load in any row,
             * if we can't read it, we won't be able to read the
             * rest of the columns.
             */
            auto const pc = columnDataNoThrow(readTypeCol);
            readType = pc.first;
            numreads = pc.second;
            return readType != nullptr;
        }
        void computeReadStart() {
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
        int64_t firstRowWithUnalignedRead(VTable const *tbl) const {
            // NB. best effort using metadata, i.e. no table scan
            KMetadata const *meta;
            auto rc = VTableOpenMetadataRead(tbl, &meta);
            if (rc == 0) {
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
                    return h - rowRange.first;
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
        SeqReader(VTable const *tbl, bool is_Aligned, bool unalignedOnly, bool filtered = false)
        : VdbBaseReader(tbl, "SEQUENCE")
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
                row = firstRowWithUnalignedRead(tbl);

            readNo = 0;
        }
        ~SeqReader() {
            if (hReadStart)
                delete [] hReadStart;
        }
        float progress() const {
            auto const total = rowRange.second;
            return total ? double(row) / total : 1.0;
        }

        bool read(Fragment *output, VdbBaseReader *other) {
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
                    if (other)
                        return other->read(output, nullptr);
                    continue;
                }
                if ((readType[thisRead] & 1) != 1 || isFiltered(thisRead))
                    continue;

                output->bases = std::string(bases + readStart[thisRead], readLen[thisRead]);
                output->spotid = std::to_string(thisRow);
                return true;
            }
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
    /* Gets aligned and unaligned reads
     * Never accesses quality scores
     */
    FastVdbReader(std::string const &acc, bool unaligned_only)
    : accession(acc)
    , alg(nullptr)
    , seq(nullptr)
    , unalignedOnly(unaligned_only)
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

            /* it is okay for this to fail */
            (void)VDatabaseOpenTableRead(db, &algtbl, "PRIMARY_ALIGNMENT");

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
        seq = new SeqReader(seqtbl, algtbl != NULL, unaligned_only);
        if (algtbl && !unaligned_only)
            alg = new AlignReader(algtbl);

        VTableRelease(algtbl);
        VTableRelease(seqtbl);
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
        return seq->progress();
    }

    bool read(Fragment* output) override
    {
        if (seq->read(output, alg)) {
#if NO_BASES
            output->bases = std::string();
#endif
            return true;
        }
        return false;
    }
};
