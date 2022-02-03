#ifndef __TAX_COLLATOR_HPP__
#define __TAX_COLLATOR_HPP__
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


#include "bm/bm.h"
#include "bm/bmstrsparsevec.h"
#include "bm/bmsparsevec_compr.h"
#include "bm/bmsparsevec_algo.h"
#include "bm/bmtimer.h"
#include "bm/bmsparsevec_serial.h"
#include "bm/bmdbg.h"
#include "bm/bmtimer.h"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>
#include "spdlog/stopwatch.h"
#include "taskflow/taskflow.hpp"
#include <taskflow/algorithm/sort.hpp>
 
#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <list>
#include <algorithm>
#include <sstream>

using namespace std;
using namespace std::chrono;

#define BEGIN_TC_NAMESPACE namespace tc {
#define END_TC_NAMESPACE }

BEGIN_TC_NAMESPACE




/**
 * @brief options to control execution
 * 
 */
static constexpr bool COMPACT_OPT = true;
static constexpr bool COLLATE_OPT = false;
static constexpr bool INCLUDE_COUNTS_OPT = true;
static constexpr bool EXCLUDE_COUNTS_OPT = false;

template <bool compact = COLLATE_OPT, bool include_counts = INCLUDE_COUNTS_OPT>
struct tax_hits_options
{
    static constexpr bool is_compact() { return compact; }          ///< Flag for compact mode: counts of unique tax_id sets
    static constexpr bool has_counts() { return include_counts; }   ///< Flag for presence of counts
};

/**
 * Spot data structure represents STATS analysis data: 
 * spot_name and associated list of tax_id with hits counters for each tax_id
 * Example: 
 * 12801   2759x3  9526    131567x2    207598  314293x4
 * 
 */
template<class Options = tax_hits_options<>>
struct Spot {
    string name;                ///< Spot name
    vector<uint32_t> tax_id;    ///< List pof tax_id 
    vector<uint32_t> counts;    ///< List of hits counter per tax_id (counts.size() == tax_id.size())
    
    /**
     * Initializes spot from vector of fields parsed out from STAT analysis output
     * 12801   2759x3  9526    131567x2    207598  314293x4
     * 
     * @param fields, first field is expected to be spot_name followed by tax_id with optional count) (9606 or 9606x3)
     * 
     */
    void init(vector<string>& fields);
    
    /**
     * @brief Copies tax_id + counts from another Spot
     * 
     * @param spot 
     */
    void add_taxa(const Spot& spot); 

    /**
     * @brief Merge with another Spot by copying tax_id+counts and perfroms normalize
     * 
     * @param spot 
     */
    void merge(const Spot& spot); 

    /**
     * @brief Makes sure tax_id+counters are sorted and unique
     * 
     */
    void normalize(); 
};


typedef bm::bvector<> bvector_type;
typedef bm::str_sparse_vector<char, bvector_type, 3> str_sv_type;

typedef bm::sparse_vector<uint32_t, bvector_type> sparse_vector_u32;
typedef bm::rsc_sparse_vector<uint32_t, sparse_vector_u32>  rsc_sparse_vector_u32;

/**
 * @brief Matrix of uint32 values where each columns is represented by rsc_sparse_vector
 * Used to store tax_id and counters
 * 
 */
struct U32_rsc_matrix 
{

    string name;              ///< Matrix name (for info and error reporting)
    uint32_t num_rows = 0;  
    uint32_t num_cols = 0;
    uint32_t curr_col = 0;    ///< Current column when matrix is being populated 
    uint32_t num_values = 0;  ///< Number of non-null values (telemetry)

    using vector_type = rsc_sparse_vector_u32;
    using size_type = vector_type::size_type;
    using value_type = vector_type::value_type;

    vector<unique_ptr<vector_type>> data;                           ///< list of columns 
    vector<unique_ptr<vector_type::back_insert_iterator>> data_bi;  ///< insert iterator valid when matrix is being populated

    /**
     * @brief Construct a new u32 rsc matrix object
     * 
     * @param _name 
     */
    U32_rsc_matrix(const string& _name);

    /**
     * @brief Serialization
     * 
     * @param file_name 
     */
    void save(const string& file_name); 

    /**
     * @brief De-serialization
     * 
     * @param file_name 
     */
    void load(const string& file_name); 

    /**
     * @brief Adds new columns 
     * 
     * @param c - number of columns to add
     */
    void add_column(int c = 1); 

    /**
     * @brief Adds new value and increments curr_column 
     * 
     * @param value 
     */
    void add_value(uint32_t value); 

    /**
     * @brief Add new null-value and increments curr_columns
     * 
     */
    void add_null(); 

    /**
     * @brief Complets the current row by 
     * padding the remaining columns with null values and rests curr_column
     * 
     */
    void end_row(); 

    /**
     * @brief Finalize insertions by flashing insert iterator and syncing rs cvector 
     * After finalize  insertions operation are not valid 
     * 
     */
    void finalize(); 

    /**
     * @brief Optmize data memory and optionally reports memory usage
     * 
     */
    void optimize(bool print_stats = true);

    /**
     * @brief Clears all columns
     * 
     */
    void clear();

    /**
     * @brief Get indexes of all rows that have values only in columns that <= col_index
     * 
     * @param col_index 
     * @param index vector row indexes/number
     */
    void get_rows(int col_index, vector<uint32_t>& index);

    /**
     * @brief Get values for row_index up to param cols size
     * 
     * @param row_index index of the row to get values for
     * @param cols how many cols to retrieve 
     */
    void get_row(uint32_t row_index, vector<uint32_t>& cols);
};

/**
 * @brief Main TaxCollator structure
 * holds list of spot name and matrices of tax_id and optional counters
 * provides collate(), group(), print() methods
 * 
 * 
 * @tparam Options 
 */
template<class Options = tax_hits_options<>>
struct Tax_hits
{
    unique_ptr<str_sv_type> spot_names;    ///< spot names list 
    unique_ptr<str_sv_type::back_insert_iterator> spot_names_bi; ///< spot names insert iterator 
 
    U32_rsc_matrix tax_ids{"Tax_id"};      ///< tax_id matrix
    U32_rsc_matrix counts{"Counts"};       ///< tax_id's counts matrix
     
    atomic<size_t> num_merges {0};         ///< Number of occurred merges (telemetry)
    //size_t num_threads{16};                ///< Max number of allowed threads 

    mutex m_output_mutex;                  ///< protects output

    vector<uint32_t> spot_index;           ///< temporary index used by collate (sort and merge)

    /**
     * @brief Construct a new Tax_hits object
     * 
     * @param _num_threads 
     * @param use_null - flag to indicate if spot_names cane be nullable (and thus pruned more efficiently)
     */
    Tax_hits(bool use_null = false);

    /**
     * @brief Adds new spot stats
     * 
     * @param spot 
     */
    void add_row(const Spot<Options>& spot);

    /**
     * @brief Finalize insertions by flashing spot name insert iterator and finalizing rsc_matrices 
     * After finalize  insertions operation are not valid 
     * 
     */
    void finalize();

    /**
     * @brief The succinct structures are optimized (re-compressed) and memory usage is reported (optionally)
     * 
     */
    void optimize(bool print_stats = true); 

    /**
     * @brief Serialize the structure
     * uses file_prefix + '.names|.taxa|.counts' filename
     * 
     * @param file_prefix 
     */
    void save(const string& file_prefix);

    /**
     * @brief Used by load to initialize self using previously saved vectors
     * 
     * @param file_prefix 
     * @return true 
     * @return false 
     */
    bool init(const string& file_prefix); 

    /**
     * @brief De-serializes from filename
     * if current directory has files with filename prefix and '.names' and '.taxa' and '.counts'
     * load will assume that those are serialized vectors from the previous run and will try to use them
     * instead of parsing the filename (saves time while debugging)
     * 
     * @param filename 
     */
    void load(const string& filename);

    /**
     * Clears all data
     * 
     */
    void clear();

    /**
     * @brief Print Tax_hits data into a stream
     * 
     * @param executor 
     * @param os 
     * @param count - number of rows ro print (-1 prints all)
     */
    void print(tf::Executor& executor, ostream& os, size_t count = -1); 
    
    /**
     * @brief Get the spot object for a given row 
     * 
     * @tparam SpotOptions 
     * @param row_index 
     * @param spot 
     * @param get_name flag to indicate if spot names has to be retrieved as well
     */
    template<class SpotOptions = Options>
    void get_spot(uint32_t row_index, Spot<SpotOptions>& spot, bool get_name = true);

    /**
     * @brief sorts spot spots names using spot_index 
     * upon completion spot_index will be sorted by spot_names
     * 
     * @param executor 
     */
    void sort(tf::Executor& executor);

    /**
     * @brief used by collate() method. Merges sorted (via spot_index) spots and saves them in the new Tax_hits structure
     * peridocally prunes itself using set_null method
     * clears spot_index upon completion
     *
     * @tparam CollateOptions 
     * @param executor 
     * @param tax_hits 
     */
    template<class CollateOptions>
    void merge(tf::Executor& executor, Tax_hits<CollateOptions>& tax_hits); 

    /**
     * @brief Self pruning using bit mask 
     * All rows with set bit will be pruned (spot_names nullified and matrices zeroed)
     * The succinct structure are optimized (re-compressed) after pruning 
     * 
     * @param executor 
     * @param bv 
     */
    void set_null(tf::Executor& executor, const bvector_type& bv); 

    /**
     * @brief Implements spot collation (sort and merge) 
     * empties itself upon completeion 
     * 
     * @tparam CollateOptions new collated TaxHits Options
     * @param executor
     * @return unique_ptr<Tax_hits<CollateOptions>> new collated Tax_hits
     */
    template<class CollateOptions>
    unique_ptr<Tax_hits<CollateOptions>> collate(tf::Executor& executor); 

    /**
     * Used by group()
     * Count the same rows and prints the count and group 
     * 
     * @param rows 
     * @param num_columns number of columns in the group
     */
    void group_columns(const vector<uint32_t>& rows, int num_columns, ostream& os);

    /**
     * Used by group()
     * Count the same rows and prints the count and group 
     * Faster alternative to group_columns()
     * 
     * @param rows 
     * @param num_columns number of columns in the group
     */
    void group_columns_bulk(const vector<uint32_t>& rows, int num_cols, ostream& os);

    /**
     * Implements compact mode: counting of unique tax_id sets
     * Uses tax_id matrix to identify rows for each set of columns 
     * (set 1- column 1; set 2 - column 1, 2; set 3 - column 1,2,3; etc) 
     * sorts each set and iterates through the results to count the same rows and
     * print compact mode output 
     * For more details See stat-compact-compute.pdf
     * 
     */
    void group(tf::Executor& executor, ostream& os); 

    /**
     * @brief Temporary buffer structure to perform asynchronous merge
     * 
     * @tparam SpotOptions 
     */
    template<class SpotOptions = Options>
    struct Merge_data 
    {
        Merge_data(const string& spot_name, uint32_t i) 
            : idx(1, i)
        {
            spot.name = spot_name;
            assert(spot_name.empty() == false);
        }
        vector<uint32_t> idx;
        Spot<SpotOptions> spot;
    };

    /**
     * @brief Run asynchronous that marges collated spots and writes them to the new Tax_hits structure
     * 
     * @tparam CollateOptions 
     * @param executor 
     * @param merge_data 
     * @param tax_hits 
     * @return tf::Future<void> 
     */
    template<class CollateOptions>
    tf::Future<void> run_merge(tf::Executor& executor, vector<Merge_data<CollateOptions>>& merge_data, Tax_hits<CollateOptions>& tax_hits); 

};

/* -----------------------------------------------------
 * Spot method implementation 
 * ----------------------------------------------------- */

template<class Options>
void Spot<Options>::init(vector<string>& fields) 
{
    assert(!fields.empty());
    name = move(fields[0]);
    tax_id.clear();
    if constexpr (Options::has_counts()) 
        counts.clear();
    for (int i = 1; i < fields.size(); ++i) {
        string& str{fields[i]};
        auto pos = str.find("x");
        try {
            if (pos == string::npos) {
                tax_id.push_back(stoi(str));
                if constexpr (Options::has_counts())
                    counts.push_back(1);
            } else {    
                tax_id.push_back(stoi(str.substr(0, pos)));
                if constexpr (Options::has_counts())
                    counts.push_back(stoi(str.substr(pos + 1)));
            }

        } catch (exception& e) {
            spdlog::error("{}", e.what());
        }
    }
}

template<class Options>
void Spot<Options>::add_taxa(const Spot<Options>& spot) 
{
    tax_id.insert(tax_id.end(), spot.tax_id.begin(), spot.tax_id.end());
    if constexpr (Options::has_counts())
        counts.insert(counts.end(), spot.counts.begin(), spot.counts.end());
}

template<class Options>
void Spot<Options>::normalize() 
{
    if constexpr (Options::has_counts()) {
        static thread_local vector<uint32_t> new_tax_id;
        static thread_local vector<uint32_t> new_counts;
        static thread_local vector<int> index;
        new_tax_id.clear();
        new_counts.clear();
        index.resize(tax_id.size());
        generate(index.begin(), index.end(), [n = 0] () mutable { return n++; });
        sort(index.begin(), index.end(), [&](const auto& l, const auto& r) {
            return tax_id[l] < tax_id[r];
        });
        int i = 0;
        int idx = index[i];
        new_tax_id.push_back(move(tax_id[idx]));
        new_counts.push_back(move(counts[idx]));
        for (i = 1; i < index.size(); ++i) {
            auto idx = index[i];
            if (tax_id[idx] == new_tax_id.back()) {
                new_counts.back() += counts[idx];
            } else {
                new_tax_id.push_back(move(tax_id[idx]));
                new_counts.push_back(move(counts[index[i]]));
            }
        }
        swap(tax_id, new_tax_id);
        swap(counts, new_counts);
    } else {
        sort( tax_id.begin(), tax_id.end() );
        tax_id.erase( unique( tax_id.begin(), tax_id.end() ), tax_id.end() );
    } 
}

template<class Options>
void Spot<Options>::merge(const Spot<Options>& spot) 
{
    add_taxa(spot);
    normalize();
}

/* -----------------------------------------------------
 * U32_rsc_matrix methods implementation 
 * ----------------------------------------------------- */

BM_DECLARE_TEMP_BLOCK(TB)

U32_rsc_matrix::U32_rsc_matrix(const string& _name) : name(_name) 
{
}

void U32_rsc_matrix::save(const string& file_name) 
{
    ofstream ofs(file_name.c_str(), ofstream::out | ofstream::binary);
    ofs.write((char*)&num_cols, sizeof(num_cols));
    bm::sparse_vector_serializer<vector_type > serializer;
    bm::sparse_vector_serial_layout<vector_type > sv_lay;
    for (int i = 0; i < num_cols; ++i) {
        serializer.serialize(*data[i], sv_lay);
        size_t sz = sv_lay.size();
        ofs.write((char*)&sz, sizeof(sz));
        const unsigned char* buf = sv_lay.data();
        ofs.write((char*)&buf[0], sz);
    }
    ofs.close();
}

void U32_rsc_matrix::load(const string& file_name) 
{
    ifstream ifs(file_name.c_str(), ofstream::in | ofstream::binary);
    ifs.read((char*)&num_cols, sizeof(num_cols));
    spdlog::info("Loading {} columns", num_cols);
    vector<char> buffer;
    for (int i = 0; i < num_cols; ++i) {
        data.push_back(std::make_unique<vector_type>(bm::use_null));
        size_t sz = 0;
        ifs.read((char*)&sz, sizeof(sz));
        buffer.resize(sz);
        ifs.read((char*)&buffer[0], sz);
        bm::sparse_vector_deserialize(*data[i], (const unsigned char*)&buffer[0], TB);
        ++num_rows;
    }
    ifs.close();
}

void U32_rsc_matrix::add_column(int c) 
{
    while (c > 0) {
        data.push_back(std::make_unique<vector_type>(bm::use_null));
        data_bi.push_back(std::make_unique<vector_type::back_insert_iterator>(data.back().get()));
        if (num_rows > 0) 
            data_bi.back()->add_null(num_rows);
        ++num_cols;
        --c;
    }
}

void U32_rsc_matrix::add_value(uint32_t value) 
{
    assert(value > 0);
    data_bi[curr_col]->add(value);
    ++num_values;
    ++curr_col;
}

void U32_rsc_matrix::add_null() 
{
    data_bi[curr_col]->add_null();
    ++curr_col;
}

void U32_rsc_matrix::end_row() 
{
    for (auto i = curr_col; i < num_cols; ++i) {
        data_bi[i]->add_null();
    }
    ++num_rows;
    curr_col = 0;
}

void U32_rsc_matrix::finalize() 
{
    for(auto& v : data_bi) 
        v->flush();
    data_bi.clear();    
    for(auto& v : data) {
        v->sync();
        ///v->optimize(TB);
    }
}

void U32_rsc_matrix::optimize(bool print_stats) 
{
    size_t mem = 0;
    U32_rsc_matrix::vector_type::statistics st;
    for(auto& v : data) {
        v->optimize(TB, bvector_type::opt_compress, &st);
        mem += st.memory_used;
    }
    if (print_stats)
        spdlog::info("{} memory: {:L}", name, mem);
}

void U32_rsc_matrix::clear() 
{
    for (auto i = 0; i < num_cols; ++i) {
        data[i]->clear(true);
    }
    num_cols = num_rows = num_values = curr_col = 0;
}

void U32_rsc_matrix::get_rows(int col_index, vector<uint32_t>& index)
{
    index.clear();
    index.shrink_to_fit();
    bvector_type bv{*data[col_index]->get_null_bvector()};    
    int j = col_index + 1;
    if (j < num_cols)
        bv.bit_sub(*data[j]->get_null_bvector()); 
    for (auto en = bv.first(); en.valid(); ++en){
        index.push_back(*en);
    }
}


void U32_rsc_matrix::get_row(uint32_t row_index, vector<uint32_t>& cols)
{
    auto sz = cols.size();
    for (int i = 0; i < sz; ++i) {
        auto& v = *data[i];
        if (v.is_null(row_index))
            throw runtime_error("IS NULL!");
        cols[i] = v.get(row_index);
    }
}


/* -----------------------------------------------------
 * Tax_hits methods implementation 
 * ----------------------------------------------------- */

template<class Options>
Tax_hits<Options>::Tax_hits(bool use_null) 
{
    spot_names = std::make_unique<str_sv_type>(use_null ? bm::use_null : bm::no_null);
    if constexpr (Options::is_compact() == false)
        spot_names_bi = std::make_unique<str_sv_type::back_insert_iterator>(spot_names.get());        
}

template<class Options>
void Tax_hits<Options>::set_null(tf::Executor& executor, const bvector_type& bv) 
{

    tf::Taskflow taskflow;
    taskflow.emplace([&]() { spot_names->set_null(bv);});
    
    taskflow.for_each(tax_ids.data.begin(), tax_ids.data.end(), [&](unique_ptr<U32_rsc_matrix::vector_type>& data) { 
        try {
            data->clear(bv);
            data->sync(false);
            data->optimize();
        } catch(exception& e) {
            cerr << e.what() << endl;
        }
    });
    if constexpr (Options::has_counts()) {
        taskflow.for_each(counts.data.begin(), counts.data.end(), [&](unique_ptr<U32_rsc_matrix::vector_type>& data) { 
            try {
                data->clear(bv);
                data->sync(false);
                data->optimize();
            } catch(exception& e) {
                cerr << e.what() << endl;
            }
        });
    }
    
    executor.run(taskflow).wait();
}

template<class Options>
void Tax_hits<Options>::clear() 
{
    spot_names->clear_all(true);
    tax_ids.clear();
    if constexpr (Options::has_counts()) 
        counts.clear();
}

static 
bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

    
template<class Options>
void Tax_hits<Options>::finalize() 
{
    if constexpr (Options::is_compact() == false) {
        spot_names_bi->flush();
        spot_names_bi.reset(0);
        spot_names->remap();
        //spot_names->optimize(TB);
    }
    tax_ids.finalize();
    if constexpr (Options::has_counts())
        counts.finalize();
    optimize();    
}


template<class Options>
void Tax_hits<Options>::optimize(bool print_stats) 
{
    if constexpr (Options::is_compact() == false) {
        spot_names->optimize(TB);
        str_sv_type::statistics st;    
//        spot_names->optimize(TB, bvector_type::opt_compress, &st); // doesn't work, returns incorrect stat
        spot_names->calc_stat(&st);

        if (print_stats) {
            spdlog::info("Spot name size: {:L}", spot_names->size());    
            spdlog::info("Spot name memory: {:L}", st.memory_used);
        }
    }
    tax_ids.optimize(print_stats);
    if constexpr (Options::has_counts())
        counts.optimize(print_stats);
}


template<class Options>
bool Tax_hits<Options>::init(const string& file_prefix) 
{
    string spot_name_file = file_prefix + ".names";
    if (!file_exists(spot_name_file))
        return false;
    tf::Executor executor{4};
    tf::Taskflow taskflow;
    taskflow.emplace([this, spot_name_file]() { 
        spdlog::info("Loading '{}'", spot_name_file);
        file_load_svector(*spot_names, spot_name_file);
        str_sv_type::statistics st;
        spot_names->calc_stat(&st);
        spdlog::info("Spot name size: {:L}", spot_names->size());
        spdlog::info("Spot name memory: {:L}", st.memory_used);
    });
    string tax_file = file_prefix + ".taxa";
    if (!file_exists(tax_file))
        return false;
    taskflow.emplace([this, tax_file]() { 
        spdlog::info("Loading '{}'", tax_file);
        tax_ids.load(tax_file);
        U32_rsc_matrix::vector_type::statistics t_st;
        size_t tax_id_mem = 0;
        size_t tax_id_num = 0;
        for (const auto& v : tax_ids.data) {
            v->sync();
            v->calc_stat(&t_st);
            tax_id_mem += t_st.memory_used;
        }
        spdlog::info("Tax_id memory {:L}", tax_id_mem);

    });
    if constexpr (Options::has_counts()) {
        string hits_file = file_prefix + ".counts";
        if (!file_exists(hits_file))
            return false;
        taskflow.emplace([this, hits_file]() { 
            spdlog::info("Loading '{}'", hits_file);
            counts.load(hits_file);
            U32_rsc_matrix::vector_type::statistics t_st;
            size_t hits_mem = 0;
            size_t hits_num = 0;
            for (const auto& v : counts.data) {
                v->sync();
                v->calc_stat(&t_st);
                hits_mem += t_st.memory_used;
            }
            spdlog::info("Counters memory {:L}", hits_mem);
        });
    }
    executor.run(taskflow).wait();
    return true;
}

template<class Options>
void Tax_hits<Options>::add_row(const Spot<Options>& spot)
{
    if constexpr (Options::is_compact()  == false) {
        if (spot.name.empty())
            return;
        *spot_names_bi = spot.name;
    }
    int spot_sz = spot.tax_id.size();

    if (tax_ids.num_cols < spot_sz) {
        int n = spot_sz - tax_ids.num_cols;
        tax_ids.add_column(n);
        if constexpr (Options::has_counts())
            counts.add_column(n);
    }
    assert(tax_ids.num_cols >= spot_sz) ;
    for (int i = 0; i < spot_sz; ++i) {
        //all_taxa.insert(spot.tax_id[i]);
        tax_ids.add_value(spot.tax_id[i]);
        if constexpr (Options::has_counts()) {
            if (spot.counts[i] > 1) 
                counts.add_value(spot.counts[i]);
            else     
                counts.add_null();
        }
    }
    tax_ids.end_row();
    if constexpr (Options::has_counts()) 
        counts.end_row();
}

static 
void serialize_vec(const string& file_name, str_sv_type& vec)
{
    bm::sparse_vector_serializer<str_sv_type> serializer;
    bm::sparse_vector_serial_layout<str_sv_type> sv_lay;
    vec.optimize(TB);
    serializer.serialize(vec, sv_lay);
    ofstream ofs(file_name.c_str(), ofstream::out | ofstream::binary);
    const unsigned char* buf = sv_lay.data();
    ofs.write((char*)&buf[0], sv_lay.size());
    ofs.close();
}

template<class Options>
void Tax_hits<Options>::save(const string& file_prefix) 
{
    if constexpr (Options::is_compact() == false) {
        serialize_vec(file_prefix + ".names", *spot_names);
        str_sv_type::statistics st;
        spot_names->calc_stat(&st);
        spdlog::info("Spot name size: {:L}, spot name memory {:L}", spot_names->size(), st.memory_used);
    }

    tax_ids.save(file_prefix + ".taxa");
    if constexpr (Options::has_counts())
        counts.save(file_prefix + ".counts");
    /*
    ofstream os(file_prefix + ".tax.stats", ios::out);
    for_each(tax_ids.data.begin(), tax_ids.data.end(), [&](unique_ptr<U32_rsc_matrix::vector_type>& data) {
        bm::print_svector_stat(os, *data);
    });
    */


    {
        U32_rsc_matrix::vector_type::statistics t_st;
        size_t tax_id_mem = 0;
        size_t tax_id_num = 0;
        for (const auto& v : tax_ids.data) {
            v->calc_stat(&t_st);
            tax_id_mem += t_st.memory_used;
        }
        spdlog::info("number of tax_id {:L}, tax_id memory {:L}", tax_ids.num_values, tax_id_mem);
    }

    if constexpr (Options::has_counts()) {
        U32_rsc_matrix::vector_type::statistics t_st;
        size_t hits_mem = 0;
        size_t hits_num = 0;
        for (const auto& v : counts.data) {
            v->calc_stat(&t_st);
            hits_mem += t_st.memory_used;
        }
        spdlog::info("Number of counter {:L}, counter memory {:L}", counts.num_values, hits_mem);
    }
}

template<class Options>
void Tax_hits<Options>::print(tf::Executor& executor, ostream& os, size_t count) 
{
    int num_cols = tax_ids.num_cols;//size();

    tf::Taskflow taskflow;

    auto sz = count == - 1 ? spot_names->size() : count;
    spdlog::info("printing {:L}", sz);
    int page_size = 500000;
    auto num_pages = sz/page_size + 1;
    if (num_pages == 1)
        page_size = sz;
    vector<int> page_index;
    page_index.resize(num_pages, 0);
    generate(page_index.begin(), page_index.end(), [n = 0] () mutable { return n++; });
//    std::locale::global(std::locale("C")); // disable comma as thousand separator
    taskflow.for_each(page_index.begin(), page_index.end(), [&](int index) { 
        
        size_t start_pos = index * page_size;
        auto it = spot_names->get_const_iterator(start_pos);

        vector<unique_ptr<U32_rsc_matrix::vector_type::const_iterator>> tax_id_b;
        for (const auto& v : tax_ids.data) {
            //v->sync();
            tax_id_b.push_back(std::make_unique<U32_rsc_matrix::vector_type::const_iterator>(v.get(), start_pos));
        }
        vector<unique_ptr<U32_rsc_matrix::vector_type::const_iterator>> hits_b;
        if constexpr (Options::has_counts()) {
            for (const auto& v : counts.data) {
                //v->sync();
                hits_b.push_back(std::make_unique<U32_rsc_matrix::vector_type::const_iterator>(v.get(), start_pos));
            }
        }
        int c = 0;
        vector<string> buffer;
        std::stringstream ss;

        while (c < page_size && it.valid()) {
            if constexpr (Options::is_compact() == false) {
                if (it.is_null()) {
                    it.advance();
                    for (int i = 0; i < num_cols; ++i) {
                        tax_id_b[i].get()->advance();
                        if constexpr (Options::has_counts())
                            hits_b[i].get()->advance();
                    }
                    ++c;
                    continue;
                }
            }
            ss << it.value();
            int i = 0;
            for (; i < num_cols; ++i) {
                auto t_it = tax_id_b[i].get();
                if (t_it->is_null())
                    break;
                ss << '\t';
                assert(t_it->value() > 0);
                ss << t_it->value(); 
                t_it->advance();
                if constexpr (Options::has_counts()) {
                    auto hits_it = hits_b[i].get();
                    if (hits_it->is_null() == false) {
                        ss << 'x' << hits_it->value(); 
                    }
                    hits_it->advance();
                }
            }
            for (; i < num_cols; ++i) {
                tax_id_b[i].get()->advance();
                if constexpr (Options::has_counts())
                    hits_b[i].get()->advance();
            }
            {        
                buffer.push_back(move(ss.str()));
                ss.str({});
                if (buffer.size() == 5000) {
                    {
                        const lock_guard<std::mutex> lock(m_output_mutex);
                        copy(buffer.begin(), buffer.end(), ostream_iterator<string>(os, "\n"));
                    }
                    buffer.clear();
                }
            }
            it.advance();
            ++c;
        }
        if (!buffer.empty()) {
            const lock_guard<std::mutex> lock(m_output_mutex);
            copy(buffer.begin(), buffer.end(), ostream_iterator<string>(os, "\n"));
        }

    }); 
    executor.run(taskflow).wait();
    os.flush();
//    std::locale::global(std::locale("en_US.UTF-8")); // enable comma as thousand separator
    return;
}


template<class Options>
template<class SpotOptions>
void Tax_hits<Options>::get_spot(uint32_t index, Spot<SpotOptions>& spot, bool get_name)
{
    if (get_name)
        spot_names->get(index, spot.name);
    spot.tax_id.clear();
    if constexpr (SpotOptions::has_counts()) 
        spot.counts.clear();
    auto sz = tax_ids.data.size();
    uint32_t value = 0;
    for (int i = 0; i < sz; ++i) {
        auto& v = *tax_ids.data[i];
        if (!v.try_get_sync(index, value))
            break;
        assert(value > 0);    
        spot.tax_id.push_back(value);
        if constexpr (SpotOptions::has_counts()) {
            auto& h = *counts.data[i];
            if (!h.try_get_sync(index, value))
                value = 1;
            assert(value > 0);    
            spot.counts.push_back(value);
        }
    }
}

static
void s_getrow(const vector<vector<uint32_t>>& columns, int row_index, vector<uint32_t>& row_data) 
{
    auto sz = columns.size();
    for (int i = 0; i < sz; ++i) {
        row_data[i] = columns[i][row_index];
    }
}


static 
void s_bulk_get_rows(U32_rsc_matrix& m, const U32_rsc_matrix::size_type* idx, U32_rsc_matrix::size_type num_rows, int cardinality, vector<vector<U32_rsc_matrix::value_type>>& values)
{
    vector<U32_rsc_matrix::size_type> tmp_buf;
    tmp_buf.resize(num_rows);
    for (int i = 0; i < cardinality; ++i) {
        values[i].resize(num_rows);
        auto &d = *m.data[i];
        d.gather(values[i].data(), idx, tmp_buf.data(), num_rows, bm::BM_UNSORTED);
    }
}

template<class Options>
void Tax_hits<Options>::group_columns_bulk(const vector<uint32_t>& rows, int cardinality, ostream& os) 
{
    spdlog::stopwatch sw; 
    static const int cCACHE_SIZE = 100000;
    static const int cPRINT_BUF_SIZE = 10000;
    auto num_rows = rows.size();

    // column values 
    vector<vector<U32_rsc_matrix::value_type>> col_values;
    col_values.resize(cardinality);

    // output buffer
    vector<string> buffer;
    buffer.reserve(cPRINT_BUF_SIZE);

    // tmp buffer for gather
    vector<U32_rsc_matrix::size_type> tmp_buf;
    tmp_buf.resize(cCACHE_SIZE);

    vector<uint32_t> last_row(cardinality), next_row(cardinality);

    stringstream ss;
    int group_count = 0;
    int batch_size = 0;
    int row = 0; 
    while (row < num_rows) {
        batch_size = min<int>(cCACHE_SIZE, num_rows - row);
        if (batch_size <= 0)
            break;
        for (int i = 0; i < cardinality; ++i) {
            col_values[i].resize(batch_size);
            auto &d = *tax_ids.data[i];
            d.gather(col_values[i].data(), &rows[row], tmp_buf.data(), batch_size, bm::BM_UNSORTED);
        }
        int idx = 0;
        if (group_count == 0) {
            for (int i = 0; i < cardinality; ++i) 
                last_row[i] = col_values[i][0];
            group_count = 1;
            idx = 1;
        }
        //s_bulk_get_rows(tax_ids, (const U32_rsc_matrix::size_type*) &rows[row], num_rows, num_cols, col_values);
        //auto num_rows = col_values.front().size();
        for (; idx < batch_size; ++idx) {
            for (int i = 0; i < cardinality; ++i) 
                next_row[i] = col_values[i][idx];
            if (last_row == next_row) {
                ++group_count;
            } else {
                ss << group_count;
                for (auto r : last_row) 
                    ss << '\t' << r;
                buffer.push_back(move(ss.str()));
                ss.str({});
                if (buffer.size() == cPRINT_BUF_SIZE) {
                    {
                        const lock_guard<std::mutex> lock(m_output_mutex);
                        copy(buffer.begin(), buffer.end(), ostream_iterator<string>(os, "\n"));
                    }
                    buffer.clear();
                }
                swap(last_row, next_row);
                group_count = 1;
            }
        }
        row += batch_size;
    }
    if (group_count > 0) {
        ss << group_count;
        for (auto r : last_row) 
            ss << '\t' << r;
        buffer.push_back(move(ss.str()));
        const lock_guard<std::mutex> lock(m_output_mutex);
        copy(buffer.begin(), buffer.end(), ostream_iterator<string>(os, "\n"));
    }
    //spdlog::info("Cardinality {} done: {:.3}", cardinality, sw);
}        



template<class Options>
void Tax_hits<Options>::group_columns(const vector<uint32_t>& rows, int num_columns, ostream& os) 
{
    //spdlog::stopwatch sw; 
    size_t idx = 0;
    auto sz = rows.size();
    vector<uint32_t> prev_row(num_columns), next_row(num_columns);
    vector<string> buffer;
    tax_ids.get_row(rows[idx], prev_row);
    int count = 1;
    std::stringstream ss;

    for(idx = 1; idx < sz; ++idx) {
        tax_ids.get_row(rows[idx], next_row);
        if (prev_row == next_row) {
            ++count;
        } else {
            ss << count;
            for (auto r : prev_row) {
                ss << '\t' << r;
            }
            buffer.push_back(move(ss.str()));
            ss.str({});

            if (buffer.size() == 5000) {
                {
                    const lock_guard<std::mutex> lock(m_output_mutex);
                    copy(buffer.begin(), buffer.end(), ostream_iterator<string>(os, "\n"));
                }
                buffer.clear();
            }
            
            swap(prev_row, next_row);
            count = 1;
        }
    };
    ss << count;
    for (auto r : prev_row) {
        ss << '\t' << r;
    }
    buffer.push_back(move(ss.str()));
    ss.str({});
    {
        const lock_guard<std::mutex> lock(m_output_mutex);
        copy(buffer.begin(), buffer.end(), ostream_iterator<string>(os, "\n"));
    }
}


template<class Options>
void Tax_hits<Options>::group(tf::Executor& executor, ostream& os) 
{
    spdlog::stopwatch sw; 
    tf::Taskflow taskflow;
    vector<vector<uint32_t>> col_groups(tax_ids.num_cols);

    typedef struct sort_cache {
        uint32_t r = -1;
        uint32_t sz = 0;
        vector<uint32_t> val;
    } sort_cache;
    vector<sort_cache> s_caches(executor.num_workers());
//    std::locale::global(std::locale("C")); // disable comma as thousand separator


    for (int col_index = 0; col_index < tax_ids.num_cols; ++col_index) {

        auto& index = col_groups[col_index];
        tax_ids.get_rows(col_index, index); 
//        spdlog::info("Cardinality {}, rows {:L}", col_index + 1, index.size());
        if (index.empty()) 
            continue;
            
        auto print_task = taskflow.emplace([&, this, col_index]() {
            try {
                group_columns_bulk(col_groups[col_index], col_index + 1, os);
            } catch (exception& e) {
                spdlog::error("{}", e.what());
            }

        });
        print_task.name(fmt::format("Cardinality {} print", col_index + 1));

        if (col_index == 0) {
            auto sort_task = taskflow.sort(index.begin(), index.end(), [this, &executor, &s_caches, cardinality = col_index + 1](uint32_t l, uint32_t r)->bool {
                try {
                    auto& v = *tax_ids.data.front();
                    auto& ch = s_caches[executor.this_worker_id()];
                    if (ch.val.size() != cardinality) {
                        ch.val.resize(cardinality, 0);
                        ch.r = r;
                        ch.val[0] = tax_ids.data[0]->get(r);
                    } else if (ch.r != r) {
                        ch.r = r; 
                        ch.val[0] = v.get(r);
                    }
                    return v.get(l) < ch.val[0];
                } catch (exception& e) {
                    spdlog::error("{}", e.what());
                }
                return false;
            });
            sort_task.name(fmt::format("Cardinality {}, rows {:L}, sort", col_index + 1, index.size()));
            sort_task.precede(print_task);
        } else {

            auto sort_task = taskflow.sort(index.begin(), index.end(), [this, &executor, &s_caches, cardinality = col_index + 1](uint32_t l, uint32_t r)->bool {
                try {
                uint32_t value_l;
                uint32_t value_r;

                //static thread_local sort_cache ch;//(col_index + 1);
                auto& ch = s_caches[executor.this_worker_id()];
                if (ch.val.size() != cardinality) {
                    ch.val.resize(cardinality, 0);
                    ch.r = r;
                    ch.sz = 1;
                    ch.val[0] = tax_ids.data[0]->get(r);
                } else if (ch.r != r) {
                    ch.r = r;
                    ch.sz = 1;
                    assert(ch.val.size() == cardinality);
                    ch.val[0] = tax_ids.data[0]->get(r);
                } 

                for (int i = 0; i < cardinality; ++i) {
                    auto& v = *tax_ids.data[i];
                    if (i < ch.sz) {
                        value_r = ch.val[i];
                    } else {
                        value_r = v.get(r);
                        ch.val[ch.sz] = value_r;
                        ++ch.sz;
                    }
                    value_l = v.get(l);
                    if (value_l == value_r)
                        continue; 
                    return (value_l < value_r);
                }
                } catch (exception& e) {
                    spdlog::error("{}", e.what());
                }
                return false;    
            });
            sort_task.name(fmt::format("Cardinality {}, rows {:L}, sort", col_index + 1, index.size()));
            sort_task.precede(print_task);
        }
    }
    struct MyObserver : public tf::ObserverInterface {

        void set_up(size_t num_workers) override final {
        }

        void on_entry(tf::WorkerView w, tf::TaskView tv) override final {
            if (!tv.name().empty()) {
                const lock_guard<std::mutex> lock(m_mutex);
                timing.emplace(tv.name(), spdlog::stopwatch());
            }
        }

        void on_exit(tf::WorkerView w, tf::TaskView tv) override final {
            if (!tv.name().empty()) {
                const lock_guard<std::mutex> lock(m_mutex);
                spdlog::info("{} took {:.3}", tv.name(), timing[tv.name()]);
            }
        }
        map<string, spdlog::stopwatch> timing;
        mutex m_mutex;                          ///< protects map

    };

    std::shared_ptr<MyObserver> observer = executor.make_observer<MyObserver>();
    executor.run(taskflow).wait();
    executor.remove_observer(std::move(observer));
    os.flush();
//    std::locale::global(std::locale("en_US.UTF-8")); // enable comma as thousand separator
    spdlog::info("Grouping took {:.3}", sw);       
} 


template<class Options>
template<class CollateOptions>
tf::Future<void> Tax_hits<Options>::run_merge(tf::Executor& executor, vector<Merge_data<CollateOptions>>& merge_data, Tax_hits<CollateOptions>& tax_hits) 
{
    tf::Taskflow taskflow;
    auto merge_task = taskflow.
    for_each(merge_data.begin(), merge_data.end(), [this] (Merge_data<CollateOptions>& data)  { 
        static thread_local Spot<CollateOptions> next_spot;
        auto it = data.idx.begin();
        auto& first_spot = const_cast<Spot<CollateOptions>&>(data.spot);
        get_spot(*it, first_spot, false);
        if (data.idx.size() > 1) {
            num_merges.fetch_add(1);
            auto it_end = data.idx.end();
            while (++it != it_end) {
                get_spot(*it, next_spot, false);
                first_spot.add_taxa(next_spot);
            }
            first_spot.normalize();
        }
        data.idx.resize(0);
        data.idx.shrink_to_fit();
    });
    auto insert_task = taskflow.emplace([&]() {
        for_each(merge_data.begin(), merge_data.end(), [&](const Merge_data<CollateOptions>& data) {
            tax_hits.add_row(data.spot); 
        });
    });
    merge_task.precede(insert_task);
    
    return executor.run(move(taskflow));
}

template<class Options>
void Tax_hits<Options>::sort(tf::Executor& executor) 
{
    auto num_rows = spot_names->size();
    spot_index.resize(num_rows);
    generate(spot_index.begin(), spot_index.end(), [n = 0] () mutable { return n++; });
    spdlog::stopwatch sw; 
    tf::Taskflow taskflow;
    taskflow.sort(spot_index.begin(), spot_index.end(), [&](uint32_t  l, uint32_t  r) {
        static thread_local string last_right_str;
        static thread_local uint32_t last_right = -1;

        if (last_right != r) {
            last_right = r; 
            spot_names->get(last_right, last_right_str);
        }
        return spot_names->compare_remap(l, last_right_str.c_str()) < 0;
    });
    executor.run(taskflow).wait();
    spdlog::info("Sorting took {:.3}", sw);       

}



template<class Options>
template<class CollateOptions>
void Tax_hits<Options>::merge(tf::Executor& executor, Tax_hits<CollateOptions>& tax_hits) 
{ 
    auto num_rows = spot_names->size();
    if (num_rows == 0)
        return;
    spdlog::stopwatch sw; 
    string prev;
    size_t pos_idx = 0;
    int total_null = 0;
    static const int BATCH_SIZE = 30000;
    int max_size = BATCH_SIZE;
    int batch_size = BATCH_SIZE + 5000;

    tf::Future<void> merge_ft;
    std::future_status status{std::future_status::ready};

    vector<Merge_data<CollateOptions>> merge_data_curr;
    merge_data_curr.reserve(max_size);

    vector<Merge_data<CollateOptions>> merge_data_batch;
    merge_data_batch.reserve(max_size);

    bvector_type null_spots_curr;
    bvector_type null_spots_batch;
    null_spots_curr.init();

    uint32_t idx = 0;
    size_t first_index = idx;
    size_t last_index = idx;

    spot_names->get(spot_index[idx], prev);
    for (idx = 1; idx < num_rows; ++idx) {
        const auto& i = spot_index[idx];
        if (spot_names->compare_remap(i, prev.c_str()) == 0) {
            last_index = idx;
        } else {

            pos_idx = spot_index[first_index];
            Merge_data<CollateOptions> data(prev, pos_idx);
            null_spots_curr.set_bit_no_check(pos_idx);
            while (++first_index <= last_index) {
                pos_idx = spot_index[first_index];
                data.idx.push_back(pos_idx);
                null_spots_curr.set_bit_no_check(pos_idx);
            }
            merge_data_curr.push_back(move(data));

            if (merge_data_curr.size() >= max_size) {
                status = std::future_status::ready;
                if (merge_ft.valid()) 
                    status = merge_ft.wait_for(std::chrono::seconds(0));
                if (status == std::future_status::ready) {
                    max_size = BATCH_SIZE;
                    total_null += merge_data_curr.size();
                    swap(merge_data_curr, merge_data_batch);
                    merge_data_curr.clear();
                    if (total_null >= 10e6) {
                        set_null(executor, null_spots_batch);
                        null_spots_batch.clear();
                        total_null = 0;
                    } 
                    null_spots_batch.bit_or(null_spots_curr);
                    null_spots_curr.clear();
                    merge_ft = run_merge(executor, merge_data_batch, tax_hits);
                } else {
                    max_size = min<int>(batch_size, max_size + 100);
                }
            }
            first_index = idx;
            spot_names->get(i, prev);
        }
    }
    pos_idx = spot_index[first_index];
    Merge_data<CollateOptions> data(prev, pos_idx);
    while (++first_index <= last_index) {
        pos_idx = spot_index[first_index];
        data.idx.push_back(pos_idx);
    }
    merge_data_curr.push_back(move(data));

    spot_index.resize(0);
    spot_index.shrink_to_fit();

    if (merge_ft.valid()) {
        merge_ft.wait();
        merge_data_batch.clear();
    }

    //null_spots_curr.clear(true);
    ///set_null(executor, null_spots_batch);
    //null_spots_batch.clear(true);

    run_merge(executor, merge_data_curr, tax_hits).wait();
    spdlog::info("Collation took {:.3}", sw);
}


template<class Options>
template<class CollateOptions>
unique_ptr<Tax_hits<CollateOptions>> Tax_hits<Options>::collate(tf::Executor& executor) 
{ 
    sort(executor);
    auto merged_tax_hits = std::make_unique<Tax_hits<CollateOptions>>(false);
    merge(executor, *merged_tax_hits);
    clear();
    merged_tax_hits->finalize();
    if constexpr (CollateOptions::is_compact()) {
        merged_tax_hits->spot_names.reset(0);
        //clear_all(true);
        if constexpr (CollateOptions::has_counts()) {
            merged_tax_hits->counts.clear();
        }
    }
    spdlog::info("Number of merges {:L}", num_merges);
    return merged_tax_hits;
}

static 
void s_split(const string& str, vector<string>& out, char c = ',')
{
    if (str.empty())
        return;
    auto begin = str.begin();
    auto end  = str.end();
    while (begin != end) {
        auto it = begin;
        while (it != end && *it != c) ++it;
        out.emplace_back(begin, it);
        if (it == end)
            break;
        begin = ++it;
    } 
}


template<class Options>
void Tax_hits<Options>::load(const string& filename)
{
    spdlog::stopwatch sw;    
    bool loaded = false;
    try {
        loaded = init(filename);
    } catch (exception& e) {

    }
    if (loaded) {
        spdlog::info("Load read took {:.3}", sw);
        return;
    }
    ifstream f(filename);
    if (f.fail())
        throw std::runtime_error(std::string("cannot open list file ") + filename);

    string s;
    vector<string> out;
    Spot<Options> last_spot;
    Spot<Options> next_spot;
    size_t num_lines = 0;
    while (!f.eof())
    {
        getline(f, s);
        if (f.fail())
            break;
        out.clear();    
        s_split(s, out, '\t');
        if (out.size() < 2)
            continue;
        next_spot.init(out);
        if (next_spot.name == last_spot.name) {
            last_spot.merge(next_spot);
        } else {
            add_row(last_spot);
            swap(last_spot, next_spot);
        }
        ++num_lines;
    }
    
    add_row(last_spot);
    spdlog::info("{:L} lines", num_lines);
    finalize();
    //serialize(filename, true);
    spdlog::info("Load read took {:.3}", sw);
}


END_TC_NAMESPACE



#endif
