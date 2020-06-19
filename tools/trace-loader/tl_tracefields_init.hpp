/*======================================================================
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
 * ==========================================================================
 *
 */

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/
/*_ WARNING: GENERATED CODE, DON'T CORRECT! _*/
/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

#ifndef _tl_tracefields_init_
#define _tl_tracefields_init_


#include <vector>

#include "tl_owp.hpp"

#define _TL_ACCESSION                   "accession"
#define _TL_AMPLIFICATION_FORWARD       "amplification_forward"
#define _TL_AMPLIFICATION_REVERSE       "amplification_reverse"
#define _TL_AMPLIFICATION_SIZE          "amplification_size"
#define _TL_ANONYMIZED_ID               "anonymized_id"
#define _TL_ASSEMBLY_ID                 "assembly_id"
#define _TL_ATTEMPT                     "attempt"
#define _TL_BASECALL_LEN                "basecall_len"
#define _TL_BASES_20                    "bases_20"
#define _TL_BASES_40                    "bases_40"
#define _TL_BASES_60                    "bases_60"
#define _TL_BASE_FILE                   "base_file"
#define _TL_CENTER_NAME                 "center_name"
#define _TL_CENTER_PROJECT              "center_project"
#define _TL_CHEMISTRY                   "chemistry"
#define _TL_CHEMISTRY_TYPE              "chemistry_type"
#define _TL_CHIP_DESIGN_ID              "chip_design_id"
#define _TL_CHROMOSOME                  "chromosome"
#define _TL_CHROMOSOME_REGION           "chromosome_region"
#define _TL_CLIP_QUALITY_LEFT           "clip_quality_left"
#define _TL_CLIP_QUALITY_RIGHT          "clip_quality_right"
#define _TL_CLIP_VECTOR_LEFT            "clip_vector_left"
#define _TL_CLIP_VECTOR_RIGHT           "clip_vector_right"
#define _TL_CLONE_ID                    "clone_id"
#define _TL_CLONE_ID_LIST               "clone_id_list"
#define _TL_COLLECTION_DATE             "collection_date"
#define _TL_COMPRESS_TYPE               "compress_type"
#define _TL_CVECTOR_ACCESSION           "cvector_accession"
#define _TL_CVECTOR_CODE                "cvector_code"
#define _TL_DEPTH                       "depth"
#define _TL_ELEVATION                   "elevation"
#define _TL_ENVIRONMENT_TYPE            "environment_type"
#define _TL_EXTENDED_DATA               "extended_data"
#define _TL_FEATURE_ID_FILE_NAME        "feature_id_file_name"
#define _TL_FEATURE_SIGNAL_FILE_NAME    "feature_signal_file_name"
#define _TL_GENE_NAME                   "gene_name"
#define _TL_HAS_ALT_BASECALL            "has_alt_basecall"
#define _TL_HAS_ALT_BASECALL_CHANGE     "has_alt_basecall_change"
#define _TL_HAS_ALT_PEAK_INDEX          "has_alt_peak_index"
#define _TL_HAS_ALT_PEAK_INDEX_CHANGE     "has_alt_peak_index_change"
#define _TL_HAS_ALT_QUALITY_SCORE       "has_alt_quality_score"
#define _TL_HAS_ALT_QUALITY_SCORE_CHANGE    "has_alt_quality_score_change"
#define _TL_HAS_PEAK_INDEX              "has_peak_index"
#define _TL_HAS_QUALITY_SCORE           "has_quality_score"
#define _TL_HI_FILTER_SIZE              "hi_filter_size"
#define _TL_HOST_CONDITION              "host_condition"
#define _TL_HOST_ID                     "host_id"
#define _TL_HOST_LOCATION               "host_location"
#define _TL_HOST_SPECIES                "host_species"
#define _TL_INDIVIDUAL_ID               "individual_id"
#define _TL_INSERT_FLANK_LEFT           "insert_flank_left"
#define _TL_INSERT_FLANK_RIGHT          "insert_flank_right"
#define _TL_INSERT_SIZE                 "insert_size"
#define _TL_INSERT_STDEV                "insert_stdev"
#define _TL_ITERATION                   "iteration"
#define _TL_LATITUDE                    "latitude"
#define _TL_LIBRARY_ID                  "library_id"
#define _TL_LOAD_DATE                   "load_date"
#define _TL_LONGITUDE                   "longitude"
#define _TL_LO_FILTER_SIZE              "lo_filter_size"
#define _TL_NCBI_PROJECT_ID             "ncbi_project_id"
#define _TL_ORGANISM_NAME               "organism_name"
#define _TL_PEAK_FILE                   "peak_file"
#define _TL_PH                          "ph"
#define _TL_PICK_GROUP_ID               "pick_group_id"
#define _TL_PLACE_NAME                  "place_name"
#define _TL_PLATE_ID                    "plate_id"
#define _TL_PMID                        "pmid"
#define _TL_POPULATION_ID               "population_id"
#define _TL_PREP_GROUP_ID               "prep_group_id"
#define _TL_PRIMER                      "primer"
#define _TL_PRIMER_CODE                 "primer_code"
#define _TL_PROGRAM_ID                  "program_id"
#define _TL_PROJECT_NAME                "project_name"
#define _TL_QUAL_FILE                   "qual_file"
#define _TL_REFERENCE_ACCESSION         "reference_accession"
#define _TL_REFERENCE_ACC_MAX           "reference_acc_max"
#define _TL_REFERENCE_ACC_MIN           "reference_acc_min"
#define _TL_REFERENCE_OFFSET            "reference_offset"
#define _TL_REFERENCE_SET_MAX           "reference_set_max"
#define _TL_REFERENCE_SET_MIN           "reference_set_min"
#define _TL_REPLACED_BY                 "replaced_by"
#define _TL_RUN_DATE                    "run_date"
#define _TL_RUN_GROUP_ID                "run_group_id"
#define _TL_RUN_LANE                    "run_lane"
#define _TL_RUN_MACHINE_ID              "run_machine_id"
#define _TL_RUN_MACHINE_TYPE            "run_machine_type"
#define _TL_SALINITY                    "salinity"
#define _TL_SEQ_LIB_ID                  "seq_lib_id"
#define _TL_SOURCE_TYPE                 "source_type"
#define _TL_SPECIES_CODE                "species_code"
#define _TL_STRAIN                      "strain"
#define _TL_STRATEGY                    "strategy"
#define _TL_SUBMISSION_TYPE             "submission_type"
#define _TL_SVECTOR_ACCESSION           "svector_accession"
#define _TL_SVECTOR_CODE                "svector_code"
#define _TL_TAXID                       "taxid"
#define _TL_TEMPERATURE                 "temperature"
#define _TL_TEMPLATE_ID                 "template_id"
#define _TL_TI                          "ti"
#define _TL_TRACE_DIRECTION             "trace_direction"
#define _TL_TRACE_END                   "trace_end"
#define _TL_TRACE_FILE                  "trace_file"
#define _TL_TRACE_FORMAT                "trace_format"
#define _TL_TRACE_LEN                   "trace_len"
#define _TL_TRACE_NAME                  "trace_name"
#define _TL_TRACE_TYPE_CODE             "trace_type_code"
#define _TL_TRANSPOSON_ACC              "transposon_acc"
#define _TL_TRANSPOSON_CODE             "transposon_code"
#define _TL_WELL_ID                     "well_id"
#define _TL_SUBMISSION_ID               "submission_id"
#define _TL_TRACE_MAX_VALUE             "trace_max_value"
#define _TL_ORGANISM_NAME_OLD           "organism_name_old"
#define _TL_CONTROL_FLAGS               "control_flags"
#define _TL_STATUS                      "status"
#define _TL_UPDATE_DATE                 "update_date"
#define _TL_COMMENTS                    "comments"
#define _TL_CLIP_ADAPTER_LEFT           "clip_adapter_left"
#define _TL_CLIP_ADAPTER_RIGHT          "clip_adapter_right"
#define _TL_HAS_EXTENDED_DATA           "has_extended_data"
#define _TL_ALT_BASECALL                "alt_basecall"
#define _TL_ALT_BASECALL_CHANGE         "alt_basecall_change"
#define _TL_ALT_PEAK_INDEX              "alt_peak_index"
#define _TL_ALT_PEAK_INDEX_CHANGE       "alt_peak_index_change"
#define _TL_ALT_QUALITY_SCORE           "alt_quality_score"
#define _TL_ALT_QUALITY_SCORE_CHANGE    "alt_quality_score_change"
#define _TL_TRACE_TABLE_TYPE            "trace_table_type"
#define _TL_POSITION_OFFSET             "position_offset"
#define _TL_FLOW_CHARS                  "flow_chars"
#define _TL_KEY_SEQUENCE                "key_sequence"


namespace _tl_ {

typedef std::vector < TL_Owp > TL_OVec;

void TL_GetFieldDescriptions ( TL_OVec & Descriptions );

}   /* namespace _tl_ */


#endif /* _tl_tracefields_init_ */
