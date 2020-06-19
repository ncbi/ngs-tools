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


#include "tl_tracefields.hpp"
#include "tl_tracefields_init.hpp"

using namespace std;
using namespace _tl_;





static
void
_add_description_of_ACCESSION ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "accession" );

    Owp . SetValue ( "comment", "\n// *\n// accession" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Genbank/EMBL/DDBJ accession number" );
    Owp . SetValue ( "field_prim", "accession" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_name", "accession" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ACCESSION () */



static
void
_add_description_of_AMPLIFICATION_FORWARD ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "amplification_forward" );

    Owp . SetValue ( "comment", "\n// *\n// amplification_forward\n// Field is dictionary from AmplificationForward\n// U32 -> ASCII. Has 400K values like: GACGACGAAGACGAAGAAGG, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The forward amplification primer sequence" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "amplification_forward" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "AmplificationForward" );
    Owp . SetValue ( "vdb_name", "amplification_forward" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_AMPLIFICATION_FORWARD () */



static
void
_add_description_of_AMPLIFICATION_REVERSE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "amplification_reverse" );

    Owp . SetValue ( "comment", "\n// *\n// amplification_reverse\n// Field is dictionary from AmplificationReverse\n// U32 -> ASCII. Has 400K values like: AACTCGCTCCAATGAGGAAA, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The reverse amplification primer sequence" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "amplification_reverse" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "AmplificationReverse" );
    Owp . SetValue ( "vdb_name", "amplification_reverse" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_AMPLIFICATION_REVERSE () */



static
void
_add_description_of_AMPLIFICATION_SIZE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "amplification_size" );

    Owp . SetValue ( "comment", "\n// *\n// amplification_size" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The expected amplification size for a pair of primers" );
    Owp . SetValue ( "field_prim", "amplification_size" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "amplification_size" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_AMPLIFICATION_SIZE () */



static
void
_add_description_of_ANONYMIZED_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "anonymized_id" );

    Owp . SetValue ( "comment", "\n// *\n// anonymized_id\n// Field is dictionary from AnonymizedId\n// U32 -> ASCII. Has 15M values like: 10685146475, 10690251308, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Anonymous ID for an individual" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "anonymized_id" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "AnonymizedId" );
    Owp . SetValue ( "vdb_name", "anonymized_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ANONYMIZED_ID () */



static
void
_add_description_of_ASSEMBLY_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "assembly_id" );

    Owp . SetValue ( "comment", "\n// *\n// assembly_id\n// Field is dictionary from AssemblyId\n// U32 -> ASCII. Has 2 values: \"NCBI BUILD 34\", and\n// \"HUMAN BUILD 33, APRIL 2003 FREEZE\"." );
    Owp . SetValue ( "description", "Public identifier for a given version of a genome assembly" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "assembly_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "table_name", "AssemblyId" );
    Owp . SetValue ( "vdb_name", "assembly_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ASSEMBLY_ID () */



static
void
_add_description_of_ATTEMPT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "attempt" );

    Owp . SetValue ( "comment", "\n// *\n// attempt" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Number of times the sequencing project has been attempted by the center and/or submitted to the Trace Archive." );
    Owp . SetValue ( "field_prim", "attempt" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "attempt" );
    Owp . SetValue ( "vdb_type", "U8" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ATTEMPT () */



static
void
_add_description_of_BASECALL_LEN ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "basecall_len" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Length of base call" );
    Owp . SetValue ( "field_prim", "basecall_len" );
    Owp . SetValue ( "vdb_name", "basecall_len" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_BASECALL_LEN () */



static
void
_add_description_of_BASES_20 ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "bases_20" );

    Owp . SetValue ( "comment", "\n// *\n// bases_20" );
    Owp . SetValue ( "description", "Number of quality scores which exceed 20" );
    Owp . SetValue ( "field_prim", "bases_20" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "bases_20" );
    Owp . SetValue ( "vdb_type", "U16" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_BASES_20 () */



static
void
_add_description_of_BASES_40 ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "bases_40" );

    Owp . SetValue ( "comment", "\n// *\n// bases_40" );
    Owp . SetValue ( "description", "Number of quality scores which exceed 40" );
    Owp . SetValue ( "field_prim", "bases_40" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "bases_40" );
    Owp . SetValue ( "vdb_type", "U16" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_BASES_40 () */



static
void
_add_description_of_BASES_60 ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "bases_60" );

    Owp . SetValue ( "comment", "\n// *\n// bases_60" );
    Owp . SetValue ( "description", "Number of quality scores which exceed 60" );
    Owp . SetValue ( "field_prim", "bases_60" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "bases_60" );
    Owp . SetValue ( "vdb_type", "U16" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_BASES_60 () */



static
void
_add_description_of_BASE_FILE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "base_file" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "File name with base calls" );
    Owp . SetValue ( "field_prim", "base_file" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_name", "base_file" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_BASE_FILE () */



static
void
_add_description_of_CENTER_NAME ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "center_name" );

    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// center_name\n// Field is dictionary from Center\n// U16 -> ASCII. Has 184 names like: 454MSC, ABC, ABI, AGI, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Name of the sequencing center" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "center_name" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "mandatory", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "Center" );
    Owp . SetValue ( "vdb_name", "center_name" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CENTER_NAME () */



static
void
_add_description_of_CENTER_PROJECT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "center_project" );

    Owp . SetValue ( "comment", "\n// *\n// center_project\n// Field is dictionary from Project\n// U32 -> ASCII. Has 2M values like: 30153, 31491, 31489, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Center defined project name" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "center_project" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "Project" );
    Owp . SetValue ( "vdb_name", "center_project" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CENTER_PROJECT () */



static
void
_add_description_of_CHEMISTRY ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "chemistry" );

    Owp . SetValue ( "comment", "\n// *\n// chemistry\n// Field is dictionary from Chemistry\n// U16 -> ASCII. Has 60 values, and one of them with id = 4 is \"\"" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Description of the chemistry used in the sequencing reaction" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "chemistry" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "Chemistry" );
    Owp . SetValue ( "vdb_name", "chemistry" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CHEMISTRY () */



static
void
_add_description_of_CHEMISTRY_TYPE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "chemistry_type" );

    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// chemistry_type\n// Field is dictionary from ChemistryType\n// U16 -> ASCII\n// Has 3 chemistry types, and 6 names:\n// 1 { T, TERM, TERMINATOR }, 2 { UNKNOWN }, 3 { P, PRIMER }" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Type of chemistry used in the sequencing reaction" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "chemistry_type" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "ChemistryType" );
    Owp . SetValue ( "vdb_name", "chemistry_type" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CHEMISTRY_TYPE () */



static
void
_add_description_of_CHIP_DESIGN_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "chip_design_id" );

    Owp . SetValue ( "comment", "\n// *\n// chip_design_id\n// Was added, cuz that value could be added\n// Mirrored feature_id_file_name" );
    Owp . SetValue ( "description", "Chip Design ID" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "feature_id_file_name" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "FeatureIdFile" );
    Owp . SetValue ( "vdb_name", "chip_design_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CHIP_DESIGN_ID () */



static
void
_add_description_of_CHROMOSOME ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "chromosome" );

    Owp . SetValue ( "comment", "\n// *\n// chromosome\n// Field is dictionary from Chromosome\n// U8 -> ASCII. Has 171 value like: NA, 11, 7, 12, 14, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Chromosome to which the trace is assigned" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "chromosome" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "Chromosome" );
    Owp . SetValue ( "vdb_name", "chromosome" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CHROMOSOME () */



static
void
_add_description_of_CHROMOSOME_REGION ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "chromosome_region" );

    Owp . SetValue ( "comment", "\n// *\n// chromosome_region\n// Field is dictionary from ChromosomeRegion\n// U32 -> ASCII. Has 5K values like: \"7:26665793-2716579226\"" );
    Owp . SetValue ( "description", "This field is a required field if the ASSEMBLY_ID field is not NULL. This field can be used to describe regions that are being re-sequenced for a given genome." );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "chromosome_region" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "ChromosomeRegion" );
    Owp . SetValue ( "vdb_name", "chromosome_region" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CHROMOSOME_REGION () */



static
void
_add_description_of_CLIP_QUALITY_LEFT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "clip_quality_left" );

    Owp . SetValue ( "description", "Left clip of the read, in base pairs, based on quality analysis" );
    Owp . SetValue ( "field_prim", "clip_quality_left" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_name", "clip_quality_left" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CLIP_QUALITY_LEFT () */



static
void
_add_description_of_CLIP_QUALITY_RIGHT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "clip_quality_right" );

    Owp . SetValue ( "description", "Right clip of the read, in base pairs, based on quality analysis" );
    Owp . SetValue ( "field_prim", "clip_quality_right" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_name", "clip_quality_right" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CLIP_QUALITY_RIGHT () */



static
void
_add_description_of_CLIP_VECTOR_LEFT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "clip_vector_left" );

    Owp . SetValue ( "comment", "\n// *\n// clip_vector_left" );
    Owp . SetValue ( "description", "Left clip of the read, in base pairs, based on vector sequence" );
    Owp . SetValue ( "field_prim", "clip_vector_left" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "clip_vector_left" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CLIP_VECTOR_LEFT () */



static
void
_add_description_of_CLIP_VECTOR_RIGHT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "clip_vector_right" );

    Owp . SetValue ( "comment", "\n// *\n// clip_vector_right" );
    Owp . SetValue ( "description", "Right clip of the read, in base pairs, based on vector sequence" );
    Owp . SetValue ( "field_prim", "clip_vector_right" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "clip_vector_right" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CLIP_VECTOR_RIGHT () */



static
void
_add_description_of_CLONE_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "clone_id" );

    Owp . SetValue ( "comment", "\n// *\n// clone_id" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The name of the clone from which the trace was derived" );
    Owp . SetValue ( "field_prim", "clone_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_name", "clone_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CLONE_ID () */



static
void
_add_description_of_CLONE_ID_LIST ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "clone_id_list" );

    Owp . SetValue ( "comment", "\n// *\n// clone_id_list" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Semi-colon delimited list of clones if the Strategy is PoolClone" );
    Owp . SetValue ( "field_prim", "dummy" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "CloneList" );
    Owp . SetValue ( "vdb_name", "clone_id_list" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CLONE_ID_LIST () */



static
void
_add_description_of_COLLECTION_DATE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "collection_date" );

    Owp . SetValue ( "comment", "\n// *\n// collection_date\n// smalldatetime" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The full date, in 'Mar 2 2006 12:00AM' format, on which an environmental sample was collected" );
    Owp . SetValue ( "field_prim", "collection_date" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "collection_date" );
    Owp . SetValue ( "vdb_type", "loctm" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_COLLECTION_DATE () */



static
void
_add_description_of_COMPRESS_TYPE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "compress_type" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Compression method of the original file" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "compress_type" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "CompressType" );
    Owp . SetValue ( "vdb_name", "compress_type" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_COMPRESS_TYPE () */



static
void
_add_description_of_CVECTOR_ACCESSION ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "cvector_accession" );

    Owp . SetValue ( "comment", "\n// *\n// cvector_accession\n// Field is dictionary from CloneVectorAcc\n// U16 -> ASCII. Has 16 values like: U80929, U13871, U75992, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Repository ( GenBank/EMBL/DDBJ) accession identifier for the cloning vector" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "cvector_accession" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "CloneVectorAcc" );
    Owp . SetValue ( "vdb_name", "cvector_accession" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CVECTOR_ACCESSION () */



static
void
_add_description_of_CVECTOR_CODE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "cvector_code" );

    Owp . SetValue ( "comment", "\n// *\n// cvector_code\n// Field is dictionary from CloneVector\n// U16 -> ASCII. Has 125 values like: PBACE3.6, PT7T3PAC, PPAC4, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Center defined code for the cloning vector" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "cvector_code" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "CloneVector" );
    Owp . SetValue ( "vdb_name", "cvector_code" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CVECTOR_CODE () */



static
void
_add_description_of_DEPTH ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "depth" );

    Owp . SetValue ( "comment", "\n// *\n// depth" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Depth (in meters) at which an environmental sample was collected" );
    Owp . SetValue ( "field_prim", "depth" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "depth" );
    Owp . SetValue ( "vdb_type", "F64" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_DEPTH () */



static
void
_add_description_of_ELEVATION ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "elevation" );

    Owp . SetValue ( "comment", "\n// *\n// elevation" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Elevation (in meters) at which an environmental sample was collected" );
    Owp . SetValue ( "field_prim", "elevation" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "elevation" );
    Owp . SetValue ( "vdb_type", "F64" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ELEVATION () */



static
void
_add_description_of_ENVIRONMENT_TYPE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "environment_type" );

    Owp . SetValue ( "comment", "\n// *\n// environment_type\n// Field is dictionary from EnvironmentType\n// U32 -> ASCII. Has 52 values like: \"BEACH SAND\",\n// \"HUMAN STOOL\", \"NORTH ATLANTIC DEEP WATER\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Type of environment from which an environmental sample was collected" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "environment_type" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "EnvironmentType" );
    Owp . SetValue ( "vdb_name", "environment_type" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ENVIRONMENT_TYPE () */



static
void
_add_description_of_EXTENDED_DATA ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "extended_data" );

    Owp . SetValue ( "comment", "\n// *\n// extended_data" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Extra ancillary information wrapped around in a EXTENDED_DATA block, where actual values are provided with a special 'field' tag'" );
    Owp . SetValue ( "field_norm", "data" );
    Owp . SetValue ( "field_prim", "ti" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "TraceExtended..ExtendedData" );
    Owp . SetValue ( "vdb_name", "extended_data" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_EXTENDED_DATA () */



static
void
_add_description_of_FEATURE_ID_FILE_NAME ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "feature_id_file_name" );

    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Reference to a common FEATURE_ID_FILE" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "feature_id_file_name" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "FeatureIdFile" );
    Owp . SetValue ( "vdb_name", "feature_id_file_name" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_FEATURE_ID_FILE_NAME () */



static
void
_add_description_of_FEATURE_SIGNAL_FILE_NAME ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "feature_signal_file_name" );

    Owp . SetValue ( "comment", "\n// *\n// feature_id_file_name\n// Field is dictionary from FeatureIdFile\n// U32 -> ASCII. Has 4k values like: PSC.CHRN_077.DES33, etc\n// *\n// feature_signal_file_name\n// Field is dictionary from FeatureSignalFile\n// U32 -> ASCII. Has only one value: N/A ... lol" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Reference to a common FEATURE_SIGNAL_FILE" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "feature_signal_file_name" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "FeatureSignalFile" );
    Owp . SetValue ( "vdb_name", "feature_signal_file_name" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_FEATURE_SIGNAL_FILE_NAME () */



static
void
_add_description_of_GENE_NAME ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "gene_name" );

    Owp . SetValue ( "comment", "\n// *\n// gene_name\n// Field is dictionary from GeneName\n// U32 -> ASCII. Has 30K values like: COI, HIF1A, PRKCB1, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Gene name or some other common identifier" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "gene_name" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "GeneName" );
    Owp . SetValue ( "vdb_name", "gene_name" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_GENE_NAME () */



static
void
_add_description_of_HAS_ALT_BASECALL ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "has_alt_basecall" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Description is MISSING" );
    Owp . SetValue ( "field_prim", "alt_basecall" );
    Owp . SetValue ( "vdb_name", "has_alt_basecall" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HAS_ALT_BASECALL () */



static
void
_add_description_of_HAS_ALT_BASECALL_CHANGE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "has_alt_basecall_change" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Description is MISSING" );
    Owp . SetValue ( "field_prim", "alt_basecall_change" );
    Owp . SetValue ( "vdb_name", "has_alt_basecall_change" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HAS_ALT_BASECALL_CHANGE () */



static
void
_add_description_of_HAS_ALT_PEAK_INDEX ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "has_alt_peak_index" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Description is MISSING" );
    Owp . SetValue ( "field_prim", "alt_peak_index" );
    Owp . SetValue ( "vdb_name", "has_alt_peak_index" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HAS_ALT_PEAK_INDEX () */



static
void
_add_description_of_HAS_ALT_PEAK_INDEX_CHANGE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "has_alt_peak_index_change" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Description is MISSING" );
    Owp . SetValue ( "field_prim", "alt_peak_index_change" );
    Owp . SetValue ( "vdb_name", "has_alt_peak_index_change" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HAS_ALT_PEAK_INDEX_CHANGE () */



static
void
_add_description_of_HAS_ALT_QUALITY_SCORE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "has_alt_quality_score" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Description is MISSING" );
    Owp . SetValue ( "field_prim", "alt_quality_score" );
    Owp . SetValue ( "vdb_name", "has_alt_quality_score" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HAS_ALT_QUALITY_SCORE () */



static
void
_add_description_of_HAS_ALT_QUALITY_SCORE_CHANGE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "has_alt_quality_score_change" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Description is MISSING" );
    Owp . SetValue ( "field_prim", "alt_quality_score_change" );
    Owp . SetValue ( "vdb_name", "has_alt_quality_score_change" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HAS_ALT_QUALITY_SCORE_CHANGE () */



static
void
_add_description_of_HAS_PEAK_INDEX ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "has_peak_index" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "When peak indecies are present for the record this filed is set to 1" );
    Owp . SetValue ( "field_prim", "peak_index" );
    Owp . SetValue ( "vdb_name", "has_peak_index" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HAS_PEAK_INDEX () */



static
void
_add_description_of_HAS_QUALITY_SCORE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "has_quality_score" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "When quality scores are present for the record this filed is set to 1" );
    Owp . SetValue ( "field_prim", "quality_score" );
    Owp . SetValue ( "vdb_name", "has_quality_score" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HAS_QUALITY_SCORE () */



static
void
_add_description_of_HI_FILTER_SIZE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "hi_filter_size" );

    Owp . SetValue ( "comment", "\n// *\n// hi_filter_size\n// Field is dictionary from HiFilterSize\n// U32 -> ASCII. Has 6 values like: \"0.1 KD-MICRONS\", \"0.22\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The largest filter used to stratify an environmental sample" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "hi_filter_size" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "HiFilterSize" );
    Owp . SetValue ( "vdb_name", "hi_filter_size" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HI_FILTER_SIZE () */



static
void
_add_description_of_HOST_CONDITION ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "host_condition" );

    Owp . SetValue ( "comment", "\n// *\n// host_condition\n// Field is dictionary from HostCondition\n// U32 -> ASCII. Has only one value: \"HEALTHY\"" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The condition of the host from which an environmental sample was obtained" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "host_condition" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "HostCondition" );
    Owp . SetValue ( "vdb_name", "host_condition" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HOST_CONDITION () */



static
void
_add_description_of_HOST_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "host_id" );

    Owp . SetValue ( "comment", "\n// *\n// host_id\n// Field is dictionary from HostIdentifier\n// U32 -> ASCII. Has 69 values like: 01DBA0010, 01DBA0011, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Unique identifier for the specific host from which an environmental sample was taken" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "host_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "HostIdentifier" );
    Owp . SetValue ( "vdb_name", "host_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HOST_ID () */



static
void
_add_description_of_HOST_LOCATION ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "host_location" );

    Owp . SetValue ( "comment", "\n// *\n// host_location\n// Field is dictionary from HostLocation\n// U32 -> ASCII. Contains 30 values like: \"BELOW CUTICLE\", \"BONE\",\n// \"CECUM\", \"CORPUS\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Specific location on the host from which an environmental sample was collected" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "host_location" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "HostLocation" );
    Owp . SetValue ( "vdb_name", "host_location" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HOST_LOCATION () */



static
void
_add_description_of_HOST_SPECIES ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "host_species" );

    Owp . SetValue ( "comment", "\n// *\n// host_species\n// Field is dictionary from HostSpecies\n// U32 -> ASCII. Has 11 values like: \"ALVINELA POMPEJANA EPIBIONT\"\n// \"ESCHRICHTIUS ROBUSTUS\", \"HOMO NEANDERTHALENSIS\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The host from which an environmental sample was obtained" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "host_species" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "HostSpecies" );
    Owp . SetValue ( "vdb_name", "host_species" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HOST_SPECIES () */



static
void
_add_description_of_INDIVIDUAL_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "individual_id" );

    Owp . SetValue ( "comment", "\n// *\n// individual_id\n// Field is dictionary from IndividualId\n// U32 -> ASCII. Has 3M values like: \"JDH LADY MANSO 62/4\",\n// \"WA-DEL RC BLCKSTR MARTHA-ET (HO-USA 13907649)\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Publicly available identifier to denote a specific individual or sample from which a trace was derived" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "individual_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "IndividualId" );
    Owp . SetValue ( "vdb_name", "individual_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_INDIVIDUAL_ID () */



static
void
_add_description_of_INSERT_FLANK_LEFT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "insert_flank_left" );

    Owp . SetValue ( "comment", "\n// *\n// insert_flank_left\n// Field is dictionary from InsertFlankLeft\n// U32 -> ASCII. Has 266 values like: AATACGACTCACTATAGGGCGAATTCGAGCTCGGTACCCGGGGATCCCAC\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Left flanking sequence at the cloning junction" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "insert_flank_left" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "InsertFlankLeft" );
    Owp . SetValue ( "vdb_name", "insert_flank_left" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_INSERT_FLANK_LEFT () */



static
void
_add_description_of_INSERT_FLANK_RIGHT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "insert_flank_right" );

    Owp . SetValue ( "comment", "\n// *\n// insert_flank_right\n// Field is dictionary from InsertFlankRight\n// U32 -> ASCII. Has 221 values like: GTGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGAGTATTCTAT, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Right flanking sequence at the cloning junction" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "insert_flank_right" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "InsertFlankRight" );
    Owp . SetValue ( "vdb_name", "insert_flank_right" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_INSERT_FLANK_RIGHT () */



static
void
_add_description_of_INSERT_SIZE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "insert_size" );

    Owp . SetValue ( "comment", "\n// *\n// insert_size" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Expected size of the insert (referred to by the value in the TEMPLATE_ID field) in base pairs" );
    Owp . SetValue ( "field_prim", "insert_size" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "insert_size" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_INSERT_SIZE () */



static
void
_add_description_of_INSERT_STDEV ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "insert_stdev" );

    Owp . SetValue ( "comment", "\n// *\n// insert_stdev" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Approximate standard deviation of value in INSERT_SIZE field" );
    Owp . SetValue ( "field_prim", "insert_stdev" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "insert_stdev" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_INSERT_STDEV () */



static
void
_add_description_of_ITERATION ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "iteration" );

    Owp . SetValue ( "comment", "\n// *\n// iteration - it did exists some tiome ago" );
    Owp . SetValue ( "description", "Attempt, Redo" );
    Owp . SetValue ( "field_prim", "attempt" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "iteration" );
    Owp . SetValue ( "vdb_type", "U16" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ITERATION () */



static
void
_add_description_of_LATITUDE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "latitude" );

    Owp . SetValue ( "comment", "\n// *\n// latitude" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The latitude measurement (using standard GPS notation) from which a sample was collected" );
    Owp . SetValue ( "field_prim", "latitude" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "latitude" );
    Owp . SetValue ( "vdb_type", "F64" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_LATITUDE () */



static
void
_add_description_of_LIBRARY_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "library_id" );

    Owp . SetValue ( "comment", "\n// *\n// library_id\n// Field is dictionary from Library\n// U16 -> ASCII. Has 11K values like: CH230, RP1, RP11, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The source of the clone identified in the CLONE_ID field" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "library_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "Library" );
    Owp . SetValue ( "vdb_name", "library_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_LIBRARY_ID () */



static
void
_add_description_of_LOAD_DATE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "load_date" );

    Owp . SetValue ( "comment", "\n// *\n// load_date\n// smalldatetime" );
    Owp . SetValue ( "description", "Date of load this data" );
    Owp . SetValue ( "field_prim", "load_date" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "load_date" );
    Owp . SetValue ( "vdb_type", "loctm" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_LOAD_DATE () */



static
void
_add_description_of_LONGITUDE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "longitude" );

    Owp . SetValue ( "comment", "\n// *\n// longitude" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The longitude measurement (using standard GPS notation) from which a sample was collected" );
    Owp . SetValue ( "field_prim", "longitude" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "longitude" );
    Owp . SetValue ( "vdb_type", "F64" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_LONGITUDE () */



static
void
_add_description_of_LO_FILTER_SIZE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "lo_filter_size" );

    Owp . SetValue ( "comment", "\n// *\n// lo_filter_size\n// Field is dictionary from LoFilterSize\n// U32 -> ASCII. Has 7 values like: \"0.0 MICRONS\", \"0.002\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The smallest filter size used to stratify an environmental sample" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "lo_filter_size" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "LoFilterSize" );
    Owp . SetValue ( "vdb_name", "lo_filter_size" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_LO_FILTER_SIZE () */



static
void
_add_description_of_NCBI_PROJECT_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "ncbi_project_id" );

    Owp . SetValue ( "comment", "\n// *\n// ncbi_project_id\n// Field is dictionary from ProjectContent\n// U32 - as is ... there are 188 values" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Project ID generated by the Genome Project database at NCBI/NLM/NIH" );
    Owp . SetValue ( "field_prim", "ncbi_project_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "ncbi_project_id" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_NCBI_PROJECT_ID () */



static
void
_add_description_of_ORGANISM_NAME ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "organism_name" );

    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// organism_name\n// Field is dictionary from OrganismName\n// U32 -> ASCII. Has 170K names like: \"THUNNUS ATLANTICUS\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Description of species for BARCODE project from which trace is derived" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "organism_name" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "OrganismName" );
    Owp . SetValue ( "vdb_name", "organism_name" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ORGANISM_NAME () */



static
void
_add_description_of_PEAK_FILE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "peak_file" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Name of file that contains the list of peak values" );
    Owp . SetValue ( "field_prim", "peak_file" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_name", "peak_file" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PEAK_FILE () */



static
void
_add_description_of_PH ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "ph" );

    Owp . SetValue ( "comment", "\n// *\n// ph" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The pH at which an environmental sample was collected" );
    Owp . SetValue ( "field_prim", "ph" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "ph" );
    Owp . SetValue ( "vdb_type", "F64" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PH () */



static
void
_add_description_of_PICK_GROUP_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "pick_group_id" );

    Owp . SetValue ( "comment", "\n// *\n// pick_group_id\n// Field is dictionary from PickGroup\n// U32 -> ASCII. Has 1.5M values like: 893877, 895133, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Id to group traces picked at the same time" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "pick_group_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "PickGroup" );
    Owp . SetValue ( "vdb_name", "pick_group_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PICK_GROUP_ID () */



static
void
_add_description_of_PLACE_NAME ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "place_name" );

    Owp . SetValue ( "comment", "\n// *\n// place_name\n// Field is dictionary from PlaceName\n// U32 -> ASCII. Has 82 values like: \"HYDROSTATION S, BERMUDA (UK)\"" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Country in which the biological sample was collected and/or common name for a given location" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "place_name" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "PlaceName" );
    Owp . SetValue ( "vdb_name", "place_name" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PLACE_NAME () */



static
void
_add_description_of_PLATE_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "plate_id" );

    Owp . SetValue ( "comment", "\n// *\n// plate_id\n// Field is dictionary from Plate\n// U32 -> ASCII. Has 7M values like: ML1B-A1156, ML1B-A1158, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Submitter defined plate id" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "plate_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "Plate" );
    Owp . SetValue ( "vdb_name", "plate_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PLATE_ID () */



static
void
_add_description_of_PMID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "pmid" );

    Owp . SetValue ( "comment", "\n// *\n// pmid" );
    Owp . SetValue ( "description", "PubMed ID, could be linked" );
    Owp . SetValue ( "field_prim", "pmid" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "pmid" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PMID () */



static
void
_add_description_of_POPULATION_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "population_id" );

    Owp . SetValue ( "comment", "\n// *\n// population_id\n// Field is dictionary from Population\n// U32 -> ASCII. Has 10 values like: \"AFRICAN\", \"CAUCASIAN\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Center provided id to designate a population from which a trace (or group of traces) was derived" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "population_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "Population" );
    Owp . SetValue ( "vdb_name", "population_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_POPULATION_ID () */



static
void
_add_description_of_PREP_GROUP_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "prep_group_id" );

    Owp . SetValue ( "comment", "\n// *\n// prep_group_id\n// Field is dictionary from PrepGroup\n// U32 -> ASCII. Has 3M values like: 6AUG03.842PMAB1, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Id to group traces prepared at the same time" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "prep_group_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "PrepGroup" );
    Owp . SetValue ( "vdb_name", "prep_group_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PREP_GROUP_ID () */



static
void
_add_description_of_PRIMER ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "primer" );

    Owp . SetValue ( "comment", "\n// *\n// primer" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The primer sequence (used in the sequencing reaction)" );
    Owp . SetValue ( "field_prim", "primer" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_name", "primer" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PRIMER () */



static
void
_add_description_of_PRIMER_CODE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "primer_code" );

    Owp . SetValue ( "comment", "\n// *\n// primer_code\n// Field is dictionary from PrimerCode\n// U32 -> ASCII. Has 95K values like: \"M13 FORWARD\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Identifier for the sequencing primer used" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "primer_code" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "PrimerCode" );
    Owp . SetValue ( "vdb_name", "primer_code" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PRIMER_CODE () */



static
void
_add_description_of_PROGRAM_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "program_id" );

    Owp . SetValue ( "comment", "\n// *\n// program_id\n// Field is dictionary from Program\n// U16 -> ASCII\n// There are 1193 values, which looks quite unique, no decriptions" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The program used to create the trace file" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "program_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "Program" );
    Owp . SetValue ( "vdb_name", "program_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PROGRAM_ID () */



static
void
_add_description_of_PROJECT_NAME ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "project_name" );

    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// project_name\n// Field is dictionary from ProjectName\n// U32 -> ASCII. Has 29 values like: \"ARALL\", \"BARCODE\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Term by which to group traces based on a common project." );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "project_name" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "ProjectName" );
    Owp . SetValue ( "vdb_name", "project_name" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_PROJECT_NAME () */



static
void
_add_description_of_QUAL_FILE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "qual_file" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Name of file containing the quality scores" );
    Owp . SetValue ( "field_prim", "qual_file" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_name", "qual_file" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_QUAL_FILE () */



static
void
_add_description_of_REFERENCE_ACCESSION ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "reference_accession" );

    Owp . SetValue ( "comment", "\n// *\n// reference_accession\n// Field is dictionary from ReferenceAccession\n// U32 -> ASCII. Has 6K names like: NT_039520.1, NT_039268.1, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Reference accession (use accession and version to specify a particular instance of a sequence) used as the basis for a re-sequencing project. In case of Comparative strategy show the basis for primers design" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "reference_accession" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "ReferenceAccession" );
    Owp . SetValue ( "vdb_name", "reference_accession" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_REFERENCE_ACCESSION () */



static
void
_add_description_of_REFERENCE_ACC_MAX ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "reference_acc_max" );

    Owp . SetValue ( "comment", "\n// *\n// reference_acc_max" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Finish position for a particular amplicon in re-sequencing or comparative projects" );
    Owp . SetValue ( "field_prim", "reference_acc_max" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "reference_acc_max" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_REFERENCE_ACC_MAX () */



static
void
_add_description_of_REFERENCE_ACC_MIN ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "reference_acc_min" );

    Owp . SetValue ( "comment", "\n// *\n// reference_acc_min" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Start position for a particular amplicon in re-sequencing or comparative projects" );
    Owp . SetValue ( "field_prim", "reference_acc_min" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "reference_acc_min" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_REFERENCE_ACC_MIN () */



static
void
_add_description_of_REFERENCE_OFFSET ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "reference_offset" );

    Owp . SetValue ( "comment", "\n// *\n// reference_offset" );
    Owp . SetValue ( "description", "Sequence offset of accession specified in REFERENCE_ACCESSION field to define the coordinate start position used as the basis for a re-sequencing project" );
    Owp . SetValue ( "field_prim", "reference_offset" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "reference_offset" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_REFERENCE_OFFSET () */



static
void
_add_description_of_REFERENCE_SET_MAX ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "reference_set_max" );

    Owp . SetValue ( "comment", "\n// *\n// reference_set_max" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Finish position for a entire re-sequencing region. This region may include several amplicons" );
    Owp . SetValue ( "field_prim", "reference_set_max" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "reference_set_max" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_REFERENCE_SET_MAX () */



static
void
_add_description_of_REFERENCE_SET_MIN ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "reference_set_min" );

    Owp . SetValue ( "comment", "\n// *\n// reference_set_min" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Start position for a entire re-sequencing region. This region may include several amplicons" );
    Owp . SetValue ( "field_prim", "reference_set_min" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "reference_set_min" );
    Owp . SetValue ( "vdb_type", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_REFERENCE_SET_MIN () */



static
void
_add_description_of_REPLACED_BY ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "replaced_by" );

    Owp . SetValue ( "comment", "\n// *\n// replaced_by" );
    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Indicates which took over" );
    Owp . SetValue ( "field_prim", "replaced_by" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "replaced_by" );
    Owp . SetValue ( "vdb_type", "U64" );
    Owp . SetValue ( "vdb_type_cast", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_REPLACED_BY () */



static
void
_add_description_of_RUN_DATE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "run_date" );

    Owp . SetValue ( "comment", "\n// *\n// run_date\n// smalldatetime" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Date the sequencing reaction was run" );
    Owp . SetValue ( "field_prim", "run_date" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "run_date" );
    Owp . SetValue ( "vdb_type", "loctm" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_RUN_DATE () */



static
void
_add_description_of_RUN_GROUP_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "run_group_id" );

    Owp . SetValue ( "comment", "\n// *\n// run_group_id\n// Field is dictionary from RunGroup\n// U32 -> ASCII. Has 5M values like: 2000-01-24-GEA-D, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Id to group traces run on the same machine" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "run_group_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "RunGroup" );
    Owp . SetValue ( "vdb_name", "run_group_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_RUN_GROUP_ID () */



static
void
_add_description_of_RUN_LANE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "run_lane" );

    Owp . SetValue ( "comment", "\n// *\n// run_lane" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Lane or capillary of the trace" );
    Owp . SetValue ( "field_prim", "run_lane" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "run_lane" );
    Owp . SetValue ( "vdb_type", "U16" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_RUN_LANE () */



static
void
_add_description_of_RUN_MACHINE_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "run_machine_id" );

    Owp . SetValue ( "comment", "\n// *\n// run_machine_id\n// Field is dictionary from Machine\n// U16 -> ASCII. Has 2K values like: 194, 137, 130, 300, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "ID of the specific sequencing machine on which a trace was obtained" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "run_machine_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "Machine" );
    Owp . SetValue ( "vdb_name", "run_machine_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_RUN_MACHINE_ID () */



static
void
_add_description_of_RUN_MACHINE_TYPE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "run_machine_type" );

    Owp . SetValue ( "comment", "\n// *\n// run_machine_type\n// Field is dictionary from MachineType\n// U16 -> ASCII. Has 160 values like: 194, 137, 130, 300, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Type or model of machine on which a trace was obtained" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "run_machine_type" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "MachineType" );
    Owp . SetValue ( "vdb_name", "run_machine_type" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_RUN_MACHINE_TYPE () */



static
void
_add_description_of_SALINITY ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "salinity" );

    Owp . SetValue ( "comment", "\n// *\n// salinity" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The salinity at which an environmental sample was collected measured in parts per thousand units (promille)" );
    Owp . SetValue ( "field_prim", "salinity" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "salinity" );
    Owp . SetValue ( "vdb_type", "F64" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_SALINITY () */



static
void
_add_description_of_SEQ_LIB_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "seq_lib_id" );

    Owp . SetValue ( "comment", "\n// *\n// seq_lib_id\n// Field is dictionary from SeqLibrary\n// U32 -> ASCII. Has 700K values like: 33802, RATBN2.5.2L, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Center specified M13/PUC library that is actually sequenced" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "seq_lib_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "SeqLibrary" );
    Owp . SetValue ( "vdb_name", "seq_lib_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_SEQ_LIB_ID () */



static
void
_add_description_of_SOURCE_TYPE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "source_type" );

    Owp . SetValue ( "can_be_ignored", "true" );
    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// source_type\n// Field is dictionary from Source\n// U8 -> ASCII. Has 7 values: G, GENOMIC, N, NON GENOMIC, S,\n//                            SYNTHETIC, VIRAL RNA" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Source of the DNA" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "source_type" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "mandatory", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "Source" );
    Owp . SetValue ( "vdb_name", "source_type" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_SOURCE_TYPE () */



static
void
_add_description_of_SPECIES_CODE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "species_code" );

    Owp . SetValue ( "can_be_ignored", "true" );
    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// species_code\n// Field is dictionary from Species\n// U32 -> ASCII. Has 12K values like: \"ARTIFICIAL SEQUENCE\",\n// \"SYNTHETIC\", \"SYNTHETIC CONSTRUCT\", etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Description of species from which trace is derived" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "species_code" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "mandatory", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "Species" );
    Owp . SetValue ( "vdb_name", "species_code" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_SPECIES_CODE () */



static
void
_add_description_of_STRAIN ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "strain" );

    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// strain\n// Field is dictionary from Strain\n// U32 -> ASCII. Has 1.2K values like: KT, SRS30216, LB400, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Strain from which a trace is derived" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "strain" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "Strain" );
    Owp . SetValue ( "vdb_name", "strain" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_STRAIN () */



static
void
_add_description_of_STRATEGY ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "strategy" );

    Owp . SetValue ( "can_be_ignored", "true" );
    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// strategy\n// Field is dictionary from Strategy (?) ( wrong size )\n// U8 -> ASCII. Has 38 values like: AFLP, BARCODE, CCS, CDNA, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Experimental STRATEGY" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "strategy" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "mandatory", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "Strategy" );
    Owp . SetValue ( "vdb_name", "strategy" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_STRATEGY () */



static
void
_add_description_of_SUBMISSION_TYPE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "submission_type" );

    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// submission_type" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Type of submission: new, update..." );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "submission_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "mandatory", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "SubmissionType" );
    Owp . SetValue ( "vdb_name", "submission_type" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_SUBMISSION_TYPE () */



static
void
_add_description_of_SVECTOR_ACCESSION ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "svector_accession" );

    Owp . SetValue ( "comment", "\n// *\n// svector_accession\n// Field is dictionary from SeqVectorAcc\n// U16 -> ASCII. Has 8 values like: AF399742, L08752, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "GenBank/EMBL/DDBJ accession of the sequencing vector" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "svector_accession" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "SeqVectorAcc" );
    Owp . SetValue ( "vdb_name", "svector_accession" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_SVECTOR_ACCESSION () */



static
void
_add_description_of_SVECTOR_CODE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "svector_code" );

    Owp . SetValue ( "comment", "\n// *\n// svector_code\n// Field is dictionary from SeqVector\n// U16 -> ASCII. Has 230 values like: POT, M13, PUC, POTW13, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Center defined code for the sequencing vector" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "svector_code" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "SeqVector" );
    Owp . SetValue ( "vdb_name", "svector_code" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_SVECTOR_CODE () */



static
void
_add_description_of_TAXID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "taxid" );

    Owp . SetValue ( "comment", "\n// *\n// taxid\n// Field is dictionary from TaxonomyTree\n// U32 -> ASCII. Has 4K values like: \"Bacteria\", \"Pelobacter\", etc" );
    Owp . SetValue ( "description", "NCBI Taxonomy ID" );
    Owp . SetValue ( "field_norm", "taxid" );
    Owp . SetValue ( "field_prim", "species_code" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "Species" );
    Owp . SetValue ( "vdb_name", "taxid" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TAXID () */



static
void
_add_description_of_TEMPERATURE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "temperature" );

    Owp . SetValue ( "comment", "\n// *\n// temperature" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "The temperature (in oC) at which an environmental sample was collected" );
    Owp . SetValue ( "field_prim", "temperature" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "temperature" );
    Owp . SetValue ( "vdb_type", "F64" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TEMPERATURE () */



static
void
_add_description_of_TEMPLATE_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "template_id" );

    Owp . SetValue ( "can_be_ignored", "true" );
    Owp . SetValue ( "comment", "\n// *\n// template_id" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Submitter defined identifier for the sequencing template" );
    Owp . SetValue ( "field_prim", "template_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_name", "template_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TEMPLATE_ID () */



static
void
_add_description_of_TI ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "ti" );

    Owp . SetValue ( "comment", "\n// *\n// ti" );
    Owp . SetValue ( "description", "Trace internal Identifier" );
    Owp . SetValue ( "field_prim", "ti" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_name", "ti" );
    Owp . SetValue ( "vdb_type", "U64" );
    Owp . SetValue ( "vdb_type_cast", "U32" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TI () */



static
void
_add_description_of_TRACE_DIRECTION ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "trace_direction" );

    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// trace_direction\n// Field is dictionary from Direction\n// U8 -> ASCII. Has 6 values: F, FORWARD, N, R, REVERSE, UNKNOWN" );
    Owp . SetValue ( "description", "Direction of the read" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "trace_direction" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "Direction" );
    Owp . SetValue ( "vdb_name", "trace_direction" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRACE_DIRECTION () */



static
void
_add_description_of_TRACE_END ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "trace_end" );

    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// trace_end\n// Field is dictionary from EndId\n// U8 -> ASCII. Has 6 values: F, FORWARD, N, R, REVERSE, UNKNOWN" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Defines the end of the template contained in the read" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "trace_end" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "EndId" );
    Owp . SetValue ( "vdb_name", "trace_end" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRACE_END () */



static
void
_add_description_of_TRACE_FILE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "trace_file" );

    Owp . SetValue ( "can_be_ignored", "true" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Filename with the trace" );
    Owp . SetValue ( "field_prim", "trace_file" );
    Owp . SetValue ( "mandatory", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "vdb_name", "trace_file" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRACE_FILE () */



static
void
_add_description_of_TRACE_FORMAT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "trace_format" );

    Owp . SetValue ( "can_be_ignored", "true" );
    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// trace_format\n// Field is dictionary from Format\n// U16 -> ASCII. Has 5 values: SCF, ZTR, AB1, ABI, SFF" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Format of the trace file (scf, abi...)" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "trace_format" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "mandatory", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "Format" );
    Owp . SetValue ( "vdb_name", "trace_format" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRACE_FORMAT () */



static
void
_add_description_of_TRACE_LEN ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "trace_len" );

    Owp . SetValue ( "deprecated", "true" );
    Owp . SetValue ( "description", "Description is MISSING" );
    Owp . SetValue ( "field_prim", "trace_len" );
    Owp . SetValue ( "vdb_name", "trace_len" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRACE_LEN () */



static
void
_add_description_of_TRACE_NAME ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "trace_name" );

    Owp . SetValue ( "comment", "\n// *\n// trace_name\n//  name of the submission." );
    Owp . SetValue ( "description", "Center defined trace identifier" );
    Owp . SetValue ( "field_prim", "trace_name" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "mandatory", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "vdb_name", "trace_name" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRACE_NAME () */



static
void
_add_description_of_TRACE_TYPE_CODE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "trace_type_code" );

    Owp . SetValue ( "can_be_ignored", "true" );
    Owp . SetValue ( "can_be_synonym", "true" );
    Owp . SetValue ( "comment", "\n// *\n// trace_type_code\n// Field is dictionary from  BlastDump (?)\n// There are 17 tinyints like: FINISHING, RANDOM, WGS, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Sequencing strategy by which the trace was obtained" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "trace_type_code" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "mandatory", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "static", "true" );
    Owp . SetValue ( "table_name", "TypeId" );
    Owp . SetValue ( "vdb_name", "trace_type_code" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRACE_TYPE_CODE () */



static
void
_add_description_of_TRANSPOSON_ACC ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "transposon_acc" );

    Owp . SetValue ( "comment", "\n// *\n// transposon_acc\n// Field is dictionary from TransposonAcc\n// U32 -> ASCII. Has 4 values: \"N/A\", \"NGB00024.2\", \"\", \"BD291493\"" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "GenBank/EMBL/DDBJ accession for transposon used in generating sequencing template" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "transposon_acc" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "TransposonAcc" );
    Owp . SetValue ( "vdb_name", "transposon_acc" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRANSPOSON_ACC () */



static
void
_add_description_of_TRANSPOSON_CODE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "transposon_code" );

    Owp . SetValue ( "comment", "\n// *\n// transposon_code\n// Field is dictionary from TransposonCode\n// U32 -> ASCII. Has 8 values like: AT-2, EZ-TN, EZ-TN5, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Center defined code for transposon used in generating sequencing template" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "transposon_code" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "searchable", "true" );
    Owp . SetValue ( "table_name", "TransposonCode" );
    Owp . SetValue ( "vdb_name", "transposon_code" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRANSPOSON_CODE () */



static
void
_add_description_of_WELL_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "well_id" );

    Owp . SetValue ( "comment", "\n// *\n// well_id\n// Field is dictionary from Well\n// U32 -> ASCII. Has 300K values like: C6, G9, H2, B5, C11, etc" );
    Owp . SetValue ( "common", "true" );
    Owp . SetValue ( "description", "Center defined well identifier for the sequencing reaction" );
    Owp . SetValue ( "field_norm", "name" );
    Owp . SetValue ( "field_prim", "well_id" );
    Owp . SetValue ( "for_dump", "true" );
    Owp . SetValue ( "rfc", "true" );
    Owp . SetValue ( "table_name", "Well" );
    Owp . SetValue ( "vdb_name", "well_id" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_WELL_ID () */



static
void
_add_description_of_SUBMISSION_ID ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "submission_id" );

    Owp . SetValue ( "comment", "\n// *\n// submission_id\n// Field is dictionary from SubmissionType\n// U16 -> ASCII. Has 15 names like: NEWTRACELESS, TRACELESS, etc\n// Should be canonized" );
    Owp . SetValue ( "field_prim", "submission_type" );
    Owp . SetValue ( "table_name", "SubmissionType" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_SUBMISSION_ID () */



static
void
_add_description_of_TRACE_MAX_VALUE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "trace_max_value" );

    Owp . SetValue ( "comment", "\n// *\n// trace_max_value" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_type", "U16" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRACE_MAX_VALUE () */



static
void
_add_description_of_ORGANISM_NAME_OLD ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "organism_name_old" );

    Owp . SetValue ( "comment", "\n// *\n// organism_name_OLD" );
    Owp . SetValue ( "field_prim", "organism_name" );
    Owp . SetValue ( "table_name", "OrganismName" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ORGANISM_NAME_OLD () */



static
void
_add_description_of_CONTROL_FLAGS ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "control_flags" );

    Owp . SetValue ( "comment", "\n// *\n// control_flags" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_type", "U16" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CONTROL_FLAGS () */



static
void
_add_description_of_STATUS ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "status" );

    Owp . SetValue ( "comment", "\n// *\n// status" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_type", "U8" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_STATUS () */



static
void
_add_description_of_UPDATE_DATE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "update_date" );

    Owp . SetValue ( "comment", "\n// *\n// update_date\n// smalldatetime" );
    Owp . SetValue ( "vdb_compression", "izip_encoding" );
    Owp . SetValue ( "vdb_type", "loctm" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_UPDATE_DATE () */



static
void
_add_description_of_COMMENTS ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "comments" );

    Owp . SetValue ( "comment", "\n// *\n// comments" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_COMMENTS () */



static
void
_add_description_of_CLIP_ADAPTER_LEFT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "clip_adapter_left" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CLIP_ADAPTER_LEFT () */



static
void
_add_description_of_CLIP_ADAPTER_RIGHT ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "clip_adapter_right" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_CLIP_ADAPTER_RIGHT () */



static
void
_add_description_of_HAS_EXTENDED_DATA ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "has_extended_data" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_HAS_EXTENDED_DATA () */



static
void
_add_description_of_ALT_BASECALL ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "alt_basecall" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ALT_BASECALL () */



static
void
_add_description_of_ALT_BASECALL_CHANGE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "alt_basecall_change" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ALT_BASECALL_CHANGE () */



static
void
_add_description_of_ALT_PEAK_INDEX ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "alt_peak_index" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ALT_PEAK_INDEX () */



static
void
_add_description_of_ALT_PEAK_INDEX_CHANGE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "alt_peak_index_change" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ALT_PEAK_INDEX_CHANGE () */



static
void
_add_description_of_ALT_QUALITY_SCORE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "alt_quality_score" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ALT_QUALITY_SCORE () */



static
void
_add_description_of_ALT_QUALITY_SCORE_CHANGE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "alt_quality_score_change" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_ALT_QUALITY_SCORE_CHANGE () */



static
void
_add_description_of_TRACE_TABLE_TYPE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "trace_table_type" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_TRACE_TABLE_TYPE () */



static
void
_add_description_of_POSITION_OFFSET ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "position_offset" );

    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_POSITION_OFFSET () */



static
void
_add_description_of_FLOW_CHARS ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "flow_chars" );

    Owp . SetValue ( "comment", "\n// *\n// flow_chars" );
    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_FLOW_CHARS () */



static
void
_add_description_of_KEY_SEQUENCE ( TL_OVec & Descriptions )
{
    TL_Owp Owp ( "key_sequence" );

    Owp . SetValue ( "comment", "\n// *\n// key_sequence" );
    Owp . SetValue ( "deprecated", "true" );

    Descriptions . insert ( Descriptions . end (), Owp ); 

}   /* _add_description_of_KEY_SEQUENCE () */

namespace _tl_ {

void
TL_GetFieldDescriptions ( TL_OVec & Descriptions )
{
    Descriptions . clear ();

    _add_description_of_ACCESSION ( Descriptions );
    _add_description_of_AMPLIFICATION_FORWARD ( Descriptions );
    _add_description_of_AMPLIFICATION_REVERSE ( Descriptions );
    _add_description_of_AMPLIFICATION_SIZE ( Descriptions );
    _add_description_of_ANONYMIZED_ID ( Descriptions );
    _add_description_of_ASSEMBLY_ID ( Descriptions );
    _add_description_of_ATTEMPT ( Descriptions );
    _add_description_of_BASECALL_LEN ( Descriptions );
    _add_description_of_BASES_20 ( Descriptions );
    _add_description_of_BASES_40 ( Descriptions );
    _add_description_of_BASES_60 ( Descriptions );
    _add_description_of_BASE_FILE ( Descriptions );
    _add_description_of_CENTER_NAME ( Descriptions );
    _add_description_of_CENTER_PROJECT ( Descriptions );
    _add_description_of_CHEMISTRY ( Descriptions );
    _add_description_of_CHEMISTRY_TYPE ( Descriptions );
    _add_description_of_CHIP_DESIGN_ID ( Descriptions );
    _add_description_of_CHROMOSOME ( Descriptions );
    _add_description_of_CHROMOSOME_REGION ( Descriptions );
    _add_description_of_CLIP_QUALITY_LEFT ( Descriptions );
    _add_description_of_CLIP_QUALITY_RIGHT ( Descriptions );
    _add_description_of_CLIP_VECTOR_LEFT ( Descriptions );
    _add_description_of_CLIP_VECTOR_RIGHT ( Descriptions );
    _add_description_of_CLONE_ID ( Descriptions );
    _add_description_of_CLONE_ID_LIST ( Descriptions );
    _add_description_of_COLLECTION_DATE ( Descriptions );
    _add_description_of_COMPRESS_TYPE ( Descriptions );
    _add_description_of_CVECTOR_ACCESSION ( Descriptions );
    _add_description_of_CVECTOR_CODE ( Descriptions );
    _add_description_of_DEPTH ( Descriptions );
    _add_description_of_ELEVATION ( Descriptions );
    _add_description_of_ENVIRONMENT_TYPE ( Descriptions );
    _add_description_of_EXTENDED_DATA ( Descriptions );
    _add_description_of_FEATURE_ID_FILE_NAME ( Descriptions );
    _add_description_of_FEATURE_SIGNAL_FILE_NAME ( Descriptions );
    _add_description_of_GENE_NAME ( Descriptions );
    _add_description_of_HAS_ALT_BASECALL ( Descriptions );
    _add_description_of_HAS_ALT_BASECALL_CHANGE ( Descriptions );
    _add_description_of_HAS_ALT_PEAK_INDEX ( Descriptions );
    _add_description_of_HAS_ALT_PEAK_INDEX_CHANGE ( Descriptions );
    _add_description_of_HAS_ALT_QUALITY_SCORE ( Descriptions );
    _add_description_of_HAS_ALT_QUALITY_SCORE_CHANGE ( Descriptions );
    _add_description_of_HAS_PEAK_INDEX ( Descriptions );
    _add_description_of_HAS_QUALITY_SCORE ( Descriptions );
    _add_description_of_HI_FILTER_SIZE ( Descriptions );
    _add_description_of_HOST_CONDITION ( Descriptions );
    _add_description_of_HOST_ID ( Descriptions );
    _add_description_of_HOST_LOCATION ( Descriptions );
    _add_description_of_HOST_SPECIES ( Descriptions );
    _add_description_of_INDIVIDUAL_ID ( Descriptions );
    _add_description_of_INSERT_FLANK_LEFT ( Descriptions );
    _add_description_of_INSERT_FLANK_RIGHT ( Descriptions );
    _add_description_of_INSERT_SIZE ( Descriptions );
    _add_description_of_INSERT_STDEV ( Descriptions );
    _add_description_of_ITERATION ( Descriptions );
    _add_description_of_LATITUDE ( Descriptions );
    _add_description_of_LIBRARY_ID ( Descriptions );
    _add_description_of_LOAD_DATE ( Descriptions );
    _add_description_of_LONGITUDE ( Descriptions );
    _add_description_of_LO_FILTER_SIZE ( Descriptions );
    _add_description_of_NCBI_PROJECT_ID ( Descriptions );
    _add_description_of_ORGANISM_NAME ( Descriptions );
    _add_description_of_PEAK_FILE ( Descriptions );
    _add_description_of_PH ( Descriptions );
    _add_description_of_PICK_GROUP_ID ( Descriptions );
    _add_description_of_PLACE_NAME ( Descriptions );
    _add_description_of_PLATE_ID ( Descriptions );
    _add_description_of_PMID ( Descriptions );
    _add_description_of_POPULATION_ID ( Descriptions );
    _add_description_of_PREP_GROUP_ID ( Descriptions );
    _add_description_of_PRIMER ( Descriptions );
    _add_description_of_PRIMER_CODE ( Descriptions );
    _add_description_of_PROGRAM_ID ( Descriptions );
    _add_description_of_PROJECT_NAME ( Descriptions );
    _add_description_of_QUAL_FILE ( Descriptions );
    _add_description_of_REFERENCE_ACCESSION ( Descriptions );
    _add_description_of_REFERENCE_ACC_MAX ( Descriptions );
    _add_description_of_REFERENCE_ACC_MIN ( Descriptions );
    _add_description_of_REFERENCE_OFFSET ( Descriptions );
    _add_description_of_REFERENCE_SET_MAX ( Descriptions );
    _add_description_of_REFERENCE_SET_MIN ( Descriptions );
    _add_description_of_REPLACED_BY ( Descriptions );
    _add_description_of_RUN_DATE ( Descriptions );
    _add_description_of_RUN_GROUP_ID ( Descriptions );
    _add_description_of_RUN_LANE ( Descriptions );
    _add_description_of_RUN_MACHINE_ID ( Descriptions );
    _add_description_of_RUN_MACHINE_TYPE ( Descriptions );
    _add_description_of_SALINITY ( Descriptions );
    _add_description_of_SEQ_LIB_ID ( Descriptions );
    _add_description_of_SOURCE_TYPE ( Descriptions );
    _add_description_of_SPECIES_CODE ( Descriptions );
    _add_description_of_STRAIN ( Descriptions );
    _add_description_of_STRATEGY ( Descriptions );
    _add_description_of_SUBMISSION_TYPE ( Descriptions );
    _add_description_of_SVECTOR_ACCESSION ( Descriptions );
    _add_description_of_SVECTOR_CODE ( Descriptions );
    _add_description_of_TAXID ( Descriptions );
    _add_description_of_TEMPERATURE ( Descriptions );
    _add_description_of_TEMPLATE_ID ( Descriptions );
    _add_description_of_TI ( Descriptions );
    _add_description_of_TRACE_DIRECTION ( Descriptions );
    _add_description_of_TRACE_END ( Descriptions );
    _add_description_of_TRACE_FILE ( Descriptions );
    _add_description_of_TRACE_FORMAT ( Descriptions );
    _add_description_of_TRACE_LEN ( Descriptions );
    _add_description_of_TRACE_NAME ( Descriptions );
    _add_description_of_TRACE_TYPE_CODE ( Descriptions );
    _add_description_of_TRANSPOSON_ACC ( Descriptions );
    _add_description_of_TRANSPOSON_CODE ( Descriptions );
    _add_description_of_WELL_ID ( Descriptions );
    _add_description_of_SUBMISSION_ID ( Descriptions );
    _add_description_of_TRACE_MAX_VALUE ( Descriptions );
    _add_description_of_ORGANISM_NAME_OLD ( Descriptions );
    _add_description_of_CONTROL_FLAGS ( Descriptions );
    _add_description_of_STATUS ( Descriptions );
    _add_description_of_UPDATE_DATE ( Descriptions );
    _add_description_of_COMMENTS ( Descriptions );
    _add_description_of_CLIP_ADAPTER_LEFT ( Descriptions );
    _add_description_of_CLIP_ADAPTER_RIGHT ( Descriptions );
    _add_description_of_HAS_EXTENDED_DATA ( Descriptions );
    _add_description_of_ALT_BASECALL ( Descriptions );
    _add_description_of_ALT_BASECALL_CHANGE ( Descriptions );
    _add_description_of_ALT_PEAK_INDEX ( Descriptions );
    _add_description_of_ALT_PEAK_INDEX_CHANGE ( Descriptions );
    _add_description_of_ALT_QUALITY_SCORE ( Descriptions );
    _add_description_of_ALT_QUALITY_SCORE_CHANGE ( Descriptions );
    _add_description_of_TRACE_TABLE_TYPE ( Descriptions );
    _add_description_of_POSITION_OFFSET ( Descriptions );
    _add_description_of_FLOW_CHARS ( Descriptions );
    _add_description_of_KEY_SEQUENCE ( Descriptions );

}   /* TL_GetFieldDescriptions () */

} /* namespace _tl_ */
