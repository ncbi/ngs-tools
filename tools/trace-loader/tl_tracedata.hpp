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
* ====================================================================/======
*
*/

#ifndef _tl_tracedata_
#define _tl_tracedata_

#include <string>
#include <map>
#include <list>

#include "tl_exception.hpp"
#include "tl_ows.hpp"
#include "tl_types.hpp"
#include "tl_util.hpp"
#include "tl_traceconfig.hpp"
#include "tl_trace_dp.hpp"

struct KDirectory;
struct KFile;
struct KMMap;

namespace _tl_ {

class TL_TraceFile;

class TL_TraceData {
public:
    TL_TraceData ();
    ~TL_TraceData ();

    void Init ( TL_TraceConfig * Config );
    void Dispose ();
    inline bool IsGood () const { return _config != NULL; }

    void LoadFile ( TL_TraceFile & File, const std::string & Path ) const;
    bool FileExists ( const std::string & Path ) const;
    uint64_t FileSize ( const std::string & Path ) const;

private:
    TL_TraceConfig * _config;
};

class TL_TraceFile {
public :
    typedef enum {
                UNK,
                SCF,        /* +++ We support it */
                ZTR,        /* +++ We support it */
                ABI,        /* +++ We support it */
                RCF,
                ALF,
                PLN,
                EXP,
                CTF,
                TTFF,
                SFF         /* +++ We do not support, but ... */
    } FType;

public :
    TL_TraceFile ();
    TL_TraceFile ( const struct KFile * File );
    ~TL_TraceFile ();

    void Init ( const struct KFile * File );
    void Dispose ();

    inline const void * Content () const { return _content; };
    inline size_t ContentSize () const { return _content_size; };

    inline FType Type () const { return _type; };
    inline const std :: string & Description () const { return _description; };

private :
    void _detectFileType ();

    void * _content;
    size_t _content_size;
    bool _allocated_content;

    const struct KMMap * _map;

    FType _type;
    std :: string _description;
};

/*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*
_*
*_ Lyrics: There are 4 files: base, qual, peak and trace
_*
*_ Base : contains actual fasta data , but with crammed first line
_* EXAMPLE
*_ >1135760112141
_* CTGCGGACGTTTAGATGATATAAATCCTTATTTTCTCTTCATAGATGTACCCATACAGGC
*_ AGCAATTTCAACAACATTCCCATACACCGGTGTTCCCCCTTACTCCCATGGAACGGGAAC
_* AGGCCACACAATAGACACCGTGATCAGAACACATGAGTACTCAAATAAGGGAAAACAATA
*_
_* Qual : contains actual quality scores, not like in FASTQ file
*_ EXAMPLE
_* >1135760112141
*_ 04 05 04 03 04 06 13 07 05 27 46 16 24 12 16 12 33
_* 57 20 27 36 34 47 27 41 57 45 53 42 45 53 57 57 57
*_ 53 61 61 61 61 15 13 33 38 61 32 43 61 61 61 61 61
_*
*_ Peak : actual peak index
_* EXAMPLE
*_ >1135760112141
_* 0002 0014 0022 0033 0053 0060 0081 0102 0116 0128 0141 0156
*_ 0167 0182 0196 0205 0219 0235 0251 0260 0271 0282 0293 0304
_* 0317 0331 0344 0356 0371 0382 0393 0407 0420 0433 0445 0457
*_
 *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*/

    /*  Basecalls, Qualities and Peak index
     */
class TL_Bases {
public :
    TL_Bases ();
    ~TL_Bases ();


    void Reset ();

        /*  That method will prepare data to store to VDB
         */
    void Prepare ();

        /*  Basecall editing.
         */
    inline uint32_t NumBases () const { return _bases . size (); };
    inline const char_a_t & Bases () const { return _bases; };
    inline char_a_t & Bases () { return _bases; };

        /*  Peak index editing.
         */
    inline const uint32_a_t & PeakIndex () const
                                            { return _peak_index; };
    inline uint32_a_t & PeakIndex () { return _peak_index; };

        /*  QualScores editing:
         */
    inline uchar_a_t & ProbA () { return _prob_A; };
    inline uchar_a_t & ProbC () { return _prob_C; };
    inline uchar_a_t & ProbG () { return _prob_G; };
    inline uchar_a_t & ProbT () { return _prob_T; };
    inline void SetValidScores ( bool Val ) { _valid_scores = Val; };

    inline const uchar_a_t & ProbSub () const { return _prob_sub; };
    inline uchar_a_t & ProbSub () { return _prob_sub; };
    inline const uchar_a_t & ProbIns () const { return _prob_ins; };
    inline uchar_a_t & ProbIns () { return _prob_ins; };
    inline const uchar_a_t & ProbDel () const { return _prob_del; };
    inline uchar_a_t & ProbDel () { return _prob_del; };
    inline void SetExtraProbs ( bool Val ) { _extra_probs = Val; };


    inline const uchar_a_t & ProbScores () const { return _prob; };
    inline uint16_t Bases20 () const { return _bases_20; };
    inline uint16_t Bases40 () const { return _bases_40; };
    inline uint16_t Bases60 () const { return _bases_60; };

        /*  Reading auxiraly files
         */
    void ReadBaseFile ( const TL_TraceFile & File );
    void ReadQualFile ( const TL_TraceFile & File );
    void ReadPeakFile ( const TL_TraceFile & File );

private:

        /*  Bases
         */
    char_a_t _bases;

        /*  Peak Index
         */
    uint32_a_t _peak_index;

        /*  Quality scores
         */
    uchar_a_t _prob_A;
    uchar_a_t _prob_C;
    uchar_a_t _prob_G;
    uchar_a_t _prob_T;
    bool _valid_scores;

    uchar_a_t _prob_sub;
    uchar_a_t _prob_ins;
    uchar_a_t _prob_del;
    bool _extra_probs;  /* for scf version 3.10 */

    uchar_a_t _prob;
    uint16_t _bases_20;
    uint16_t _bases_40;
    uint16_t _bases_60;
};

class TL_Samples {
public :
    TL_Samples ();
    TL_Samples ( const TL_Samples & Samples );
    ~TL_Samples ();

    TL_Samples & operator = ( const TL_Samples & Samples );

    void Reset ();

        /*  That method will prepare data to store to VDB
         */
    void Prepare ();

    inline size_t Size () const { return _sample_A . size (); };

    inline const uint16_a_t & SampleA () const { return _sample_A; };
    inline uint16_a_t & SampleA () { return _sample_A; };

    inline const uint16_a_t & SampleC () const { return _sample_C; };
    inline uint16_a_t & SampleC () { return _sample_C; };

    inline const uint16_a_t & SampleG () const { return _sample_G; };
    inline uint16_a_t & SampleG () { return _sample_G; };

    inline const uint16_a_t & SampleT () const { return _sample_T; };
    inline uint16_a_t & SampleT () { return _sample_T; };

    inline uint16_t MaxTraceVal () const { return _max_trace_val; };

    inline const uint64_a_t & SampleCombined () const
                                        { return _sample_combined; };

    inline const char_a_t & FlowChars () const { return _flow_chars; };
    inline char_a_t & FlowChars () { return _flow_chars; };

    inline const char_a_t & KeySequence () const
                            { return _key_sequence; };
    inline char_a_t & KeySequence ()
                            { return _key_sequence; };

private :

    uint16_a_t _sample_A;
    uint16_a_t _sample_C;
    uint16_a_t _sample_G;
    uint16_a_t _sample_T;

    uint64_a_t _sample_combined;

    uint16_t _max_trace_val;

    char_a_t _flow_chars;
    char_a_t _key_sequence;
};

class  TL_Traces {
public :
    TL_Traces ();
    ~TL_Traces ();

    void ReadTraces (
                    const std::string &Name,
                    const std::string & ProgramId,
                    const TL_TraceFile & File
                    );

    void Reset ();

        /*  That method will prepare data to store to VDB
         */
    void Prepare ();

    inline const std::string & Name () const { return _name; };
    inline const std::string & ProgramID () const { return _program_id; };

    inline const TL_Bases & Bases () const { return _bases; };
    inline TL_Bases & Bases () { return _bases; };

    inline const TL_Samples & Samples () const { return _samples; };
    inline TL_Samples & Samples () { return _samples; };

    inline const std::string & Comments () const
                                        { return _comments; };
    inline void SetComments ( const std::string & S )
                                        { _comments = S; };

    inline const std::string & ExtendedData () const
                                        { return _extended_data; };
    inline void SetExtendedData ( const std::string & S )
                                        { _extended_data = S; };

        /*  If a clipping value is not computed,
         *  the field should be set to 0
         */
    inline uint32_t ClipQualityRight () const 
                                        { return _clip_quality_right; };
    inline void SetClipQualityRight ( uint32_t Clip )
                                        { _clip_quality_right = Clip; };
    inline uint32_t ClipQualityLeft () const 
                                        { return _clip_quality_left; };
    inline void SetClipQualityLeft ( uint32_t Clip )
                                        { _clip_quality_left = Clip; };

        /*  If a clipping value is not computed,
         *  the field should be set to 0
         */
        /* The Same as a CLIP_VECTOR_(RIGHT/LEFT) but poppier */
    inline uint32_t ClipAdapterRight () const 
                                        { return _clip_adapter_right; };
    inline void SetClipAdapterRight ( uint32_t Clip )
                                        { _clip_adapter_right = Clip; };
    inline uint32_t ClipAdapterLeft () const 
                                        { return _clip_adapter_left; };
    inline void SetClipAdapterLeft ( uint32_t Left )
                                        { _clip_adapter_right = Left; };

    inline const uchar_a_t & PrivateData () const
                                        { return _private_data; };

    inline void SetPrivateData ( const uchar_a_t & D )
                                        { _private_data = D; };

    inline bool IsOverrideQuals () const { return _override_quals; };
    inline void OverrideQuals ( bool Over ) { _override_quals = Over; };

    inline bool IsOverridePeaks () const { return _override_peaks; };
    inline void OverridePeaks ( bool Over ) { _override_peaks = Over; };

private :
    std::string _name;
    std::string _program_id;

    uint32_t _clip_quality_right;
    uint32_t _clip_quality_left;

    uint32_t _clip_adapter_right;
    uint32_t _clip_adapter_left;


    bool _override_peaks;
    bool _override_quals;

    int _format;

    TL_Bases _bases;
    TL_Samples _samples;

    std::string _comments;
    std::string _extended_data;
    uchar_a_t _private_data;
};

    /*  Simple abstract class with pure virtual methods
     *  It will be feeded to General Writer to retrieve data to write
     */
class TL_TraceInfo;

class TL_Trace_DP : public TL_TraceDataProvider {
public:
    TL_Trace_DP ();
    virtual ~TL_Trace_DP ();

    void Set (
            const TL_TraceData & TraceData,
            const TL_TraceInfo & TraceInfo,
            size_t RowIdx
            );
    void Reset ();

    const std::string & Name () const;
    const std::string & ProgramID () const;
    const char_a_t & Bases () const;
    const uint32_a_t & PeakIndex () const;
    const uchar_a_t & ProbScores () const;
    uint16_t Bases20 () const;
    uint16_t Bases40 () const;
    uint16_t Bases60 () const;
    const uint64_a_t & SampleCombined () const;
    uint16_t MaxTraceVal () const;
    uint32_t ClipQualityRight () const;
    uint32_t ClipQualityLeft () const;
    const char_a_t & FlowChars () const;
    const char_a_t & KeySequence () const;
    const std::string & Comments () const;
    const std::string & ExtendedData () const;

private:
    TL_Bases _bases;
    TL_Traces _traces;

    std::string _name;
    std::string _program_id;
};

}   /* namespace _tl_ */

#endif /* _tl_tracedata_ */
