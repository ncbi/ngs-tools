/* 
 * ===========================================================================
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
 * Author:  Alex Kotliarov
 *
 * File Description:
 */
#include <string>
#include <utility>

namespace libwgmlst
{

//
// CBlastHit represents `BLASTn` alignment.
// Initialize it with values that conform 
// to BLASTn tabular report:
// - 1-based offsets
// - subject's (start, stop) pair must be reversed if 
// hit is on '-' strand.
//
class CBlastHit
{
public:
    enum class CStrand {plus=1, minus};

private:
    std::string m_Query;
    std::string m_Subject;
    int m_Length; // length of alignment
    std::pair<int, int> m_QueryRange;
    std::pair<int, int> m_SubjectRange;
    int m_SubjectEnd;
    double m_IdentityPercent;

    std::string m_Btop;
    CStrand m_Strand;

public:
    CBlastHit() = delete;
    CBlastHit(std::string const& query,
              std::string const& subject,
              int align_len,
              int qry_start,
              int qry_end,
              int subj_start,
              int subj_end,
              double identity_pct,
              std::string const& btop,
              CStrand strand)
    : m_Query(query),
      m_Subject(subject),
      m_Length(align_len),
      m_QueryRange(std::make_pair(qry_start, qry_end)),
      m_SubjectRange(std::make_pair(subj_start, subj_end)),
      m_IdentityPercent(identity_pct),
      m_Btop(btop),
      m_Strand(strand)
    {
    }
       
    std::string const& Query() const { return m_Query; }
    std::string const& Subject() const { return m_Subject; }
    int Length() const { return m_Length; }
    std::pair<int, int> const& QueryRange() const { return m_QueryRange; }
    std::pair<int, int> const& SubjectRange() const { return m_SubjectRange; }
    double IdentityPct() const { return m_IdentityPercent; }
    std::string const& Btop() const { return m_Btop; }
    CStrand Strand() const { return m_Strand; }

};

class CBlastHitsSource
{
public:
    virtual ~CBlastHitsSource() {}

    virtual void First() = 0;
    virtual void Next() = 0;
    explicit virtual operator bool () const = 0;    
    virtual CBlastHit Get() const = 0; 
};

} // namespace libwgmlst
