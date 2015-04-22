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

public class total_count
{
    private int total;
    private int features;
    private int no_features;
    private int ambiguous;
    private int low_mapq;
    private int other_strand;

	public total_count()
	{
        total = features = no_features = ambiguous = low_mapq = other_strand = 0;
	}

    public void count( count_result cr )
    {
        total++;
        switch( cr )
        {
            case FEATURE        : features++; break;
            case NO_FEATURE     : no_features++; break;
            case AMBIGUOUS      : ambiguous++; break;
        }
    }

    public void rejected_mapq()
    {
        total++;
        low_mapq++;
    }

    public void on_other_strand()
    {
        total++;
        other_strand++;
    }

    public void add( total_count other )
    {
        total += other.total;
        features += other.features;
        no_features += other.no_features;
        ambiguous += other.ambiguous;
        low_mapq += other.low_mapq;
        other_strand += other.other_strand;
    }

    @Override public String toString()
    {
        StringBuffer sb = new StringBuffer();
        sb.append( "\ntotal = " );
        sb.append( total );
        sb.append( "\nfeatures = " );
        sb.append( features );
        sb.append( "\nno_features = " );
        sb.append( no_features );
        sb.append( "\nambiguous = " );
        sb.append( ambiguous );
        sb.append( "\nlow mapq = " );
        sb.append( low_mapq );
        sb.append( "\nother strand = " );
        sb.append( other_strand );

        return sb.toString();
    }
}
