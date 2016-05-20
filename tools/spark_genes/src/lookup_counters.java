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

import java.io.*;

public class lookup_counters
{
    private long n_alignments;
    private long n_alig_counted;
    private long n_alig_not_counted;
    private long n_alig_ambiguous;
    private long n_alig_no_features;
    private long n_alig_mapq_low;
    private long n_unaligned;

    public lookup_counters()
    {
        clear();
    }

    public lookup_counters( final lookup_counters other )
    {
        n_alignments = other.n_alignments;
        n_alig_counted = other.n_alig_counted;
        n_alig_not_counted = other.n_alig_not_counted;
        n_alig_ambiguous = other.n_alig_ambiguous;
        n_alig_no_features = other.n_alig_no_features;
        n_alig_mapq_low = other.n_alig_mapq_low;
        n_unaligned = other.n_unaligned;
    }

    public void add( final lookup_counters other )
    {
        n_alignments += other.n_alignments;
        n_alig_counted += other.n_alig_counted;
        n_alig_not_counted += other.n_alig_not_counted;
        n_alig_ambiguous += other.n_alig_ambiguous;
        n_alig_no_features += other.n_alig_no_features;
        n_alig_mapq_low += other.n_alig_mapq_low;
        n_unaligned += other.n_unaligned;
    }

    public void clear()
    {
        n_alignments = 0;
        n_alig_counted = 0;
        n_alig_not_counted = 0;
        n_alig_ambiguous = 0;
        n_alig_no_features = 0;
        n_alig_mapq_low = 0;
        n_unaligned = 0;
    }

    public boolean equals( final lookup_counters other )
    {
        return (
            n_alignments == other.n_alignments &&
            n_alig_counted == other.n_alig_counted &&
            n_alig_not_counted == other.n_alig_not_counted &&
            n_alig_ambiguous == other.n_alig_ambiguous &&
            n_alig_no_features == other.n_alig_no_features &&
            n_alig_mapq_low == other.n_alig_mapq_low &&
            n_unaligned == other.n_unaligned );
    }

    String to_String()
    {
        long total = n_alig_counted + n_alig_not_counted + n_alig_ambiguous + n_alig_no_features + n_alig_mapq_low;
        StringBuffer sb = new StringBuffer();
        sb.append( String.format( "alignments = %,d\n", n_alignments ) );
        sb.append( String.format( "alignment counts = %,d\n", n_alig_counted ) );
        sb.append( String.format( "alignment does not count = %,d\n", n_alig_not_counted ) );
        sb.append( String.format( "alignment ambiguous = %,d\n", n_alig_ambiguous ) );
        sb.append( String.format( "alignment has not features = %,d\n", n_alig_no_features ) );
        sb.append( String.format( "alignment with mapq too low = %,d\n", n_alig_mapq_low ) );
        if ( n_unaligned > 0 )
            sb.append( String.format( "unaligned reads = %,d\n", n_unaligned ) );
        sb.append( String.format( "total processed = %,d", total ) );
        return sb.toString();
    }

    public void report()
    {
        System.out.println( to_String() );
    }

    public void save_to_file( PrintWriter out )
    {
        out.println( to_String() );
    }

    public void inc_alignments( long value ) { n_alignments += value; }
    public boolean inc_alignments_and_check( int mod_value ) { return ( ++n_alignments % mod_value == 0 ); }

    public void counted( long value ) { n_alig_counted += value; }
    public void not_counted( long value ) { n_alig_not_counted += value; }
    public void ambiguous( long value ) { n_alig_ambiguous += value; }
    public void no_feature( long value ) { n_alig_no_features += value; }
    public void mapq_too_low( long value ) { n_alig_mapq_low += value; }
    public void unaligned( long value ) { n_unaligned += value; }
}
