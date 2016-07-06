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

public class interval implements Comparable<interval>, java.io.Serializable
{
    private interval_type iv_type;
    private String ref;
    private String id;
    private long start;
    private int len;
    private boolean reverse;

    public interval( final String ref, final String id, final long start, final int len, final boolean reverse )
    {
        this.ref = ref;
        this.id = id;

        if ( id.equals( "gap" ) )
            this.iv_type = interval_type.GAP;
        else if ( id.equals( "amb" ) )
            this.iv_type = interval_type.AMBIGUOUS;
        else
            this.iv_type = interval_type.FEATURE;

        this.start = start;
        this.len = len;
        this.reverse = reverse;
    }

    public String get_ref() { return ref; }
    public boolean has_ref( final String ref ) { return this.ref.equals( ref ); }
    public String get_id() { return id; }
    public long get_start() { return start; }
    public int get_len() { return len; }
    public long get_end() { return ( start + len ) - 1; }
    public boolean get_reverse() { return reverse; }
    public interval_type get_type() { return iv_type; }

    @Override public int compareTo( final interval other )
    {
        return (int)( start - other.get_start() ); /* For Ascending order*/
    }

    @Override public String toString()
    {
        StringBuffer sb = new StringBuffer();
        sb.append( ref );
        sb.append( "\t" );
        sb.append( start );
        sb.append( "." );
        sb.append( len );
        sb.append( "\t" );
        sb.append( get_end() );
        sb.append( "\t" );
        sb.append( id );
        sb.append( "\t" );
        sb.append( reverse );
        return sb.toString();
    }
}