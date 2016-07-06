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
import java.util.*;

public class gtf_feature implements Comparable<gtf_feature>, java.io.Serializable
{
	private String id, ref;
	private boolean reverse;
    private int ranges, allocated, counter;
    long [] start;
    int [] len;
    boolean [] amb;

    public gtf_feature( int to_alloc, final String id, final String ref, boolean reverse )
    {
        this.id = id;
        this.ref = ref;
        this.reverse = reverse;
        counter = 0;
        ranges = 0;
        allocated = to_alloc;
        start = new long[ to_alloc ];
        len = new int[ to_alloc ];
        amb = new boolean[ to_alloc ];
    }

    public void add_seg( long start, int len )
    {
        if ( len > 0 && ranges < allocated )
        {
            this.start[ ranges ] = start;
            this.len[ ranges++ ] = len;
        }
    }

    public void add_seg_1( long start, int len, boolean amb )
    {
        if ( len > 0 && ranges < allocated )
        {
            this.start[ ranges ] = start;
            this.len[ ranges ] = len;
            this.amb[ ranges ++ ] = amb;
        }
    }

    public String get_id() { return id; }
    public String get_ref() { return ref; }
    public boolean has_ref( final String ref ) { return this.ref.equals( ref ); }
    public boolean is_reverse() { return reverse; }
    public long get_start() { return start[ 0 ]; }
    public long get_end() { return start[ ranges - 1 ] + len[ ranges - 1 ] - 1; }

    public int get_len()
    {
        long e = get_end();
        e -= get_start();
        e++;
        return (int)e;
    }

    public int get_ranges() { return ranges; }

    public void update_segments( final interval_list l )
    {
        int n = l.get_count();
        if ( n > allocated )
        {
            allocated = n;
            start = new long[ n ];
            len = new int[ n ];
            amb = new boolean[ n ];
        }
        ranges = 0;
        for ( int i = 0; i < n; ++i )
            add_seg_1( l.get_start( i ), l.get_len( i ), l.is_amb( i ) );
    }

    public void inc_counter() { counter++; };
    public int get_counter() { return counter; };

    public long get_start_at( final int idx )
    {
        if ( idx < ranges ) return start[ idx ]; else return 0;
    }

    public int get_len_at( final int idx )
    {
        if ( idx < ranges ) return len[ idx ]; else return 0;
    }

    public boolean get_amb_at( final int idx )
    {
        if ( idx < ranges ) return amb[ idx ]; else return false;
    }

    public void translate_ref( translater tr )
    {
        if ( tr != null )
            ref = tr.translate( ref );
    }

    public void head( StringBuffer sb, translater tr )
    {
        sb.append( id );
        sb.append( "\t" );
        if ( reverse ) sb.append( "-\t" ); else sb.append( "+\t" );
        if ( tr != null )
            sb.append( tr.translate( ref ) );
        else
            sb.append( ref );
        sb.append( "\t" );
    }

    @Override public int compareTo( final gtf_feature other )
    {
        return (int)( get_start() - other.get_start() ); /* For Ascending order*/
    }

    @Override public String toString()
    {
        StringBuffer sb = new StringBuffer();
        head( sb, null );
        for ( int i = 0; i < ranges; ++ i )
        {
            if ( i > 0 ) sb.append( ";" );
            if ( amb[ i ] ) sb.append( "A" ); else sb.append( "N" );
            sb.append( start[ i ] );
            sb.append( "." );
            sb.append( len[ i ] );
        }
        return sb.toString();
    }
}