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
package Bio;

import data.CLogger;
import job.BioRunState;
import ngs.*;

public class BioDumper implements Runnable
{
    protected final BioStats stats;
    protected final BioFormatter fmt;
    protected final BioFilterSet filter_set;
    protected final String name;
    
    protected BioRunState rs;
    private long read_count_or_size;
    
    public void set_read_count_or_size( final long value ) { read_count_or_size = value; }
    public long get_read_count_or_size() { return read_count_or_size; }
    public void set_run_state( final BioRunState rs ) { this.rs = rs; }
    
    private void sleep_for( int time_to_sleep )
    {
        try
        {
            Thread.sleep( time_to_sleep );
        }
        catch ( InterruptedException iex )
        {
            CLogger.logfmt( "job >%s< error: %s" , name , iex.toString() );
        }
    }

    /* used by all derived dumpers to get a BioRecord to work with */
    public BioRead get_empty_bio_read()
    {
        BioRead res = null;
        if ( rs != null )
        {
            while ( res == null )
            {
                res = rs.get_from_stock();
                if ( res == null ) sleep_for( 200 );
            }
        }
        else
            res = new BioRead();
        return res;
    }

    public BioRecord get_empty_bio_record()
    {
        BioRecord res = null;
        if ( rs != null )
        {
            while ( res == null )
            {
                res = rs.get_empty_record();
                if ( res == null ) sleep_for( 200 );
            }
        }
        else
            res = new BioRecord();
        return res;
    }
    
    /* called by run, overwritten by derived dumpers */    
    public boolean produce( BioRead read ) { return false; }
    
    /* called by run, overwritten by derived dumpers */    
    public boolean next() { return false; }
            
    public void add_filter( BioFilter filter ) { filter_set.addFilter( filter ); }
    public BioStats get_stats() { return stats; }
    
    public boolean pass( ReadIterator iter, String bases )
    {
        boolean res = true;
        if ( filter_set != null )
            res = filter_set.pass( stats, iter, bases );
        return res;
    }

    public boolean pass( AlignmentIterator iter, String bases )
    {
        boolean res = true;
        if ( filter_set != null )
            res = filter_set.pass( stats, iter, bases );
        return res;
    }

    public boolean pass( ReadIterator rd_iter, FragmentIterator frag_iter, String bases )
    {
        boolean res = true;
        if ( filter_set != null )
        {
            res = filter_set.pass( stats, rd_iter, frag_iter, bases );
        }
        return res;
    }
    
    @Override public void run()
    {
        if ( rs != null )
        {
            rs.set_running( true );
            while ( rs.is_running() )
            {
                if ( next() )
                {
                    BioRead read = get_empty_bio_read();
                    if ( produce( read ) )
                        rs.put( read );
                    else
                    {
                        rs.put_to_stock( read );
                        rs.set_running( false );
                    }
                }
                else
                    rs.set_running( false );
            }
        }
    }
    
    public BioDumper( final BioFormatter fmt, final String name )
    {
        this.rs = null;
        this.name = name;
        this.stats = new BioStats();
        this.filter_set = new BioFilterSet();
        this.fmt = fmt;
        this.fmt.set_stats( this.stats );
        this.read_count_or_size = 0;
    }
}
