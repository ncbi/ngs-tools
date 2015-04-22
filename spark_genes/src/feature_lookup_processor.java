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
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.*;
import org.apache.spark.api.java.function.*;

import java.util.*;
import java.io.*;
import ngs.*;
import scala.Tuple2;

class node_job implements java.io.Serializable
{
    long start;
    long count;

    node_job( long start, long count )
    {
        this.start = start;
        this.count = count;
    }
}

class node_job_divider
{
    private final long start;
    private final long count;

    node_job_divider( final long start, final long count )
    {
        this.start = start;
        this.count = count;
    }

    ArrayList< node_job > run( long sections )
    {
        ArrayList< node_job > res = new ArrayList< node_job >();
        long n = count;
        long s = 1;
        long per_section = ( n / sections ) + sections;
        while ( n > 0 )
        {
            if ( n < per_section ) per_section = n;
            res.add( new node_job( s, per_section ) );
            s += per_section;
            n -= per_section;
        }
        return res;
    }
}


class alignment_status
{
    String id;
    int found_features;

    int counts;
    int ambiguous;
    int mapq_too_low;
    int no_feature;

    void clear()
    {
        id = null;
        found_features = counts = ambiguous = mapq_too_low = no_feature = 0;
    }
}

class alignment_checker
{
    private final gtf_lookup lookup;
    private final count_mode cm;
    private final int min_mapq;
    private final boolean ignore_orientation;
	
    interval_list feature_intervals;
    interval_list alignment_intervals;
    interval_list_cmp comparator;
    List< gtf_feature > found;

    alignment_checker( final gtf_lookup lookup, final count_mode cm, final int min_mapq,
                       final boolean ignore_orientation )
    {
        this.lookup = lookup;
        this.cm = cm;
        this.min_mapq = min_mapq;
        this.ignore_orientation = ignore_orientation;
		
        feature_intervals = new interval_list( 512 );
        alignment_intervals = new interval_list( 512 );
        comparator = new interval_list_cmp( false );
        found = new ArrayList< gtf_feature >();
    }

    void check( AlignmentIterator alig_iter, alignment_status status )
    {
        status.clear();
        try
        {
            if ( min_mapq < 1 || alig_iter.getMappingQuality() >= min_mapq )
            {
                long start = alig_iter.getAlignmentPosition();
                long len = alig_iter.getAlignmentLength();
				String ref_name = alig_iter.getReferenceSpec(); /* how to ask for canonical? */
				
                status.found_features = lookup.lookup( ref_name, start, len, found );
                if ( status.found_features > 0 )
                {
                    long n_feature_counts = 0;
                    long n_feature_ambiguous = 0;
                    boolean alignment_reverse = ignore_orientation ? false : alig_iter.getIsReversedOrientation();

                    alignment_intervals.clear();
                    alignment_intervals.from_cigar( start, alig_iter.getShortCigar( false ) );

                    for ( gtf_feature f : found )
                    {
                        boolean feature_reverse = ignore_orientation ? false : f.is_reverse();
                        if ( alignment_reverse == feature_reverse )
                        {
                            int n = f.get_ranges();
                            feature_intervals.clear();
                            for ( int i = 0; i < n; ++i )
                                feature_intervals.add( f.get_start_at( i ), f.get_len_at( i ), f.get_amb_at( i ) );

                            count_result cr = comparator.walk( feature_intervals, alignment_intervals, cm );
                            switch( cr )
                            {
                                case FEATURE        :   status.id = f.get_id(); status.counts++; break;
                                case NO_FEATURE     :   break;
                                case AMBIGUOUS      :   status.ambiguous++; break;
                            }
                        }
                    }
                }
            }
            else
                status.mapq_too_low++;
        }
        catch ( ErrorMsg err )
        {
            System.err.println( err.toString() );
        }
    }
}


class alignment_loop
{
    private final boolean show_progress;
    private alignment_checker checker;
    private alignment_status status;

    alignment_loop( final gtf_lookup lookup, final count_mode cm, final int min_mapq,
                    final boolean ignore_orientation, final boolean show_progress )
    {
        
        this.show_progress = show_progress;
        status = new alignment_status();
        checker = new alignment_checker( lookup, cm, min_mapq, ignore_orientation );
    }


    lookup_counters run( AlignmentIterator alig_iter, feature_counts ftc )
    {
        lookup_counters res = new lookup_counters();
        try
        {
            while ( alig_iter.nextAlignment() )
            {
                if ( checker == null )
                {
                    System.err.println( "checker is null" );
                    return res;
                }

                if ( alig_iter == null )
                {
                    System.err.println( "alig_iter is null" );
                    return res;
                }

                if ( status == null )
                {
                    System.err.println( "status is null" );
                    return res;
                }

                checker.check( alig_iter, status );

                if ( status.mapq_too_low > 0 )
                    res.mapq_too_low( 1 );
                else if ( status.ambiguous > 0 )
                    res.ambiguous( 1 );
                else if ( status.found_features == 0 )
                    res.no_feature( 1 );
                else if ( status.counts == 1 )
                {
                    if ( status.id != null && ftc != null )
                        ftc.count( status.id, 1 );
                    res.counted( 1 );
                }
                else
                    res.not_counted( 1 );

                if ( res.inc_alignments_and_check( 10000 ) )
                {
                    if ( show_progress ) System.out.print( "." );
                }
            }
        }
        catch ( ErrorMsg err )
        {
            System.err.println( err.toString() );
        }

        return res;
    }
}


class subset_thread implements Runnable
{
	private Thread thread;
    private final SparkGenesOpt opt;
    private final gtf_lookup lookup;
    private final long first_row;
    private final long rows;
    lookup_counters result;
    feature_counts ftc;

    subset_thread( String thread_id, final gtf_lookup lookup, feature_counts ftc,
                   final SparkGenesOpt opt, final long first_row, final long rows )
    {
        thread = new Thread( this, thread_id );
        this.lookup = lookup;
        this.ftc = ftc;
        this.opt = opt;
        this.first_row = first_row;
        this.rows = rows;
        this.result = null;
    }

    public void start() { thread.start(); }

	public void run()
    {
        try
        {
            ReadCollection ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( opt.get_Accession() );
            AlignmentIterator ngs_iter = ngs_run.getAlignmentRange ( first_row, rows, Alignment.all );

            alignment_loop loop = new alignment_loop( lookup, opt.get_CountMode(), opt.get_MinMapq(),
                                                       opt.get_IgnoreOrientation(), opt.get_ShowProgress() );
            result = loop.run( ngs_iter, ftc );
        }
        catch ( ErrorMsg err )
        {
            System.out.println( err.toString() );
        }
	}

    public void join( lookup_counters global )
    {
        try
        {
            thread.join();
            if ( result != null ) global.add( result );
        }
        catch ( InterruptedException e )
        {
            System.out.println( e.toString() );
        }
    }
}


// node_job ---> [ Feature-ID / AMB ]
class node_job_mapper implements FlatMapFunction < node_job, String >
{
	private final String accession;
    private final gtf_lookup lookup;
    private final count_mode cm;
    private final int min_mapq;
    private final boolean ignore_orientation;

    private alignment_checker checker;
    private alignment_status status;
	private ReadCollection ngs_run;
    private boolean initialized;
    private LinkedList< String > iter_res;

	public node_job_mapper( final String accession, final gtf_lookup lookup, final int min_mapq,
                            final count_mode cm, final boolean ignore_orientation )
	{
		this.accession = accession;
        this.lookup = lookup;
		this.min_mapq = min_mapq;
		this.cm = cm;
        this.ignore_orientation = ignore_orientation;

        this.initialized = false;
		this.ngs_run = null;
        this.checker = null;
        this.status = null;
        this.iter_res  =null;
	}

    private boolean initialize()
    {
        boolean res = false;
        try
        {
            ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( accession );
            status = new alignment_status();
            checker = new alignment_checker( lookup, cm, min_mapq, ignore_orientation );
            iter_res = new LinkedList< String >();
            res = true;
        }
        catch ( ErrorMsg err )
        {
            System.out.println( err.toString() );
        }
        return res;
    }

	public Iterable< String > call( node_job job )
	{
        if ( !initialized ) initialized = initialize();

        if ( initialized )
        {
            iter_res.clear();
            try
            {
                AlignmentIterator ngs_iter = ngs_run.getAlignmentRange ( job.start, job.count, Alignment.all );
                while ( ngs_iter.nextAlignment() )
                {
                    checker.check( ngs_iter, status );

                    if ( status.mapq_too_low > 0 )
                        iter_res.add( "LOW" );
                    else if ( status.ambiguous > 0 )
                        iter_res.add( "AMB" );
                    else if ( status.found_features == 0 )
                        iter_res.add( "NOFT" );
                    else if ( status.counts == 1 )
                        iter_res.add( status.id );
                    else
                        iter_res.add( "NC" );
                }
            }
            catch ( ErrorMsg err )
            {
                System.out.println( err.toString() );
            }
        }
        return iter_res;
    }
}

public class feature_lookup_processor
{
    private final gtf_lookup lookup;
	
    feature_lookup_processor( Iterable< gtf_feature > features, long bin_size, ref_cmp_res refs )
    {
        lookup = new gtf_lookup( features, bin_size, refs );
    }

    void report() { lookup.report(); }

    long get_alignment_count( final String acc )
    {
        long res = 0;
        try
        {
            ReadCollection ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( acc );
            res = ngs_run.getAlignmentCount();
        }
        catch ( ErrorMsg err )
        {
            System.err.println( err.toString() );
        }
        return res;
    }

    lookup_counters run_sequentially( feature_counts ftc, SparkGenesOpt opt )
    {
		if ( opt.get_ShowProgress() )
			System.out.println( "Counting features in sequentialls on local machine" );

        lookup_counters res = null;
        try
        {
            ReadCollection ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( opt.get_Accession() );
            AlignmentIterator alig_iter = ngs_run.getAlignments( Alignment.all );
            alignment_loop loop = new alignment_loop( lookup, opt.get_CountMode(), opt.get_MinMapq(),
                                                       opt.get_IgnoreOrientation(), opt.get_ShowProgress() );
            res = loop.run( alig_iter, ftc );
        }
        catch ( ErrorMsg err )
        {
            System.err.println( err.toString() );
        }
        if ( opt.get_ShowProgress() ) System.out.println( "-" );
        return res;
    }

    lookup_counters run_in_threads( feature_counts ftc, SparkGenesOpt opt )
    {
		int n_threads = opt.get_NumSlices();

		if ( opt.get_ShowProgress() )
			System.out.println( "Counting features in parallel threads on local machine in " + n_threads + " threads" );
	
        lookup_counters res = new lookup_counters();
        ArrayList< subset_thread > threads = new ArrayList< subset_thread >();
        node_job_divider d = new node_job_divider( 1, get_alignment_count( opt.get_Accession() ) );
        ArrayList< node_job > jobs = d.run( n_threads );
        for ( node_job j : jobs )
        {
            subset_thread t = new subset_thread( String.format( "t%d", j.start ), lookup, ftc, opt, j.start, j.count );
            threads.add( t );
            t.start();
        }
        for ( subset_thread r : threads ) r.join( res );
        if ( opt.get_ShowProgress() ) System.out.println( "-" );
        return res;
    }


    lookup_counters run_sparked( feature_counts ftc, SparkGenesOpt opt )
    {
        lookup_counters res = new lookup_counters();

        int slices = opt.get_NumSlices();
        long alig_count = get_alignment_count( opt.get_Accession() );
		long n_jobs = ( alig_count / 1500 );
		
        node_job_divider d = new node_job_divider( 1, alig_count );
        ArrayList< node_job > jobs = d.run( n_jobs );

        if ( opt.get_ShowProgress() )
            System.out.println( jobs.size() + " jobs scheduled" );

		SparkConf sc = new SparkConf();
        sc.setAppName( "SparkGenes" );

        if ( opt.get_SparkMaster() != null )
            sc.setMaster( opt.get_SparkMaster() );

		JavaSparkContext jsc = new JavaSparkContext( sc );

        JavaRDD<node_job> rdd_jobs = jsc.parallelize( jobs, slices );

        node_job_mapper njm = new node_job_mapper( opt.get_Accession(), lookup,
                                   opt.get_MinMapq(), opt.get_CountMode(),
								   opt.get_IgnoreOrientation() );

        JavaRDD<String> rdd_counts = rdd_jobs.flatMap( njm ); // <--- here is the work done!

		Map< String, Long > counters = rdd_counts.countByValue(); // let spark count the values

        res.inc_alignments( alig_count );
        for ( Map.Entry< String, Long > entry : counters.entrySet() )
        {
            String key = entry.getKey();
            long value = entry.getValue();
            if ( key.equals( "AMB" ) )
                res.ambiguous( value );
            else if ( key.equals( "NOFT" ) )
                res.no_feature( value );
            else if ( key.equals( "LOW" ) )
                res.mapq_too_low( value );
            else if ( key.equals( "NC" ) )
                res.not_counted( value );
            else
            {
                ftc.count( key, value );
                res.counted( value );
            }
        }

		sc.stop();
		
        return res;
    }
}