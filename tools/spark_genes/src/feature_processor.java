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
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.api.java.function.VoidFunction;

import java.util.*;
import java.io.*;
import ngs.*;
import scala.Tuple2;

class feature_counter
{
    private Reference ngs_ref;
    private SparkGenesOpt opt;
    private interval_list_cmp cmp;
    private interval_list read;
    private total_count tc;

    public feature_counter( Reference ngs_ref, SparkGenesOpt opt, total_count tc )
    {
        this.ngs_ref = ngs_ref;
        this.opt = opt;
        this.cmp = new interval_list_cmp( false );
        this.read = new interval_list( 512 );
        this.tc = tc;
    }

    public void process( gtf_feature feature )
    {
        interval_list ft_intervals = new interval_list( feature.get_ranges() + 5 );
        int n = feature.get_ranges();
        for ( int i = 0; i < n; ++i )
            ft_intervals.add( feature.get_start_at( i ), feature.get_len_at( i ), feature.get_amb_at( i ) );

        try
        {
            AlignmentIterator ngs_iter = ngs_ref.getAlignmentSlice ( feature.get_start(), feature.get_len() );
            while ( ngs_iter.nextAlignment() )
            {
                boolean alignment_reverse = ngs_iter.getIsReversedOrientation();
                if ( alignment_reverse != feature.is_reverse() )
                {
                    tc.on_other_strand();
                }
                else
                {
                    int alignment_mapq = ngs_iter.getMappingQuality();
                    if ( alignment_mapq < opt.get_MinMapq() )
                    {
                        tc.rejected_mapq();
                    }
                    else
                    {
                        read.clear();
                        read.from_cigar( ngs_iter.getAlignmentPosition(), ngs_iter.getShortCigar( false ) );
                        count_result cr = cmp.walk( ft_intervals, read, opt.get_CountMode() );
                        tc.count( cr );
                        if ( cr == count_result.FEATURE )
                            feature.inc_counter();    
                    }
                }
            }
        }
        catch ( ErrorMsg err )
        {
            System.out.println( err.toString() );
        }
    }
}

class Ref_Thread implements Runnable
{
	private Thread thread;
    private SparkGenesOpt opt;
    private String ref;
    private total_count tc;
    private ArrayList<gtf_feature> features;

	public Ref_Thread( String ref, SparkGenesOpt opt, ArrayList<gtf_feature> features )
    {
        this.features = features;
        this.ref = ref;
        this.opt = opt;
        this.tc = new total_count();
        thread = new Thread( this, ref );
	}

    public void start() { thread.start(); }

    public void join( PrintWriter writer, total_count gloabl_counter )
    {
        if ( writer != null )
        {
            try { thread.join(); }
            catch ( InterruptedException e ) { System.out.println( e.toString() ); }

            for ( gtf_feature f : features )
                writer.println( f.get_id() + "\t" + f.get_counter() );
            gloabl_counter.add( tc );
        }
        else
        {
            try { thread.join(); }
            catch ( InterruptedException e ) { System.out.println( e.toString() ); }
        }
    }

	public void run()
    {
        try
        {
            ReadCollection ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( opt.get_Accession() );
            try
            {
                Reference ngs_ref = ngs_run.getReference( ref );
                feature_counter fc = new feature_counter( ngs_ref, opt, tc );
                for ( gtf_feature f : features )
                    fc.process( f );
            }
            catch ( ErrorMsg err ) { System.out.println( ">" + ref + "< not found" ); }
        }
        catch ( ErrorMsg err )
        {
            System.out.println( err.toString() );
        }
	}
}


public class feature_processor
{
    List<gtf_feature> features;

    feature_processor( List<gtf_feature> features )
    {
        this.features = features;
    }

    private long run_local_parallel( SparkGenesOpt opt )
    {
        ArrayList<gtf_feature> res_features = new ArrayList<gtf_feature>();
        ArrayList<Ref_Thread> threads = new ArrayList<Ref_Thread>();
        long processed = 0;
        String curr_ref = "";

        for ( gtf_feature f : features )
        {
            if ( !f.has_ref( curr_ref ) )
            {
                if ( res_features.size() > 0 )
                {
                    Ref_Thread r = new Ref_Thread( curr_ref, opt, res_features );
                    threads.add( r );
                    res_features = new ArrayList<gtf_feature>();
                    r.start();
                    
                }
                /* we have entered a new chromosome! */
                curr_ref = f.get_ref();
                if ( opt.get_ShowProgress() )
                    System.out.println( "\nenter: " + curr_ref );
            }
            res_features.add( f );
            processed++;
        }
        if ( res_features.size() > 0 )
        {
            Ref_Thread r = new Ref_Thread( curr_ref, opt, res_features );
            threads.add( r );
            r.start();
        }

        try
        {
            total_count tc = new total_count();
            PrintWriter writer = new PrintWriter( new BufferedWriter( new FileWriter( opt.get_Outputfile() ) ) );
            for ( Ref_Thread r : threads ) r.join( writer, tc );
            writer.println( tc );
            writer.close();
        }
        catch ( IOException ioe )
        {
            System.out.println( ioe.getMessage() );
        }

        return processed;
    }


    private long run_seq_read_coll_writer( SparkGenesOpt opt, ReadCollection ngs_run, PrintWriter writer )
    {
        long processed = 0;
        String curr_ref = "";
        feature_counter fc = null;
        total_count tc = new total_count();

        for ( gtf_feature f : features )
        {
            if ( !f.has_ref( curr_ref ) )
            {
                curr_ref = f.get_ref();
                try
                {
                    fc = new feature_counter( ngs_run.getReference( curr_ref ), opt, tc );
                }
                catch ( ErrorMsg err )
                {
                    System.out.println( "\n" + err.toString() + " >" + curr_ref + "< not found" );
                    fc = null;
                }

                /* we have entered a new chromosome! */
                if ( opt.get_ShowProgress() && fc != null )
                    System.out.println( "\nenter: " + curr_ref );
            }

            if ( fc != null ) fc.process( f );
            writer.println( f.get_id() + "\t" + f.get_counter() );
            processed++;
        }
        writer.println( tc );
        return processed;
    }

    private long run_seq_read_coll( SparkGenesOpt opt, ReadCollection ngs_run )
    {
        long processed = 0;
        try
        {
            PrintWriter writer = new PrintWriter( new BufferedWriter( new FileWriter( opt.get_Outputfile() ) ) );
            processed = run_seq_read_coll_writer( opt, ngs_run, writer );
            writer.close();
        }
        catch ( IOException ioe )
        {
            System.out.println( ioe.getMessage() );
        }
        return processed;
    }

    private long run_local_seq( SparkGenesOpt opt )
    {
        long processed = 0;
        try
        {
            ReadCollection ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( opt.get_Accession() );
            processed = run_seq_read_coll( opt, ngs_run );
        }
        catch ( ErrorMsg err )
        {
            System.out.println( err.toString() );
        }
        return processed;
    }


    private long run_sparked( SparkGenesOpt opt )
    {
        long processed = 0;

		SparkConf sc = new SparkConf();
        sc.setAppName( "SparkGenes" );
        sc.setMaster( opt.get_SparkMaster() );
		JavaSparkContext jsc = new JavaSparkContext( sc );

        JavaRDD<gtf_feature> rdd_features = jsc.parallelize( features, opt.get_NumSlices() );

        alignment_counter counter = new alignment_counter( opt.get_Accession(),
                                                            opt.get_MinMapq(),
                                                            opt.get_CountMode() );

		JavaPairRDD<String, Integer> feature_counts = rdd_features.mapToPair( counter );

		List< Tuple2< String, Integer > > counters = feature_counts.collect();

		try
		{
			PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter( opt.get_Outputfile() ) ) );
			for( Tuple2 t : counters ) out.println( t._1 + " : " + t._2 );
			out.close();
		}
		catch ( IOException ioe )
		{
			System.out.println( "Trouble writing to file: " + ioe.getMessage() );
		}

        processed = counters.size();

        jsc.stop();
        return processed;
    }

    public long run( SparkGenesOpt opt )
    {
        long res = 0;
        switch ( opt.get_RunMode() )
        {
            case LOCAL_P : res = run_local_parallel( opt ); break;
            case LOCAL_S : res = run_local_seq( opt ); break;
            case SPARKED : res = run_sparked( opt ); break;
        }
        if ( opt.get_ShowProgress() )
            System.out.println( "===> done. / processed = " + res );
        return res;
    }
}

