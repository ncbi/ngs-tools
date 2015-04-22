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

class alignment_attrib implements java.io.Serializable
{
    public String id;
    public interval_type iv_type;
    public long count;

    private void set( final String id, final interval_type iv_type, final long count )
    {
        this.id = id;
        this.iv_type = iv_type;
        this.count = count;
    }
    
    alignment_attrib() { set( "", interval_type.NONE, 0 ); }
    alignment_attrib( final String id, final interval_type iv_type ) { set( id, iv_type, 1 ); }
    alignment_attrib( final String id, final interval_type iv_type, final long count ) { set( id, iv_type, count ); }

    alignment_attrib( final alignment_attrib a, alignment_attrib b )
    {
        if ( a == null )
        {
            if ( b == null )
                set( "", interval_type.NONE, 0 );
            else
                set( b.id, b.iv_type, b.count );
        }
        else
        {
            if ( b == null )
                set( a.id, a.iv_type, a.count );
            else
            {
                switch( a.iv_type )
                {
                    case FEATURE    :  switch( b.iv_type )
                                        {
                                            case FEATURE    : if ( a.id.equals( b.id ) )
                                                                    set( a.id, interval_type.FEATURE, a.count + b.count );
                                                               else
                                                                    set( a.id, interval_type.AMBIGUOUS, a.count + b.count );
                                                               break;

                                            case AMBIGUOUS  : set( a.id, interval_type.AMBIGUOUS, a.count + b.count ); break;
                                            case GAP        : set( a.id, a.iv_type, a.count ); break;
                                            case NONE       : set( a.id, a.iv_type, a.count ); break;
                                        }
                                        break;

                    case AMBIGUOUS  :  set( a.id, interval_type.AMBIGUOUS, a.count + b.count );
                                        break;

                    case GAP        :  switch( b.iv_type )
                                        {
                                            case FEATURE    : set( b.id, b.iv_type, b.count ); break;
                                            case AMBIGUOUS  : set( a.id, interval_type.AMBIGUOUS, a.count + b.count ); break;
                                            case GAP        : set( a.id, interval_type.GAP, a.count + b.count ); break;
                                            case NONE       : set( a.id, a.iv_type, a.count ); break;
                                        }
                                        break;

                    case NONE       :  set( b.id, b.iv_type, b.count );
                                        break;
                }
            }
        }
    }

    @Override public String toString()
    {
        StringBuffer sb = new StringBuffer();
        sb.append( iv_type );
        sb.append( ":" );
        sb.append( id );
        sb.append( ":" );
        sb.append( count );
        return sb.toString();
    }
}


class amb_filter implements Function< alignment_attrib, Boolean >
{
    public Boolean call ( alignment_attrib t )
    {
        return ( t.iv_type == interval_type.AMBIGUOUS );
    }
}

class no_feature_filter implements Function< alignment_attrib, Boolean >
{
    public Boolean call ( alignment_attrib t )
    {
        return ( t.iv_type == interval_type.GAP );
    }
}


class alignment_attrib_pair_maker implements PairFlatMapFunction < interval, String, alignment_attrib >
{
	private String accession, reference;
    private boolean initialized;
	private ReadCollection ngs_run;
	private Reference ngs_ref;

	public alignment_attrib_pair_maker( String accession )
	{
        this.initialized = false;
		this.accession = accession;
		this.reference = "";
		this.ngs_run = null;
        this.ngs_ref = null;
	}

    private boolean open_accession()
    {
        boolean res = false;
        try
        {
            ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( accession );
            res = true;
        }
        catch ( ErrorMsg err )
        {
            System.out.println( err.toString() );
            ngs_run = null;
        }
        return res;
    }

    private void open_reference( String a_refname )
    {
        reference = a_refname;
        try
        {
            ngs_ref = ngs_run.getReference( a_refname );
        }
        catch ( ErrorMsg err )
        {
            System.out.println( err.toString() );
            ngs_ref = null;
        }
    }

	public Iterable< scala.Tuple2< String, alignment_attrib > > call( interval intv )
	{
        ArrayList< scala.Tuple2< String, alignment_attrib > > res = new ArrayList< scala.Tuple2< String, alignment_attrib > >();

        if ( !initialized )
            initialized = open_accession();

		if ( initialized )
		{
			if ( !intv.has_ref( reference ) )
                open_reference( intv.get_ref() );

            String feature_id = intv.get_id();
            interval_type iv_type = intv.get_type();
            int count = 0;

			if ( ngs_ref != null )
			{
                try
                {
                    AlignmentIterator ngs_iter = ngs_ref.getAlignmentSlice( intv.get_start(), intv.get_len() );
                    while ( ngs_iter.nextAlignment() )
                    {
                        count++;
                        res.add( new scala.Tuple2<>( ngs_iter.getAlignmentId(), new alignment_attrib( feature_id, iv_type ) ) );
                    }
                }
                catch ( ErrorMsg err )
                {
                    System.out.println( err.toString() );
                }
            }

            // we have to record that we have no alignments for this feature, later we want to know
            // which features have zero alignments in this run...
            // we use the feature-id as a key ...
            if ( count == 0 && ( iv_type == interval_type.FEATURE ) )
            {
                res.add( new scala.Tuple2<>( feature_id, new alignment_attrib( feature_id, interval_type.NONE, 0 ) ) );
            }

        }
        return res;
    }
}


class attrib_reducer implements Function2 < alignment_attrib, alignment_attrib, alignment_attrib >
{
	public alignment_attrib call( alignment_attrib a, alignment_attrib b )
    {
        return new alignment_attrib( a, b );
    }
}


class feature_id_pair_maker implements PairFlatMapFunction < alignment_attrib, String, Integer >
{
	public Iterable< scala.Tuple2< String, Integer > > call( alignment_attrib attr )
    {
        ArrayList< scala.Tuple2< String, Integer > > res = new ArrayList< scala.Tuple2< String, Integer > >();
        switch( attr.iv_type )
        {
            case FEATURE    : res.add( new scala.Tuple2<>( attr.id, 1 ) ); break;
            case AMBIGUOUS  : res.add( new scala.Tuple2<>( attr.id, 0 ) ); break;
            case GAP        : res.add( new scala.Tuple2<>( attr.id, 0 ) ); break;
            case NONE       : res.add( new scala.Tuple2<>( attr.id, 0 ) ); break;
        }
        return res;
    }
}


class feature_count_reducer implements Function2 < Integer, Integer, Integer >
{
	public Integer call( Integer a, Integer b ) {  return a + b; }
}


public class interval_processor
{
    ArrayList<interval> intervals;

    interval_processor( ArrayList<interval> intervals )
    {
        this.intervals = intervals;
    }

    private boolean process_interval( Reference ngs_ref, PrintWriter writer, interval intv )
    {
        boolean res = false;
        try
        {
            AlignmentIterator ngs_iter = ngs_ref.getAlignmentSlice ( intv.get_start(), intv.get_len() );
            while ( ngs_iter.nextAlignment() )
            {
                writer.println( ngs_iter.getAlignmentId() + "\t" + intv.get_id() );
            }
            res = true;
        }
        catch ( ErrorMsg err )
        {
            System.out.println( err.toString() );
        }
        return res;
    }

    private long run_local_seq( SparkGenesOpt opt )
    {
        long processed = 0;
        try
        {
            PrintWriter writer = new PrintWriter( new BufferedWriter( new FileWriter( opt.get_Outputfile() ) ) );
            try
            {
                ReadCollection ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( opt.get_Accession() );
                Reference ngs_ref = null;
                String reference = "";
                for ( interval intv : intervals )
                {
                    if ( !reference.equals( intv.get_ref() ) )
                    {
                        reference = intv.get_ref();
                        System.out.println( "\nenter : " + reference );
                        try
                        {
                            ngs_ref = ngs_run.getReference( reference );
                        }
                        catch ( ErrorMsg err )
                        {
                            System.out.println( ">" + reference + "< not found" );
                            ngs_ref = null;
                        }
                    }

                    if ( ngs_ref != null )
                    {
                        if ( process_interval( ngs_ref, writer, intv ) )
                            processed++;
                        if ( processed % 1000 == 0 )
                            System.out.print( "." );        
                    }
                }
            }
            catch ( ErrorMsg err )
            {
                System.out.println( err.toString() );
            }
            writer.close();
        }
        catch ( IOException ioe )
        {
            System.out.println( ioe.getMessage() );
        }
        return processed;
    }


    private long run_sparked( SparkGenesOpt opt )
    {
		SparkConf sc = new SparkConf();
        sc.setAppName( "SparkGenes" );
        sc.setMaster( opt.get_SparkMaster() );
		JavaSparkContext jsc = new JavaSparkContext( sc );

        // (1) make a RDD containing the whole GTF as sequence of intervals with the gaps artifically inserted
        JavaRDD< interval > RDD_of_intervals = jsc.parallelize( intervals, opt.get_NumSlices() );

        // (2) make a PairRDD containing pairs of alignment-ID's vs a attrib object ( FEATURE-ID, matches, gap, amb )
        alignment_attrib_pair_maker pair_maker1 = new alignment_attrib_pair_maker( opt.get_Accession() );
        JavaPairRDD< String, alignment_attrib > PairRDD_of_align_vs_attrib = RDD_of_intervals.flatMapToPair( pair_maker1 );

        // (3) reduce by the alignment-ID, accumulating status in the attrib object, stripping off the alignment-ID
        JavaRDD< alignment_attrib > RDD_of_attribs = PairRDD_of_align_vs_attrib.reduceByKey( new attrib_reducer() ).values();

        // (4) counting the not ambiguous and not matching alignments
        long reads_counted_as_ambiguous = RDD_of_attribs.filter( new amb_filter() ).count();
        long reads_not_counted = RDD_of_attribs.filter( new no_feature_filter() ).count();

        // (5) make a RDD containing FEATURE-ID's vs counters
        JavaPairRDD< String, Integer > PairRDD_of_Features = RDD_of_attribs.flatMapToPair( new feature_id_pair_maker() );
        JavaPairRDD< String, Integer > PairRDD_of_reduces_Features = PairRDD_of_Features.reduceByKey( new feature_count_reducer() );

        //PairRDD_of_reduces_Features.saveAsTextFile( opt.get_Outputfile() );
    
        java.util.List< scala.Tuple2< String, Integer > > feature_counts = PairRDD_of_reduces_Features.collect();

        long reads_counted = 0;
        try
        {
            PrintWriter writer = new PrintWriter( new BufferedWriter( new FileWriter( opt.get_Outputfile() ) ) );
            for ( scala.Tuple2< String, Integer > t : feature_counts )
            {
                writer.println( t._1 + "\t" + t._2 );
                reads_counted += t._2;
            }
            writer.println( "\n\nambiguous reads = " + reads_counted_as_ambiguous );
            writer.println( "no feature = " + reads_not_counted );
            writer.println( "as feature = " + reads_counted );
            writer.close();

        }
        catch ( IOException ioe )
        {
            System.out.println( ioe.getMessage() );
        }

        long reads_total = reads_counted_as_ambiguous + reads_not_counted + reads_counted;

        System.out.println( "ambiguous reads = " + reads_counted_as_ambiguous );
        System.out.println( "no feature = " + reads_not_counted );
        System.out.println( "as feature = " + reads_counted );
        System.out.println( "reads total = " + reads_total );

        return reads_total;
    }

    public long run( SparkGenesOpt opt )
    {
        long res = 0;
        switch ( opt.get_RunMode() )
        {
            case LOCAL_S : res = run_local_seq( opt ); break;
            case SPARKED : res = run_sparked( opt ); break;
        }
        if ( opt.get_ShowProgress() )
            System.out.println( "===> done. / processed = " + res );
        return res;
    }

}