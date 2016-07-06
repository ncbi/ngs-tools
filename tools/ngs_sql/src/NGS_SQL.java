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
import java.io.*;

import org.apache.spark.*;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.*;
import org.apache.spark.api.java.function.*;
import org.apache.spark.sql.api.java.*;


public final class NGS_SQL
{
	private static boolean make_table( NGS_TAB_Opt opt, int slices,
									JavaSparkContext spark_context,
									JavaSQLContext sql_context )
	{
		// make the jobs, the node-factory forks inside depending on the table-type...
		List< job > jobs = new LinkedList< job >();
		job_maker job_factory = new job_maker( opt );
		if ( !job_factory.make_jobs( jobs ) )
			return false;
		
		// parallelize the jobs...
		JavaRDD< job > rdd_jobs = spark_context.parallelize( jobs, slices );
		
		EnumSet< ngs_field > fields = opt.get_fields();
		
		// make the functor...
		ngs_functor functor = null;
		switch( opt.get_tab_type() )
		{
			case ALIG 	: functor = new ngs_alig_functor( opt.get_src(), fields ); break;
			case SPOT 	: functor = new ngs_spot_functor( opt.get_src(), fields ); break;
			case FRAG 	: functor = new ngs_frag_functor( opt.get_src(), fields ); break;
			case REF  	: functor = new ngs_ref_functor( opt.get_src(), fields ); break;
			case RDGRP	: functor = new ngs_rdgrp_functor( opt.get_src(), fields ); break;
			case PILEUP	: functor = new ngs_pileup_functor( opt.get_src(), fields ); break;
		}

		if ( functor == null ) return false;
		
		// let the functor transform the jobs into a RDD of sql-rows
		JavaRDD< Row > rows = rdd_jobs.flatMap( functor );

		// make a schema matching the requested table-type
		List< StructField > sfl = new ArrayList< StructField >();
		for ( ngs_field f : fields )
			sfl.add( DataType.createStructField( f.toString().toLowerCase(), f.sql_type(), true ) );
		StructType schema = DataType.createStructType( sfl );
		
		// transform the RDD of rows into a RDD with rows and a schema
		JavaSchemaRDD src_tab = sql_context.applySchema( rows, schema );

		// give it the name that can be used in the SQL-statement later
		src_tab.registerTempTable( opt.get_tab_name() );
		
		return true;
	}

	
	public static void save_as_textfiles( NGS_SQL_Opt opt, JavaRDD< String > res_txt, String header )
	{
		// remove the content of the temp. directory if it exists...
		post_processor.delete_recursive( new File( opt.get_Temp() ) );

		// let the nodes write into a temp. output folder
		res_txt.saveAsTextFile( opt.get_Temp() );
		
		// concatenate all the text-files written...
		long lines = post_processor.concat( opt.get_Temp(), opt.get_Out(), header );

		// print summary
		System.out.println( "\n" + lines + " lines written." );
	}

	public static void save_via_collectAsync( NGS_SQL_Opt opt, JavaRDD< String > res_txt, String header )
	{
        try
        {
            PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter( opt.get_Out() ) ) );

			if ( header.length() > 0 )
				out.println( header );
				
			JavaFutureAction< List< String > > future_strings = res_txt.collectAsync();
			while ( !future_strings.isDone() )
			{
				long received = 0;
				try
				{
					List< String > sl = future_strings.get();
					for ( String s : sl )
					{
						out.println( s );
						received++;
					}
				}
				catch ( Exception e )
				{
					System.err.println( e );
				}
				
				System.out.println( received + " lines received"  );
				out.flush();
			}
			out.close();
        }
        catch ( IOException ioe )
        {
            System.err.println( "Trouble writing to file: " + opt.get_Out() );
        }
	}
	
	public static void save_via_LocalIterator( NGS_SQL_Opt opt, JavaRDD< String > res_txt, String header )
	{
        
		Iterator< String > iter = res_txt.toLocalIterator();
		long received = post_processor.write_iter( iter, opt.get_Out(), header );
		System.out.println( received + " lines received"  );
	}

	public static void save_via_collectPartitions( NGS_SQL_Opt opt, JavaRDD< String > res_txt, String header )
	{
        try
        {
            PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter( opt.get_Out() ) ) );

			if ( header.length() > 0 )
				out.println( header );

			List< Partition > partition_list = res_txt.partitions();
			Iterator< Partition > partition_iter = partition_list.iterator();
			
			System.out.println( partition_list.size() + " partitions in result-rdd" );
			
			List< Integer > partition_id_list = new ArrayList< Integer >();
			List< String >[] result;

			long received = 0;
			boolean done = false;
			while ( !done )
			{
				partition_id_list.clear();
				for ( int i = 0; i < opt.get_ResRequests(); ++i )
				{
					if ( partition_iter.hasNext() )
					{
						Partition p = partition_iter.next();
						partition_id_list.add( p.index() );
					}
				}
				done = ( partition_id_list.size() < 1 );
				if ( !done )
				{
					int[] partition_id_array = new int[ partition_id_list.size() ];
					
					for ( int i = 0; i < partition_id_list.size(); i++ )
						partition_id_array[ i ] = partition_id_list.get( i );
						
					result = res_txt.collectPartitions( partition_id_array );
					for ( List< String > part : result )
					{
						for ( String s : part )
						{
							out.println( s );
							received++;
						}
					}
				}
			}

			System.out.println( received + " lines received"  );
			out.close();
        }
        catch ( IOException ioe )
        {
            System.err.println( "Trouble writing to file: " + opt.get_Out() );
        }
	}
	
	public static void main( String[] args ) throws Exception
	{
		long startTime = System.currentTimeMillis();
		
        NGS_SQL_Opt opt = new NGS_SQL_Opt( args[ 0 ] );

		if ( !opt.valid() )
		{
            System.out.println( "invalid options!" );
			opt.report();
            return;
		}

		SparkConf spark_conf = new SparkConf();
		spark_conf.setAppName( NGS_SQL.class.getSimpleName() );
        spark_conf.setMaster( opt.get_Master() );

		String exec_mem = opt.get_Exec_Mem();
		if ( exec_mem.length() > 0 )
			spark_conf.set( "spark.executor.memory", exec_mem );
		
		String drv_mem = opt.get_Drv_Mem();
		if ( drv_mem.length() > 0 )
			spark_conf.set( "spark.driver.memory", drv_mem );

		String max_res = opt.get_Max_Res();
		if ( max_res.length() > 0 )
			spark_conf.set( "spark.driver.maxResultSize", max_res );
		
		JavaSparkContext spark_context = new JavaSparkContext( spark_conf );
		JavaSQLContext sql_context = new JavaSQLContext( spark_context );

		int n = opt.get_table_count();
		boolean table_ok = true;
		for ( int i = 0; i < n && table_ok; ++i )
			table_ok = make_table( opt.get( i ), opt.get_Slices(), spark_context, sql_context );

		if ( !table_ok ) return;
		
		// SQL can be run over RDDs that have been registered as tables.
		JavaSchemaRDD sql_result = sql_context.sql( opt.get_Sql() );
		
		// write the result-set of the SQL-statement into the requested directory
		JavaRDD< Row > row_result = sql_result.wrapRDD( sql_result.rdd() );

		out_fmt format = opt.get_format();
		String header = "";
		
		try
		{
			JavaRDD< String > res_txt;
			
			if ( format == out_fmt.FASTQ )
			{
				fastq_functor fqf = new fastq_functor( sql_result.schema().getFields() );
				res_txt = row_result.flatMap( fqf );
			}
			else
			{
				String delim = "\t";
				if ( format == out_fmt.CSV ) delim = ",";

				res_line_functor rf = new res_line_functor( sql_result.schema().getFields(), delim );
				res_txt = row_result.map( rf );
				
				header = post_processor.header( sql_result.schema().getFields(), delim );
			}
			
			//save_via_LocalIterator( opt, res_txt, header );
			//save_via_collectAsync( opt, res_txt, header );
			save_via_collectPartitions( opt, res_txt, header );
			
		}
		catch ( Exception e )
		{
            System.err.println( "\n>>>>>> " + e + "\n" );
		}
		
		spark_context.stop();
		
		long stopTime = System.currentTimeMillis();
		long elapsedTime = stopTime - startTime;
		System.out.println( "time = " + ( elapsedTime / 1000 ) + " seconds." );
	}
}
