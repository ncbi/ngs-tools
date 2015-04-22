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
import ngs.*;

public final class SparkGenes
{
	private static String first_line_of( final String filename )
	{
		String res = "";
        try
        {
            BufferedReader r = new BufferedReader( new InputStreamReader( new FileInputStream( filename ) ) );
			res = r.readLine();
        }
        catch ( IOException e )
        {
            e.printStackTrace();
        }
		return res;
	}

	private static List< gtf_feature > read_features_from_gtf( SparkGenesOpt opt )
    {
        gtf_processor p = new gtf_processor( opt.get_GTFfile(),
                                              opt.get_FeatureType(), 
                                              opt.get_FeatureID(),
                                              opt.get_ShowProgress() );

        if ( opt.has_RefTranslation() )
            p.set_translater( new translater( opt.get_RefTranslation() ) );

        ArrayList<gtf_feature> features = new ArrayList<gtf_feature>();
        p.to_list( features, opt.get_PreScanGtf() );
        return features;
    }

	private static List< gtf_feature > read_features_from_pre( SparkGenesOpt opt )
    {
        feature_reader p = new feature_reader( opt.get_GTFfile(),
												opt.get_ShowProgress() );
        ArrayList<gtf_feature> features = new ArrayList<gtf_feature>();
        p.to_list( features );
        return features;
    }

    private static long count_alignments_of_these_refs( final String acc, Iterable<String> refs )
    {
        long res = 0;
        try
        {
            ReadCollection ngs_run = gov.nih.nlm.ncbi.ngs.NGS.openReadCollection( acc );
            for ( String spec : refs )
            {
                Reference ngs_ref = ngs_run.getReference( spec );
                res += ngs_ref.getAlignmentCount();
            }
        }
        catch ( ErrorMsg err )
        {
            System.err.println( err.toString() );
        }
        return res;
    }

    private static List<gtf_feature> filter_features_by_refs( Iterable<gtf_feature> src, Set<String> refs )
    {
        List<gtf_feature> res = new LinkedList<gtf_feature>();
        for ( gtf_feature f : src )
        {
            if ( refs.contains( f.get_ref() ) ) res.add( f );
        }
        return res;
    }

    // -----------------------------------------------------------------------------------------------------------

	private static void process_gtf_and_save( final List<gtf_feature> features, final String output_file )
	{
        try
        {
            int max_ranges_per_feature = 0;
            PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter( output_file ) ) );
			out.println( "#preprocessed" );
            for ( gtf_feature f : features )
            {
                int n = f.get_ranges();
                if ( n > max_ranges_per_feature ) max_ranges_per_feature = n;
                out.println( f );
            }
            out.close();
            System.out.println( "max intervals : " + max_ranges_per_feature );
            System.out.println( "written to file : " + output_file );
        }
        catch ( IOException ioe )
        {
            System.out.println( "Trouble writing to file: " + ioe.getMessage() );
        }
    }


    private static void save_counts_and_conters( feature_counts ftc,
                                                 lookup_counters c,
                                                 String path )
    {
        try
        {
            PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter( path ) ) );
            if ( ftc != null ) ftc.save_to_file( out );
            if ( c != null ) c.save_to_file( out );
            out.close();
        }
        catch ( IOException ioe )
        {
            System.out.println( "Trouble writing to file: " + ioe.getMessage() );
        }
    }

	private static void print_elapsed_time( long startTime )
	{
		long stopTime = System.currentTimeMillis();
		long elapsedTime = stopTime - startTime;
		System.out.println( "time = " + ( elapsedTime / 1000 ) + " seconds." );
	}
	
    private static lookup_counters run_lookup_seq( feature_lookup_processor p,
                                                   feature_counts ftc,
                                                   SparkGenesOpt opt )
    {
		long startTime = System.currentTimeMillis();
		
        ftc.clear();
        lookup_counters res = p.run_sequentially( ftc, opt );
        if ( opt.get_ShowProgress() ) res.report();

        if ( opt.get_MeasureTime() )
			print_elapsed_time( startTime );

		return res;
    }

    private static lookup_counters run_lookup_parallel( feature_lookup_processor p,
                                                        feature_counts ftc,
                                                        SparkGenesOpt opt )
    {
		long startTime = System.currentTimeMillis();

        ftc.clear();
        lookup_counters res = p.run_in_threads( ftc, opt );
        if ( opt.get_ShowProgress() ) res.report();

        if ( opt.get_MeasureTime() )
			print_elapsed_time( startTime );

		return res;
    }

    private static lookup_counters run_lookup_sparked( feature_lookup_processor p,
                                                       feature_counts ftc,
                                                       SparkGenesOpt opt )
    {
		long startTime = System.currentTimeMillis();
		
        ftc.clear();
        lookup_counters res = p.run_sparked( ftc, opt );
        if ( opt.get_ShowProgress() ) res.report();
		
        if ( opt.get_MeasureTime() )
			print_elapsed_time( startTime );

		return res;
    }


	public static void main( String[] args ) throws Exception
	{
        SparkGenesOpt opt = new SparkGenesOpt( args );
        
        if ( opt.get_HelpRequested() )
        {
            opt.print_help();
            return;
        }

        if ( opt.get_GTFfile() == null )
        {
            System.out.println( "please specify a GTF- or preprocessed file! example: -f file.gtf" );
            return;
        }

        if ( opt.get_RunMode() == run_mode.PRE_GTF && opt.get_Raw_Outputfile() == null )
		{
			// we do need a raw (not computed) ouput-file in case of pre-processing
            System.out.println( "please specify an output-file! example: -o features.txt" );
            return;
		}

		List< gtf_feature > features;
		String line1 = first_line_of( opt.get_GTFfile() );
		if ( line1.equals( "#preprocessed" ) )
		{
			if ( opt.get_ShowProgress() )
				System.out.println( "reading a preprocessed gtf-file" );
			features = read_features_from_pre( opt );
		}
		else
		{
			if ( opt.get_ShowProgress() )
				System.out.println( "reading a raw gtf-file" );
		    features = read_features_from_gtf( opt );
		}
		
        if ( opt.get_ShowProgress() )
            System.out.println( features.size() + " features read from " + opt.get_GTFfile() );

        if ( opt.get_RunMode() == run_mode.PRE_GTF )
        {
            process_gtf_and_save( features, opt.get_Outputfile() );
            return;
        }

        if ( opt.get_Accession() == null )
        {
            System.out.println( "please specify an accession: -a SRR123456" );
            return;
        }

		ref_cmp_res refs = ref_cmp.compare( features, opt.get_Accession() );		
		
		if ( opt.get_RunMode() == run_mode.REFS )
		{
			int l = refs.longest_string();
			
			System.out.println( "\n\nComparison of reference-names between" );
			System.out.println( "GTF-file  : " + opt.get_GTFfile() );
			System.out.println( "Accession : " + opt.get_Accession() );
			String fmt = "%nGTF%" + l + "s%n";
			System.out.printf( fmt, "ACC" );
			System.out.println( "-------------------------------------------------------------" );
		
			for ( StringPair sp : refs.in_gtf_and_acc ) System.out.println( "     ..... " + sp + " ....." );
			System.out.println( "-------------------------------------------------------------" );
			
			fmt = "%" + ( l + l ) + "s%n";
			for ( StringPair sp : refs.only_in_acc ) System.out.printf( fmt, sp );		
			System.out.println( "-------------------------------------------------------------" );
			
			for ( String s : refs.only_in_gtf ) System.out.println( s );

			return;
		}
		
        feature_lookup_processor processor = new feature_lookup_processor( features,
																			opt.get_LookupBinSize(),
																			refs );
        if ( opt.get_ShowProgress() )
            processor.report();

        feature_counts ftc = new feature_counts( features );

        lookup_counters c = null;
        switch( opt.get_RunMode() )
        {
            case LOCAL_S : c = run_lookup_seq( processor, ftc, opt ); break;
            case LOCAL_P : c = run_lookup_parallel( processor, ftc, opt ); break;
            case SPARKED : c = run_lookup_sparked( processor, ftc, opt ); break;
        }

		
        if ( opt.get_CountUnaligned() && c != null )
        {
            if ( opt.get_ShowProgress() )
                System.out.println( "counting unaligned reads:" );
            long n_unaligned = unaligned.count( opt );
            c.unaligned( n_unaligned );
            if ( opt.get_ShowProgress() )
                System.out.println( n_unaligned + " unaligned reads" );
        }

        save_counts_and_conters( ftc, c, opt.get_Outputfile() );
	}
}
