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

class feature_writer
{
    translater translater;

    feature_writer() { init(); }
    void set_translater( translater dict ) { translater = dict; }    
    void init() { translater = null; }
    void close() { ; }
    void write( gtf_feature f, interval_list l ) { ; }
}

class feature_txt_writer extends feature_writer
{
    private PrintWriter out;
    private boolean just_name;

    feature_txt_writer( final String filename, final boolean just_name )
    {
        super();
        this.just_name = just_name;
        try
        {
            this.out = new PrintWriter( new BufferedWriter( new FileWriter( filename ) ) );
        }
        catch ( IOException ioe )
        {
            System.out.println( "Trouble writing to file: " + ioe.getMessage() );
        }
    }

    @Override void close() { out.close(); }

    //void write_string( final String s ) { out.println( s ); }

    @Override void write( gtf_feature f, interval_list l )
    {
        StringBuffer sb = new StringBuffer();
        if ( just_name )
        {
            sb.append( f.get_id() );
        }
        else
        {
            f.head( sb, translater );
            l.toStartLen( sb );
        }
        out.println( sb.toString() );
    }
}

class feature_list_writer extends feature_writer
{
    private ArrayList<gtf_feature> out_list;

    feature_list_writer( ArrayList<gtf_feature> out_list )
    {
        super();
        this.out_list = out_list;
    }

    @Override void write( gtf_feature f, interval_list l )
    {
        f.translate_ref( translater );
        f.update_segments( l );
        out_list.add( f );
    }
}

/* private class for gtf_processor */
class feature_window
{
    private ArrayList<gtf_feature> to_report;
    private ArrayList<gtf_feature> reported;
    private feature_writer writer;
    private int start, end, max_in_window;
    private interval_list l1;
    private interval_list l2;

    feature_window( feature_writer writer )
    {
        to_report = new ArrayList<gtf_feature>();
        reported = new ArrayList<gtf_feature>();
        this.writer = writer;
        l1 = new interval_list( 2048 );
        l2 = new interval_list( 2048 );
        start = end = max_in_window = 0;
    }

    int get_max_in_window() { return max_in_window; }

    private boolean first_in_to_report_ends_before( long pos )
    {
        if ( to_report.size() < 1 ) return false;
        long e = to_report.get( 0 ).get_end();
        return ( e < pos );
    }

    private boolean first_in_reported_ends_before( long pos )
    {
        if ( reported.size() < 1 ) return false;
        long e = reported.get( 0 ).get_end();
        return ( e < pos );
    }

    private void feature_to_interval_list( gtf_feature f, interval_list l )
    {
        l.clear();
        int n = f.get_ranges();
        for ( int i = 0; i < n; ++i )
            l.add( f.get_start_at( i ), f.get_len_at( i ), false );
    }

    private void report( gtf_feature f, interval_list l )
    {
        while ( first_in_reported_ends_before( f.get_start() ) )
            reported.remove( 0 );
        writer.write( f, l );
        reported.add( f );
    }

    private void overlap_with_reported( interval_list l )
    {
        for ( gtf_feature f : reported )
        {
            feature_to_interval_list( f, l2 );
            l.detect_amb( l2 );
        }
    }

    private void overlap_with_to_report( interval_list l )
    {
        for ( gtf_feature f : to_report )
        {
            feature_to_interval_list( f, l2 );
            l.detect_amb( l2 );
        }
    }

    private void overlap_loop( long next_start )
    {
        while ( first_in_to_report_ends_before( next_start ) )
        {
            gtf_feature f = to_report.get( 0 );    // get the first feature
            to_report.remove( 0 );                  // take it out of the to_report-list

            feature_to_interval_list( f, l1 );      // put its ranges into l1
            overlap_with_reported( l1 );            // overlap it with the already reported features
            overlap_with_to_report( l1 );           // overlap it with the still to report features

            report( f, l1 );
        }
    }

    public void put( gtf_feature f )
    {
        /* this is just to give an idea how big my window is... */
        int stored = to_report.size() + reported.size();
        if ( stored > max_in_window ) max_in_window = stored;
        overlap_loop( f.get_start() );
        to_report.add( f );
    }

    public void clear()
    {
        overlap_loop( Long.MAX_VALUE );
        reported.clear();
    }

    public void close()
    {
        clear();
        writer.close();
    }
}


public class gtf_processor
{
    private String gtf_line_type;   // "exon" or "gene" or "transcript" 
    private String gtf_id_key;      // "gene_id" or "transcript_id" 
    private boolean show_progress;
    private feature_asm asm;
    private ArrayList<String> ref_filter;
    private translater translater;

    gtf_processor( final String gtf_file, final String gtf_line_type, final String gtf_id_key, final boolean show_progress )
    {
        this.gtf_line_type = gtf_line_type;
        this.gtf_id_key = gtf_id_key;
        this.show_progress = show_progress;
        translater = null;
        ref_filter = new ArrayList<String>();
        asm = new feature_asm( new gtf_reader( gtf_file, gtf_line_type, gtf_id_key, show_progress ) );
    }

    void clear_ref_filter() { ref_filter.clear(); }
    void add_ref_filter( final String a_ref ) { ref_filter.add( a_ref ); }
    void set_translater( translater dict ) { translater = dict; }

    private void common_loop_without_prescan( feature_writer writer )
    {
        ArrayList<gtf_feature> temp = new ArrayList<gtf_feature>();
        feature_window window = new feature_window( writer );
        String curr_ref = "";
        boolean skip_ref = false;
        gtf_feature f;
        while ( ( f = asm.next() ) != null )
        {
            if ( !f.has_ref( curr_ref ) )
            {
                /* we have entered a new chromosome :
                ( 1 ) handle what we have accumulated in temp */
                Collections.sort( temp );
                for ( gtf_feature obj : temp ) window.put( obj );
                temp.clear();
                window.clear();

                curr_ref = f.get_ref();
                if ( ref_filter.size() > 0 )
                    skip_ref = (!ref_filter.contains( curr_ref ) );

                if ( show_progress )
                {
                    String translated_ref = curr_ref;
                    if ( translater != null ) translated_ref = translater.translate( translated_ref );
                    System.out.println( "\nenter: " + translated_ref + " ( skip : " + skip_ref + " )" );
                }
            }

            if ( !skip_ref ) temp.add( f );
        }
        window.close();

        if ( show_progress )
            System.out.println( "===> done. / max in window = " + window.get_max_in_window() );
    }


    private void common_loop_with_prescan( feature_writer writer )
    {
        gtf_prescan prescan = new gtf_prescan( asm, show_progress );
        prescan.run();
        asm.reset();

        ArrayList<gtf_feature> unsorted = new ArrayList<gtf_feature>();
        feature_window window = new feature_window( writer );
        String curr_ref = "";
        boolean is_sorted = true;
        boolean skip_ref = false;
        gtf_feature f;
        while ( ( f = asm.next() ) != null )
        {
            if ( !f.has_ref( curr_ref ) )
            {
                /* we have entered a new chromosome :
                ( 1 ) what we have accumulated in unsorted */
                if ( !is_sorted )
                {
                    Collections.sort( unsorted );
                    for ( gtf_feature obj : unsorted ) window.put( obj );
                    unsorted.clear();
                }
                window.clear();

                curr_ref = f.get_ref();
                is_sorted = prescan.is_chromosome_sorted( curr_ref );
                if ( ref_filter.size() > 0 )
                    skip_ref = (!ref_filter.contains( curr_ref ) );

                if ( show_progress )
                {
                    String translated_ref = curr_ref;
                    if ( translater != null ) translated_ref = translater.translate( translated_ref );
                    System.out.println( "\nenter: " + translated_ref + " ( sorted : " + is_sorted + " skip : " + skip_ref + " )" );
                }
            }

            if ( !skip_ref )
            {
                if ( is_sorted )
                    window.put( f );
                else
                    unsorted.add( f );
            }
        }
        window.close();

        if ( show_progress ) System.out.println( "===> done. / max in window = " + window.get_max_in_window() );
    }

    private void common_loop( feature_writer writer, final boolean with_prescan )
    {
        if ( translater != null ) writer.set_translater( translater );
        if ( with_prescan )
            common_loop_with_prescan( writer );
        else
            common_loop_without_prescan( writer );
    }

    void to_file( final String outfile, final boolean with_prescan, final boolean just_name )
    {
        common_loop( new feature_txt_writer( outfile, just_name ), with_prescan );
    }

    void to_list( ArrayList<gtf_feature> features, final boolean with_prescan )
    {
        common_loop( new feature_list_writer( features ), with_prescan );
    }


}

