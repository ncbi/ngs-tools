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

public class interval_parser
{
    private List<interval> window;

    public interval_parser()
    {
        window = new LinkedList<interval>();
    }

    private void put( interval i )
    {
        window.add( i );
    }

    private void sort_window_and_write_out( String curr_ref, List<interval> intervals )
    {
        if ( window.size() > 0 )
        {
            Collections.sort( window );
            long last_end = 0;
            long last_start = 0;

            for ( interval i : window )
            {
                long this_start = i.get_start();
                long this_end = i.get_end();

                if ( this_start > ( last_end + 1 ) )
                {
                    // we have to insert a gap!
                    intervals.add( new interval( curr_ref, "gap", last_end + 1, (int)( this_start - last_end ) - 1, false ) );
                }

                if ( this_start > last_start && this_end != last_end && this_start > last_end )
                {
                    intervals.add( i );
                    last_start = this_start;
                    last_end = this_end;
                }
            }
            // and a final gap to the end of the reference and beyond!
            intervals.add( new interval( curr_ref, "gap", last_end + 1, Integer.MAX_VALUE, false ) );
            window.clear();
        }
    }

    public void parse( Iterable<gtf_feature> features, List<interval> intervals )
    {
        String curr_ref = "";

        for ( gtf_feature f : features )
        {
            if ( !f.has_ref( curr_ref ) )
            {
                // enter new reference
                sort_window_and_write_out( curr_ref, intervals );
                curr_ref = f.get_ref();
            }

            /* loop through the intervals of this feature */
            int n = f.get_ranges();
            boolean f_reverse = f.is_reverse();
            String f_id = f.get_id();
            for ( int idx = 0; idx < n; ++idx )
            {
                long feature_start = f.get_start_at( idx );
                int len = f.get_len_at( idx );

                if ( f.get_amb_at( idx ) )
                    window.add( new interval( curr_ref, "amb", feature_start, len, f_reverse ) );
                else
                    window.add( new interval( curr_ref, f_id, feature_start, len, f_reverse ) );
            }
        }
        sort_window_and_write_out( curr_ref, intervals );
    }

    public void write_intervals( Iterable<interval> intervals, final String filename )
    {
        try
        {
            PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter( filename ) ) );
            for ( interval i : intervals ) out.println( i );
            out.close();
        }
        catch ( IOException ioe )
        {
            System.out.println( "Trouble writing to file: " + ioe.getMessage() );
        }
    }
}