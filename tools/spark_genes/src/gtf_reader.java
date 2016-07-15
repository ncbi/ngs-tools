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
import java.io.*;

/*
    - the gtf_reader wraps a gtf-file ( given as filename ) into a BufferedReader
    - it can be reset by calling reset() ... throws away the BufferedReader and makes an new one

    - it is an iterator, the user calls next() to get the next line from the gtf-file
    - next() returns false when the file is processed
    - the user has to provide an initialized array of 5 String's when calling next()
    - the iterator will put into it:
        feature[ 0 ] ... chromosome         gtf[ 0 ]
        feature[ 1 ] ... start              gtf[ 3 ]
        feature[ 2 ] ... end                gtf[ 4 ]
        feature[ 3 ] ... strand             gtf[ 6 ]
        feature[ 4 ] ... id                 gtf[ 8 ].id
    - the iterator skips gtf-lines starting with '#' ( comments )
    - the iterator skips gtf-lines which do not contain the col_2_filter-String in the 2nd column
    - the iterator extracts an ID from the 8th column based on the id-filter-String
    - the iterator can print progess ( a dot after every 1000 lines ) if show_progress == true
*/

public class gtf_reader
{
    private BufferedReader buf_reader;
    private String line, comment_char, split_char, id_delim_char, filename, col_2_filter, id_filter;
    private String[] line_tokens;
    private int lines, features, id_filter_ofs;
    private boolean show_progress;

	gtf_reader ( final String filename, final String col_2_filter, final String id_filter, final boolean show_progress )
	{
        comment_char = "#";
        split_char = "\t";
        id_delim_char = ";";
        this.filename = filename;
        this.col_2_filter = col_2_filter;
        this.id_filter = id_filter;
        id_filter_ofs = id_filter.length() + 2;
        this.show_progress = show_progress;
        reset();
	}

    void reset()
    {
        try
        {
            if ( buf_reader != null )
                buf_reader.close();
            buf_reader = new BufferedReader( new InputStreamReader( new FileInputStream( filename ) ) );
            lines = 0;
            features = 0;
        }
        catch ( IOException e )
        {
            e.printStackTrace();
            buf_reader = null;
        }
    }

    int get_lines() { return lines; }
    int get_features() { return features; }

    // 0...no more features, 1...skip, 2...valid feature
    private int read_feature()
    {
        if ( buf_reader != null )
        {
            try
            {
                line = buf_reader.readLine();
                if ( line != null )
                {
                    lines++;
                    if ( show_progress && ( lines % 1000 == 0 ) ) System.out.print( "." );

                    if ( line.startsWith( comment_char ) )
                        return 1;

                    line_tokens = line.split( split_char );
                    if ( line_tokens.length < 9 )
                        return 1;
                    if ( !line_tokens[ 2 ].equals( col_2_filter ) )
                        return 1;

                    return 2;
                }
            }
            catch ( IOException e )
            {
                e.printStackTrace();
            }
        }
        return 0;
    }

    /*
        feature[ 0 ] ... chromosome         gtf[ 0 ]
        feature[ 1 ] ... start              gtf[ 3 ]
        feature[ 2 ] ... end                gtf[ 4 ]
        feature[ 3 ] ... strand             gtf[ 6 ]
        feature[ 4 ] ... id                 gtf[ 8 ].id
    */
    boolean next( String[] feature )
    {
        int rf;
        while ( ( rf = read_feature() ) != 0 )
        {
            if ( rf == 2 )
            {
                int id_idx = line_tokens[ 8 ].indexOf( id_filter );
                if ( id_idx >= 0 )
                {
                    int delim_idx = line_tokens[ 8 ].indexOf( id_delim_char, id_idx );
                    if ( delim_idx >= 0 )
                    {
                        feature[ 0 ] = line_tokens[ 0 ];    // chromosome
                        feature[ 1 ] = line_tokens[ 3 ];    // start
                        feature[ 2 ] = line_tokens[ 4 ];    // end
                        feature[ 3 ] = line_tokens[ 6 ];    // strand
                        feature[ 4 ] = line_tokens[ 8 ].substring( id_filter_ofs, delim_idx - 1 ); // id

                        features++;
                        return true;
                    }
                }
            }
        }
        return false;
    }
}
