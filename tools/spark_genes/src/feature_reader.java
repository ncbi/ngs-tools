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
import java.util.ArrayList;

public class feature_reader
{
    private BufferedReader buf_reader;
    private String line;
    private final String comment_char, split_char, delim_char;
    private String filename;
    private String[] line_tokens;
    private String[] sections;
    private String[] sec_parts;
    private int lines;
    private final boolean show_progress;

	feature_reader ( final String filename, final boolean show_progress )
	{
        this.filename = filename;
        comment_char = "#";
        split_char = "\t";
        delim_char = ";";
        lines = 0;
        this.show_progress = show_progress;
        reset();
	}

    public void reset()
    {
        try
        {
            buf_reader = new BufferedReader( new InputStreamReader( new FileInputStream( filename ) ) );
        }
        catch ( IOException e )
        {
            e.printStackTrace();
            buf_reader = null;
        }
    }

    // 0...no more features, 1...skip, 2...valid feature
    private int read()
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
                    if ( line_tokens.length < 4 )
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
        line_tokens[ 0 ] ... id
        line_tokens[ 1 ] ... strand
        line_tokens[ 2 ] ... ref
        line_tokens[ 3 ] ... sections
    */
    public gtf_feature next()
    {
        int rf;
        while ( ( rf = read() ) != 0 )
        {
            if ( rf == 2 )
            {
                sections = line_tokens[ 3 ].split( delim_char );
                int n = sections.length;

                gtf_feature res = new gtf_feature( n, line_tokens[ 0 ], line_tokens[ 2 ], line_tokens[ 1 ].equals( "-" ) );

                for ( int i = 0; i < n; ++ i )
                {
                    sec_parts = sections[ i ].substring( 1 ).split( "\\." );
                    if ( sec_parts.length > 1 )
                    {
                        res.add_seg_1( Integer.parseInt( sec_parts[ 0 ] ),
                                       Integer.parseInt( sec_parts[ 1 ] ),
                                       sections[ i ].startsWith( "A" ) );
                    }
                }
                return res;
            }
        }
        if ( show_progress ) System.out.println( "." );
        return null;
    }

    public void to_list( ArrayList<gtf_feature> features )
    {
        gtf_feature f;
        while ( ( f = next() ) != null ) features.add( f );
    }
}
