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
import java.util.*;

public class translater
{
    private BufferedReader buf_reader;
    private String line, comment_char, split_char;
    private String[] line_tokens;
    private Map<String, String> dict;

	translater ( final String filename )
	{
        comment_char = "#";
        split_char = "=";
        try
        {
            buf_reader = new BufferedReader( new InputStreamReader( new FileInputStream( filename ) ) );
        }
        catch ( IOException e )
        {
            buf_reader = null;
        }
        dict = new HashMap<String, String>();
        while ( next() )
        {
            dict.put( line_tokens[ 0 ].trim(), line_tokens[ 1 ].trim() );
        }
	}

    // 0...no more features, 1...skip, 2...valid feature
    private int read_line()
    {
        if ( buf_reader != null )
        {
            try
            {
                line = buf_reader.readLine();
                if ( line != null )
                {
                    if ( line.startsWith( comment_char ) )
                        return 1;

                    line_tokens = line.split( split_char );
                    if ( line_tokens.length < 2 )
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

    private boolean next()
    {
        int rf;
        while ( ( rf = read_line() ) != 0 )
        {
            if ( rf == 2 )
                return true;
        }
        return false;
    }

    public String translate( final String key )
    {
        String res = dict.get( key );
        if ( res == null ) res = key;
        return res;
    }
}
