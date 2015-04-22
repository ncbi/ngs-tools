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
import org.apache.spark.sql.api.java.*;

public class post_processor
{
	public static void delete_recursive( File path )
	{
		File[] currList;
		Stack< File > stack = new Stack< File >();
		stack.push( path );
		while ( !stack.isEmpty() )
		{
			if ( stack.lastElement().isDirectory() )
			{
				currList = stack.lastElement().listFiles();
				if ( currList.length > 0 )
				{
					for ( File curr: currList )
					{
						stack.push( curr );
					}
				}
				else
				{
					stack.pop().delete();
				}
			}
			else
			{
				stack.pop().delete();
			}
		}		
	}
	
	public static long concat( final String src, final String dst, final String header )
	{
		long res = 0;
        try
        {
            PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter( dst ) ) );

			if ( header != null && header.length() > 0 )
				out.println( header );
				
			File dir = new File( src );
			if ( dir.exists() )
			{
				String line;
				String[] entries = dir.list();
				for ( String s : entries )
				{
					if ( s.startsWith( "part-" ) )
					{
						File f = new File( dir.getPath(), s );
						BufferedReader buf_reader = new BufferedReader( new InputStreamReader( new FileInputStream( f.getPath() ) ) );
						do
						{
							line = buf_reader.readLine();
							if ( line != null )
							{
								out.println( line );
								if ( ( ++res % 1000 ) == 0 ) System.out.print( "." );
							}
						} while ( line != null );
						f.delete();
					}
				}
			}
			out.close();
        }
        catch ( IOException ioe )
        {
            System.out.println( "Trouble writing to file: " + dst );
        }
		return res;
	}
	
	public static String header( StructField schema_fields[], final String delim )
	{
		StringBuffer sb = new StringBuffer();
		boolean c = false;
		for ( StructField sf : schema_fields )
		{
			if ( c )
				sb.append( delim );
			else
				c = true;
			sb.append( sf.getName() );
		}
		return sb.toString();
	}
	
	public static long write_iter( Iterator< String > iter, final String dst, final String header )
	{
		long res = 0;
        try
        {
            PrintWriter out = new PrintWriter( new BufferedWriter( new FileWriter( dst ) ) );

			if ( header.length() > 0 )
				out.println( header );

			while( iter.hasNext() )
			{
				out.println( iter.next() );
				if ( ( ++res % 1000 ) == 0 ) System.out.print( "." );
			}
			out.close();
        }
        catch ( IOException ioe )
        {
            System.out.println( "Trouble writing to file: " + dst );
        }
		return res;
	}
	
}