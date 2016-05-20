/* ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
=========================================================================== */
package HttpWrap;

import java.io.*;
import java.net.*;
import data.CLogger;

public class HttpDnldRangeReader extends HttpDnldUrlWrapper
{
    public long get_size()
    {
        long res = 0;
        if ( valid() )
        {
            HttpURLConnection conn = make_connection( "HEAD" );
            if ( conn != null )
            {
                try
                {
                    conn.connect();
                    if ( conn.getResponseCode() == HttpURLConnection.HTTP_OK )
                        res = conn.getContentLengthLong();
                    conn.disconnect();
                }
                catch ( IOException ex ) { CLogger.log( ex.toString() ); }
            }
        }
        return res;
    }
    
    public void read( long start, int count, HttpDnldContent c )
    {
        c.available = 0;
        if ( valid() )
        {
            HttpURLConnection conn = make_connection( "GET" );
            if ( conn != null )
            {
                try
                {
                    String range = String.format( "bytes=%d-%d", start, start + count - 1 );
                    conn.addRequestProperty( "Range", range );
                    conn.connect();
                    int resp_code = conn.getResponseCode();
                    if ( resp_code == HttpURLConnection.HTTP_OK ||
                         resp_code == HttpURLConnection.HTTP_PARTIAL )
                    {
                        c.available = Integer.parseInt( conn.getHeaderField( "Content-Length" ) );
                        if ( c.available > 0 )
                        {
                            DataInputStream dis = new DataInputStream( conn.getInputStream() );
                            dis.readFully( c.data, 0, c.available );
                        }
                    }
                    conn.disconnect();
                }
                catch ( IOException ex ) { CLogger.log( ex.toString() ); }
            }
        }
    }
    
    public HttpDnldRangeReader( String url )
    {
        super( url );
    }
}
