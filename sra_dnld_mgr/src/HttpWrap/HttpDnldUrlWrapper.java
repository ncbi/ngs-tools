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

public class HttpDnldUrlWrapper
{
    private URL url;
    private String agent;
    private int conn_timeout;
    private int read_timeout;
    
    public final String dflt_agent = "dnld-mgr";
            
    public boolean valid() { return url != null; }

    public HttpURLConnection make_connection( String RequestMethod )
    {
        HttpURLConnection conn = null;
        try
        {
            conn = ( HttpURLConnection )url.openConnection();
            conn.setConnectTimeout( conn_timeout );
            conn.setReadTimeout( read_timeout );
            conn.setRequestMethod( RequestMethod );
            conn.addRequestProperty( "User-Agent", agent );
        }
        catch ( IOException ex ) { CLogger.log( ex.toString() ); }
        return conn;
    }
    
    public void set_agent( String agent ) { this.agent = agent; }
    public void set_conn_timeout( int value ) { conn_timeout = value; }
    public void set_read_timeout( int value ) { read_timeout = value; }
    
    private void set_defaults()
    {
        url = null;        
        agent = dflt_agent;
        conn_timeout = 5000;
        read_timeout = 5000;
    }
    
    public HttpDnldUrlWrapper( String url )
    {
        set_defaults();
        
        try { this.url = new URL( url ); }
        catch ( MalformedURLException ex ) { CLogger.log( ex.toString() ); }
    }

}
