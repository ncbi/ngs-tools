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
import data.Settings;

public class HttpDnldCGIResolver extends HttpDnldUrlWrapper
{
    private final int v_maj, v_min;

    private class Lines
    {
        String[] lines;
        int count;

        public Lines( int n )
        {
            lines = new String[ n ];
            count = 0;
        }
    }
    
    private Lines post( String parameters )
    {
        Lines res  = new Lines( 5 );
        if ( valid() )
        {
            HttpURLConnection conn = make_connection( "POST" );
            conn.setUseCaches( false );
            conn.setDoOutput( true );
            try ( DataOutputStream wr = new DataOutputStream( conn.getOutputStream() ) )
            {
                wr.writeBytes( parameters );

                InputStream is = conn.getInputStream();
                try ( BufferedReader rd = new BufferedReader( new InputStreamReader( is ) ) )
                {
                    boolean done = false;                   
                    while ( !done )
                    {
                        done = ( res.count >= res.lines.length );
                        if ( !done )
                        {
                            res.lines[ res.count ] = rd.readLine();
                            done = ( res.lines[ res.count ] == null );
                            if ( !done ) res.count++;
                        }
                    }
                }
            }
            catch ( IOException ex ) { CLogger.log( ex.toString() ); }
        }
        return res;
    }
    
    private boolean split_1_0( String s, HttpDnldResolved r )
    {
        String[] parts = s.split( "\\|" );        
        if ( parts.length > 0 ) r.accession = parts[ 0 ];
        if ( parts.length > 1 ) r.dnld_ticket = parts[ 1 ];
        if ( parts.length > 2 ) r.url = parts[ 2 ];
        if ( parts.length > 3 ) r.result_code = Integer.parseInt( parts[ 3 ] );
        if ( parts.length > 4 ) r.msg = parts[ 4 ];
        return ( r.result_code == 200 );
    }
    
    private boolean split_1_1( String s, HttpDnldResolved r )
    {
        String[] parts = s.split( "\\|" );
        if ( parts.length > 0 ) r.accession = parts[ 0 ];
        if ( parts.length > 1 ) r.obj_id = parts[ 1 ];
        if ( parts.length > 2 ) r.name = parts[ 2 ];
        if ( parts.length > 3 )
        {
            try { r.size = Long.parseLong( parts[ 3 ] ); }
            catch( NumberFormatException e ) { r.size = 0; }
        }
        if ( parts.length > 4 ) r.mod_date = parts[ 4 ];
        if ( parts.length > 5 ) r.md5 = parts[ 5 ];
        if ( parts.length > 6 ) r.dnld_ticket = parts[ 6 ];
        if ( parts.length > 7 ) r.url = parts[ 7 ];
        if ( parts.length > 8 ) r.result_code = Integer.parseInt( parts[ 8 ] );
        if ( parts.length > 9 ) r.msg = parts[ 9 ];
        return ( r.result_code == 200 );
    }
    
    public boolean resolve( String accession, HttpDnldResolved r )
    {
        boolean found = false;
        if ( valid() )
        {
            StringBuilder sb = new StringBuilder();
            
            sb.append( String.format( "version=%d.%d", v_maj, v_min ) );
            sb.append( '&' );
            
            sb.append( String.format( "acc=%s", accession ) );
            sb.append( '&' );

            sb.append( String.format( "accept-proto=http" ) );
            
            Lines lines = post( sb.toString() );
            if ( lines.count > 1 )
            {
                switch ( lines.lines[ 0 ] )
                {
                    case "#1.0":
                        found = split_1_0( lines.lines[ 1 ], r );
                        break;
                    case "#1.1":
                        found = split_1_1( lines.lines[ 1 ], r );
                        break;
                }
            }
        }
        return found;
    }
    
    public HttpDnldCGIResolver( String url )
    {
        super( url );
        this.v_maj = 1;
        this.v_min = 1;
    }
    
    public HttpDnldCGIResolver()
    {
        super( Settings.getInstance().get_resolver() );
        this.v_maj = 1;
        this.v_min = 1;
        Settings settings = Settings.getInstance();
        set_conn_timeout( settings.get_conn_timeout() );
        set_read_timeout( settings.get_read_timeout() );
        set_agent( settings.get_user_agent() );
    }
}
