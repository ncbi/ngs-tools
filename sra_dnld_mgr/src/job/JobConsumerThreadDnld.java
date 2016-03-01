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

package job;

import HttpWrap.*;
import data.CLogger;
import data.Settings;
import java.io.*;

public class JobConsumerThreadDnld extends JobConsumerThread
{
    private HttpDnldRangeReader range_reader;
    private FileOutputStream out_stream;
    private final HttpDnldCGIResolver resolver;
    private final HttpDnldResolved remote_location;
    private final HttpDnldContent content;
    private final int blocksize;

    private boolean make_output_stream()
    {
        boolean res = false;
        try
        {
            boolean append_mode = ( data.job.get_progress() > 0 );
            if ( !data.job.does_downloadpath_exist() )
            {
                res = data.job.create_downloadpath();
                CLogger.logfmt( "creating path: %s = %s" ,
                        data.job.get_downloadpath(),
                        res ? "OK" : "ERROR" );
            }
            else
                res = true;
            if ( res )
            {
                String output_filename = data.job.get_output_filename();
                out_stream = new FileOutputStream( output_filename, append_mode );
                res = ( out_stream != null );
                CLogger.logfmt( "creating/opening file: %s = %s" ,
                        output_filename,
                        res ? "OK" : "ERROR" );
            }
        }
        catch ( IOException ex ) { CLogger.log( ex.toString() ); }
        return res;
    }

    private void update_md5( String value )
    {
        if ( !value.equals( data.job.get_md5() ) )
            data.job.set_md5( value, true );
    }
    
    @Override public boolean perform_start()
    {
        boolean res = super.perform_start();
        if ( res )
        {
            if ( resolver.resolve( data.job.get_full_source(), remote_location ) )
            {
                Settings settings = Settings.getInstance();
                
                long remote_size = remote_location.get_size();
                
                CLogger.logfmt( "job >%s< resolved to >%s< ( size: %d )" ,
                        data.job.get_short_source(), remote_location.get_url(), remote_size );
                
                // update the display about the md5-sum ( not all download-urls have one )
                update_md5( remote_location.get_md5() );
                
                range_reader = new HttpDnldRangeReader( remote_location.get_url() );
                range_reader.set_conn_timeout( settings.get_conn_timeout() );
                range_reader.set_read_timeout( settings.get_read_timeout() );
                range_reader.set_agent( settings.get_user_agent() );

                if ( remote_size == 0 )
                {
                    remote_size = range_reader.get_size();
                    CLogger.logfmt( "job >%s< discovered size: %d" ,
                            data.job.get_short_source(), remote_size );
                }
                
                if ( remote_size != data.job.get_max() )
                    data.job.set_max( remote_size, true );
                
                // update the display about the size of the download
                update_maximum( remote_size );
                
                make_output_stream();
            }
            else
                CLogger.logfmt( "error resolving >%s<" , data.job.get_full_source() );
        }
        return res;
    }

    @Override public boolean perform_step()
    {
        boolean res = true;
        if ( range_reader != null )
        {
            range_reader.read( current, blocksize, content );
            if ( out_stream != null )
            {
                try
                {
                    if ( content.available == 0 )
                    {
                        CLogger.logfmt( "perform_step( %s ).content = %d" ,
                                data.job.get_full_source(), content.available );
                        res = false;
                    }
                    else
                    {
                        out_stream.write( content.data, 0, content.available );
                        current += content.available;
                    }
                }
                catch ( IOException ex )
                {
                    CLogger.log( ex.toString() );
                    res = false;
                }
            }
            else
                CLogger.logfmt( "perform_step( %s ).out_stream is null" , data.job.get_full_source() );
        }
        res &= super.perform_step();
        return res;
    }

    @Override public void perform_done()
    {
        if ( out_stream != null )
        {
            try { out_stream.close(); }
            catch ( IOException ex ) { CLogger.log( ex.toString() ); }
        }
        super.perform_done();
   }
    
    public JobConsumerThreadDnld( final JobThreadData data )
    {
        super( data );
        blocksize = Settings.getInstance().get_dnld_block_size();
        remote_location = new HttpDnldResolved();
        content = new HttpDnldContent( blocksize );
        resolver = new HttpDnldCGIResolver();
        range_reader = null;
        out_stream = null;
        CLogger.logfmt( "HttpDownload created for >%s<" ,
                data.job.get_short_source() );
    }
}