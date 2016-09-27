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

import Bio.*;
import data.CLogger;
import java.io.*;

public class JobConsumerThreadBio extends JobConsumerThread
{
    private final static long START_VALUE = 1;
    private final static int MAX_BIO_READS_IN_Q = 200;
    
    private BioDumper producer;
    private BufferedWriter writer0;
    private BufferedWriter writer1;    
    private final BioRunState rs;
    private final JobSubFormat subformat;
    private long step_time;
    private long rejected;
    
    private BufferedWriter make_writer( String filename, boolean append )
    {
        BufferedWriter res = null;
        try
        {
            if ( append )
                CLogger.logfmt( "job >%s< re-opening outputfile: %s " ,
                        data.job.get_short_source(), filename );
            else
                CLogger.logfmt( "job >%s< creating outputfile: %s " ,
                        data.job.get_short_source(), filename );
            res = new BufferedWriter( new FileWriter( filename, append ), 10000 );
        }
        catch ( FileNotFoundException ex )
        {
            if ( append )
                CLogger.logfmt( "job >%s< error re-opening outputfile: %s " ,
                        data.job.get_short_source(), ex.toString() );
            else
                CLogger.logfmt( "job >%s< error creating outputfile: %s " ,
                        data.job.get_short_source(), ex.toString() );
        }
        catch ( IOException ex )
        {
            if ( append )
                CLogger.logfmt( "job >%s< error re-opening outputfile: %s " ,
                        data.job.get_short_source(), ex.toString() );
            else
                CLogger.logfmt( "job >%s< error creating outputfile: %s " ,
                        data.job.get_short_source(), ex.toString() );
        }
        return res;
    }

    private boolean write( BioRecord rec )
    {
        boolean res = true;
        try
        {
            if ( !rec.is_empty() )
                writer0.write( rec.get() );
        }
        catch( IOException ex )
        {
            CLogger.logfmt( "job >%s< error performing step: %s",
                    data.job.get_short_source(), ex.toString() );
            res = false;
        }
        return res;
    }
    
    private boolean write_splitted( BioRecord rec, int frag_nr )
    {
        boolean res = true;
        try
        {
            if ( !rec.is_empty() )
            {
                switch( frag_nr )
                {
                    case 0  : writer0.write( rec.get() ); break;
                    default : writer1.write( rec.get() ); break;
                }
            }
        }
        catch( IOException ex )
        {
            CLogger.logfmt( "job >%s< error performing step: %s",
                    data.job.get_short_source(), ex.toString() );
            res = false;
        }
        return res;
    }
    
    @Override public boolean can_start()
    {
        boolean res = super.can_start();
        String source = data.job.get_short_source();

        if ( res )
        {
            res = gov.nih.nlm.ncbi.ngs.NGS.isSupported();
            if ( !res )
                CLogger.logfmt( "job >%s< NGS not supported" , source );
        }
        if ( res )
        {
            res = gov.nih.nlm.ncbi.ngs.NGS.isValid( data.job.get_full_source());
            if ( !res )
                CLogger.logfmt( "job >%s< has an invalid NGS accession" , source );
        }
        
        if ( res )
        {
            /* the correct start-value is needed by the dumper-creation! */
            current = ( current <= START_VALUE ) ? START_VALUE : current + 1;
            producer = BioDumperFactory.make_dumper( data.job, current );
            res = ( producer != null );
        }
        
        if ( res )
        {
            update_maximum( producer.get_read_count_or_size() );
            producer.set_run_state( rs );
            
            boolean append_mode = ( current > START_VALUE );
            if ( data.job.get_subformat().equals( JobSubFormat.FRAGMENTS_SPLITED ) )
            {
                writer0 = make_writer( data.job.get_output_filename( 1 ), append_mode );
                writer1 = make_writer( data.job.get_output_filename( 2 ), append_mode );
            }
            else
            {
                writer0 = make_writer( data.job.get_output_filename(), append_mode );
            }
        }
        
        if ( res )
        {
            rs.clear();
            Thread t = new Thread( producer ); /* because the dumper is a runnable! */
            t.start();
        }

        return res;
    }

    /* otherwise the super-routine will not be called!!! */
    @Override public boolean perform_start()
    {
        return super.perform_start();
    }

    private boolean write_bio_read( final BioRead read )
    {
        boolean res = true;
        int i = 0;
        while ( res && !read.is_empty() )
        {
            BioRecord rec = read.get_record();
            if ( rec != null )
            {
                switch ( subformat )
                {
                    case SPOT               : res &= write( rec ); break;
                    case FRAGMENTS          : res &= write( rec ); break;
                    case FRAGMENTS_SPLITED  : res &= write_splitted( rec, i ); break;
                }
                rs.put_record_back( rec );
            }
            else
                rejected++;
            i++;
        }
        return res;
    }
    
    private boolean run_for( final long mili_secs )
    {
        boolean res = true;
        boolean done = false;
        long start = System.currentTimeMillis();
        long elapsed = 0;
        long reads_handled = 0;
        while ( res && !done && ( elapsed < mili_secs ) )
        {
            BioRead read = rs.get();
            if ( read == null )
            {
                // the producer thread has not yet/in the mean time produced spot's
                done = !rs.is_running();
                if ( done )
                {
                    set_finished();
                }
                else
                {
                    // maybe sleep a little bit, if we are not done ...
                    sleep_for_ms( 1 );
                }
            }
            else
            {
                res = write_bio_read( read );
                if ( res ) current = read.read_id;
                rs.put_to_stock( read );
                reads_handled++;
            }
            elapsed = ( System.currentTimeMillis() - start );
        }

        if ( reads_handled > 2000 )
        {
            // there are too many reads per time-slot, we can reduce
            // the time-slot but not below 125 ms
            if ( step_time > 175 )
                step_time -= 50;
        }
        else if ( reads_handled < 1000 )
        {
            // there are not anoth reads per time-slot, we can increase
            // the time-slot but not above 1000 ms
            if ( step_time < 900 )
                step_time += 50;
        }
        return res;
    }

    /* one step means up to PER_STEP spots before we update and save */
    @Override public boolean perform_step()
    {
        boolean res = rs.is_running();
        if ( !res ) res = rs.has_entries(); // notice the !res
        if ( res ) res = run_for( step_time );
        if ( res ) res = super.perform_step(); // updates the current in the job...
        return res;
    }

    private void close_writer( final BufferedWriter writer, final int nr )
    {
        try
        {
            writer.close();
            System.out.printf( "closing writer for '%s'\n", data.job.get_short_source() );
        }
        catch ( IOException ex )
        {
            CLogger.log( String.format( "job >%s< error closing stream #%d: %s",
                    data.job.get_short_source(), nr, ex.toString() ) );
        }

    }
    
    @Override public void perform_done()
    {
        super.perform_done();
        data.job.set_rejected( rejected );
        data.job.store();
        rs.set_running( false );
        if ( writer0 != null ) close_writer( writer0, 0 );
        if ( writer1 != null ) close_writer( writer1, 1 );
        
        /* to force de-allocation of the internal read-collection! */
        if ( producer != null )
        {
            producer = null;
            System.gc();
            System.runFinalization();
        }
        rs.clear();
    }
    
    public JobConsumerThreadBio( final JobThreadData data )
    {
        super( data );
        rs = new BioRunState( MAX_BIO_READS_IN_Q );
        subformat = data.job.get_subformat();
        step_time = 250;
        rejected = 0;
        producer = null;
        writer0 = null;
        writer1 = null;        
    }
}
