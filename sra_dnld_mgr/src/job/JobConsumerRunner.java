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

import Bio.BioSpec;
import Bio.BioVerifier;
import data.*;

public class JobConsumerRunner
{
    private final JobThreadData data;
    private JobConsumerThread thread;
    
    public boolean start()
    {
        boolean res = false;
        if ( thread == null )
        {
            BioVerifier v = new BioVerifier();
            BioSpec spec = v.verify( data.job.get_full_source() );
            if ( spec.is_valid() )
            {
                JobFormat job_fmt = data.job.get_format();
                switch( job_fmt )
                {
                    case DOWNLOAD   : thread = new JobConsumerThreadDnld( data ); break;
                    case FASTA      : thread = new JobConsumerThreadBio( data ); break;
                    case FASTQ      : thread = new JobConsumerThreadBio( data ); break;
                    default : CLogger.logfmt( "job( %s ).start: invalid output-format: %s",
                        data.job.get_short_source(), job_fmt.toString() );
                }
                res = ( thread != null );
                if ( res ) thread.start();
            }
            else
            {
                CLogger.logfmt( "job( %s ).start: verify source failed",
                        data.job.get_short_source() );
            }
        }
        else
        {
            CLogger.logfmt( "job( %s ).start: thread did already exist",
                    data.job.get_short_source() );
        }
        return res;
    }

    private boolean stop()
    {
        boolean res = false;
        if ( thread != null )
        {
            thread.request_exit();
            try
            {
                thread.join();
                res = true;
                thread = null;
            }
            catch( InterruptedException ex )
            {
                CLogger.logfmt( "job >%s< error joining thread: %s " ,
                    data.job.get_short_source(), ex.toString() );
            }
        }
        return res;
    }
    
    public boolean pause()
    {
        boolean res = stop();
        if ( res )
        {
            /*
            if ( data.job.change_job_state( JobState.READY ) )
                data.notifier.set_state_changed( JobState.READY );
            */
        }
        return res;
    }
    
    public void reset()
    {
        stop();
        data.notifier.put_progress( 0, 0 );
        data.job.set_runtime( 0 );
        data.job.set_rejected( 0 );
        data.job.set_progress( 0, true );
        // no if (), because we want the notifier fired, even if the state has not changed...
        // that is necessary if we are already in READY, but now the curr-row is zero
        
        JobState prev_state = data.job.get_state();
        data.job.change_job_state( JobState.READY );
        data.notifier.put_state( prev_state, data.job.get_state() );
    }
    
    public JobConsumerRunner( JobData job,
                      StateAndProgressNotifier notifier )
    {
        data = new JobThreadData( job, notifier );
        thread = null;
    }
}
