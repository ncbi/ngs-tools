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

import data.CLogger;

/*----------------------------------------------------
    This is the common Thread-Class for
    - JobThreadBio  ( writes FASTA/FASTQ )
    - JobThreadDnld ( downloads the Accession )
-----------------------------------------------------*/
class JobConsumerThread extends Thread
{
    protected final JobThreadData data;
    protected long current, maximum, run_time, total_elapsed_time;
    protected volatile boolean exit_flag;
    private boolean finished;
    
    public boolean is_finished() { return finished; }
    public void set_finished() { finished = true; }
    
    public void sleep_for_ms( long time_to_sleep )
    {
        try
        {
            Thread.sleep( time_to_sleep );
        }
        catch ( InterruptedException iex )
        {
        }
    }

    /* can be called from the GUI to stop the job
        is called by JobRunner.stop()  */
    public void request_exit() { exit_flag = true; }

    public void update_maximum( long value )
    {
        if ( maximum != value )
        {
            maximum = value;
            data.job.set_max( maximum, true );
            data.notifier.put_maximum( maximum, total_elapsed_time );
        }
    }
    
    public boolean can_start()
    {
        boolean res = data.job.is_runnable();
        if ( res )
        {
            current = data.job.get_progress();
            maximum = data.job.get_max();
            res = ( maximum == 0 || current < maximum );
            if ( !res )
                CLogger.logfmt( "job( %s ).can_start: failed",
                        data.job.get_short_source() );
        }
        else
        {
            CLogger.logfmt( "job( %s ).can_start: not runnable",
                    data.job.get_short_source() );
        }
        return res;
    }

    private void change_state_and_notify( final JobState new_state )
    {
        JobState prev_state = data.job.get_state();
        if ( data.job.change_job_state( new_state ) )
            data.notifier.put_state( prev_state, new_state );
    }
    
    public boolean perform_start()
    {
        JobState prev_state = data.job.get_state();
        data.job.change_job_state( JobState.RUNNING );
        data.notifier.put_state( prev_state, data.job.get_state() );
        if ( maximum == 0 )
            data.notifier.put_start();
        run_time = System.currentTimeMillis();
        return true;
    }    
    
    public boolean perform_step()
    {
        long now = System.currentTimeMillis();
        long elapsed = ( now - run_time );
        total_elapsed_time += elapsed;
        data.job.set_runtime( total_elapsed_time );
        run_time= now;
    
        if ( ( maximum > 0 && current > maximum ) || is_finished() ) current = maximum;

        data.job.set_progress( current, true );
        data.notifier.put_progress( current, total_elapsed_time );
        return true;
    }

    public void perform_done()
    {
        if ( data.job.get_state() == JobState.RUNNING )
        {
            change_state_and_notify( exit_flag ? JobState.PAUSED : JobState.DONE );
            if ( maximum == 0 )
                data.notifier.put_stop();
        }
    }

    @Override public void run()
    {
        if ( can_start() )
        {
            boolean running = perform_start();    // <======
            while ( running && !exit_flag )
            {
                running = ( data.job.get_state() == JobState.RUNNING );
                if ( running ) running = perform_step(); // <======
            }
            perform_done(); // <======
        }
    }

    public JobConsumerThread( final JobThreadData data )
    {
        super();
        this.data = data;
        current = 0;
        maximum = 0;
        run_time = 0;
        total_elapsed_time = data.job.get_runtime();
        exit_flag = false;
        finished = false;
    }
}
