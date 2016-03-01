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

import data.QConnect;
import java.awt.EventQueue;
import java.util.ArrayList;
import java.util.List;

public class StateAndProgressNotifier
    extends QConnect< StateAndProgressEvent >
    implements Runnable
{
    private final List<ProgressListenerInterface> progress_listeners;
    private final List<ProgressListenerInterface> state_listeners;
    
    private long progress;
    private long maximum;
    private long elapsed_time;
    
    public long get_progress() { return progress; }
    public long get_maximum() { return maximum; }
    public long get_elapsed_time() { return elapsed_time; }
    
    public void add_progress_listener( ProgressListenerInterface listener )
    {
        progress_listeners.add( listener );
    }

    public void add_state_listener( ProgressListenerInterface listener )
    {
        state_listeners.add( listener );
    }

    @Override public StateAndProgressEvent make_entry() { return new StateAndProgressEvent(); }

    private void put_ev( final StateAndProgressType type )
    {
        StateAndProgressEvent ev = get_from_stock();
        if ( ev != null )
        {
            ev.type = type;
            ev.value = 0;
            ev.elapsed_time = elapsed_time;
            put( ev );
            EventQueue.invokeLater( this ); // call run() in the context of the gui-task
        }
    }

    private void put_ev( final StateAndProgressType type, final long value )
    {
        StateAndProgressEvent ev = get_from_stock();
        if ( ev != null )
        {
            ev.type = type;
            ev.value = value;
            ev.elapsed_time = elapsed_time;
            put( ev );
            EventQueue.invokeLater( this ); // call run() in the context of the gui-task
        }
    }

    /* these are called in the context of the thread */
    final public void put_progress( long value, long elapsed )
    {
        elapsed_time = elapsed;
        put_ev( StateAndProgressType.PROGRESS, value );
    }
    
    final public void put_maximum( long value, long elapsed )
    {
        elapsed_time = elapsed;
        put_ev( StateAndProgressType.MAXIMUM, value );
    }

    final public void put_start() { put_ev( StateAndProgressType.START ); }
    final public void put_stop() { put_ev( StateAndProgressType.STOP ); }

    final public void put_state( final JobState prev_state, final JobState new_state )
    {
        put_ev( StateAndProgressType.STATE );
    }
    
    @Override public void run()
    {
        while( has_entries() )
        {
            StateAndProgressEvent ev = get();
            if ( ev != null )
            {
                // notify all listeners...
                
                if ( ev.type == StateAndProgressType.STATE )
                {
                    for ( ProgressListenerInterface li : state_listeners )
                        li.on_state_progress_event( ev );
                }
                else
                {
                    for ( ProgressListenerInterface li : progress_listeners )
                        li.on_state_progress_event( ev );
                }
                
                switch( ev.type )
                {
                    case PROGRESS : this.progress = ev.value; break;
                    case MAXIMUM  : this.maximum = ev.value; break;
                }
                put_to_stock( ev );
            }
        }
    }
   
    public StateAndProgressNotifier( final long maximum,
                                     final long progress,
                                     final long elapsed )
    {
        super( 50 );
        progress_listeners = new ArrayList<>();
        state_listeners = new ArrayList<>();
        this.progress = progress;
        this.maximum = maximum;
        this.elapsed_time = elapsed;
        put_maximum( maximum, elapsed );
        put_progress( progress, elapsed );
    }
}
