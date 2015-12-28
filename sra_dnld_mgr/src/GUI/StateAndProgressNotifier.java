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
package GUI;

import data.QConnect;
import java.awt.EventQueue;

public class StateAndProgressNotifier
    extends QConnect< StateAndProgressEvent >
    implements Runnable
{
    private final ProgressPanel progress_panel;
    private final Runnable on_state_changed;

    @Override public StateAndProgressEvent make_entry() { return new StateAndProgressEvent(); }
    
    public void put_ev( final StateAndProgressType type, final long value )
    {
        StateAndProgressEvent ev = get_from_stock();
        if ( ev != null )
        {
            ev.type = type;
            ev.value = value;
            put( ev );
            EventQueue.invokeLater( this ); // call run() in the context of the gui-task
        }
    }

    /* these are called in the context of the thread */
    final public void put_progress( long value ) { put_ev( StateAndProgressType.PROGRESS, value ); }
    final public void put_state() { put_ev( StateAndProgressType.STATE, 0 ); }
    final public void put_maximum( long value ) { put_ev( StateAndProgressType.MAXIMUM, value ); }
    final public void put_start() { put_ev( StateAndProgressType.START, 0 ); }
    final public void put_stop() { put_ev( StateAndProgressType.STOP, 0 ); }

    @Override public void run()
    {
        while( has_entries() )
        {
            StateAndProgressEvent ev = get();
            if ( ev != null )
            {
                switch( ev.type )
                {
                    case PROGRESS : progress_panel.set_progress( ev.value ); break;
                    case MAXIMUM  : progress_panel.set_maximum( ev.value ); break;
                    case STATE    : on_state_changed.run(); break;
                    case START    : progress_panel.start(); break;
                    case STOP     : progress_panel.stop(); break;
                }
                put_to_stock( ev );
            }
        }
    }
   
    public StateAndProgressNotifier( final ProgressPanel p,
                                     final Runnable on_state_changed,
                                     final long maximum,
                                     final long progress )
    {
        super( 50 );
        
        progress_panel = p;
        this.on_state_changed = on_state_changed;
        put_maximum( maximum );
        put_progress( progress );
    }
}
