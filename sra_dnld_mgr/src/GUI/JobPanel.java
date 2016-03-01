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

import data.TimeDiff;
import job.JobData;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import job.JobConsumerRunner;
import job.JobFormat;
import job.ProgressListenerInterface;
import job.StateAndProgressEvent;
import job.StateAndProgressNotifier;

public class JobPanel extends JPanel
{
    static final long serialVersionUID = 1;

    private final JobPanelMouseListener on_mouse_event = new JobPanelMouseListener();
    private final ProgressListener progress_event = new ProgressListener();
    
    private final MainWindow parent;
    private final JobButtons buttons;
    private final ProgressPanel progress_panel;
    private final StateAndProgressNotifier notifier;
    private boolean active;
    final JobData job;

    class JobPanelMouseListener extends MouseAdapter
    {
        @Override public void mouseClicked( MouseEvent e ) { clicked(); }
    }

    class ProgressListener implements ProgressListenerInterface
    {
        @Override public void on_state_progress_event( StateAndProgressEvent ev )
        { on_progress_event( ev ); }
    }
    
    // needed by MainWindow to count how many jobs are running
    public boolean is_running() { return buttons.is_running(); }
    public boolean is_done() { return buttons.is_done(); }
    
    // needed by MainWindow to start/pause all
    public boolean start() { return buttons.start(); }
    
    public void pause()
    {
        buttons.pause();
    }
    
    public void reset() { buttons.reset(); }
    
    // needeb by MainWindow to remove completed Jobs
    public boolean has_name( String name ) { return buttons.has_name( name ); }
   
    public final void set_active( boolean active )
    {
        this.active = active;
        setBackground( active ? Color.green : Color.lightGray );
        if ( active )
        {
            StatusBar.set_acc( getName() );
            
            long maximum = notifier.get_maximum();
            long progress = notifier.get_progress();
            long elapsed  = notifier.get_elapsed_time();
            StatusBar.set_is_dnld( job.get_format().equals( JobFormat.DOWNLOAD ) );
            StatusBar.set_max( maximum);
            StatusBar.set_pro( progress );
            StatusBar.set_time( elapsed );
            StatusBar.set_rpm_or_bpm( TimeDiff.calc_rpm( elapsed, progress ) );
            StatusBar.set_left( TimeDiff.calc_time_left( elapsed, progress, maximum ) );
        }
    }
    
    public boolean is_active() { return active; }
    
    public void clicked() { parent.on_panel_clicked( this ); }
    
    public void add_progress_listener( final ProgressListenerInterface listener )
    {
        notifier.add_progress_listener( listener );
    }

    public void add_state_listener( final ProgressListenerInterface listener )
    {
        notifier.add_state_listener( listener );
    }

    private void on_progress_event( final StateAndProgressEvent ev )
    {
        switch( ev.type )
        {
            case PROGRESS : progress_panel.set_progress( ev.value );
                            if ( active )
                            {
                                StatusBar.set_pro( ev.value );
                                StatusBar.set_time( ev.elapsed_time );
                                StatusBar.set_rpm_or_bpm( TimeDiff.calc_rpm( ev.elapsed_time, ev.value ) );
                                StatusBar.set_left( TimeDiff.calc_time_left( ev.elapsed_time, ev.value, notifier.get_maximum() ) );
                            }
                            break;
                
            case MAXIMUM  : progress_panel.set_maximum( ev.value );
                            if ( active )
                            {
                                StatusBar.set_max( ev.value );
                                StatusBar.set_time( ev.elapsed_time );
                            }
                            break;

            case START    : progress_panel.start(); break;
            case STOP     : progress_panel.stop(); break;
        }
    }
    
    public JobPanel( final MainWindow parent,
                     final JobData job,
                     final JobDeleteEvent jde )
    {
        super();
        
        this.parent = parent;
        this.job = job;
        this.setName( job.get_short_source() );

        setPreferredSize( new Dimension( 400, 50 ) );
        setMinimumSize( getPreferredSize() );
        setMaximumSize( new Dimension( Short.MAX_VALUE, 50 ) );
        setLayout( new BorderLayout( 5, 0 ) );

        set_active( false );
        setBorder( BorderFactory.createMatteBorder( 2, 2, 2, 2, Color.WHITE ) );
                
        /* on the left we have a lable with the name of the job */
        ALabel l = new ALabel( job.get_short_source(), job.get_bio_type().get_color() );
        add( l, BorderLayout.LINE_START );

        /* in the middle there is a panel with a label and a progress-bar */
        progress_panel = new ProgressPanel( job );
        add( progress_panel, BorderLayout.CENTER );

        /* we need this to notify the parent and the button-panel about state-changes */
        //EventRelay on_state_change = new EventRelay( done_signal, null );

        /* this class updates the progress-panel ( bar in percent and label )
           and propagates state-change-events originating from the job-thread */
        notifier = new StateAndProgressNotifier(
                job.get_max(), job.get_progress(), job.get_runtime() );
        notifier.add_progress_listener( progress_event );
        
        /* the download-task maintains a thread to run the job-loop */
        JobConsumerRunner runner = new JobConsumerRunner( job, notifier );
        
        /* a panel with the Start, Pause, Stop and Edit - buttons */
        buttons = new JobButtons( this, job, runner, jde );
        
        /* the button-panel also wants to be notified about state-changes */
        notifier.add_state_listener( buttons );
        add( buttons, BorderLayout.LINE_END );
        
        this.addMouseListener( on_mouse_event );
    }
}
