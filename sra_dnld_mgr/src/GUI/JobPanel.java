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

import data.CLogger;
import data.TimeDiff;
import job.JobState;
import job.JobData;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import job.JobConsumerRunner;
import job.ProgressListenerInterface;
import job.StateAndProgressEvent;
import job.StateAndProgressNotifier;
import job.StateAndProgressType;

class JobButtons extends JPanel 
    implements ActionListener, Runnable, ProgressListenerInterface
{
    static final long serialVersionUID = 1;
    
    private final JobPanel parent;
    private final JButton Start, Pause, Reset, Edit, Delete;
    private final JobData job;
    private final JobConsumerRunner runner;
    private final JobDeleteEvent jde;
            
    @Override public void actionPerformed( ActionEvent ae )
    {
        JButton b = (JButton) ae.getSource();
        parent.clicked();
        if ( b == Start ) start();
        else if ( b == Pause ) pause();
        else if ( b == Reset ) reset();
        else if ( b == Edit ) edit();
        else if ( b == Delete ) delete_job();
    }
    
    public final void update_buttons_states()
    {
        JobState js = job.get_state();
        switch ( js )
        {
            case INVALID    : Start.setEnabled( false );
                              Pause.setEnabled( false );
                              Reset.setEnabled( false );
                              Edit.setEnabled( true );
                              Delete.setEnabled( true );
                              break;

            case READY      : Start.setEnabled( true );
                              Pause.setEnabled( false );
                              Reset.setEnabled( job.get_progress() > 0 );
                              Edit.setEnabled( true );
                              Delete.setEnabled( true );                             
                              break;

            case RUNNING    : Start.setEnabled( false );
                              Pause.setEnabled( true );
                              Reset.setEnabled( true );
                              Edit.setEnabled( false );
                              Delete.setEnabled( false );                              
                              break;

            case PAUSED     : Start.setEnabled( true );
                              Pause.setEnabled( false );
                              Reset.setEnabled( true );
                              Edit.setEnabled( false );
                              Delete.setEnabled( false );                              
                              break;

            case ERROR      : Start.setEnabled( false );
                              Pause.setEnabled( false );
                              Reset.setEnabled( false );
                              Edit.setEnabled( true );
                              Delete.setEnabled( true );                              
                              break;

            case DONE       : Start.setEnabled( false );
                              Pause.setEnabled( false );
                              Reset.setEnabled( true );
                              Edit.setEnabled( true );
                              Delete.setEnabled( true );                              
                              break;
        }
    }
            
    public boolean start()
    {
        CLogger.logfmt( "start job %s", job.get_short_source() );
        Start.setEnabled( false );
        boolean res = runner.start();
        if ( !res )
        {
            JOptionPane.showMessageDialog( this, "this job cannot be started, the accession is invalid",
                        "Error", JOptionPane.INFORMATION_MESSAGE );
            Start.setEnabled( true );
        }
        return res;
    }
    
    public void pause()
    {
        CLogger.logfmt( "pause job %s", job.get_short_source() );
        runner.pause();
    }
    
    public void reset()
    {
        CLogger.logfmt( "reset job %s", job.get_short_source() );
        runner.reset();
    }
    
    public void edit()
    {
        CLogger.logfmt( "edit job %s", job.get_short_source() );
        boolean res = JobWindow.edit( job );
        if ( res ) res = job.store(); /* edit does not store, the caller has to */
        if ( res ) update_buttons_states();
    }

    public void delete_job()
    {
        CLogger.logfmt( "delete job %s", job.get_short_source() );        
        if ( !is_running() ) { jde.delete_job( job ); }
    }
    
    /* this runnable has been notified about a change of job-state */
    @Override public void run()
    {
        update_buttons_states();
    }

    @Override  public void on_state_progress_event( final StateAndProgressEvent ev )
    {
        if ( ev.type == StateAndProgressType.STATE )
            update_buttons_states();
    }

    public boolean has_name( String name ) { return name.equals( job.get_short_source() ); }
    public boolean is_running() { return( job.get_state() == JobState.RUNNING ); }
    public boolean is_done() { return( job.get_state() == JobState.DONE ); }
    
    private JButton add_btn( final ImageIcon icon, final String txt )
    {
        JButton b = ResourceImages.make_img_button( icon, txt, 2, this );
        b.setToolTipText( txt );
        add( b, BorderLayout.WEST );
        return b;
    }

    public JobButtons( final JobPanel parent,
                       final JobData job,
                       final JobConsumerRunner runner,
                       final JobDeleteEvent jde )
    {
        super();
        this.parent = parent;
        this.job = job;
        this.runner = runner;
        this.jde = jde;
        setOpaque( false );
        
        Start  = add_btn( ResourceImages.get_start_img(), "Start job" );
        Pause  = add_btn( ResourceImages.get_pause_img(), "Pause job" );
        Reset  = add_btn( ResourceImages.get_stop_img(), "Reset job" );
        Edit   = add_btn( ResourceImages.get_settings_img(), "Edit job" );
        Delete = add_btn( ResourceImages.get_del_img(), "Delete job" );

        update_buttons_states();
   }
}

class EventRelay implements Runnable
{
    private Runnable r1, r2;
    
    @Override public void run()
    {
        if ( r1 != null ) r1.run();
        if ( r2 != null ) r2.run();
    }

    public void set_relay1( Runnable r ) { r1 = r; }
    public void set_relay2( Runnable r ) { r2 = r; }
    
    public EventRelay( Runnable r1, Runnable r2 )
    {
        this.r1 = r1;
        this.r2 = r2;
    }
}

class JobPanelMouseListener extends MouseAdapter
{
    private final JobPanel panel;
    @Override public void mouseClicked( MouseEvent e ) { panel.clicked(); }
    public JobPanelMouseListener( final JobPanel panel ) { this.panel = panel; }
}

public class JobPanel
    extends JPanel
    implements ProgressListenerInterface
{
    static final long serialVersionUID = 1;

    private final MainWindow parent;
    private final JobButtons buttons;
    private final ProgressPanel progress_panel;
    private final StateAndProgressNotifier notifier;
    private boolean active;
    
    // needed by MainWindow to count how many jobs are running
    public boolean is_running() { return buttons.is_running(); }
    public boolean is_done() { return buttons.is_done(); }
    
    // needed by MainWindow to start/pause all
    public boolean start() { return buttons.start(); }
    public void pause() { buttons.pause(); }
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
            StatusBar.set_max( maximum);
            StatusBar.set_pro( progress );
            StatusBar.set_time( elapsed );
            StatusBar.set_rpm( TimeDiff.calc_rpm( elapsed, progress ) );
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

    @Override public void on_state_progress_event( final StateAndProgressEvent ev )
    {
        switch( ev.type )
        {
            case PROGRESS : progress_panel.set_progress( ev.value );
                            if ( active )
                            {
                                StatusBar.set_pro( ev.value );
                                StatusBar.set_time( ev.elapsed_time );
                                StatusBar.set_rpm( TimeDiff.calc_rpm( ev.elapsed_time, ev.value ) );
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
        notifier.add_progress_listener( this );
        
        /* the download-task maintains a thread to run the job-loop */
        JobConsumerRunner runner = new JobConsumerRunner( job, notifier );
        
        /* a panel with the Start, Pause, Stop and Edit - buttons */
        buttons = new JobButtons( this, job, runner, jde );
        
        /* the button-panel also wants to be notified about state-changes */
        notifier.add_state_listener( buttons );
        add( buttons, BorderLayout.LINE_END );
        
        this.addMouseListener( new JobPanelMouseListener( this ) );
    }
}
