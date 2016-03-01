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
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import job.JobConsumerRunner;
import job.JobData;
import job.JobState;
import job.ProgressListenerInterface;
import job.StateAndProgressEvent;
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
