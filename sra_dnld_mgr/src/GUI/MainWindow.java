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

import Bio.BioSpec;
import job.*;
import data.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

class Main_Window_Action implements ActionListener
{
    public final MainWindow main_window;
    @Override public void actionPerformed( ActionEvent ae ) { }
    public Main_Window_Action( MainWindow main_window ) { this.main_window = main_window; }
}

class Edit_Global_Settings extends Main_Window_Action
{
    public Edit_Global_Settings( MainWindow main_window ) { super( main_window ); }
    @Override public void actionPerformed( ActionEvent ae ) { SettingsWindow.edit(); }
}

class Exit_Action extends Main_Window_Action
{
    public Exit_Action( MainWindow main_window ) { super( main_window ); }
    @Override public void actionPerformed( ActionEvent ae ) { main_window.on_exit(); }
}

class Stop_All_Action extends Main_Window_Action
{
    public Stop_All_Action( MainWindow main_window ) { super( main_window ); }
    @Override public void actionPerformed( ActionEvent ae ) { main_window.stop_all_threads(); }
}

class Start_All_Action extends Main_Window_Action
{
    public Start_All_Action( MainWindow main_window ) { super( main_window ); }
    @Override public void actionPerformed( ActionEvent ae ) { main_window.start_download_jobs(); }
}

class New_Job_Action extends Main_Window_Action
{
    public New_Job_Action( MainWindow main_window ) { super( main_window ); }
    @Override public void actionPerformed( ActionEvent ae ) { main_window.NewJob(); }
}

class CleanUp_Action extends Main_Window_Action
{
    public CleanUp_Action( MainWindow main_window ) { super( main_window ); }
    @Override public void actionPerformed( ActionEvent ae ) { main_window.CleanUp(); }
}

class KeyDispatcher implements KeyEventDispatcher
{
    public final MainWindow main_window;
    public KeyDispatcher( MainWindow main_window ) { this.main_window = main_window; }

    @Override public boolean dispatchKeyEvent( KeyEvent ke )
    {
        boolean res = false;
        if ( ke.getID() == KeyEvent.KEY_PRESSED )
            res = main_window.on_key( ke );
        return res;
    }
}

public class MainWindow extends JFrame
{
    static final long serialVersionUID = 1;
    
    private JToolBar toolbar;
    private JPanel jobs;
    private final JobDoneEvent job_done_event;
    private final JobDeleteEvent job_delete_event;
    private final KeyDispatcher key_dispatcher;
    
    public final void enable_key_processor( boolean enable )
    {
        KeyboardFocusManager kfm = KeyboardFocusManager.getCurrentKeyboardFocusManager();
        if ( enable )
            kfm.addKeyEventDispatcher( key_dispatcher );
        else
            kfm.removeKeyEventDispatcher( key_dispatcher );
    }
    
    public JobPanel get_active_panel()
    {
        JobPanel res = null;
        int i, n_panels = jobs.getComponentCount();
        for ( i = 0; res == null && i < n_panels; ++i )
        {
            JobPanel vj = ( JobPanel ) jobs.getComponent( i );
            if ( vj != null && vj.is_active() ) res = vj;
        }
        return res;
    }
    
    public int get_active_panel_id()
    {
        int res = -1;
        int i, n_panels = jobs.getComponentCount();
        for ( i = 0; res < 0 && i < n_panels; ++i )
        {
            JobPanel vj = ( JobPanel ) jobs.getComponent( i );
            if ( vj != null && vj.is_active() ) res = i;
        }
        return res;
    }

    public int get_panel_id( final JobPanel p )
    {
        int res = -1;
        int i, n_panels = jobs.getComponentCount();
        for ( i = 0; res < 0 && i < n_panels; ++i )
        {
            JobPanel vj = ( JobPanel ) jobs.getComponent( i );
            if ( vj == p ) res = i;
        }
        return res;
    }
    
    public void set_panel_active( int idx, boolean state )
    {
        JobPanel vj = ( JobPanel ) jobs.getComponent( idx );
        if ( vj != null ) vj.set_active( state );
    }
    
    public int nxt_idx( final int count, final int idx, final int fwd )
    {
        int res = idx + fwd;
        if ( res < 0 ) res = 0;
        if ( res >= count ) res = count -1;
        return res;
    }
    
    public boolean switch_panel( final int count, final int from, final int to )
    {
        boolean res = false;
        if ( from < 0 )
        {
            if ( count > 0 )
            {
                set_panel_active( 0, true );
                res = true;
            }
        }
        else if ( count > 1 )
        {
            set_panel_active( from, false );
            set_panel_active( to, true );
            res = true;
        }
        return res;
    
    }
    
    public boolean nxt_panel( final int fwd )
    {
        int count = jobs.getComponentCount();
        int idx = get_active_panel_id();
        return switch_panel( count, idx, nxt_idx( count, idx, fwd ) );
    }

    public boolean first_panel()
    {
        return switch_panel( jobs.getComponentCount(), get_active_panel_id(), 0 );
    }

    public boolean last_panel()
    {
        int count = jobs.getComponentCount();
        return switch_panel( count, get_active_panel_id(), count - 1 );
    }

    public boolean start_stop()
    {
        CLogger.log( "start_stop()" );
        
        JobPanel p = get_active_panel();
        boolean res = ( p != null );
        if ( res )
        {
            if ( p.is_running() )
                p.pause();
            else
            {
                if ( p.is_done() )
                    p.reset();
                else
                    p.start();
            }
        }
        return res;
    }
    
    public boolean on_key( final KeyEvent ke )
    {
        boolean res = false;
        switch( ke.getKeyCode() )
        {
            case KeyEvent.VK_UP        : res = nxt_panel( -1 ); break;
            case KeyEvent.VK_DOWN      : res = nxt_panel( +1 ); break; 
            case KeyEvent.VK_PAGE_UP   : res = nxt_panel( -7 ); break;
            case KeyEvent.VK_PAGE_DOWN : res = nxt_panel( +7 ); break;
            case KeyEvent.VK_HOME      : res = first_panel(); break;
            case KeyEvent.VK_END       : res = last_panel(); break;
            case KeyEvent.VK_SPACE     : res = start_stop(); break;
        }
        return res;
    }
    
    public void on_panel_clicked( final JobPanel p )
    {
        CLogger.log( "on_panel_clicked" );
        int active_id = get_active_panel_id();
        int new_panel_id = get_panel_id( p );
        if ( active_id != new_panel_id )
            switch_panel( jobs.getComponentCount(), active_id, new_panel_id );
    }
    
    private class JobDoneEvent
        implements Runnable, ProgressListenerInterface
    {
        @Override public void run()
        {
            if ( Settings.getInstance().get_autostart() )
                start_download_jobs();
        }

        @Override public void on_state_progress_event( final StateAndProgressEvent ev )
        {

        }
    }

    public void re_arrange()
    {
        Rectangle R = getBounds();
        try
        {
            setVisible( false );
            pack();
            setBounds( R );
        }
        finally
        {
            setVisible( true );
        }
    }
    
    public boolean JobExists( String src )
    {
        JobList jl = new JobList( Settings.getInstance().get_jobpath() );
        return jl.contains_job( String.format( "%s.JOB", src ) );
    }

    private boolean new_visible_job( JobData job )
    {
        boolean res = job.is_valid();
        if ( res )
        {
            /*
            CLogger.logfmt( "job( %s ) loaded / state = %s",
                    job.get_short_source(), job.get_state().to_string() );
            */
            
            JobPanel p = new JobPanel( this, job, job_delete_event );
            p.add_state_listener( job_done_event );
            p.set_active( jobs.getComponentCount() == 0 );
            jobs.add( p );
        }
        return res;
    }
    
    public void NewJob()
    {
        BioSpec spec = new BioSpec();
        if ( AccessionWindow.select( spec ) )
        {
            JobData job = new JobData( spec );
                    
            /* this does not save the job, because we want to verify it! */
            if ( JobWindow.edit( job ) )
            {
                if ( JobExists( job.get_short_source() ) )
                {
                    JOptionPane.showMessageDialog( this, "this job already exists",
                        "Error", JOptionPane.INFORMATION_MESSAGE );
                }
                else
                {
                    job.set_state( JobState.READY, false );
                    job.set_valid( true );
                    if ( job.save( Settings.getInstance().get_jobpath() ) )
                    {
                        if ( new_visible_job( job ) )
                        {
                            re_arrange();
                            if ( Settings.getInstance().get_autostart() )
                                start_download_jobs();
                        }
                    }
                }
            }
        }
    }
    
    public int running_downloads()
    {
        int res = 0;
        int i, n_panels = jobs.getComponentCount();
        for ( i = 0; i < n_panels; ++i )
        {
            JobPanel vj = ( JobPanel ) jobs.getComponent( i );
            if ( vj != null && vj.is_running() )
                res ++;
        }
        return res;
    }
    
    public void start_loop( int to_start )
    {
        int i, n_panels = jobs.getComponentCount();
        for ( i = 0; i < n_panels && to_start > 0; ++i )
        {
            JobPanel vj = ( JobPanel ) jobs.getComponent( i );
            if ( vj != null && !vj.is_running() )
            {
                if ( vj.start() ) to_start--;
            }
        }
    }

    public void start_download_jobs()
    {
        start_loop( Settings.getInstance().get_maxdownloads() - running_downloads() );
    }
    
    public void stop_all_threads()
    {
        int i, n_panels = jobs.getComponentCount();
        for ( i = 0; i < n_panels; ++i )
        {
            JobPanel vj = ( JobPanel ) jobs.getComponent( i );
            if ( vj != null ) vj.pause();
        }
    }

    private JToolBar make_toolbar()
    {
        JToolBar res = new JToolBar();
        res.add( ResourceImages.make_img_button( ResourceImages.get_exit_img(), "exit", 2, new Exit_Action( this ) ) );
        res.add( ResourceImages.make_img_button( ResourceImages.get_settings_img(), "settings", 2, new Edit_Global_Settings( this ) ) );
        res.addSeparator();
        res.add( ResourceImages.make_img_button( ResourceImages.get_start_img(), "start downloads", 2, new Start_All_Action( this) ) );
        res.add( ResourceImages.make_img_button( ResourceImages.get_pause_img(), "pause all downloads", 2, new Stop_All_Action( this ) ) );
        res.addSeparator();
        res.add( ResourceImages.make_img_button( ResourceImages.get_add_img(), "add new job", 2, new New_Job_Action( this ) ) );
        res.add( ResourceImages.make_img_button( ResourceImages.get_clear_img(), "remove finished jobs", 2, new CleanUp_Action( this ) ) );
        res.setFloatable( false );
        return res;
    }

    public JobPanel find_job_panel( String name )
    {
        JobPanel res = null;
        int i, n_panels = jobs.getComponentCount();
        for ( i = 0; i < n_panels && res == null; ++i )
        {
            JobPanel vj = ( JobPanel ) jobs.getComponent( i );
            if ( vj != null )
            {
                if ( vj.has_name( name ) )
                    res = vj;
            }
        }
        return res;
    }
    
    private boolean RemoveJob( JobData job )
    {
        String jobname = job.get_short_source();
        boolean res = job.delete_file();
        if ( res )
        {
            JobPanel vj = find_job_panel( jobname );
            if ( vj != null )
            {
                CLogger.logfmt( "removing job >%s<" , jobname );
                jobs.remove( vj );
            }
        }
        return res;
    }
    
    public void RemoveJobAndReArrange( JobData job )
    {
        if ( RemoveJob( job ) )
            re_arrange();
    }
    
    public void CleanUp()
    {
        CLogger.log( "cleaning up..." );
        int removed = 0;
        JobList jl = new JobList( Settings.getInstance().get_jobpath() );
        for ( String filename : jl.make_list() )
        {
            JobData job = new JobData();
            if ( job.load_from( filename ) )
            {
                if ( job.is_complete() )
                {
                    if ( RemoveJob( job ) )
                        removed++;
                }
            }
        }
        if ( removed > 0 ) re_arrange();
    }
    
    public void on_exit()
    {
        Settings.getInstance().set_position( getBounds(), true );
        stop_all_threads();
        CLogger.stop( "End Application" );
        dispose();
        System.exit( 0 );
    }
    
    private void populateContentPane( Container pane )
    {
        toolbar = make_toolbar();
        pane.add( toolbar, BorderLayout.PAGE_START );
        
        jobs = new JPanel();
        jobs.setLayout( new BoxLayout( jobs, BoxLayout.PAGE_AXIS ) );
        
        JScrollPane scroll = new JScrollPane( jobs,
            JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
        
        pane.add( scroll, BorderLayout.CENTER );

        JobList jl = new JobList( Settings.getInstance().get_jobpath() );
        for ( String filename : jl.make_list() )
        {
            new_visible_job( new JobData( filename ) );
        }

        add( StatusBar.getInstance(), BorderLayout.SOUTH );
        
        if ( Settings.getInstance().get_autostart() )
            start_download_jobs();
    }
    
    public MainWindow()
    {
        super( "SRA download" );
        this.setIconImage( ResourceImages.get_logo_img().getImage() );
        
        setJMenuBar( new TheMenuBar( this ) );
        StatusBar.make_status_bar( this.getWidth(), 64 );
        
        job_done_event = new JobDoneEvent();
        job_delete_event = new JobDeleteEvent( this );
        
        populateContentPane( getContentPane() );
        pack();
        
        Dimension MinAndPrefDim = new Dimension( jobs.getWidth() + 50, 400 );
        setPreferredSize( MinAndPrefDim );
        setMinimumSize( MinAndPrefDim );
        pack();
        
        setLocationRelativeTo( null );
        setBounds( Settings.getInstance().get_position() );

        /* to save the position of this frame on exit */
        setDefaultCloseOperation( JFrame.DO_NOTHING_ON_CLOSE );
        addWindowListener( new WindowAdapter()
            {
                @Override public void windowClosing( WindowEvent e )
                {
                    on_exit();
                }
            }
        );
        
        key_dispatcher = new KeyDispatcher( this );
        enable_key_processor( true );

        setVisible( true );
    }
}
