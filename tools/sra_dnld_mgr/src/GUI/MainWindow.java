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
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import javax.swing.*;

public class MainWindow extends JFrame
{
    static final long serialVersionUID = 1;

    private final on_settings settings_event = new on_settings();
    private final on_exit exit_event = new on_exit();
    private final on_stop_all stop_all_event = new on_stop_all();
    private final on_start_all start_all_event = new on_start_all();
    private final on_new_job now_job_event = new on_new_job();
    private final on_cleanup cleanup_event = new on_cleanup();
    private final on_key_dispatch key_dispatch_event = new on_key_dispatch();
    private final on_progress progress_event = new on_progress();
    private final on_window_close window_close_event = new on_window_close();
    private final on_job_delete job_delete_event = new on_job_delete();
            
    private final JToolBar toolbar = make_toolbar();
    private final JPanel jobs = new JPanel();
    
    private String to_test = null;
    private String scratch_space = null;

    private class on_settings implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { SettingsWindow.edit(); }
    }

    private class on_exit implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { on_exit( 0 ); }
    }

    private class on_stop_all implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { stop_all_threads(); }
    }

    private class on_start_all implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { start_download_jobs(); }
    }

    private class on_new_job implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { NewJob(); }
    }

    private class on_cleanup implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { CleanUp(); }
    }

    private class on_key_dispatch implements KeyEventDispatcher
    {
        @Override public boolean dispatchKeyEvent( KeyEvent ke )
        {
            boolean res = false;
            if ( ke.getID() == KeyEvent.KEY_PRESSED ) res = on_key( ke );
            return res;
        }
    }

    private class on_progress implements ProgressListenerInterface
    {
        @Override public void on_state_progress_event( StateAndProgressEvent ev )
        {
            if ( ev.type == StateAndProgressType.STATE &&
                 ev.prev_state == JobState.RUNNING &&
                 ev.new_state == JobState.DONE )
                job_has_finished();
        }
    }

    private class on_window_close extends WindowAdapter
    {
        @Override public void windowClosing( WindowEvent we ) { on_exit( 0 ); }
    }

    private class on_job_delete implements JobDeleteEvent
    {
        @Override public void delete_job( JobData job ) { ask_and_delete_job( job ); }
    }
    
    /* -----------------------------------------------------------------------*/

    public final void enable_key_processor( boolean enable )
    {
        KeyboardFocusManager kfm = KeyboardFocusManager.getCurrentKeyboardFocusManager();
        if ( enable )
            kfm.addKeyEventDispatcher( key_dispatch_event );
        else
            kfm.removeKeyEventDispatcher( key_dispatch_event );
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
        if ( p != null )
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
        return ( p != null );
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
        int active_id = get_active_panel_id();
        int new_panel_id = get_panel_id( p );
        if ( active_id != new_panel_id )
            switch_panel( jobs.getComponentCount(), active_id, new_panel_id );
    }
    
    public void re_arrange( final boolean change_visiblity )
    {
        Rectangle R = getBounds();
        try
        {
            if ( change_visiblity ) setVisible( false );
            pack();
            setBounds( R );
        }
        finally
        {
            if ( change_visiblity ) setVisible( true );
        }
    }
    
    public boolean JobExists( String src )
    {
        JobList jl = new JobList( Settings.getInstance().get_jobpath() );
        return jl.contains_job( String.format( "%s.JOB", src ) );
    }

    public void ask_and_delete_job( JobData job )
    {
        Boolean perform_delete = true;
        
        if ( Settings.getInstance().get_confirm_delete() )
        {
            int response = JOptionPane.showConfirmDialog(
                this,
                String.format( "delete job %s", job.get_short_source() ),
                "question",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.INFORMATION_MESSAGE );

            perform_delete = ( response == JOptionPane.YES_OPTION );
        }
        
        if ( perform_delete )
            RemoveJobAndReArrange( job, true );
    }
    
    public boolean new_visible_job( JobData job )
    {
        boolean res = job.is_valid();
        if ( res )
        {
            JobPanel p = new JobPanel( this, job, job_delete_event );
            p.add_state_listener( progress_event );
            p.set_active( jobs.getComponentCount() == 0 );
            jobs.add( p );
        }
        return res;
    }
    
    public void re_arrange_and_autostart( final Boolean change_visibility )
    {
        re_arrange( change_visibility );
        if ( Settings.getInstance().get_autostart() )
            start_download_jobs();
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
                            re_arrange_and_autostart( true );
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

    public final void start_download_jobs()
    {
        start_loop( Settings.getInstance().get_maxdownloads() - running_downloads() );
    }

    private int test_has_finished()
    {
        /* take the job's output, calculate md5, compare with expected md5 */
        int exit_code = 0;
        JobPanel p = get_active_panel();
        JobData job = p.job;
        String output_filename = job.get_output_filename();
        String created_md5 = md5_extractor.md5_of_file( output_filename);
        String expected_md5 = job.get_expected_md5();

        try
        {
            Files.deleteIfExists( Paths.get( output_filename ) );
        }
        catch ( IOException ex )
        {
            System.out.printf( "delete '%s' : %s\n", output_filename, ex.toString() );
        }

        if ( expected_md5.length() == 0 )
        {
            System.out.printf( "no expected md5-sum found in job-file, writing calculated md5 '%s' into job-file\n", created_md5 );
            job.set_expected_md5( created_md5 );
            job.store();
            exit_code = 3;
        }
        else
        {
            if ( created_md5.equalsIgnoreCase( expected_md5 ) )
            {
                System.out.printf( "md5-sum of '%s' matches expected md5-sum from job-file\n", created_md5 );
            }
            else
            {
                System.out.printf( "calculated md5-sump of '%s' does not match expected md5-sum of '%s'\n", created_md5, expected_md5 );
                exit_code = 3;
            }
        }
        return exit_code;
    }
    
    public final void job_has_finished()
    {
        if ( to_test != null )
            on_exit( test_has_finished() );
        else if ( Settings.getInstance().get_autostart() )
            start_download_jobs();
    }
    
    public final void stop_all_threads()
    {
        int i, n_panels = jobs.getComponentCount();
        for ( i = 0; i < n_panels; ++i )
        {
            JobPanel vj = ( JobPanel ) jobs.getComponent( i );
            if ( vj != null )
            {
                if ( vj.is_running() ) vj.pause();
            }
        }
    }

    private JToolBar make_toolbar()
    {
        JToolBar res = new JToolBar();
        res.add( ResourceImages.make_img_button( ResourceImages.get_exit_img(), "exit", 2, exit_event ) );
        res.add( ResourceImages.make_img_button( ResourceImages.get_settings_img(), "settings", 2, settings_event ) );
        res.addSeparator();
        res.add( ResourceImages.make_img_button( ResourceImages.get_start_img(), "start downloads", 2, start_all_event ) );
        res.add( ResourceImages.make_img_button( ResourceImages.get_pause_img(), "pause all downloads", 2, stop_all_event ) );
        res.addSeparator();
        res.add( ResourceImages.make_img_button( ResourceImages.get_add_img(), "add new job", 2, now_job_event ) );
        res.add( ResourceImages.make_img_button( ResourceImages.get_clear_img(), "remove finished jobs", 2, cleanup_event ) );
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
    
    public void RemoveJobAndReArrange( final JobData job, final Boolean change_visibility )
    {
        if ( RemoveJob( job ) )
            re_arrange( change_visibility );
    }
    
    public void CleanUp()
    {
        CLogger.log( "cleaning up..." );
        int removed = 0;
        JobList jl = new JobList( Settings.getInstance().get_jobpath() );
        for ( String filename : jl.make_list() )
        {
            JobData job = new JobData( filename );
            if ( job.is_valid() && job.is_complete() )
            {
                if ( RemoveJob( job ) ) removed++;
            }
        }
        if ( removed > 0 ) re_arrange( true );
    }
    
    public void on_exit( int exit_code )
    {
        Settings.getInstance().set_position( getBounds(), true );
        stop_all_threads();
        CLogger.stop( "End Application" );
        dispose();
        System.exit( exit_code );
    }
    
    private void populateContentPane( Container pane )
    {
        pane.add( toolbar, BorderLayout.PAGE_START );
        
        jobs.setLayout( new BoxLayout( jobs, BoxLayout.PAGE_AXIS ) );
        
        JScrollPane scroll = new JScrollPane( jobs,
            JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
        
        pane.add( scroll, BorderLayout.CENTER );

        if ( to_test == null )
        {
            /* no job to execute from cmd-line: read jobs from path */
            JobList jl = new JobList( Settings.getInstance().get_jobpath() );
            for ( String filename : jl.make_list() )
                new_visible_job( new JobData( filename ) );
        }
        else
        {
            /* we have a job to execute from cmd-line: */
            JobData job = new JobData( to_test );
            /* reset the job in any case... */
            job.reset_as_test( scratch_space );
            new_visible_job( job );
        }

        add( StatusBar.getInstance(), BorderLayout.SOUTH );
        
        if ( to_test == null && Settings.getInstance().get_autostart() )
            start_download_jobs();
    }
    
    public MainWindow( final String to_test, final String scratch_space )
    {
        super( "SRA download" );
        
        if ( to_test != null )
        {
            this.to_test = to_test;
            this.setTitle( "SRA download testing: " + to_test );
        }
        
        if ( scratch_space != null )
            this.scratch_space = scratch_space;
        
        this.setIconImage( ResourceImages.get_logo_img().getImage() );
        
        setJMenuBar( new TheMenuBar( this ) );
        StatusBar.make_status_bar( this.getWidth(), 64 );
        populateContentPane( getContentPane() );
       
        Dimension MinAndPrefDim = new Dimension( jobs.getWidth() + 50, 400 );
        setPreferredSize( MinAndPrefDim );
        setMinimumSize( MinAndPrefDim );
        pack();
        
        setLocationRelativeTo( null );
        setBounds( Settings.getInstance().get_position() );

        setDefaultCloseOperation( JFrame.DO_NOTHING_ON_CLOSE );
        addWindowListener( window_close_event );
        enable_key_processor( true );
        
        if ( to_test == null )
            setVisible( true );
    }
}
