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

import Bio.BioRefSpec;
import Bio.BioReferenceEnumerator;
import data.Settings;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import javax.swing.*;
import job.JobData;
import job.JobList;

public class ReferenceWindow extends JDialog
{
    static final long serialVersionUID = 1;

    private static ReferenceWindow INSTANCE = null;
    public static ReferenceWindow getInstance() { return INSTANCE; }
    
    private final KeyEventReceiver KeyEventReceiverImpl = new KeyEventReceiver();
    private final RefSpecReceiverImpl receiver = new RefSpecReceiverImpl();
    private final PropertyChaneReceiver PropertyChaneReceiverImpl = new PropertyChaneReceiver();
    private final DoneEventReceiver DoneEventImpl = new DoneEventReceiver();
    private final CreateEventReceiver CreateEventImpl = new CreateEventReceiver();
    
    private BioReferenceEnumerator ref_enum = null;
    private final MainWindow parent;
    private final JPanel ref_panel;
    private final JPanel bottom_panel;
    private final JProgressBar progress;
    
    private class RefSpecReceiverImpl implements BioRefSpecReceiver
    {
        @Override public void on_received_spec( final BioRefSpec spec )
        {
            ref_panel.add( new ReferencePanel( spec ) );
            pack();
        }
    }

    private class PropertyChaneReceiver implements PropertyChangeListener
    {
        @Override public void propertyChange( final PropertyChangeEvent event )
        {
            final String propname = event.getPropertyName();
            if ( propname.equals( "progress" ) )
            {
                progress.setValue( ( Integer )event.getNewValue() );
            }
            else if ( propname.equals( "state" ) )
            {
                if ( (SwingWorker.StateValue)event.getNewValue() == SwingWorker.StateValue.DONE )
                {
                    indicate_working( false );
                }
            }
        }
    }

    private class KeyEventReceiver implements MyKeyEventReceiver
    {
        @Override public boolean on_key( final KeyEvent e )
        {
            boolean res = false;
            switch( e.getKeyCode() )
            {
                case KeyEvent.VK_ESCAPE  : setVisible( false );
                                           res = true;
                                           break;
            }
            return res;
        }
    }
    
    private class DoneEventReceiver implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { done(); }
    }

    private class CreateEventReceiver implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae ) { create_jobs(); }
    }

    public static void make_instance( final MainWindow parent )
    {
        if ( INSTANCE == null )
            INSTANCE = new ReferenceWindow( parent );
    }

    public static void show_refs( final JobData job )
    {
        if ( INSTANCE != null && job.is_csra() )
            INSTANCE.show_refs_inst( job );        
    }

    private void indicate_working( Boolean working )
    {
        progress.setIndeterminate( working );
        progress.setVisible( working );
        bottom_panel.setVisible( !working );
    }
    
    private void show_refs_inst( final JobData job )
    {
        setLocationRelativeTo( parent );
        setTitle( String.format( "references of %s", job.get_full_source() ) );
        ref_panel.removeAll();
        indicate_working( true );
        
        MyKeyEventListener.set_receiver( KeyEventReceiverImpl );
        
        ref_enum = new BioReferenceEnumerator( receiver, job.get_full_source() );
        ref_enum.addPropertyChangeListener( PropertyChaneReceiverImpl );
        ref_enum.execute();
        
        setVisible( true );
    }

    private void done()
    {
        setVisible( false );
    }
    
    private void create_jobs()
    {
        JobList jl = new JobList( Settings.getInstance().get_jobpath() );
        
        indicate_working( true );
        int jobs_made = 0;
        int i, n_refs = ref_panel.getComponentCount();
        for ( i = 0; i < n_refs; ++i )
        {
            ReferencePanel ref = ( ReferencePanel ) ref_panel.getComponent( i );
            if ( ref != null )
            {
                JobData job = ref.make_job( jl );
                if ( job != null )
                {
                    if ( parent.new_visible_job( job ) )
                        jobs_made++;
                }
            }
        }
        if ( jobs_made > 0 )
            parent.re_arrange_and_autostart( false );
        
        indicate_working( false );
        done();
    }
    
    private JPanel make_ref_panel()
    {
        JPanel res = new JPanel();
        res.setLayout( new BoxLayout( res, BoxLayout.PAGE_AXIS ) );
        return res;
    }
    
    private JScrollPane make_scroll_pane( final JPanel c )
    {
        JScrollPane res = new JScrollPane( c,
            JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
            JScrollPane.HORIZONTAL_SCROLLBAR_NEVER );
        return res;
    }

    private JProgressBar make_progressbar()
    {
        JProgressBar res = new JProgressBar();
        res.setIndeterminate( false );
        res.setVisible( false );
        return res;
    }
    
    final JButton make_btn( final String caption, final ActionListener listener )
    {
        JButton b;
        b = new JButton( caption );
        b.setPreferredSize( new Dimension( 150, 35 ) );
        if ( listener != null )  b.addActionListener( listener );
        return b;
    }
    
    private JPanel make_bottom_panel()
    {
        JPanel res = new JPanel();
        Dimension dim = new Dimension( 500, 40 );
        res.setPreferredSize( dim );
        res.setMinimumSize( dim );
        res.setMaximumSize( new Dimension( Short.MAX_VALUE, 40 ) );

        res.setLayout( new BoxLayout( res, BoxLayout.LINE_AXIS ) );
        res.add( Box.createHorizontalGlue() );
        res.add( make_btn( "Done", DoneEventImpl ) );
        res.add( Box.createHorizontalGlue() );
        res.add( make_btn( "Create download-jobs", CreateEventImpl ) );
        res.add( Box.createHorizontalGlue() );
        
        res.setVisible( false );
        return res;
    }
    
    public ReferenceWindow( final MainWindow parent )
    {
        super( parent, "references", JDialog.DEFAULT_MODALITY_TYPE );        
        this.parent = parent;
    
        Dimension dim = new Dimension( 500, 300 );
        setPreferredSize( dim );
        setMinimumSize( dim );
        
        Container pane = getContentPane();

        ref_panel = make_ref_panel();
        pane.add( make_scroll_pane( ref_panel ), BorderLayout.CENTER );
        
        progress = make_progressbar();
        pane.add( progress, BorderLayout.NORTH );

        bottom_panel = make_bottom_panel();
        pane.add( bottom_panel, BorderLayout.SOUTH );
    }
}
