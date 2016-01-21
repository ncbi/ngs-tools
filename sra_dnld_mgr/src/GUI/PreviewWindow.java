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

import Bio.BioDumperPreview;
import data.Settings;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.KeyEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import javax.swing.BoxLayout;
import javax.swing.JDialog;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;
import job.JobData;

public class PreviewWindow extends JDialog
{
    static final long serialVersionUID = 1;
    
    private static PreviewWindow INSTANCE = null;
    public static PreviewWindow getInstance() { return INSTANCE; }
    
    public static void make_instance( final MainWindow parent )
    {
        if ( INSTANCE == null )
            INSTANCE = new PreviewWindow( parent );
    }
    
    public static void show_preview( final JobData job )
    {
        if ( INSTANCE != null )
            INSTANCE.preview( job );
    }
    
    private final KeyEventReceiver KeyEventReceiverImpl = new KeyEventReceiver();
    private final PropertyChaneReceiver PropertyChaneReceiverImpl = new PropertyChaneReceiver();
    
    private final MainWindow parent;
    private final JTextArea txt = new JTextArea();
    private final JScrollPane scroll = new JScrollPane( txt );  
    private final JProgressBar progress = new JProgressBar();
    
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
                switch ( (SwingWorker.StateValue) event.getNewValue() )
                {
                    case DONE    : progress.setVisible( false ); break;

                    case STARTED : 
                    case PENDING : progress.setVisible( true ); break;
                }
            }
        }
    }
    
    private void preview( final JobData job )
    {
        setLocationRelativeTo( parent );
        setTitle( String.format( "preview %s", job.get_short_source() ) );
        txt.setText( "" );
        
        BioDumperPreview preview = new BioDumperPreview( job, txt,
                Settings.getInstance().get_preview_rows() );
        preview.addPropertyChangeListener( PropertyChaneReceiverImpl );
        preview.execute();
        
        MyKeyEventListener.set_receiver( KeyEventReceiverImpl );
        setVisible( true );
    }

    public PreviewWindow( MainWindow parent )
    {
        super( parent, "preview", JDialog.DEFAULT_MODALITY_TYPE );
        this.parent = parent;

        Dimension dim = new Dimension( 500, 300 );
        setPreferredSize( dim );
        setMinimumSize( dim );

        Container pane = getContentPane();
        pane.setLayout( new BoxLayout( pane, BoxLayout.PAGE_AXIS ) );
        
        txt.setEditable( false );
        pane.add( scroll );

        progress.setIndeterminate( true );
        pane.add( progress, BoxLayout.Y_AXIS );
    }
}
