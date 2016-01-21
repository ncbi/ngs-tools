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

import java.awt.*;
import java.awt.event.KeyEvent;
import javax.swing.JDialog;

public class DlgWithMaxSize extends JDialog
    implements MyKeyEventReceiver, SaveCancelEventHandler
{
    static final long serialVersionUID = 1;
    
    public final MainWindow parent;
    public boolean show_result;
    
    @Override public boolean on_key( final KeyEvent e )
    {
        boolean res = false;
        switch( e.getKeyCode() )
        {
            case KeyEvent.VK_ESCAPE  : cancel_dlg(); res = true; break;
        }
        return res;
    }

    @Override public void on_save_cancel_filter_event( final SaveCancelEventType event_type )
    {
        switch( event_type )
        {
            case SAVE   : save_and_close_dlg(); break;
            case CANCEL : cancel_dlg(); break;
        }
    }

    /*  this is unfortunately neccessary,
        because setMaximumSize( ... ) is ignored by JFrame */
    @Override public void paint( Graphics g )
    {
        Dimension d = getSize();
        Dimension m = getMaximumSize();
        if ( d.width > m.width || d.height > m.height )
        {
            d.width = Math.min( m.width, d.width );
            d.height = Math.min( m.height, d.height );
            
            Point p = getLocation();
            setVisible( false );
            setSize( d );
            setLocation( p );
            setVisible( true );
        }
        else
        {
            super.paint( g );
        }
    }
    
    public Save_Cancel_Panel add_save_cancel_panel( Container pane )
    {
        Save_Cancel_Panel res = new Save_Cancel_Panel();
        res.add_btn_handler( this );
        pane.add( res );
        return res;
    }

    public void cancel_dlg()
    {
        show_result = false;
        setVisible( false );
    }

    public void save_and_close_dlg()
    {
        show_result = true;
        setVisible( false );
    }
    
    public boolean show_dialog()
    {
        setLocationRelativeTo( parent );
        show_result = false;
        parent.enable_key_processor( false );
        MyKeyEventListener.set_receiver( this );
        setVisible( true );
        MyKeyEventListener.set_receiver( null );
        parent.enable_key_processor( true );
        return show_result;
    }

    private void setHeight( int value )
    {
        Dimension dim = this.getMinimumSize();
        dim.height = value;
        this.setMinimumSize( dim );

        dim = this.getMaximumSize();
        dim.height = value;
        this.setMaximumSize( dim );

        dim = this.getSize();
        dim.height = value;
        setSize( dim );
    }

    public void adjust_height( final int offset )
    {
        pack();
        int h = offset;
        Container pane = getContentPane();
        for ( Component c : pane.getComponents() ) h += c.getHeight();
        setHeight( h );
    }
    
    public int get_max_label_text_width( Container pane )
    {
        int res = 0;
        int count = pane.getComponentCount();
        for ( int i = 0; i < count; ++i )
        {
            Component c = pane.getComponent( i );
            if ( c != null && ( c instanceof DlgPanel ) )
            {
                int w = ( ( DlgPanel )c ).get_label_text_width();
                if ( w > res ) res = w;
            }
        }
        return res;
    }
    
    public void resize_labels( Container pane )
    {
        int w = get_max_label_text_width( pane ) + 20;
        int count = pane.getComponentCount();
        for ( int i = 0; i < count; ++i )
        {
            Component c = pane.getComponent( i );
            if ( c != null && ( c instanceof DlgPanel ) )
                ( ( DlgPanel ) c ).set_label_width( w );
        }
    }
    
    public DlgWithMaxSize( MainWindow parent,
            final String caption, final Dimension dim )
    {
        super( parent, caption, JDialog.DEFAULT_MODALITY_TYPE );
        this.parent = parent;

        setPreferredSize( dim );
        setMinimumSize( dim );
        setMaximumSize( dim );
        show_result = false;
    }
}
