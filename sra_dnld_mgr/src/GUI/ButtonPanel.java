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
import java.awt.event.*;
import javax.swing.*;

public class ButtonPanel extends DlgPanel
{
    private final JTextField tf = make_input( false );
    private final JButton b;
    private final on_button button_event = new on_button();
    private ActionListener relay;
    
    private class on_button implements ActionListener
    {
        @Override public void actionPerformed( ActionEvent ae )
        {
            if ( relay != null ) relay.actionPerformed( ae );
        }
    }

    public String get_text() { return tf.getText(); }
    public void set_text( String value ) { tf.setText( value ); }
    public void set_action_listener( final ActionListener al ) { relay = al; }

    public void set_enabled( final boolean enabled ) { b.setEnabled( enabled ); }
    
    public ButtonPanel( final String caption, final String hint,
                        final ActionListener al, ImageIcon icon )
    {
        super( caption, DFLT_PANEL_WIDTH, icon == null ? 0 : icon.getIconHeight() + 4 );
        relay = al;
        
        tf.setText( hint );
        add( tf, BorderLayout.CENTER );

        if ( icon == null )
            b = make_btn( "...", 75, 32, button_event, null );
        else
            b = make_btn( "", 75, 32, button_event, icon );
        add( b, BorderLayout.LINE_END );
    }

}
