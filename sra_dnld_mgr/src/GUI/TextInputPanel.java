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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

public class TextInputPanel extends DlgPanel
    implements DocumentListener, ActionListener
{
    static final long serialVersionUID = 1;
    
    private final JTextField tf;
    private final Timer timer;
    private ActionListener listener;
    private boolean blink_state;
            
    public String get_text() { return tf.getText(); }
    public void set_text( String value ) { tf.setText( value ); }
    public void set_editable( boolean value ) { tf.setEditable( value ); }

    private void has_changed()
    {
        if ( listener != null ) listener.actionPerformed( null );
    }
    
    @Override public void actionPerformed( ActionEvent ae )
    {
        tf.setBackground( blink_state ? Color.lightGray : Color.DARK_GRAY );
        blink_state = !blink_state;
    }
    
    public void blink( final boolean on )
    {
        if ( on )
        {
            blink_state = false;
            tf.setBackground( Color.DARK_GRAY );
            timer.start();
        }
        else
        {
            tf.setBackground( Color.white );
            timer.stop();
        }
    }
    
    @Override public void insertUpdate( DocumentEvent de ) { has_changed(); }
    @Override public void removeUpdate( DocumentEvent de ) { has_changed(); }
    @Override public void changedUpdate( DocumentEvent de ) { has_changed(); }
    public void add_listener( ActionListener listener ) { this.listener = listener; }
    
    public TextInputPanel( final String caption, final boolean editable )
    {
        super( caption, 75, 0 );
        listener = null;
        
        tf = make_input( editable );
        tf.getDocument().addDocumentListener( this );
        add( tf, BorderLayout.CENTER );

        blink_state = false;
        timer = new Timer( 300, this );
        timer.setInitialDelay( 300 );
    }
}
