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
import javax.swing.text.*;

public final class IntInputPanel extends DlgPanel implements ActionListener
{
    static final long serialVersionUID = 1;
    private final JTextField tf;
    private final JCheckBox checkbox;
    private final JLabel unit_label;
 
    @Override public void actionPerformed( ActionEvent event )
    {
        tf.setEditable( checkbox.isSelected() );
        tf.setEnabled( checkbox.isSelected() );
    }

    @Override public JCheckBox make_checkbox( boolean enabled )
    {
        JCheckBox res = super.make_checkbox( enabled );
        res.addActionListener( this );
        return res;
    }
    
    @Override public JTextField make_input( boolean editable )
    {
        JTextField res = super.make_input( editable );
        
        PlainDocument doc = (PlainDocument) res.getDocument();
        doc.setDocumentFilter( new MyIntFilter() );
        
        return res;
    }

    public int get_value() { return Integer.parseInt( tf.getText() ); }
    public void set_value( int value ) { tf.setText( Integer.toString( value ) ); }
    
    public boolean is_editable()
    {
        if ( checkbox != null )
            return checkbox.isSelected();
        return false;
    }

    public void set_editable( boolean value )
    {
        if ( checkbox != null )
            checkbox.setSelected( value );
        tf.setEditable( value );
        tf.setEnabled( value );
    }
    
    public void set_enabled( boolean value )
    {
        if ( checkbox != null )
            checkbox.setEnabled( value );
        tf.setEnabled( value );
    }

    public IntInputPanel( final String caption, final String unit,
                          boolean editable, boolean show_enabled )
    {
        super( caption, DFLT_PANEL_WIDTH );

        JPanel p = new JPanel( new BorderLayout() );
        
        if ( show_enabled )
        {
            checkbox = make_checkbox( editable );
            p.add( checkbox, BorderLayout.LINE_START );
        }
        else
        {
            checkbox = null;
        }
        
        tf = make_input( editable );
        p.add( tf, BorderLayout.CENTER );
        
        add( p, BorderLayout.CENTER );
        
        unit_label = make_label( unit, 50 );
        add( unit_label, BorderLayout.LINE_END );
    }
}


class MyIntFilter extends DocumentFilter
{
    @Override public void insertString( FilterBypass fb,
           int offset, String string,
           AttributeSet attr ) throws BadLocationException
    {

        Document doc = fb.getDocument();
        StringBuilder sb = new StringBuilder();
        sb.append( doc.getText( 0, doc.getLength() ) );
        sb.insert( offset, string );

        if ( test( sb.toString() ) )
        {
            super.insertString(fb, offset, string, attr);
        }
        else
        {
            // warn the user and don't allow the insert
        }
    }

    private boolean test( String text )
    {
        boolean res = true;
        try { Integer.parseInt( text ); }
        catch ( NumberFormatException e ) { res = false; }
        return res;
    }

    @Override public void replace( FilterBypass fb,
           int offset, int length, String text,
           AttributeSet attrs ) throws BadLocationException
    {
        Document doc = fb.getDocument();
        StringBuilder sb = new StringBuilder();
        sb.append( doc.getText( 0, doc.getLength() ) );
        sb.replace( offset, offset + length, text );

        if ( test( sb.toString() ) )
        {
            super.replace( fb, offset, length, text, attrs );
        }
        else
        {
         // warn the user and don't allow the insert
        }
    }

    @Override public void remove( FilterBypass fb,
           int offset, int length ) throws BadLocationException
    {
        Document doc = fb.getDocument();
        StringBuilder sb = new StringBuilder();
        sb.append( doc.getText( 0, doc.getLength() ) );
        sb.delete( offset, offset + length );

        if ( test( sb.toString() ) )
        {
            super.remove( fb, offset, length );
        }
        else
        {
         // warn the user and don't allow the insert
        }
    }
}