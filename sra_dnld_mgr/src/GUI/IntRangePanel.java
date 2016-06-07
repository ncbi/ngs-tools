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

public final class IntRangePanel extends DlgPanel implements ActionListener
{
    static final long serialVersionUID = 1;
    private final JCheckBox checkbox;
    private final JTextField tf_start;
    private final JTextField tf_count;
 
    @Override public void actionPerformed( ActionEvent event )
    {
        tf_start.setEditable( checkbox.isSelected() );
        tf_start.setEnabled( checkbox.isSelected() );
        tf_count.setEditable( checkbox.isSelected() );
        tf_count.setEnabled( checkbox.isSelected() );
    }

    @Override public JCheckBox make_checkbox( boolean enabled )
    {
        JCheckBox res = super.make_checkbox( enabled );
        res.addActionListener( this );
        return res;
    }
    
    @Override public JTextField make_input( boolean editable, int w )
    {
        JTextField res = super.make_input( editable, w );
        
        PlainDocument doc = (PlainDocument) res.getDocument();
        doc.setDocumentFilter( new MyRangeFilter() );
        
        return res;
    }

    public long get_start_value() { return Long.parseLong( tf_start.getText() ); }
    public void set_start_value( long value ) { tf_start.setText( Long.toString( value ) ); }
    
    public long get_count_value() { return Long.parseLong(tf_count.getText() ); }
    public void set_count_value( long value ) { tf_count.setText( Long.toString( value ) ); }

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
        tf_start.setEditable( value );
        tf_start.setEnabled( value );
        tf_count.setEditable( value );
        tf_count.setEnabled( value );        
    }
    
    public void set_enabled( boolean value )
    {
        if ( checkbox != null )
            checkbox.setEnabled( value );
        tf_start.setEnabled( value );
        tf_count.setEnabled( value );
    }

    private JPanel make_labeled( JTextField tf, String label, int w )
    {
        JPanel p = new JPanel( new BorderLayout() );
        p.add( make_label( label, w ), BorderLayout.LINE_START );
        p.add( tf, BorderLayout.LINE_END );
        return p;
    }
    
    private JPanel make_range( JTextField tfs, JTextField tfc, int w )
    {
        JPanel p = new JPanel( new BorderLayout() );
        p.add( make_labeled( tfs, "start=", w ) , BorderLayout.LINE_START );
        p.add( make_labeled( tfc, "count=", w ), BorderLayout.LINE_END );
        return p;
    }

    public IntRangePanel( final String caption, final String unit,
                          boolean editable, boolean show_enabled )
    {
        super( caption, DFLT_PANEL_WIDTH, 0 );

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
        
        tf_start = make_input( editable, 80 );
        tf_count = make_input( editable, 80 );
        
        p.add( make_range( tf_start, tf_count, 45 ), BorderLayout.CENTER );
        add( p, BorderLayout.CENTER );
        
        add( make_label( unit, 50 ), BorderLayout.LINE_END );
    }
}


class MyRangeFilter extends DocumentFilter
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
        try { Long.parseLong( text ); }
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