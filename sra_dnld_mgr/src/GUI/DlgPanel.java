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
import javax.swing.*;

class DlgPanel extends JPanel
{
    static final long serialVersionUID = 1;
    static final int DFLT_PANEL_WIDTH = 75;
            
    private JLabel L;
    
    public final JLabel make_label( String caption, int width )
    {
        JLabel res;
        res = new JLabel( caption );
        res.setPreferredSize( new Dimension( width, 1 ) );
        res.setHorizontalAlignment( SwingConstants.CENTER );
        res.setBorder( BorderFactory.createMatteBorder( 1, 1, 1, 1, Color.BLACK ) );
        return res;
    }
    
    public JCheckBox make_checkbox( boolean enabled )
    {
        JCheckBox res = new JCheckBox( "enabled" );
        res.setPreferredSize( new Dimension( 100, 5 ) );
        res.setSelected( enabled );
        return res;
    }

    public JComboBox<String> make_combo_box( boolean enabled )
    {
        JComboBox<String> res = new JComboBox<>();
        res.setPreferredSize( new Dimension( 100, 5 ) );
        ( ( JLabel )res.getRenderer() ).setHorizontalAlignment( SwingConstants.CENTER );
        res.setEnabled( enabled );
        res.setEditable( false );
        return res;
    }

    public JComboBox<String> make_combo_box( final String[] list )
    {
        JComboBox<String> res = new JComboBox<>( list );
        res.setPreferredSize( new Dimension( 100, 5 ) );
        ( ( JLabel )res.getRenderer() ).setHorizontalAlignment( SwingConstants.CENTER );
        res.setEnabled( true );
        res.setEditable( false );
        return res;
    }

    public JTextField make_input( boolean editable )
    {
        JTextField res = new JTextField();
        res.setPreferredSize( new Dimension( 100, 5 ) );
        res.setEditable( editable );
        res.setEnabled( editable );
        res.setHorizontalAlignment( JTextField.CENTER );
        return res;
    }

    public int get_label_text_width()
    {
        int res = 0;
        if ( L != null )
        {
            FontMetrics fm = L.getFontMetrics( L.getFont() );
            res = fm.stringWidth( L.getText() );
        }
        return res;
    }
    
    public void set_label_width( int w )
    {
        if ( L != null )
        {
            Dimension dim = L.getSize();
            dim.width = w;
            L.setMinimumSize( dim );
            L.setPreferredSize( dim );
            L.setSize( dim );
        }
    }
    
    public DlgPanel( final String caption, final int width )
    {
        super();
        L = null;
        
        setLayout( new BorderLayout( 5, 0 ) );
        setBorder( BorderFactory.createMatteBorder( 1, 1, 1, 1, Color.WHITE ) );

        setPreferredSize( new Dimension( Short.MAX_VALUE, 25 ) );
        setMinimumSize( getPreferredSize() );
        setMaximumSize( getPreferredSize() );

        if ( width > 0 )
        {
            L = make_label( caption, width );
            add( L , BorderLayout.LINE_START );
        }
    }
}
