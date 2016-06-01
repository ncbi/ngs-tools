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

import Bio.BioReadType;
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JPanel;


public class JobReadTypePanel extends DlgPanel implements ActionListener
{
    static final long serialVersionUID = 1;
    
    private final JCheckBox checkbox;
    private final JComboBox<String> types;

    public BioReadType get_value()
    {
        return BioReadType.from_ordinal( types.getSelectedIndex() + 1 );
    }

    public void set_value( final BioReadType format )
    {
        types.setSelectedIndex( format.to_ordinal() - 1 );
    }
    
    public void set_editable( boolean value )
    {
        checkbox.setSelected( value );
        types.setEnabled( value );
    }

    public boolean get_editable() { return checkbox.isSelected(); }

    @Override public void actionPerformed( ActionEvent e )
    {
        types.setEnabled( checkbox.isSelected() );
    }

    public JobReadTypePanel( String caption )
    {
        super( caption, DFLT_PANEL_WIDTH, 0 );
        JPanel p = new JPanel( new BorderLayout() );
        
        checkbox = make_checkbox( false );
        checkbox.addActionListener( this );
        p.add( checkbox, BorderLayout.LINE_START );
        
        types = make_combo_box( BioReadType.to_StringList( BioReadType.INVALID ) );
        p.add( types, BorderLayout.CENTER );

        add( p, BorderLayout.CENTER );
        
    }
    
}
