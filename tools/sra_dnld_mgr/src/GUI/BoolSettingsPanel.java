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

import java.awt.BorderLayout;
import java.awt.Dimension;
import javax.swing.JCheckBox;

public class BoolSettingsPanel extends DlgPanel
{
    static final long serialVersionUID = 1;
    
    private final JCheckBox checkbox;
    private boolean value;
    
    public boolean get_value()
    {
        if ( checkbox != null )
            return checkbox.isSelected();
        else
            return false;
    }

    public void set_value( boolean value )
    {
        this.value = value;
        if ( checkbox != null )
            checkbox.setSelected( value );
    }

    public boolean has_changed()
    {
        return ( value != get_value() );
    }
    
    public void set_enabled( boolean value )
    {
        if ( checkbox != null )
            checkbox.setEnabled( value );
    }

    public BoolSettingsPanel( final String caption, final boolean init_value )    
    {
        super( caption, DFLT_PANEL_WIDTH, 0 );

        value = init_value;
        checkbox = new JCheckBox( "enabled" );
        checkbox.setPreferredSize( new Dimension( 100, 5 ) );
        checkbox.setSelected( init_value );
        add( checkbox, BorderLayout.CENTER );
    }
}
