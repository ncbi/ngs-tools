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
import javax.swing.JLabel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class IntSettingsPanel extends DlgPanel implements ChangeListener
{
    static final long serialVersionUID = 1;
    
    JSlider slider;
    JLabel slider_value;

    @Override public void stateChanged( ChangeEvent e )
    {
        slider_value.setText( Integer.toString( slider.getValue() ) );
    }    
    
    private JSlider make_slider( int min, int max, int value )
    {
        JSlider res = new JSlider( JSlider.HORIZONTAL, min, max, value );
        res.setPreferredSize( new Dimension( 100, 5 ) );
        res.addChangeListener( this );
        return res;
    }

    public int get_value() { return slider.getValue(); }

    public void set_value( int value )
    {
        slider.setValue( value );
        slider_value.setText( Integer.toString( value ) );
    }
    
    public IntSettingsPanel( final String caption, final int min, final int max )
    {
        super( caption, DFLT_PANEL_WIDTH );
        
        slider = make_slider( min, max, min );
        add( slider, BorderLayout.CENTER );
        
        slider_value = make_label( Integer.toString( min ), 50 );
        add( slider_value, BorderLayout.LINE_END );
    }
}
