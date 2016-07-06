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

import Bio.BioRefSpec;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import javax.swing.BorderFactory;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import job.JobData;
import job.JobList;

public class ReferencePanel extends JPanel
{
    static final long serialVersionUID = 1;
    
    private final BioRefSpec spec;
    private JCheckBox checkbox;
    
    public JobData make_job( final JobList jl )
    {
        JobData res = null;
        if ( checkbox.isSelected() )
        {
            res = spec.make_job( jl );
            checkbox.setSelected( res == null );
        }
        return res;
    }
    
    private JCheckBox make_checkbox( Boolean enabled )
    {
        JCheckBox res = new JCheckBox( "download" );
        res.setPreferredSize( new Dimension( 75, 20 ) );
        res.setSelected( enabled );
        res.setEnabled( enabled );
        return res;
    }
    
    private JLabel make_label( String caption, int width )
    {
        JLabel res;
        res = new JLabel( caption );
        res.setPreferredSize( new Dimension( width, 1 ) );
        res.setHorizontalAlignment( SwingConstants.CENTER );
        res.setBorder( BorderFactory.createMatteBorder( 1, 1, 1, 1, Color.BLACK ) );
        return res;
    }
    
    ReferencePanel( final BioRefSpec spec )
    {
        this.spec = spec;
        setLayout( new BorderLayout( 5, 0 ) );
        setBorder( BorderFactory.createMatteBorder( 1, 1, 1, 1, Color.WHITE ) );

        setPreferredSize( new Dimension( 300, 35 ) );
        setMinimumSize( getPreferredSize() );
        setMaximumSize( new Dimension( Short.MAX_VALUE, 35 ) );
        setLayout( new BorderLayout( 5, 0 ) );
        
        add( make_label( spec.name(), 150 ), BorderLayout.LINE_START );
        add( make_label( spec.desc(), 150 ), BorderLayout.CENTER );
        checkbox = make_checkbox( !spec.get_downloaded() );
        add( checkbox, BorderLayout.LINE_END );
    }
}
