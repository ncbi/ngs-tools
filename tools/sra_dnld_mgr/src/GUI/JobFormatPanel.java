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

import job.JobFormat;
import java.awt.BorderLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import job.JobSubFormat;

public class JobFormatPanel extends DlgPanel implements ItemListener
{
    static final long serialVersionUID = 1;
    
    private final JComboBox<String> format_box;
    private final JComboBox<String> sub_format_box;
    private final ItemListener listener;
    
    @Override public void itemStateChanged( ItemEvent e )
    {
        if ( listener != null )
            listener.itemStateChanged( e );
    }
    
    public JobFormat get_format()
    {
        return JobFormat.from_ordinal( format_box.getSelectedIndex() + 1 );
    }

    public void set_format( final JobFormat format )
    {
        format_box.setSelectedIndex( format.to_ordinal() - 1 );
    }

    public JobSubFormat get_sub_format()
    {
        return JobSubFormat.from_ordinal( sub_format_box.getSelectedIndex() + 1 );
    }
    
    public void set_sub_format( final JobSubFormat sub_format )
    {
        if ( sub_format.equals( JobSubFormat.INVALID ) )
            sub_format_box.setSelectedIndex( JobSubFormat.SPOT.to_ordinal() - 1 );    
        else
            sub_format_box.setSelectedIndex( sub_format.to_ordinal() - 1 );
    }

    public void set_sub_format_enabled( final boolean enabled )
    {
        sub_format_box.setEnabled( enabled );
    }

    public JobFormatPanel( final String caption, final ItemListener listener )
    {
        super( caption, 75, 0 );
        
        this.listener = listener;
        
        JPanel p = new JPanel( new BorderLayout() );
        
        format_box = make_combo_box( JobFormat.to_StringList( JobFormat.INVALID ) );
        format_box.addItemListener( this );
        p.add( format_box, BorderLayout.LINE_START );
        
        sub_format_box = make_combo_box( JobSubFormat.to_StringList( JobSubFormat.INVALID ) );
        p.add( sub_format_box, BorderLayout.CENTER );
        
        add( p, BorderLayout.CENTER );
    }
}
