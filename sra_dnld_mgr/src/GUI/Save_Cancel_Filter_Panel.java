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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JPanel;


public class Save_Cancel_Filter_Panel extends JPanel
    implements ActionListener
{
    static final long serialVersionUID = 1;
    private final JButton b_save;
    private final JButton b_cancel;
    private final JButton b_filter;
    private final JButton b_preview;
    
    private final List< SaveCancelFilterEventHandler > btn_handlers;
   
    private JButton make_txt_btn( String caption,
                                ImageIcon icon )
    {
        JButton b;
        if ( icon != null )
            b = new JButton( caption, icon );
        else
            b = new JButton( caption );
        
        b.setPreferredSize( new Dimension( 110, 20 ) );
        b.addActionListener( this );
        return b;
    }

    @Override public void actionPerformed( ActionEvent ae )
    {
        if ( ae != null )
        {
            Object src = ae.getSource();
            if ( src != null )
            {
                if ( src instanceof JButton )
                {
                    JButton btn = ( JButton )src;
                    
                    SaveCancelFilterEventType ev_type = SaveCancelFilterEventType.INVALID;
                    if ( btn == b_save ) ev_type = SaveCancelFilterEventType.SAVE;
                    else if ( btn == b_cancel ) ev_type = SaveCancelFilterEventType.CANCEL;
                    else if ( btn == b_filter ) ev_type = SaveCancelFilterEventType.FILTER;
                    else if ( btn == b_preview ) ev_type =SaveCancelFilterEventType.PREVIEW;
                    
                    if ( ev_type != SaveCancelFilterEventType.INVALID )
                    {
                        for ( SaveCancelFilterEventHandler h : btn_handlers )
                            h.on_save_cancel_filter_event( ev_type );
                    }
                }
            }
        }
    }
    
    public void set_filter_btn_enabled( final boolean enabled )
    {
        if ( b_filter != null )
            b_filter.setEnabled( enabled );
    }
    
    public void set_preview_btn_enabled( final boolean enabled )
    {
        if ( b_preview != null )
            b_preview.setEnabled( enabled );
    }
    
    public void add_btn_handler( final SaveCancelFilterEventHandler btn_handler )
    {
        btn_handlers.add( btn_handler );
    }
    
    public Save_Cancel_Filter_Panel()
    {
        super();

        this.btn_handlers = new ArrayList<>();
        
        setMaximumSize( new Dimension( Short.MAX_VALUE, 50 ) );
        setLayout( new BoxLayout( this, BoxLayout.LINE_AXIS ) );
        setBorder( BorderFactory.createMatteBorder( 1, 1, 1, 1, Color.WHITE ) );
        
        add( Box.createHorizontalGlue() );
        
        b_save = make_txt_btn( "Save", ResourceImages.get_check_img() );
        add( b_save );
        
        add( Box.createRigidArea( new Dimension( 10, 10 ) ) );
        
        b_cancel = make_txt_btn( "Cancel", ResourceImages.get_del_img() );
        add( b_cancel );
        
        add( Box.createRigidArea( new Dimension( 10, 10 ) ) );
        
        b_filter = make_txt_btn( "Filter", ResourceImages.get_filter_img() );
        add( b_filter );

        add( Box.createRigidArea( new Dimension( 10, 10 ) ) );
        
        b_preview = make_txt_btn( "Preview", ResourceImages.get_preview_img() );
        add( b_preview );
        
        add( Box.createHorizontalGlue() );
    }

}
