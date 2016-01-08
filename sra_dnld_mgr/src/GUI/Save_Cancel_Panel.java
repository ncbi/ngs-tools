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

public class Save_Cancel_Panel extends JPanel implements ActionListener
{
    static final long serialVersionUID = 1;
    private final JButton b_save;
    private final JButton b_cancel;    
    private final List< SaveCancelFilterEventHandler > btn_handlers;
    
    final JButton make_txt_btn( String caption,
                                ImageIcon icon )
    {
        JButton b;
        b = new JButton( caption, icon );
        b.setPreferredSize( new Dimension( 150, 5 ) );
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

                    if ( ev_type != SaveCancelFilterEventType.INVALID )
                    {
                        for ( SaveCancelFilterEventHandler h : btn_handlers )
                            h.on_save_cancel_filter_event( ev_type );
                    }
                }
            }
        }
    }
    
    public void set_save_btn_status( final boolean enabled )
    {
        b_save.setEnabled( enabled );
    }
    
    public void add_btn_handler( final SaveCancelFilterEventHandler btn_handler )
    {
        btn_handlers.add( btn_handler );
    }

    public Save_Cancel_Panel()
    {
        super();
        this.btn_handlers = new ArrayList<>();
        
        setMaximumSize( new Dimension( Short.MAX_VALUE, 50 ) );
        setLayout( new BoxLayout( this, BoxLayout.LINE_AXIS ) );
        setBorder( BorderFactory.createMatteBorder( 1, 1, 1, 1, Color.WHITE ) );
        
        add( Box.createHorizontalGlue() );
        b_save = make_txt_btn( "Save", ResourceImages.get_check_img() );
        add( b_save );
        add( Box.createRigidArea( new Dimension( 50, 10 ) ) );
        b_cancel = make_txt_btn( "Cancel", ResourceImages.get_del_img() );
        add( b_cancel );
        add( Box.createHorizontalGlue() );
    }
}
