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
import java.io.File;
import javax.swing.*;

public class PathChooserPanel extends DlgPanel implements ActionListener
{
    private JTextField tf;

    @Override public void actionPerformed( ActionEvent e )
    {
        JFileChooser fc = new JFileChooser( tf.getText() );
        fc.setFileSelectionMode( JFileChooser.DIRECTORIES_ONLY );
        int res = fc.showOpenDialog( tf );
        if ( res == JFileChooser.APPROVE_OPTION )
        {
            File f = fc.getSelectedFile();
            try { tf.setText( f.getCanonicalPath() ); }
            catch ( Exception ex ) { }
        }
    }

    public String get_text() { return tf.getText(); }
    public void set_text( String value ) { tf.setText( value ); }
    
    public PathChooserPanel( final String caption )
    {
        super( caption, DFLT_PANEL_WIDTH, 0 );
        
        tf = make_input( false, 100 );
        add( tf, BorderLayout.CENTER );
        
        JButton b = make_btn( "...", 50, 5, this, null );
        add( b, BorderLayout.LINE_END );
    }
}
