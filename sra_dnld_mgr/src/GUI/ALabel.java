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
import javax.swing.border.Border;

public class ALabel extends JLabel
{
    static final long serialVersionUID = 1;
    
    @Override public Dimension getMaximumSize()
    {
        Dimension size = getPreferredSize();
        size.height = Short.MAX_VALUE;
        return size;
    }

    @Override public Dimension getMinimumSize()
    {
        return getMaximumSize();
    }

    private void decorate( final Color c )
    {
        setOpaque( true );
        setBackground( c );
        
        Border one = BorderFactory.createEtchedBorder();
        Border two = BorderFactory.createMatteBorder( 2, 2, 2, 2, Color.WHITE );
        setBorder( BorderFactory.createCompoundBorder( one, two ) );
    }
    
    private void common_init( final String caption, final Color c )
    {
        decorate( c );
        setPreferredSize( new Dimension( 100, 40 ) );
        setText( caption );
        
        setVerticalAlignment( SwingConstants.CENTER );
        setHorizontalAlignment( SwingConstants.CENTER );

    }
    
    public ALabel( final String caption )
    {
        super();
        common_init( caption, Color.LIGHT_GRAY );
    }
    
    public ALabel( final String caption, final Color c )
    {
        super();
        common_init( caption, c );
    }
    
}
