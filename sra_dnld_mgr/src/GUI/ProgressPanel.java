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

import job.JobData;
import java.awt.*;
import java.text.NumberFormat;
import javax.swing.*;

public class ProgressPanel extends JPanel
{
    static final long serialVersionUID = 1;
    
    private final JLabel info;
    private final JProgressBar pro;
    private long progress, maximum;
    
    private String get_progress_string( long value, long maximum )
    {
        return String.format( "%s of %s",
                NumberFormat.getInstance().format( value ),
                NumberFormat.getInstance().format( maximum )
                );
    }
    
    public void update()
    {
        info.setText( get_progress_string( progress, maximum ) );
        int int_max = JobData.to_blocks( maximum );
        if ( pro.getMaximum() != int_max ) pro.setMaximum( int_max );
        pro.setValue( JobData.to_blocks( progress ) );
    }
    
    public void set_progress( final long value ) { progress = value; update(); }
    public void set_maximum( long value ) { maximum = value; update(); }
    public void start() { pro.setIndeterminate( true ); }
    public void stop() { pro.setIndeterminate( false ); }

    public ProgressPanel( JobData job )
    {
        setLayout( new BorderLayout( 0, 2 ) );
        setOpaque( false );
        
        info = new JLabel( get_progress_string( job.get_progress(), job.get_max() ) );
        add( info, BorderLayout.PAGE_START );
    
        pro = new JProgressBar();

        pro.setMaximum( JobData.to_blocks( job.get_max() ) );
        pro.setMinimum( 0 );
        pro.setValue( JobData.to_blocks( job.get_progress() ) );
        pro.setBorder( BorderFactory.createMatteBorder( 2, 2, 2, 2, Color.GRAY ) );
        pro.setStringPainted( true );
        
        add( pro, BorderLayout.CENTER );
    }
}
